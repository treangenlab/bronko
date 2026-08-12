use crate::call::*;
use crate::lcb::*;

use log::*;
use dashmap::DashMap;

pub const ANCHOR_MARGIN: usize = 2; //pull anchors this far off the breakpoint, coverage at the cliff is unreliable
pub const ANCHOR_DEPTH_CONSISTENCY: f64 = 0.50; //min/max depth across an anchor window must stay above this
pub const ANCHOR_SEARCH_MAX: usize = 50; //how far an anchor may slide looking for a clean window
pub const MAX_END_EXTENSION: usize = 100; //longest end walk, so a wrong turn cannot rewrite a long stretch

#[derive(Debug, Clone)]
pub struct Breakpoint {
    pub pos: usize,      // 1-based position on the low-coverage side of the depth change
    pub direction: i8,   // +1 for a rise in depth, -1 for a drop
    pub depth_high: u64, // larger total depth of the two consecutive positions the change spans
    pub depth_low: u64,  // smaller total depth of the two consecutive positions the change spans
}

// per-sequence breakpoint detection output, fed to the planning step
#[derive(Debug, Clone)]
pub struct SeqIndelInfo {
    pub seq: String,
    pub total: Vec<u64>,                      // total depth (fwd+rev across all bases) per position
    pub breakpoints: Vec<Breakpoint>,         // all detected breakpoints, in increasing position order
    pub pairs: Vec<(Breakpoint, Breakpoint)>, // drop -> following-rise pairs (candidate indels)
}

// a planned internal indel, carrying the anchors the planner validated so reconstruction reuses them
#[derive(Debug, Clone)]
pub struct IndelEvent {
    pub drop: Breakpoint,
    pub rise: Breakpoint,
    pub left_anchor: (usize, usize),  // 0-based inclusive
    pub right_anchor: (usize, usize), // 0-based inclusive
    pub left_depth: u64,              // min total depth across the left anchor
    pub right_depth: u64,             // min total depth across the right anchor
}

// which terminus an end reconstruction rebuilds toward
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EndSide {
    Five,  // 5' start: anchor walks left toward position 0
    Three, // 3' end: anchor walks right toward the terminus
}

// where to start an end walk when an indel sits too close to a terminus to bridge
#[derive(Debug, Clone)]
pub struct EndAnchor {
    pub seq: String,
    pub side: EndSide,
    pub anchor_first_base: usize, // 0-based first base of the anchor k-mer
}

// 0-based inclusive left anchor window for a drop at 1-based `drop_pos`; the -2 also converts to 0-based
fn left_anchor_window(drop_pos: usize, k: usize) -> Option<(usize, usize)> {
    let left_end = drop_pos.checked_sub(ANCHOR_MARGIN + 2)?;
    if left_end + 1 < k {
        return None;
    }
    Some((left_end + 1 - k, left_end))
}

// 0-based inclusive right anchor window for a rise at 1-based `rise_pos`
fn right_anchor_window(rise_pos: usize, k: usize, ref_len: usize) -> Option<(usize, usize)> {
    let start = (rise_pos - 1) + ANCHOR_MARGIN;
    if start + k > ref_len {
        return None;
    }
    Some((start, start + k - 1))
}

// slide the left anchor off the drop until it lands on a clean window; the fixed offset often
// straddles an unrelated coverage step
fn find_left_anchor(
    drop_pos: usize,
    k: usize,
    total: &[u64],
    bps: &[Breakpoint],
) -> Option<(usize, usize)> {
    let (_, end) = left_anchor_window(drop_pos, k)?;
    for back in 0..=ANCHOR_SEARCH_MAX {
        let e = end.checked_sub(back)?;
        if e + 1 < k {
            return None;
        }
        let w = (e + 1 - k, e);
        if window_is_clean(w, total, bps) {
            return Some(w);
        }
    }
    None
}

// same, sliding the right anchor further past the rise
fn find_right_anchor(
    rise_pos: usize,
    k: usize,
    ref_len: usize,
    total: &[u64],
    bps: &[Breakpoint],
) -> Option<(usize, usize)> {
    let (start, _) = right_anchor_window(rise_pos, k, ref_len)?;
    for fwd in 0..=ANCHOR_SEARCH_MAX {
        let s = start + fwd;
        if s + k > ref_len {
            return None;
        }
        let w = (s, s + k - 1);
        if window_is_clean(w, total, bps) {
            return Some(w);
        }
    }
    None
}

// does the cliff at 1-based `pos` overlap the window? the cliff spans positions [pos-2, pos-1]
fn cliff_in_window(pos: usize, window: (usize, usize)) -> bool {
    let a = pos.saturating_sub(2);
    let b = pos.saturating_sub(1);
    a <= window.1 && b >= window.0
}

// clean = inside the sequence, no breakpoint cliff inside it, and consistent depth across it
fn window_is_clean(window: (usize, usize), total: &[u64], breakpoints: &[Breakpoint]) -> bool {
    let (s, e) = window;
    if e >= total.len() || s > e {
        return false;
    }
    for bp in breakpoints {
        if cliff_in_window(bp.pos, window) {
            return false;
        }
    }
    let mut lo = u64::MAX;
    let mut hi = 0u64;
    for &d in &total[s..=e] {
        lo = lo.min(d);
        hi = hi.max(d);
    }
    if hi == 0 {
        return false;
    }
    (lo as f64 / hi as f64) >= ANCHOR_DEPTH_CONSISTENCY
}

// total depth (fwd+rev over all bases) per position
fn total_depth(fwd: &OutputData, rev: &OutputData) -> Vec<u64> {
    (0..fwd.counts.len())
        .map(|i| (0..4).map(|b| fwd.counts[i][b] + rev.counts[i][b]).sum())
        .collect()
}

// detect breakpoints (sharp depth changes between neighbors) and pair each drop with the next rise,
// which is the coverage signature of an indel. Returns the pairs grouped by sequence.
pub fn identify_indel_breakpoints(
    output: &DashMap<String, OutputData>,
    output_rev: &DashMap<String, OutputData>,
    max_ratio: f64,
    min_depth: u64,           // minimum total depth on the major side of the breakpoint
    max_depth: u64,           // maximum total depth on the minor side of the breakpoint
    max_pair_distance: usize, // maximum drop -> rise distance to pair as an indel
) -> Vec<SeqIndelInfo> {
    let mut grouped: Vec<SeqIndelInfo> = Vec::new();

    // loop through each sequence pileup for the genome
    for entry in output.iter() {
        let seq = entry.key();
        let fwd = entry.value();
        let rev = output_rev.get(seq).expect("Missing reverse strand for sequence");

        let len = fwd.counts.len();
        if len < 2 {
            continue;
        }

        // breakpoints found within this sequence, in increasing position order
        let mut breakpoints: Vec<Breakpoint> = Vec::new();

        let total = total_depth(fwd, &rev); // get the per-position depth for the sequence

        // flag cliffs between consecutive positions (low < max_ratio*high, within the depth bounds)
        for i in 1..len {
            let low = total[i - 1].min(total[i]);
            let high = total[i - 1].max(total[i]);
            if high < min_depth {
                continue;
            }
            if low > max_depth {
                continue;
            }
            if (low as f64 / high as f64) < max_ratio {
                breakpoints.push(Breakpoint {
                    pos: i + 1, // 1-based; the change spans total[i-1] -> total[i]
                    direction: if total[i] >= total[i - 1] { 1 } else { -1 },
                    depth_high: high,
                    depth_low: low,
                });
            }
        }

        trace!("{} breakpoints found on sequence {}: {:?}", breakpoints.len(), seq, breakpoints);

        if breakpoints.is_empty() {
            continue;
        }

        // pair a drop (negative) with the next rise (positive) within max_pair_distance
        let mut seq_pairs: Vec<(Breakpoint, Breakpoint)> = Vec::new();
        let mut pending_drop: Option<Breakpoint> = None;
        for bp in &breakpoints {
            match bp.direction {
                -1 => pending_drop = Some(bp.clone()), // most recent drop opens a potential indel
                1 => {
                    if let Some(drop) = pending_drop.take() {
                        if bp.pos.saturating_sub(drop.pos) <= max_pair_distance {
                            seq_pairs.push((drop, bp.clone()));
                        }
                    }
                }
                _ => {}
            }
        }

        grouped.push(SeqIndelInfo {
            seq: seq.clone(),
            total,
            breakpoints,
            pairs: seq_pairs,
        });
    }

    grouped
}

// turn breakpoint pairs into a reconstruction plan. Pairs too close to seed separately are merged,
// events reaching a terminus are handed to end reconstruction, and the rest become internal events with
// clean anchors. Returns (internal events by sequence, end anchors).
pub fn plan_indel_events(
    infos: &[SeqIndelInfo],
    k: usize,
    min_depth: u64,
    gap_max_depth: u64, // a position at or under this counts as unassembled, so worth rebuilding
) -> (Vec<(String, Vec<IndelEvent>)>, Vec<EndAnchor>) {
    let mut internal: Vec<(String, Vec<IndelEvent>)> = Vec::new();
    let mut ends: Vec<EndAnchor> = Vec::new();

    for info in infos {
        let total = &info.total;
        let ref_len = total.len();
        let bps = &info.breakpoints;

        // merge pairs with no room between them for a solid seed, since they are one disruption.
        // a merge only extends the rise rightward, so one left-to-right pass chains them
        let min_separation = k + ANCHOR_MARGIN;
        let mut merged: Vec<(Breakpoint, Breakpoint)> = Vec::new();
        for (drop, rise) in &info.pairs {
            let mut collide = false;
            if let Some((_, prev_rise)) = merged.last() {
                collide = drop.pos.saturating_sub(prev_rise.pos) < min_separation;
            }
            if collide {
                merged.last_mut().unwrap().1 = rise.clone();
            } else {
                merged.push((drop.clone(), rise.clone()));
            }
        }

        // ---- classify each merged event: internal bridge vs. end handoff ----
        let first_good = total.iter().position(|&d| d >= min_depth);
        let last_good = total.iter().rposition(|&d| d >= min_depth);

        let mut seq_internal: Vec<IndelEvent> = Vec::new();
        for (drop, rise) in merged {
            // only rebuild where the reads actually failed. Every called SNV makes a depth cliff as
            // its neighbours lose their deposits, but that dip still assembles, so require a position
            // with essentially no coverage
            let lo = drop.pos.saturating_sub(1).min(ref_len);
            let hi = rise.pos.min(ref_len);
            if !(lo..hi).any(|p| total[p] <= gap_max_depth) {
                continue;
            }

            // slide each anchor out to the nearest clean window rather than demanding the fixed offset
            let lw = find_left_anchor(drop.pos, k, total, bps);
            let rw = find_right_anchor(rise.pos, k, ref_len, total, bps);

            // an anchor running off the sequence or out of covered reference means this side reaches
            // a terminus, so hand it to end reconstruction instead of bridging
            let near_five = match lw {
                None => true,
                Some((s, _)) => first_good.map_or(true, |fg| s < fg),
            };
            let near_three = match rw {
                None => true,
                Some((_, e)) => last_good.map_or(true, |lg| e > lg),
            };

            if near_five {
                // rebuild the 5' end: anchor at the event's clean right anchor, walk left to the terminus
                if let Some((s, _)) = rw {
                    ends.push(EndAnchor { seq: info.seq.clone(), side: EndSide::Five, anchor_first_base: s });
                }
            }
            if near_three {
                // rebuild the 3' end: anchor at the event's clean left anchor, walk right to the terminus
                if let Some((s, _)) = lw {
                    ends.push(EndAnchor { seq: info.seq.clone(), side: EndSide::Three, anchor_first_base: s });
                }
            }
            if near_five || near_three {
                continue; // handled as an end (or un-anchorable); never bridged
            }

            // internal event: needs clean anchors on both sides. The bridge is gated on the anchors'
            // own depth, since a searched anchor can sit in quite different coverage
            if let (Some(left_anchor), Some(right_anchor)) = (lw, rw) {
                let depth_at = |(s, e): (usize, usize)| {
                    total[s..=e].iter().copied().min().unwrap_or(0)
                };
                seq_internal.push(IndelEvent {
                    left_depth: depth_at(left_anchor),
                    right_depth: depth_at(right_anchor),
                    drop,
                    rise,
                    left_anchor,
                    right_anchor,
                });
            }
        }

        if !seq_internal.is_empty() {
            internal.push((info.seq.clone(), seq_internal));
        }
    }

    (internal, ends)
}

#[derive(Debug, Clone)]
pub struct ReconstructedIndel {
    pub seq: String,
    pub drop_pos: usize,  // drop breakpoint position (1-based)
    pub rise_pos: usize,  // rise breakpoint position (1-based)
    pub ref_start: usize, // 1-based position of the first base of the left anchor k-mer
    pub ref_end: usize,   // 1-based position of the last base of the right anchor k-mer
    pub allele: Vec<u8>,  // reconstructed ASCII bases spanning [ref_start, ref_end] inclusive, spliced over that footprint
}

// both-strand read support for extending a (k-1)-mer by each base. The map is keyed on canonical buckets,
// so support sits in either the forward or the rc bucket (complemented base) and both are summed.
fn flanking_base_counts(
    km1: u64,
    k: usize,
    flanking_base_map: &DashMap<u64, [[u64; 4]; 2]>,
    extend_right: bool,
) -> [u64; 4] {
    // build a full k-mer with a placeholder base at the free position (the bucket is invariant to it)
    let full_fwd = if extend_right {
        (km1 << 2) & ((1u64 << (2 * k)) - 1) // placeholder at the last position
    } else {
        km1 // placeholder (0) at the first position; km1 occupies the low k-1 bases
    };
    let rc_full = reverse_complement_u64(full_fwd, k);

    let fwd_buckets = assign_buckets(full_fwd, k);
    let rc_buckets = assign_buckets(rc_full, k);

    // forward orientation uses the free-position bucket; rc orientation uses the mirrored bucket
    let (b_fwd, b_rc) = if extend_right {
        (fwd_buckets[k - 1], rc_buckets[0])
    } else {
        (fwd_buckets[0], rc_buckets[k - 1])
    };

    let fwd_row = flanking_base_map.get(&b_fwd).map(|r| *r); // copy the [[u64;4];2] out of the guard
    let rc_row = flanking_base_map.get(&b_rc).map(|r| *r);

    // sum both rc rows; only one of fwd_row/rc_row is non-empty so nothing is double counted
    let mut counts = [0u64; 4];
    for b in 0..4usize {
        let cf = fwd_row.map(|x| x[0][b] + x[1][b]).unwrap_or(0);
        let cr = rc_row.map(|x| x[0][b ^ 0b11] + x[1][b ^ 0b11]).unwrap_or(0); // complement of b on the rc bucket
        counts[b] = cf + cr;
    }
    counts
}

// best-supported extending base and its count; None if there's no support (last base wins on ties)
fn pick_major_base(counts: &[u64; 4]) -> Option<(u8, u64)> {
    let (b, &c) = counts.iter().enumerate().max_by_key(|(_, c)| **c)?;
    (c > 0).then_some((b as u8, c))
}

// bridge a drop and a rise by walking the de Bruijn graph from `left_kmer` to `right_kmer`, taking the
// major-supported base each step. Every base's support must sit inside the anchors' depth band, else the
// bridge is rejected. Returns the allele spanning the whole [left_kmer, right_kmer] footprint.
fn reconstruct_indel(
    left_kmer: &[u8],
    left_depth: u64,
    right_kmer: &[u8],
    right_depth: u64,
    ref_gap: usize, // reference bases between the two anchors, so the walk budget fits the geometry
    k: usize,
    flanking_base_map: &DashMap<u64, [[u64; 4]; 2]>,
) -> Option<Vec<u8>> {
    const MAX_BRIDGE: usize = 50; // largest novel (inserted) span we attempt to bridge
    const DEPTH_LO: f64 = 0.5;
    const DEPTH_HI: f64 = 2.0;

    if left_kmer.len() != k || right_kmer.len() != k || k < 2 {
        return None;
    }

    let km1 = k - 1;
    let mask_k = (1u64 << (2 * k)) - 1;
    let mask_km1 = (1u64 << (2 * km1)) - 1;

    // the walk crosses the reference gap, then overlaps the right anchor by k-1 before its [0] bucket
    // can match; MAX_BRIDGE is the slack for a novel insertion
    let max_ext = ref_gap + km1 + MAX_BRIDGE;

    let left_val = kmer_to_u64(left_kmer);
    let right_val = kmer_to_u64(right_kmer);

    // end bucket: marks the running k-mer overlapping the right anchor by k-1. assign_buckets[0] ignores
    // the top base, so `right_val >> 2` (placeholder up top, right's first k-1 bases below) yields it.
    let end_bucket = assign_buckets(right_val >> 2, k)[0];

    // pick the major base extending the (k-1) suffix `suffix` to the right; returns (base, support)
    let pick_base = |suffix: u64| -> Option<(u8, u64)> {
        let counts = flanking_base_counts(suffix, k, flanking_base_map, true);
        let (b, support) = pick_major_base(&counts)?;
        trace!(
            "picked base {} (support {}) extending suffix after getting counts {:?}",
            nucleotide_bits_to_char(b as u64),
            support,
            counts
        );
        Some((b, support))
    };

    // walk right from the left anchor, collecting (base, support), until the running k-mer overlaps the
    // right anchor by k-1 (its [0] bucket hits end_bucket)
    let mut cur = left_val;
    let mut bridge: Vec<(u8, u64)> = Vec::new();
    if assign_buckets(cur, k)[0] != end_bucket {
        let mut met = false;
        for _ in 0..max_ext {
            let (b, support) = pick_base(cur & mask_km1)?;
            cur = ((cur << 2) | b as u64) & mask_k;
            bridge.push((b, support));
            if assign_buckets(cur, k)[0] == end_bucket {
                met = true;
                break;
            }
        }
        if !met {
            return None; // never reached the right anchor within MAX_EXT
        }
    }

    // depth band over the whole bridge, bounded by the envelope of the two anchors; a per-side band
    // cannot describe a bridge crossing a real coverage step
    let lo = DEPTH_LO * left_depth.min(right_depth) as f64;
    let hi = DEPTH_HI * left_depth.max(right_depth) as f64;
    for &(_, support) in &bridge {
        if (support as f64) < lo || (support as f64) > hi {
            trace!("bridge support {} outside depth band [{:.0}, {:.0}]", support, lo, hi);
            return None;
        }
    }

    // assemble: left_kmer + bridging bases + right_kmer's final base (its leading k-1 already end the bridge)
    let mut allele: Vec<u8> = left_kmer.to_vec();
    for &(b, _) in &bridge {
        allele.push(nucleotide_bits_to_char(b as u64) as u8);
    }
    allele.push(right_kmer[k - 1]);

    trace!(
        "reconstructed indel bridge of {} bases: {}",
        bridge.len(),
        String::from_utf8_lossy(&allele)
    );

    Some(allele)
}

// reconstruct an allele for each event from plan_indel_events, using the anchors the planner validated
pub fn reconstruct_indels(
    output: &DashMap<String, OutputData>,
    flanking_base_map: &DashMap<u64, [[u64; 4]; 2]>,
    grouped_pairs: &[(String, Vec<IndelEvent>)],
    k: usize,
) -> Vec<ReconstructedIndel> {
    let mut results: Vec<ReconstructedIndel> = Vec::new();

    for (seq, events) in grouped_pairs {
        let ref_entry = match output.get(seq) {
            Some(e) => e,
            None => continue,
        };
        let ref_bases = &ref_entry.ref_bases;

        for event in events {
            let (drop, rise) = (&event.drop, &event.rise);
            trace!(
                "attempting to reconstruct indel between drop at {} and rise at {} on sequence {}",
                drop.pos,
                rise.pos,
                seq
            );
            let (left_start, left_end) = event.left_anchor;
            let (right_start, right_end) = event.right_anchor;
            if left_end >= ref_bases.len() || right_end >= ref_bases.len() {
                continue;
            }

            let left_kmer = &ref_bases[left_start..=left_end];
            let right_kmer = &ref_bases[right_start..=right_end];

            // anchors spanning ambiguous bases cannot be encoded, so skip the event
            if !is_acgt(left_kmer) || !is_acgt(right_kmer) {
                continue;
            }

            trace!(
                "reconstructing indel between drop at {} and rise at {} on sequence {} with left k-mer {} and right k-mer {}",
                left_end + 1, //1-based
                right_start + 1, //1-based
                seq,
                String::from_utf8_lossy(left_kmer),
                String::from_utf8_lossy(right_kmer)
            );

            if let Some(allele) = reconstruct_indel(
                left_kmer,
                event.left_depth,
                right_kmer,
                event.right_depth,
                right_start.saturating_sub(left_end + 1), // reference bases between the anchors
                k,
                flanking_base_map,
            ) {
                results.push(ReconstructedIndel {
                    seq: seq.clone(),
                    drop_pos: drop.pos,
                    rise_pos: rise.pos,
                    ref_start: left_start + 1, // 1-based footprint the allele replaces
                    ref_end: right_end + 1,
                    allele,
                });
            }
        }
    }

    results
}

// walk out from a single anchor k-mer to rebuild a sequence end, taking the major-supported base each
// step. With no second anchor to meet it runs to `max_extension` or a dead end. Returns bases in
// left-to-right reference order.
fn reconstruct_end(
    anchor_kmer: &[u8],
    extend_right: bool,
    max_extension: usize,
    k: usize,
    flanking_base_map: &DashMap<u64, [[u64; 4]; 2]>,
) -> Vec<u8> {
    if anchor_kmer.len() != k || k < 2 {
        return Vec::new();
    }

    let mask_k = (1u64 << (2 * k)) - 1;
    let mask_km1 = (1u64 << (2 * (k - 1))) - 1;
    let shift_km1 = 2 * (k - 1);

    let mut cur = kmer_to_u64(anchor_kmer);
    let mut bases: Vec<u8> = Vec::new();

    for _ in 0..max_extension {
        // the (k-1)-mer the next base attaches to: trailing suffix when extending right, leading prefix
        // when extending left
        let km1 = if extend_right { cur & mask_km1 } else { cur >> 2 };
        let counts = flanking_base_counts(km1, k, flanking_base_map, extend_right);

        // pick the best-supported extending base, or stop at a dead end with no read evidence
        let b = match pick_major_base(&counts) {
            Some((b, _)) => b as u64,
            None => break,
        };
        if extend_right {
            cur = ((cur << 2) | b) & mask_k;
        } else {
            cur = (b << shift_km1) | (cur >> 2);
        }
        bases.push(nucleotide_bits_to_char(b) as u8);
    }

    if !extend_right {
        bases.reverse(); // collected nearest-anchor-first; flip to left-to-right reference order
    }
    bases
}

// rebuild each sequence's 5' and 3' ends: anchor just inside the covered region and walk outward to the
// terminus. Each end returns a length-preserving ReconstructedIndel so it splices through the same
// machinery without shifting coordinates, and it only improves the consensus -- the caller never reports
// it as a variant. `end_anchors` overrides where a walk starts when the planner handed one off.
pub fn reconstruct_ends(
    output: &DashMap<String, OutputData>,
    output_rev: &DashMap<String, OutputData>,
    flanking_base_map: &DashMap<u64, [[u64; 4]; 2]>,
    k: usize,
    min_depth: u64,
    end_anchors: &[EndAnchor],
) -> Vec<ReconstructedIndel> {
    const END_ANCHOR_MARGIN: usize = ANCHOR_MARGIN; // pull the anchor this far inside the confidently-covered region

    let mut results: Vec<ReconstructedIndel> = Vec::new();

    for entry in output.iter() {
        let seq = entry.key();
        let fwd = entry.value();
        let rev = match output_rev.get(seq) {
            Some(r) => r,
            None => continue,
        };
        let ref_bases = &fwd.ref_bases;
        let len = ref_bases.len();
        if len < 2 * k {
            continue; // too short to anchor and reconstruct meaningfully
        }

        let total = total_depth(fwd, &rev);

        let first_good = total.iter().position(|&d| d >= min_depth);
        let last_good = total.iter().rposition(|&d| d >= min_depth);

        // breakpoint-derived anchor overrides for this sequence, if the planner handed any off
        let override_five = end_anchors.iter()
            .find(|e| e.seq == *seq && e.side == EndSide::Five)
            .map(|e| e.anchor_first_base);
        let override_three = end_anchors.iter()
            .find(|e| e.seq == *seq && e.side == EndSide::Three)
            .map(|e| e.anchor_first_base);

        // anchor just inside the first/last covered position (or at the planner's override) and walk out
        let anchor_start_5 = override_five.or_else(|| {
            first_good.and_then(|fg| (fg > 0).then(|| fg + END_ANCHOR_MARGIN))
        });
        let anchor_start_3 = override_three.or_else(|| {
            last_good.and_then(|lg| {
                if lg + 1 < len && lg >= END_ANCHOR_MARGIN {
                    let anchor_end = lg - END_ANCHOR_MARGIN;
                    (anchor_end + 1 >= k).then(|| anchor_end + 1 - k)
                } else {
                    None
                }
            })
        });

        for (side, anchor_start) in [(EndSide::Five, anchor_start_5), (EndSide::Three, anchor_start_3)] {
            let anchor_start = match anchor_start {
                Some(a) => a,
                None => continue,
            };
            let anchor_end = anchor_start + k - 1; // 0-based last base
            if anchor_end >= len || !is_acgt(&ref_bases[anchor_start..=anchor_end]) {
                continue;
            }
            let anchor_kmer = &ref_bases[anchor_start..=anchor_end];

            // room to the terminus, capped so a wrong turn cannot rewrite a long stretch
            let (extend_right, room) = match side {
                EndSide::Five => (false, anchor_start),
                EndSide::Three => (true, (len - 1) - anchor_end),
            };
            let walked = reconstruct_end(anchor_kmer, extend_right, room.min(MAX_END_EXTENSION), k, flanking_base_map);
            if walked.is_empty() {
                continue;
            }

            // footprint spans exactly the rebuilt region plus the anchor
            let n = walked.len();
            let (drop_pos, rise_pos, ref_start, ref_end, allele) = match side {
                EndSide::Five => {
                    let mut allele = walked;
                    allele.extend_from_slice(anchor_kmer);
                    (1, anchor_start + 1, anchor_start - n + 1, anchor_end + 1, allele)
                }
                EndSide::Three => {
                    let mut allele = anchor_kmer.to_vec();
                    allele.extend_from_slice(&walked);
                    (anchor_end + 1, len, anchor_start + 1, anchor_end + n + 1, allele)
                }
            };

            results.push(ReconstructedIndel {
                seq: seq.clone(),
                drop_pos,
                rise_pos,
                ref_start,
                ref_end,
                allele,
            });
        }
    }

    results
}

// sorted, validated, non-overlapping indels for `seq_name` as (start_0based, end_exclusive, indel).
// malformed footprints and any overlapping an already-kept one are dropped.
fn ordered_indels<'a>(
    seq_name: &str,
    indels: &'a [ReconstructedIndel],
    len: usize,
) -> Vec<(usize, usize, &'a ReconstructedIndel)> {
    let mut seq_indels: Vec<&ReconstructedIndel> =
        indels.iter().filter(|i| i.seq == seq_name).collect();
    seq_indels.sort_by_key(|i| i.ref_start);

    let mut out: Vec<(usize, usize, &ReconstructedIndel)> = Vec::new();
    let mut cursor = 0usize; // next 0-based base after the last kept footprint
    for indel in seq_indels {
        if indel.ref_start < 1 || indel.ref_end > len || indel.ref_start > indel.ref_end {
            continue; // malformed footprint
        }
        let start = indel.ref_start - 1;
        let end = indel.ref_end; // exclusive
        if start < cursor {
            continue; // overlaps an already-kept indel; skip (non-overlapping assumed)
        }
        out.push((start, end, indel));
        cursor = end;
    }
    out
}

// same length as its footprint means a block substitution, not an indel: it patches in place and
// leaves coordinates untouched
fn is_block_substitution(indel: &ReconstructedIndel, len: usize) -> bool {
    if indel.ref_start < 1 || indel.ref_end > len || indel.ref_start > indel.ref_end {
        return false;
    }
    indel.allele.len() == indel.ref_end - (indel.ref_start - 1)
}

// splice the indels for `seq_name` into `bases` (same length and coordinates), returning the rebuilt sequence
pub fn splice_indels(bases: &[u8], seq_name: &str, indels: &[ReconstructedIndel]) -> Vec<u8> {
    // block substitutions only fill gaps, so overlapping ones can all apply; only length-changing
    // indels need the non-overlap discipline in ordered_indels
    let mut patched = bases.to_vec();
    for ind in indels.iter().filter(|i| i.seq == seq_name && is_block_substitution(i, bases.len())) {
        let start = ind.ref_start - 1;
        for (off, &b) in ind.allele.iter().enumerate() {
            if patched[start + off] == b'N' {
                patched[start + off] = b;
            }
        }
    }

    let real: Vec<ReconstructedIndel> = indels
        .iter()
        .filter(|i| i.seq == seq_name && !is_block_substitution(i, bases.len()))
        .cloned()
        .collect();

    let mut new_seq: Vec<u8> = Vec::with_capacity(patched.len());
    let mut cursor = 0usize; // 0-based index of the next base to copy
    for (start, end, indel) in ordered_indels(seq_name, &real, patched.len()) {
        new_seq.extend_from_slice(&patched[cursor..start]); // bases up to the footprint
        new_seq.extend_from_slice(&indel.allele); // a real length change: replace the footprint
        cursor = end;
    }
    new_seq.extend_from_slice(&patched[cursor..]); // remaining tail

    new_seq
}

// Splice every reconstructed indel into the selected genome's reference sequence(s), returning the
// rebuilt sequences as (name, bases). Sequences with no indels pass through unchanged.
pub fn apply_indels_to_reference(
    output: &DashMap<String, OutputData>,
    indels: &[ReconstructedIndel],
) -> Vec<(String, Vec<u8>)> {
    let mut rebuilt: Vec<(String, Vec<u8>)> = Vec::new();
    for entry in output.iter() {
        let seq = entry.key();
        let new_seq = splice_indels(&entry.value().ref_bases, seq, indels);
        rebuilt.push((seq.clone(), new_seq));
    }
    rebuilt
}

// coordinate translation for one sequence between an indel-corrected sequence (consensus + spliced
// indels) and the original reference, all 0-based.
//   - `corrected_to_orig[i]` = the original position corrected `i` came from, or None if `i` is inside a
//     spliced-in allele.
//   - `orig_to_corrected[j]` = the corrected position original base `j` maps to, or None if `j` fell
//     inside a replaced indel footprint.
#[derive(Debug, Clone)]
pub struct CoordMap {
    pub corrected_to_orig: Vec<Option<usize>>,
    pub orig_to_corrected: Vec<Option<usize>>,
}

// coordinate map for `seq_name`, mirroring splice_indels so the map matches the spliced sequence
pub fn build_coord_map(ref_len: usize, seq_name: &str, indels: &[ReconstructedIndel]) -> CoordMap {
    let mut corrected_to_orig: Vec<Option<usize>> = Vec::with_capacity(ref_len);
    let mut orig_to_corrected: Vec<Option<usize>> = vec![None; ref_len];

    // only length-changing indels shift coordinates, so filter as splice_indels does
    let real: Vec<ReconstructedIndel> = indels
        .iter()
        .filter(|i| i.seq == seq_name && !is_block_substitution(i, ref_len))
        .cloned()
        .collect();

    let mut cursor = 0usize; // next original 0-based base to copy
    for (start, end, indel) in ordered_indels(seq_name, &real, ref_len) {
        // copied run [cursor, start): 1:1 in both directions
        for orig in cursor..start {
            orig_to_corrected[orig] = Some(corrected_to_orig.len());
            corrected_to_orig.push(Some(orig));
        }
        // allele span: corrected positions with no original coordinate
        for _ in 0..indel.allele.len() {
            corrected_to_orig.push(None);
        }
        // footprint [start, end) originals already left as None in orig_to_corrected
        cursor = end;
    }
    // tail [cursor, ref_len): 1:1
    for orig in cursor..ref_len {
        orig_to_corrected[orig] = Some(corrected_to_orig.len());
        corrected_to_orig.push(Some(orig));
    }

    CoordMap { corrected_to_orig, orig_to_corrected }
}