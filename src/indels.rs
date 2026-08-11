use crate::call::*;
use crate::lcb::*;

use log::*;
use dashmap::DashMap;

pub const ANCHOR_MARGIN: usize = 2; // pull anchor k-mers this many bases off the breakpoint, since coverage right at the cliff is unreliable
pub const ANCHOR_DEPTH_CONSISTENCY: f64 = 0.50; // an anchor window is only trusted if min/max total depth across it stays >= this fraction

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

// which terminus an end reconstruction rebuilds toward
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EndSide {
    Five,  // 5' start: anchor walks left toward position 0
    Three, // 3' end: anchor walks right toward the terminus
}

// an end-reconstruction anchor for when an indel sits too close to a terminus to place a clean outer
// anchor; its clean inner anchor starts the end walk instead of the generic first_good/last_good anchor
#[derive(Debug, Clone)]
pub struct EndAnchor {
    pub seq: String,
    pub side: EndSide,
    pub anchor_first_base: usize, // 0-based first base of the anchor k-mer
}

// 0-based inclusive window of the left anchor k-mer for a drop at 1-based `drop_pos`. The anchor ends
// ANCHOR_MARGIN bases before the drop (the -2 also converts 1-based to 0-based and steps off the cliff).
// None if it runs off the start of the sequence.
fn left_anchor_window(drop_pos: usize, k: usize) -> Option<(usize, usize)> {
    let left_end = drop_pos.checked_sub(ANCHOR_MARGIN + 2)?;
    if left_end + 1 < k {
        return None;
    }
    Some((left_end + 1 - k, left_end))
}

// 0-based inclusive window of the right anchor k-mer for a rise at 1-based `rise_pos`. The anchor starts
// ANCHOR_MARGIN bases past the rise. None if it runs off the end of the sequence.
fn right_anchor_window(rise_pos: usize, k: usize, ref_len: usize) -> Option<(usize, usize)> {
    let start = (rise_pos - 1) + ANCHOR_MARGIN;
    if start + k > ref_len {
        return None;
    }
    Some((start, start + k - 1))
}

// does the cliff of the breakpoint at 1-based `pos` overlap the 0-based inclusive window? The cliff spans
// total[i-1] -> total[i], i.e. positions [pos-2, pos-1], so test that interval against the window.
fn cliff_in_window(pos: usize, window: (usize, usize)) -> bool {
    let a = pos.saturating_sub(2);
    let b = pos.saturating_sub(1);
    a <= window.1 && b >= window.0
}

// a window is clean if it lies inside the sequence, holds no other breakpoint's cliff, and keeps
// consistent depth (min/max across the window >= ANCHOR_DEPTH_CONSISTENCY)
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

// detect indel breakpoints in the genome's pileup and pair them into candidate indels. A breakpoint is a
// sharp total-depth change between neighbors (min/max < max_ratio); a drop followed by a nearby rise is
// the coverage signature of an indel. Returns the (drop, rise) pairs grouped by sequence name.
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

// turn the breakpoint pairs into a reconstruction plan, making sure every anchor k-mer sits over clean,
// consistently-covered reference settings. Neighboring pairs with colliding anchor windows are merged into one
// event; events too close to a terminus to anchor are handed off to end reconstruction, and the rest
// become internal pairs (events left unclean by a stray breakpoint or depth wobble are dropped).
// Returns (internal pairs grouped by sequence, end anchors) for reconstruct_indels / reconstruct_ends.
pub fn plan_indel_events(
    infos: &[SeqIndelInfo],
    k: usize,
    min_depth: u64,
) -> (Vec<(String, Vec<(Breakpoint, Breakpoint)>)>, Vec<EndAnchor>) {
    let mut internal: Vec<(String, Vec<(Breakpoint, Breakpoint)>)> = Vec::new();
    let mut ends: Vec<EndAnchor> = Vec::new();

    for info in infos {
        let total = &info.total;
        let ref_len = total.len();
        let bps = &info.breakpoints;

        // ---- merge pass: collapse neighboring pairs whose anchor windows collide ----
        // a merge only extends a pair's rise rightward, so one left-to-right pass chains colliding pairs
        let mut merged: Vec<(Breakpoint, Breakpoint)> = Vec::new();
        for (drop, rise) in &info.pairs {
            // copy prev rise pos out to release the immutable borrow of `merged` before the mutable one
            let mut collide = false;
            if let Some((_, prev_rise)) = merged.last() {
                let prev_rise_pos = prev_rise.pos;
                let rw = right_anchor_window(prev_rise_pos, k, ref_len);
                let lw = left_anchor_window(drop.pos, k);
                collide = match (rw, lw) {
                    (Some(rw), Some(lw)) => {
                        let overlap = rw.0 <= lw.1 && lw.0 <= rw.1;
                        overlap
                            || cliff_in_window(drop.pos, rw)    // this pair's drop sits in prev's right anchor
                            || cliff_in_window(prev_rise_pos, lw) // prev's rise sits in this pair's left anchor
                    }
                    // a missing window means no room for a clean anchor between them -> merge
                    _ => true,
                };
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

        let mut seq_internal: Vec<(Breakpoint, Breakpoint)> = Vec::new();
        for (drop, rise) in merged {
            let lw = left_anchor_window(drop.pos, k);
            let rw = right_anchor_window(rise.pos, k, ref_len);

            // an outer anchor that runs off the sequence or out of covered reference means this side
            // reaches a terminus: hand off to end reconstruction instead of bridging
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
                    if window_is_clean((s, s + k - 1), total, bps) {
                        ends.push(EndAnchor { seq: info.seq.clone(), side: EndSide::Five, anchor_first_base: s });
                    }
                }
            }
            if near_three {
                // rebuild the 3' end: anchor at the event's clean left anchor, walk right to the terminus
                if let Some((s, e)) = lw {
                    if window_is_clean((s, e), total, bps) {
                        ends.push(EndAnchor { seq: info.seq.clone(), side: EndSide::Three, anchor_first_base: s });
                    }
                }
            }
            if near_five || near_three {
                continue; // handled as an end (or un-anchorable); never bridged
            }

            // internal event: both outer anchors must be clean reference, else drop it
            match (lw, rw) {
                (Some(lw), Some(rw))
                    if window_is_clean(lw, total, bps) && window_is_clean(rw, total, bps) =>
                {
                    seq_internal.push((drop, rise));
                }
                _ => {}
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

// sum both-strand read support for extending a (k-1)-mer by each base. Since the map is keyed on canonical
// buckets, support can live in the forward bucket or the reverse-complement bucket (complemented base), so
// both are summed. `extend_right` extends to the right, else to the left. Returns per-base [A,C,G,T] counts.
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

    // both strands of a canonical k-mer share a bucket split across the two rc rows, so sum both rows;
    // only one of fwd_row/rc_row is non-empty, so there's no double counting
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

// reconstruct the allele bridging a drop and a rise by walking the de Bruijn graph left-to-right from
// `left_kmer` to `right_kmer`, taking the major-supported base each step until the running k-mer overlaps
// the right anchor. Each base's support must stay within [DEPTH_LO, DEPTH_HI] of the nearer anchor depth
// (left_depth over the first half of the bridge, right_depth over the second), else the bridge is
// rejected. Gives up after MAX_BRIDGE bases or at a dead end. The returned allele spans the full
// [left_kmer, right_kmer] footprint, ready to splice back in.
fn reconstruct_indel(
    left_kmer: &[u8],
    left_depth: u64,
    right_kmer: &[u8],
    right_depth: u64,
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

    // the walk must cover the right anchor's k-1 overlap before cur's [0] bucket can match, so cap the
    // extension at that overlap plus the largest novel bridge we allow
    let max_ext = km1 + MAX_BRIDGE;

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

    // depth band, per side: first half of the bridge gated by left_depth, second half by right_depth
    let l = bridge.len();
    for (i, &(_, support)) in bridge.iter().enumerate() {
        let depth = if i * 2 < l { left_depth } else { right_depth };
        let lo = DEPTH_LO * depth as f64;
        let hi = DEPTH_HI * depth as f64;
        if (support as f64) < lo || (support as f64) > hi {
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

// reconstruct an allele for each anchor-validated pair from plan_indel_events. The reference is looked up
// by name from the genome's output table; the geometry helpers re-derive the anchor windows here.
pub fn reconstruct_indels(
    output: &DashMap<String, OutputData>,
    flanking_base_map: &DashMap<u64, [[u64; 4]; 2]>,
    grouped_pairs: &[(String, Vec<(Breakpoint, Breakpoint)>)],
    k: usize,
) -> Vec<ReconstructedIndel> {
    let mut results: Vec<ReconstructedIndel> = Vec::new();

    for (seq, pairs) in grouped_pairs {
        let ref_entry = match output.get(seq) {
            Some(e) => e,
            None => continue,
        };
        let ref_bases = &ref_entry.ref_bases;


        for (drop, rise) in pairs {
            // left k-mer ends ANCHOR_MARGIN bases before the last covered base ahead of the drop;
            // right k-mer starts ANCHOR_MARGIN bases after the first covered base past the rise
            trace!(
                "attempting to reconstruct indel between drop at {} and rise at {} on sequence {}",
                drop.pos,
                rise.pos,
                seq
            );
            let (left_start, left_end) = match left_anchor_window(drop.pos, k) {
                Some(w) => w,
                None => continue,
            };
            let (right_start, right_end) = match right_anchor_window(rise.pos, k, ref_bases.len()) {
                Some(w) => w,
                None => continue,
            };

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
                drop.depth_high,
                right_kmer,
                rise.depth_high,
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

// walk the de Bruijn graph outward from a single anchor k-mer, taking the major-supported base, to rebuild
// a sequence end. `extend_right` walks toward the 3' terminus (appending), else toward the 5' start
// (prepending). With no second anchor to meet, it runs until `max_extension` bases (the caller caps this
// at the reference length) or a dead end. Returned bases are in left-to-right reference order.
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

// reconstruct the 5' and/or 3' ends of each sequence: anchor a k-mer just inside the covered region (the
// first/last position reaching `min_depth`, pulled in by END_ANCHOR_MARGIN) and walk outward to the
// terminus, capped at the reference boundary. Each end comes back as a length-preserving ReconstructedIndel
// so it splices through the same machinery, but it only improves the consensus -- the caller does not
// report it as a variant, and coordinates stay unchanged. Unreached terminal bases are left as the
// consensus produced them (N without coverage). `end_anchors` lets the planner override where an end walk
// begins when an indel sits too close to a terminus to bridge; otherwise the generic anchor is used.
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

        // ---- 5' start: anchor just inside the first covered position (or at the override), walk left ----
        // 0-based first base of the anchor k-mer: the planner's override wins, else the generic first_good anchor
        let anchor_start_5 = override_five.or_else(|| {
            first_good.and_then(|fg| (fg > 0).then(|| fg + END_ANCHOR_MARGIN))
        });
        if let Some(anchor_start) = anchor_start_5 {
            let anchor_end = anchor_start + k - 1; // 0-based last base
            if anchor_end < len && is_acgt(&ref_bases[anchor_start..=anchor_end]) {
                let anchor_kmer = &ref_bases[anchor_start..=anchor_end];
                // walk left at most `anchor_start` bases (reaching position 0 = the reference start)
                let head = reconstruct_end(anchor_kmer, false, anchor_start, k, flanking_base_map);
                if !head.is_empty() {
                    // footprint exactly spans the rebuilt region [anchor_start - head.len(), anchor_end]
                    let ref_start = anchor_start - head.len() + 1; // 1-based
                    let mut allele = head;
                    allele.extend_from_slice(anchor_kmer);
                    results.push(ReconstructedIndel {
                        seq: seq.clone(),
                        drop_pos: 1,
                        rise_pos: anchor_start + 1,
                        ref_start,
                        ref_end: anchor_end + 1,
                        allele,
                    });
                }
            }
        }

        // ---- 3' end: anchor just inside the last covered position (or at the override), walk right ----
        // 0-based first base of the anchor k-mer: the planner's override wins, else the generic last_good anchor
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
        if let Some(anchor_start) = anchor_start_3 {
            let anchor_end = anchor_start + k - 1; // 0-based last base
            if anchor_end < len && is_acgt(&ref_bases[anchor_start..=anchor_end]) {
                let anchor_kmer = &ref_bases[anchor_start..=anchor_end];
                let max_ext = (len - 1) - anchor_end; // reach the last reference position at most
                let tail = reconstruct_end(anchor_kmer, true, max_ext, k, flanking_base_map);
                if !tail.is_empty() {
                    // footprint exactly spans the rebuilt region [anchor_start, anchor_end + tail.len()]
                    let ref_end = anchor_end + tail.len() + 1; // 1-based
                    let mut allele = anchor_kmer.to_vec();
                    allele.extend_from_slice(&tail);
                    results.push(ReconstructedIndel {
                        seq: seq.clone(),
                        drop_pos: anchor_end + 1,
                        rise_pos: len,
                        ref_start: anchor_start + 1,
                        ref_end,
                        allele,
                    });
                }
            }
        }
    }

    results
}

// sorted, validated, non-overlapping indels for `seq_name`, each as (start_0based, end_exclusive, indel).
// 1-based inclusive [ref_start, ref_end] footprints are converted to 0-based [start, end); malformed
// footprints and any overlapping an already-kept one are dropped so callers stay consistent with each other.
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

// splice the indels for `seq_name` into a coordinate-matched buffer `bases` (reference or per-position
// consensus, same length and coordinates), returning the rebuilt sequence. Each `allele` replaces its
// footprint; see ordered_indels for the validation/overlap handling.
pub fn splice_indels(bases: &[u8], seq_name: &str, indels: &[ReconstructedIndel]) -> Vec<u8> {
    let mut new_seq: Vec<u8> = Vec::with_capacity(bases.len());
    let mut cursor = 0usize; // 0-based index of the next base to copy
    for (start, end, indel) in ordered_indels(seq_name, indels, bases.len()) {
        new_seq.extend_from_slice(&bases[cursor..start]); // bases up to the footprint
        new_seq.extend_from_slice(&indel.allele); // the reconstructed allele
        cursor = end;
    }
    new_seq.extend_from_slice(&bases[cursor..]); // remaining tail

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

// build the coordinate map for `seq_name`; `ref_len` is the original reference length. Mirrors
// `splice_indels`: copied runs map 1:1 and overlapping/malformed indels are skipped identically, so the
// map stays consistent with the spliced sequence.
pub fn build_coord_map(ref_len: usize, seq_name: &str, indels: &[ReconstructedIndel]) -> CoordMap {
    let mut corrected_to_orig: Vec<Option<usize>> = Vec::with_capacity(ref_len);
    let mut orig_to_corrected: Vec<Option<usize>> = vec![None; ref_len];

    let mut cursor = 0usize; // next original 0-based base to copy
    for (start, end, indel) in ordered_indels(seq_name, indels, ref_len) {
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