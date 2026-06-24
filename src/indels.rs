use crate::call::OutputData;
use crate::lcb::*;

use log::*;
use dashmap::DashMap;

// Pull each anchor k-mer this many bases further from the breakpoint, into more confidently-covered
// reference, since coverage (and base calls) right at the cliff can be unreliable. Shared by internal
// indel reconstruction and (as the end-anchor margin) by end reconstruction so the geometry matches.
pub const ANCHOR_MARGIN: usize = 2;

// An anchor k-mer window is only trusted if its total depth stays consistent across the whole window:
// the smaller depth in the window must be at least this fraction of the larger. Catches gradual depth
// drift that never trips a single-step breakpoint cliff.
const ANCHOR_DEPTH_CONSISTENCY: f64 = 0.5;

#[derive(Debug, Clone)]
pub struct Breakpoint {
    pub pos: usize,      // 1-based position on the low-coverage side of the depth change
    pub direction: i8,   // +1 for a rise in depth, -1 for a drop
    pub depth_high: u64, // larger total depth of the two consecutive positions the change spans
    pub depth_low: u64,  // smaller total depth of the two consecutive positions the change spans
}

// Per-sequence output of breakpoint detection: the full per-position total depth and every detected
// breakpoint (not just the paired ones), so the planning step can test candidate anchor windows for
// depth consistency and for collisions with neighboring variants.
#[derive(Debug, Clone)]
pub struct SeqIndelInfo {
    pub seq: String,
    pub total: Vec<u64>,                      // total depth (fwd+rev across all bases) per position
    pub breakpoints: Vec<Breakpoint>,         // all detected breakpoints, in increasing position order
    pub pairs: Vec<(Breakpoint, Breakpoint)>, // drop -> following-rise pairs (candidate indels)
}

// Which terminus an end reconstruction rebuilds toward.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EndSide {
    Five,  // 5' start: anchor walks left toward position 0
    Three, // 3' end: anchor walks right toward the terminus
}

// An explicit end-reconstruction anchor, produced when an indel event sits too close to a terminus to
// place a clean outer anchor. The breakpoint's clean inner anchor becomes the start of the end walk,
// replacing the generic first_good/last_good anchor on that side.
#[derive(Debug, Clone)]
pub struct EndAnchor {
    pub seq: String,
    pub side: EndSide,
    pub anchor_first_base: usize, // 0-based first base of the anchor k-mer
}

// 0-based inclusive window [start, end] of the left anchor k-mer for a drop at 1-based `drop_pos`.
// Mirrors the original geometry in reconstruct_indels: the anchor ends ANCHOR_MARGIN bases before the
// last confidently-covered base ahead of the drop (the extra -2 converts 1-based to 0-based and steps
// off the cliff's low-coverage base). None if it would run off the start of the sequence.
fn left_anchor_window(drop_pos: usize, k: usize) -> Option<(usize, usize)> {
    let left_end = drop_pos.checked_sub(ANCHOR_MARGIN + 2)?;
    if left_end + 1 < k {
        return None;
    }
    Some((left_end + 1 - k, left_end))
}

// 0-based inclusive window [start, end] of the right anchor k-mer for a rise at 1-based `rise_pos`.
// The anchor starts ANCHOR_MARGIN bases into the confidently-covered region past the rise. None if it
// would run off the end of the sequence.
fn right_anchor_window(rise_pos: usize, k: usize, ref_len: usize) -> Option<(usize, usize)> {
    let start = (rise_pos - 1) + ANCHOR_MARGIN;
    if start + k > ref_len {
        return None;
    }
    Some((start, start + k - 1))
}

// Does the depth cliff of the breakpoint at 1-based `pos` fall inside the 0-based inclusive window?
// A breakpoint with pos = i+1 spans the two consecutive positions total[i-1] -> total[i] (0-based
// i-1 and i), so treat it as occupying [pos-2, pos-1] and test for interval overlap with the window.
fn cliff_in_window(pos: usize, window: (usize, usize)) -> bool {
    let a = pos.saturating_sub(2);
    let b = pos.saturating_sub(1);
    a <= window.1 && b >= window.0
}

// An anchor window is clean if it lies inside the sequence, contains no other breakpoint's cliff, and
// holds consistent depth (smaller/larger total depth across the window >= ANCHOR_DEPTH_CONSISTENCY).
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

// Detection of potential indel breakpoints in the selected genome's pileup, then pairs them into
// candidate indels. We walk consecutive positions and flag a breakpoint wherever total depth (fwd+rev)
// changes sharply between neighbors: min(d1,d2)/max(d1,d2) < max_ratio, i.e. the smaller side falls under
// `max_ratio` of the larger. The high-coverage side must also reach `min_depth` to avoid low-coverage
// noise. Flagged breakpoints are then paired as a drop (negative) followed by a rise (positive) within a
// reasonable distance, which is the coverage signature of an indel. Returns the (drop, rise) pairs grouped
// by sequence name, so multi-sequence genomes keep each pair associated with the sequence it was found on.
pub fn identify_indel_breakpoints(
    output: &DashMap<String, OutputData>,
    output_rev: &DashMap<String, OutputData>,
    max_ratio: f64,
    min_depth: u64, // this is the minimum depth of the major side of the breakpoint (aka if this is set to 1000, then the major side must have at least 1000 total depth for the breakpoint to be considered, regardless of the depth after drop)
    max_depth: u64, // this is the maximum depth of the minor side of the breakpoint (aka if this is set to 100, then the depth after the drop must be less than 100 for the breakpoint to be considered)
    max_pair_distance: usize, // this is the maximum distance between a drop and the following rise to be paired as an indel
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

        // total depth (fwd+rev over all bases) per position
        let mut total = vec![0u64; len];
        for i in 0..len {
            let mut sum = 0u64;
            for b in 0..4 {
                sum += fwd.counts[i][b] + rev.counts[i][b];
            }
            total[i] = sum;
        }

        // flag a sharp cliff between consecutive positions: the low side is < max_ratio of the high side; high side is >= min_depth; low side is <= max_depth
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

// Turn the detected breakpoint pairs into a concrete reconstruction plan, enforcing that every anchor
// k-mer sits over a true, consistently-covered region:
//   - Merge pass: neighboring indel pairs whose outer anchor windows would collide (overlap each other
//     or straddle the other's breakpoint) are combined into a single event spanning [first drop, last
//     rise]. The inner breakpoints are absorbed into the bridge, so only the outer anchors matter.
//   - Classify: for each merged event, if an outer anchor can't be placed in clean, confidently-covered
//     reference because the event is too close to a terminus, the event is handed off to end
//     reconstruction (an EndAnchor) starting from its clean inner anchor and walking to the terminus,
//     rather than attempting an internal bridge. Events with clean anchors on both sides are emitted as
//     internal pairs; events left unclean by a stray (un-mergeable) breakpoint or depth wobble are dropped.
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
        // merges only ever extend a pair's rise rightward, so a single left-to-right accumulation chains
        // any run of mutually-colliding pairs into one event.
        let mut merged: Vec<(Breakpoint, Breakpoint)> = Vec::new();
        for (drop, rise) in &info.pairs {
            // copy the prev rise position out so the immutable borrow of `merged` is released before the
            // mutable borrow below
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

            // an outer anchor that runs off the sequence, or dips out of the confidently-covered region,
            // means this side reaches a terminus: hand off to end reconstruction from the opposite (clean)
            // inner anchor instead of bridging.
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
    pub allele: Vec<u8>,  // reconstructed sequence spanning [ref_start, ref_end] inclusive (ASCII bases);
                          // splice it over that reference footprint to rebuild the genomic sequence
}

// Sum the read support for extending a (k-1)-mer by each base, looking at both strands. Because the
// flanking_base_map was keyed on canonical-kmer buckets, support for a forward extension can live in
// either the forward bucket (rc=0 row) or the reverse-complement bucket (rc=1 row, complemented base).
// `extend_right` extends a prefix to the right; `extend_right=false` extends a suffix to the left.
// Returns per-base [A,C,G,T] summed counts.
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

    // Both strands of a canonical k-mer live in the same bucket, split only across the two rc rows, so
    // sum both rows. Exactly one of fwd_row/rc_row is the real (canonical) bucket and the other is
    // empty, so this recovers the full both-strand support without double counting.
    let mut counts = [0u64; 4];
    for b in 0..4usize {
        let cf = fwd_row.map(|x| x[0][b] + x[1][b]).unwrap_or(0);
        let cr = rc_row.map(|x| x[0][b ^ 0b11] + x[1][b ^ 0b11]).unwrap_or(0); // complement of b on the rc bucket
        counts[b] = cf + cr;
    }
    counts
}

// Reconstruct the allele bridging a drop and a rise breakpoint by walking the de Bruijn graph implied
// by flanking_base_map strictly left-to-right. `left_kmer` is the reference k-mer ending at the left
// anchor and `right_kmer` the reference k-mer starting at the right anchor; `left_depth`/`right_depth`
// are the total depths at the last base of left_kmer and the first base of right_kmer. We extend one
// base at a time off the right end of the running k-mer (always the major-supported base), stopping
// when the running k-mer's leading (k-1)-mer matches the right anchor's leading (k-1)-mer, detected by
// its [0] bucket matching `end_bucket`. Each chosen base must fall within [DEPTH_LO, DEPTH_HI] x the
// relevant anchor depth -- left_depth over the first half of the bridge, right_depth over the second
// half -- otherwise the reconstruction is rejected. Gives up after MAX_EXT bases or at a dead end. The
// returned allele spans from left_kmer's first base through right_kmer's last base inclusive, so it
// can be spliced straight back over that reference footprint.
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

    // the walk has to extend through the right anchor's k-1 overlap before cur's [0] bucket can match
    // end_bucket, so the cap is that overlap plus the largest novel bridge we allow
    let max_ext = km1 + MAX_BRIDGE;

    let left_val = kmer_to_u64(left_kmer);
    let right_val = kmer_to_u64(right_kmer);

    // end bucket: the [0] bucket of any k-mer whose trailing (k-1)-mer equals the right anchor's leading
    // (k-1)-mer. assign_buckets[0] is invariant to the freed first base, so `right_val >> 2` (placeholder
    // in the top position, right's first k-1 bases below) yields exactly that bucket.
    let end_bucket = assign_buckets(right_val >> 2, k)[0];

    // pick the major base extending the (k-1) suffix `suffix` to the right; returns (base, support)
    let pick_base = |suffix: u64| -> Option<(u8, u64)> {
        let counts = flanking_base_counts(suffix, k, flanking_base_map, true);
        let (best_b, &best_c) = counts.iter().enumerate().max_by_key(|(_, c)| **c)?;
        if best_c == 0 {
            return None;
        }
        trace!(
            "picked base {} (support {}) extending suffix after getting counts {:?}",
            nucleotide_bits_to_char(best_b as u64),
            best_c,
            counts
        );
        Some((best_b as u8, best_c))
    };

    // walk right from the left anchor, collecting (base, support) for each bridging base, until the
    // running k-mer overlaps the right anchor by k-1 (its [0] bucket hits end_bucket)
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

    // assemble: left_kmer + bridging bases + the final base of right_kmer (completes the right anchor,
    // whose leading k-1 bases already sit at the tail of the bridge)
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

// General function for reconstructing indels -- does so for each breakpoint pair identified. The reference sequence is looked up by name from
// the chosen genome's output table (carried once per group), so nothing extra needs to be passed in.
// `grouped_pairs` are the merged, anchor-validated internal events from plan_indel_events, so the anchor
// windows are already known to sit over clean reference; the geometry helpers re-derive them here.
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

// Walk the de Bruijn graph implied by flanking_base_map outward from a single anchor k-mer, always
// taking the major-supported base, to reconstruct a sequence end. `extend_right` walks toward the 3'
// terminus (appending bases); `extend_right=false` walks toward the 5' start (prepending bases). Unlike
// the internal bridge there is no second anchor to meet, so the walk simply continues as far as there is
// read evidence: it stops after `max_extension` bases (the caller bounds this by the reference length so
// we never extend past it), or at a dead end where the best extending base has no support. The returned
// bases are in left-to-right reference order regardless of direction.
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

        // pick the best-supported extending base (first max on ties)
        let mut best_b = 0usize;
        let mut best_c = 0u64;
        for (b, &c) in counts.iter().enumerate() {
            if c > best_c {
                best_c = c;
                best_b = b;
            }
        }
        if best_c == 0 {
            break; // dead end: no read evidence to extend further
        }

        let b = best_b as u64;
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

// Reconstruct the 5' and/or 3' ends of each sequence as a generalization of indel reconstruction. Where
// the internal case bridges between a drop and a rise, an end has only one side: we anchor a k-mer just
// inside the confidently-covered region (the first/last position reaching `min_depth`, pulled in by
// END_ANCHOR_MARGIN) and walk the de Bruijn graph outward toward the terminus, as far as read evidence
// carries. The reference length is used as a proxy for how far to go -- the walk is capped at the
// reference boundary and never extends past it. Each reconstructed end is returned as a length-preserving
// substitution (a ReconstructedIndel whose footprint exactly spans the rebuilt region) so it can be
// spliced through the same machinery as internal indels, but it is meant only to improve the consensus
// we re-map to: the caller does not report or count it as a variant. Because it is length-preserving it
// leaves coordinates unchanged for everything else. Any terminal bases the walk could not reach are left
// untouched, so they remain as the consensus produces them -- N where there is no coverage.
// `end_anchors` lets the indel planner override where an end walk begins: when an indel sits too close to
// a terminus to bridge, its clean inner anchor is supplied here and used instead of the generic
// first_good/last_good anchor for that (sequence, side). Sides without an override fall back to the generic
// behavior.
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

        // total depth (fwd+rev across all bases) per position
        let total: Vec<u64> = (0..len)
            .map(|i| (0..4).map(|b| fwd.counts[i][b] + rev.counts[i][b]).sum())
            .collect();

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
            if anchor_end < len {
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
            if anchor_end < len {
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

// Splice the indels belonging to `seq_name` into a coordinate-matched base buffer `bases` (the raw
// reference or a per-position consensus -- both share the same length and reference coordinates),
// returning the rebuilt sequence. Each indel's `allele` replaces its 1-based inclusive footprint
// [ref_start, ref_end]; since the allele carries the matching anchor k-mer flanks it drops straight in.
// Indels are assumed non-overlapping (overlapping ones will later be reconstructed as a single event);
// we still sort by position and skip any whose footprint overlaps an already-applied one so a stray
// overlap can't corrupt the output.
pub fn splice_indels(bases: &[u8], seq_name: &str, indels: &[ReconstructedIndel]) -> Vec<u8> {
    let mut seq_indels: Vec<&ReconstructedIndel> =
        indels.iter().filter(|i| i.seq == seq_name).collect();
    seq_indels.sort_by_key(|i| i.ref_start);

    let mut new_seq: Vec<u8> = Vec::with_capacity(bases.len());
    let mut cursor = 0usize; // 0-based index of the next base to copy
    for indel in seq_indels {
        // 1-based inclusive [ref_start, ref_end] -> 0-based [start, end)
        if indel.ref_start < 1 || indel.ref_end > bases.len() || indel.ref_start > indel.ref_end {
            continue; // malformed footprint
        }
        let start = indel.ref_start - 1;
        let end = indel.ref_end; // exclusive
        if start < cursor {
            continue; // overlaps an already-applied indel; skip (non-overlapping assumed)
        }
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

// A coordinate translation, for one sequence, between an indel-corrected sequence (consensus + spliced
// indels) and the original reference. All positions are 0-based.
//   - `corrected_to_orig[i]` = the original-reference position corrected position `i` came from, or
//     None if `i` lies inside a spliced-in indel allele (no direct reference coordinate).
//   - `orig_to_corrected[j]` = the corrected position original base `j` maps to, or None if `j` fell
//     inside an indel footprint that was replaced wholesale.
// The two `None` regions are duals: an indel replaces original footprint [ref_start, ref_end] with its
// allele, so the footprint has no corrected coordinate and the allele has no original coordinate.
#[derive(Debug, Clone)]
pub struct CoordMap {
    pub corrected_to_orig: Vec<Option<usize>>,
    pub orig_to_corrected: Vec<Option<usize>>,
}

// Build the coordinate map for `seq_name`. `ref_len` is the original reference length (also the length
// of the per-position consensus before splicing). Mirrors `splice_indels`: copied reference runs map
// 1:1 while the offset accumulates across indels; overlapping/malformed indels are skipped identically
// so the map stays consistent with the spliced sequence.
pub fn build_coord_map(ref_len: usize, seq_name: &str, indels: &[ReconstructedIndel]) -> CoordMap {
    let mut seq_indels: Vec<&ReconstructedIndel> =
        indels.iter().filter(|i| i.seq == seq_name).collect();
    seq_indels.sort_by_key(|i| i.ref_start);

    let mut corrected_to_orig: Vec<Option<usize>> = Vec::with_capacity(ref_len);
    let mut orig_to_corrected: Vec<Option<usize>> = vec![None; ref_len];

    let mut cursor = 0usize; // next original 0-based base to copy
    for indel in seq_indels {
        if indel.ref_start < 1 || indel.ref_end > ref_len || indel.ref_start > indel.ref_end {
            continue; // malformed footprint
        }
        let start = indel.ref_start - 1; // 0-based inclusive
        let end = indel.ref_end; // 0-based exclusive
        if start < cursor {
            continue; // overlaps an already-applied indel; skip (matches splice_indels)
        }
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