use crate::call::OutputData;
use crate::lcb::*;

use log::*;
use dashmap::DashMap;

#[derive(Debug, Clone)]
pub struct Breakpoint {
    pub pos: usize,      // 1-based position on the low-coverage side of the depth change
    pub direction: i8,   // +1 for a rise in depth, -1 for a drop
    pub depth_high: u64, // larger total depth of the two consecutive positions the change spans
    pub depth_low: u64,  // smaller total depth of the two consecutive positions the change spans
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
) -> Vec<(String, Vec<(Breakpoint, Breakpoint)>)> {
    let mut grouped: Vec<(String, Vec<(Breakpoint, Breakpoint)>)> = Vec::new();

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

        // pair a drop (negative) with the next rise (positive) within max_pair_distance
        let mut seq_pairs: Vec<(Breakpoint, Breakpoint)> = Vec::new();
        let mut pending_drop: Option<Breakpoint> = None;
        for bp in breakpoints {
            match bp.direction {
                -1 => pending_drop = Some(bp), // most recent drop opens a potential indel
                1 => {
                    if let Some(drop) = pending_drop.take() {
                        if bp.pos.saturating_sub(drop.pos) <= max_pair_distance {
                            seq_pairs.push((drop, bp));
                        }
                    }
                }
                _ => {}
            }
        }

        if !seq_pairs.is_empty() {
            grouped.push((seq.clone(), seq_pairs));
        }
    }

    grouped
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
    const MAX_BRIDGE: usize = 20; // largest novel (inserted) span we attempt to bridge
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
pub fn reconstruct_indels(
    output: &DashMap<String, OutputData>,
    flanking_base_map: &DashMap<u64, [[u64; 4]; 2]>,
    grouped_pairs: &[(String, Vec<(Breakpoint, Breakpoint)>)],
    k: usize,
) -> Vec<ReconstructedIndel> {
    // pull each anchor k-mer this many bases further from the breakpoint, into more confidently-covered
    // reference, since coverage (and base calls) right at the cliff can be unreliable
    const ANCHOR_MARGIN: usize = 2;

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
            let left_end = match drop.pos.checked_sub(ANCHOR_MARGIN+2) { // -2 because -1 for from 1-based to 0-based, -1 because based off second nucleotide
                Some(v) => v, // 0-based, last base of the left anchor k-mer
                None => continue,
            };
            if left_end + 1 < k {
                continue;
            }
            let left_start = left_end + 1 - k; // 0-based, first base of the left anchor k-mer
            let right_start = rise.pos - 1 + ANCHOR_MARGIN; // 0-based, first base of the right anchor k-mer
            if right_start + k > ref_bases.len() {
                continue;
            }
            let right_end = right_start + k - 1; // 0-based, last base of the right anchor k-mer

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
