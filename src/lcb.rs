use crate::consts::MAX_KMER_SIZE;

pub fn assign_buckets(kmer: u64, k: usize) -> Vec<u64> {
    let mut buckets = [0u64; 32];
    let mut num_a = [0u64; 32];
    let mut val = [0u64; 32];
    let mut mu = [0u64; 32];

    // Initialize mask and power variables
    let mut mask = 3u64 << ((k - 1) * 2);
    let mut p = 1u64 << ((k - 1) * 2);
    let mut cur = kmer & mask;

    val[0] = kmer - cur;
    mu[0] = if cur != 0 {
        p + ((cur >> 2) * (k as u64 - 1))
    } else {
        val[0]
    };
    let mut sum_mu = mu[0];

    for i in 1..k {
        num_a[i] = num_a[i - 1] + if cur == 0 { 1 } else { 0 };

        mask >>= 2;
        cur = kmer & mask;
        p >>= 2;
        val[i] = val[i - 1] - cur;
        mu[i] = if cur != 0 {
            p + ((cur >> 2) * (k as u64 - i as u64 - 1))
        } else {
            val[i]
        };
        sum_mu += mu[i];
    }

    mask = 3u64 << ((k - 1) * 2);
    for i in 0..k {
        cur = kmer & mask;
        mask >>= 2;

        let p = sum_mu - mu[i] + val[i] - num_a[i] * cur + 1 + num_a[i];
        buckets[i] = p;
    }

    buckets[..k].to_vec()
}



pub fn nt_to_bits(nt: u8) -> u8 {
    match nt {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0, // todo: handle other bases
    }
}

pub fn nucleotide_bits_to_char(bits: u64) -> char {
    match bits {
        0b00 => 'A',
        0b01 => 'C',
        0b10 => 'G',
        0b11 => 'T',
        _ => 'N',
    }
}

pub fn kmer_to_u64(kmer: &[u8]) -> u64 {
    let mut val = 0u64;
    for &base in kmer {
        val <<= 2;
        val |= nt_to_bits(base) as u64;
    }
    val
}

pub fn reverse_complement_u64(kmer_val: u64, k: usize) -> u64 {
    //reverse bits in place in 3 iterations rather than iterating k times. 
    assert!(k <= 32, "Ksize must be less than 32"); //for u64 memory. to increase k, consider using u128
    let rev = kmer_val.reverse_bits() >> (64 - 2 * k);
    let rev_fixed = ((rev >> 1) & 0x5555_5555_5555_5555u64)
        | ((rev & 0x5555_5555_5555_5555u64) << 1);
    rev_fixed ^ ((1u64 << (2 * k)).wrapping_sub(1))
}


pub fn canonical_kmer(kmer: &[u8], k: usize) -> (u64, bool) {
    let fwd = kmer_to_u64(kmer);
    let rev = reverse_complement_u64(fwd, k);
    if fwd < rev {
        (fwd, false)  // canonical is same as forward
    } else {
        (rev, true)   // canonical is different (reverse complement)
    }
}

pub fn canonical_kmer_u64(kmer: u64, k: usize) -> (u64, bool) { 
    let rev = reverse_complement_u64(kmer, k);
    if kmer < rev {
        (kmer, false)  // canonical is same as forward
    } else {
        (rev, true)   // canonical is different (reverse complement)
    }
}

pub fn seq_to_canon_kmers(seq: Vec<u8>, k: usize) -> Vec<(u64, bool)> {
    let mut results = Vec::with_capacity(seq.len().saturating_sub(k) + 1);
    if seq.len() < k {
        return results;
    }

    let mut current_kmer = 0u64;
    let mut valid_bases = 0;

    for i in 0..seq.len() {
        let base = seq[i];
        let bits = match base {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => {
                // Non-ACGT base, reset
                valid_bases = 0;
                current_kmer = 0;
                continue;
            }
        };

        current_kmer = ((current_kmer << 2) | bits) & ((1 << (2 * k)) - 1);
        valid_bases += 1;

        if valid_bases >= k {
            results.push(canonical_kmer_u64(current_kmer, k))
        }
    }

    results
}

pub fn parse_specify_pattern(pattern: &str, k: usize) -> Result<Vec<bool>, String> {
    if pattern.len() > MAX_KMER_SIZE {
        return Err(format!("--specify-pattern is too long ({} chars) - must be {} characters or fewer", pattern.len(), MAX_KMER_SIZE));
    }
    if pattern.len() != k {
        return Err(format!("--specify-pattern length ({}) must exactly match --kmer-size ({})", pattern.len(), k));
    }
    if !pattern.chars().all(|c| c == '0' || c == '1') {
        return Err("--specify-pattern must contain only '0' and '1'".to_string());
    }
    if !pattern.contains('1') {
        return Err("--specify-pattern must keep at least one position".to_string());
    }
    Ok(pattern.chars().map(|c| c == '1').collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn astring(){
        assert_eq!(vec![1, 2, 3, 4], assign_buckets(0, 4))
    }

    #[test]
    fn kstring(){
        assert_eq!(vec![238258108556, 47877379752, 215381104296, 227729135272, 235782198952, 237342480040, 238258108557, 238236915369, 238248449705, 238254544553, 238258108558, 238257944234, 238258089642, 238258095018, 238258106282, 238258108559, 238258108483, 238258108525, 238258108547], assign_buckets(41547505179, 19))
    }
}
