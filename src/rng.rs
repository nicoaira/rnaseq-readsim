//! Pseudo-random numbers: xoshiro256++ generators seeded through SplitMix64.
//!
//! The generator is implemented here (about thirty lines) instead of taken
//! from a crate so that the output of the simulator for a given seed does not
//! change when a dependency is updated.
//!
//! References: D. Blackman and S. Vigna, "Scrambled linear pseudorandom number
//! generators", ACM Trans. Math. Softw. 47 (2021), doi:10.1145/3460772
//! (xoshiro256++); G. L. Steele, D. Lea and C. H. Flood, "Fast splittable
//! pseudorandom number generators", OOPSLA 2014, doi:10.1145/2660193.2660195
//! (SplitMix64).

const GOLDEN_GAMMA: u64 = 0x9E37_79B9_7F4A_7C15;

/// The SplitMix64 output function (a bijection on 64-bit integers).
#[inline]
pub fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// SplitMix64, used only to expand a 64-bit seed into a xoshiro state.
pub struct SplitMix64(u64);

impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        SplitMix64(seed)
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(GOLDEN_GAMMA);
        mix64(self.0)
    }
}

/// xoshiro256++ generator.
#[derive(Clone, Debug)]
pub struct Rng {
    s: [u64; 4],
}

impl Rng {
    /// Creates a generator from an explicit state (must not be all zero).
    pub fn from_state(s: [u64; 4]) -> Self {
        assert!(s != [0; 4], "xoshiro256++ state must not be all zero");
        Rng { s }
    }

    /// Creates the generator of one independent stream.
    ///
    /// The simulator gives every chunk of fragments its own stream, identified
    /// by the chunk index, so that the output does not depend on which thread
    /// produced which chunk. The seed and the stream index are hashed
    /// separately and combined before the state is expanded with SplitMix64.
    pub fn for_stream(seed: u64, stream: u64) -> Self {
        let key = mix64(seed).wrapping_add(mix64(!stream));
        let mut sm = SplitMix64::new(key);
        let mut s = [0u64; 4];
        loop {
            for word in &mut s {
                *word = sm.next_u64();
            }
            if s != [0; 4] {
                return Rng { s };
            }
        }
    }

    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let s = &mut self.s;
        let result = s[0].wrapping_add(s[3]).rotate_left(23).wrapping_add(s[0]);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = s[3].rotate_left(45);
        result
    }

    /// A uniform integer in `0..n` (`n` must be at least 1).
    ///
    /// Uses the multiply-and-shift mapping without a rejection step; the
    /// deviation from uniformity is below `n / 2^64`, which is far below
    /// anything a sequencing experiment could resolve.
    #[inline]
    pub fn below(&mut self, n: u64) -> u64 {
        ((self.next_u64() as u128 * n as u128) >> 64) as u64
    }

    /// A uniform number in `[0, 1)` with 53 random bits.
    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
}

/// Converts a probability into a threshold `t` such that `next_u64() < t`
/// is true with that probability.
pub fn threshold_u64(p: f64) -> u64 {
    if p <= 0.0 || p.is_nan() {
        0
    } else if p >= 1.0 {
        u64::MAX
    } else {
        // 2^64 as f64 is exact; the cast saturates.
        (p * 18_446_744_073_709_551_616.0) as u64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn xoshiro256plusplus_reference_values() {
        // First outputs for the state [1, 2, 3, 4]; the same vector is used by
        // the test-suite of the `rand_xoshiro` crate, and the first two values
        // are easy to check by hand: rotl(1 + 4, 23) + 1 = 41943041.
        let mut rng = Rng::from_state([1, 2, 3, 4]);
        let expected: [u64; 10] = [
            41943041,
            58720359,
            3588806011781223,
            3591011842654386,
            9228616714210784205,
            9973669472204895162,
            14011001112246962877,
            12406186145184390807,
            15849039046786891736,
            10450023813501588000,
        ];
        for e in expected {
            assert_eq!(rng.next_u64(), e);
        }
    }

    #[test]
    fn splitmix64_reference_values() {
        // Reference outputs of SplitMix64 for seed 1234567.
        let mut sm = SplitMix64::new(1234567);
        assert_eq!(sm.next_u64(), 6457827717110365317);
        assert_eq!(sm.next_u64(), 3203168211198807973);
        assert_eq!(sm.next_u64(), 9817491932198370423);
    }

    #[test]
    fn streams_differ_and_are_reproducible() {
        let mut a = Rng::for_stream(42, 0);
        let mut b = Rng::for_stream(42, 1);
        let mut a2 = Rng::for_stream(42, 0);
        let mut c = Rng::for_stream(43, 0);
        let va: Vec<u64> = (0..8).map(|_| a.next_u64()).collect();
        let vb: Vec<u64> = (0..8).map(|_| b.next_u64()).collect();
        let va2: Vec<u64> = (0..8).map(|_| a2.next_u64()).collect();
        let vc: Vec<u64> = (0..8).map(|_| c.next_u64()).collect();
        assert_eq!(va, va2);
        assert_ne!(va, vb);
        assert_ne!(va, vc);
    }

    #[test]
    fn below_is_in_range_and_roughly_uniform() {
        let mut rng = Rng::for_stream(7, 0);
        let mut counts = [0u32; 5];
        for _ in 0..50_000 {
            counts[rng.below(5) as usize] += 1;
        }
        for c in counts {
            assert!((9_500..10_500).contains(&c), "count {c}");
        }
        for _ in 0..1000 {
            assert_eq!(rng.below(1), 0);
        }
    }

    #[test]
    fn thresholds() {
        assert_eq!(threshold_u64(0.0), 0);
        assert_eq!(threshold_u64(1.0), u64::MAX);
        assert_eq!(threshold_u64(0.5), 1u64 << 63);
        assert_eq!(threshold_u64(-1.0), 0);
    }
}
