//! The fragment-length distribution and the effective transcript length.

use crate::alias::AliasTable;
use crate::rng::Rng;

/// A normal distribution restricted to the integers `min..=max`.
///
/// `P(L) ∝ exp(-(L - mean)² / (2 sd²))` for integer `L` in `min..=max`
/// (the normal density evaluated at the integers, renormalised). With
/// `sd = 0` every fragment has length `round(mean)`.
#[derive(Clone, Debug)]
pub struct FragmentModel {
    min: usize,
    max: usize,
    /// `pmf[i]` is `P(L = min + i)`.
    pmf: Vec<f64>,
    /// `cum_p[i]` is `P(L <= min + i)`.
    cum_p: Vec<f64>,
    /// `cum_lp[i]` is the sum of `L * P(L)` for `L <= min + i`.
    cum_lp: Vec<f64>,
    alias: AliasTable,
}

impl FragmentModel {
    pub fn new(mean: f64, sd: f64, min: usize, max: usize) -> Result<Self, String> {
        if !mean.is_finite() || mean <= 0.0 {
            return Err(format!(
                "the mean fragment length must be a positive number (got {mean})"
            ));
        }
        if !sd.is_finite() || sd < 0.0 {
            return Err(format!(
                "the fragment length standard deviation must be zero or positive (got {sd})"
            ));
        }
        if min < 1 {
            return Err("the minimum fragment length must be at least 1".into());
        }
        if max < min {
            return Err(format!(
                "the maximum fragment length ({max}) is below the minimum ({min})"
            ));
        }
        if max > 100_000_000 {
            return Err(format!("the maximum fragment length ({max}) is too large"));
        }
        let n = max - min + 1;
        let mut pmf = vec![0.0f64; n];
        if sd == 0.0 {
            let l = mean.round();
            if l < min as f64 || l > max as f64 {
                return Err(format!(
                    "with a standard deviation of 0 the mean fragment length ({mean}) must lie \
                     between the minimum ({min}) and the maximum ({max})"
                ));
            }
            pmf[l as usize - min] = 1.0;
        } else {
            // Exponents relative to the largest one, so that a range far in
            // the tail of the normal still has a usable distribution.
            let exponent = |l: usize| -> f64 {
                let z = (l as f64 - mean) / sd;
                -0.5 * z * z
            };
            let best = (min..=max).map(exponent).fold(f64::NEG_INFINITY, f64::max);
            for (i, p) in pmf.iter_mut().enumerate() {
                *p = (exponent(min + i) - best).exp();
            }
        }
        let total: f64 = pmf.iter().sum();
        for p in &mut pmf {
            *p /= total;
        }
        let mut cum_p = Vec::with_capacity(n);
        let mut cum_lp = Vec::with_capacity(n);
        let (mut a, mut b) = (0.0f64, 0.0f64);
        for (i, &p) in pmf.iter().enumerate() {
            a += p;
            b += p * (min + i) as f64;
            cum_p.push(a);
            cum_lp.push(b);
        }
        let alias = AliasTable::new(&pmf)?;
        Ok(FragmentModel {
            min,
            max,
            pmf,
            cum_p,
            cum_lp,
            alias,
        })
    }

    pub fn min(&self) -> usize {
        self.min
    }

    pub fn max(&self) -> usize {
        self.max
    }

    /// `P(L = len)`.
    pub fn probability(&self, len: usize) -> f64 {
        if len < self.min || len > self.max {
            0.0
        } else {
            self.pmf[len - self.min]
        }
    }

    /// The mean of the distribution itself (before any transcript is involved).
    pub fn mean(&self) -> f64 {
        *self.cum_lp.last().expect("non-empty distribution")
    }

    /// Draws a fragment length.
    #[inline]
    pub fn sample(&self, rng: &mut Rng) -> usize {
        self.min + self.alias.sample(rng)
    }

    /// The effective length of a transcript of `len` bases:
    /// `sum over L of P(L) * max(len - L + 1, 1)`.
    ///
    /// `len - L + 1` is the number of start positions of a fragment of `L`
    /// bases; when the drawn length exceeds the transcript, the fragment is
    /// the whole transcript, which has exactly one start position.
    pub fn effective_length(&self, len: usize) -> f64 {
        if len < self.min {
            return 1.0;
        }
        let i = (len - self.min).min(self.pmf.len() - 1);
        let p_fit = self.cum_p[i].min(1.0);
        let lp_fit = self.cum_lp[i];
        (len as f64 + 1.0) * p_fit - lp_fit + (1.0 - p_fit)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distribution_is_normalised_and_truncated() {
        let m = FragmentModel::new(180.0, 55.0, 40, 600).unwrap();
        let total: f64 = (0..1000).map(|l| m.probability(l)).sum();
        assert!((total - 1.0).abs() < 1e-12);
        assert_eq!(m.probability(39), 0.0);
        assert_eq!(m.probability(601), 0.0);
        assert!(m.probability(180) > m.probability(100));
        // symmetric around the mean where both sides are inside the range
        assert!((m.probability(150) - m.probability(210)).abs() < 1e-15);
        let mut rng = Rng::for_stream(5, 0);
        for _ in 0..10_000 {
            let l = m.sample(&mut rng);
            assert!((40..=600).contains(&l));
        }
    }

    #[test]
    fn sampled_mean_matches() {
        let m = FragmentModel::new(180.0, 55.0, 40, 600).unwrap();
        let mut rng = Rng::for_stream(9, 0);
        let n = 200_000;
        let sum: usize = (0..n).map(|_| m.sample(&mut rng)).sum();
        let mean = sum as f64 / n as f64;
        assert!((mean - m.mean()).abs() < 0.6, "{mean} vs {}", m.mean());
    }

    #[test]
    fn fixed_length() {
        let m = FragmentModel::new(250.0, 0.0, 40, 600).unwrap();
        let mut rng = Rng::for_stream(1, 0);
        for _ in 0..100 {
            assert_eq!(m.sample(&mut rng), 250);
        }
        assert_eq!(m.effective_length(1000), 751.0);
        assert_eq!(m.effective_length(250), 1.0);
        assert_eq!(m.effective_length(100), 1.0);
        assert!(FragmentModel::new(700.0, 0.0, 40, 600).is_err());
    }

    #[test]
    fn effective_length_matches_direct_sum() {
        let m = FragmentModel::new(180.0, 55.0, 40, 600).unwrap();
        for len in [1usize, 39, 40, 41, 99, 180, 181, 599, 600, 601, 5000] {
            let direct: f64 = (40..=600)
                .map(|l| {
                    let starts = if l <= len { (len - l + 1) as f64 } else { 1.0 };
                    m.probability(l) * starts
                })
                .sum();
            let fast = m.effective_length(len);
            assert!(
                (direct - fast).abs() < 1e-8 * direct.max(1.0),
                "len {len}: {direct} vs {fast}"
            );
        }
        // long transcripts: length - mean + 1
        assert!((m.effective_length(5000) - (5001.0 - m.mean())).abs() < 1e-8);
    }

    #[test]
    fn rejects_bad_parameters() {
        assert!(FragmentModel::new(180.0, 55.0, 0, 600).is_err());
        assert!(FragmentModel::new(180.0, 55.0, 700, 600).is_err());
        assert!(FragmentModel::new(180.0, -1.0, 40, 600).is_err());
        assert!(FragmentModel::new(f64::NAN, 55.0, 40, 600).is_err());
        // far in the tail: still a valid distribution
        let m = FragmentModel::new(100_000.0, 1.0, 40, 600).unwrap();
        assert!((m.probability(600) - 1.0).abs() < 1e-9);
    }
}
