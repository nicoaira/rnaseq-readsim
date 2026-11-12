//! Walker/Vose alias tables: sampling from a discrete distribution in O(1).
//!
//! Reference: M. D. Vose, "A linear algorithm for generating random numbers
//! with a given distribution", IEEE Trans. Softw. Eng. 17 (1991),
//! doi:10.1109/32.92917.

use crate::rng::{Rng, threshold_u64};

/// Builds the alias method's tables for non-negative weights.
///
/// Returns, for every column `i`, the probability of keeping `i` and the
/// index to use otherwise. A column whose weight is exactly zero is never
/// returned by sampling.
pub fn vose(weights: &[f64]) -> Result<(Vec<f64>, Vec<u32>), String> {
    let n = weights.len();
    if n == 0 {
        return Err("no categories to sample from".into());
    }
    if n > u32::MAX as usize {
        return Err("too many categories for an alias table".into());
    }
    let mut total = 0.0f64;
    let mut heaviest = 0usize;
    for (i, &w) in weights.iter().enumerate() {
        if !w.is_finite() || w < 0.0 {
            return Err(format!("weight {w} is not a finite, non-negative number"));
        }
        total += w;
        if w > weights[heaviest] {
            heaviest = i;
        }
    }
    if total <= 0.0 || !total.is_finite() {
        return Err("the weights do not have a positive, finite sum".into());
    }

    let scale = n as f64 / total;
    let mut prob: Vec<f64> = weights.iter().map(|w| w * scale).collect();
    let mut alias: Vec<u32> = (0..n as u32).collect();
    let mut small: Vec<usize> = Vec::new();
    let mut large: Vec<usize> = Vec::new();
    for (i, &p) in prob.iter().enumerate() {
        if p < 1.0 {
            small.push(i);
        } else {
            large.push(i);
        }
    }
    while let (Some(&s), Some(&l)) = (small.last(), large.last()) {
        small.pop();
        alias[s] = l as u32;
        prob[l] = (prob[l] + prob[s]) - 1.0;
        if prob[l] < 1.0 {
            large.pop();
            small.push(l);
        }
    }
    // Whatever is left over is there because of rounding and has probability
    // one, except that a zero-weight column must stay unreachable.
    for &l in &large {
        prob[l] = 1.0;
    }
    for &s in &small {
        if weights[s] == 0.0 {
            prob[s] = 0.0;
            alias[s] = heaviest as u32;
        } else {
            prob[s] = 1.0;
        }
    }
    Ok((prob, alias))
}

/// An alias table with 64-bit thresholds, for large or very uneven
/// distributions (transcripts, fragment lengths).
#[derive(Clone, Debug)]
pub struct AliasTable {
    threshold: Vec<u64>,
    alias: Vec<u32>,
}

impl AliasTable {
    pub fn new(weights: &[f64]) -> Result<Self, String> {
        let (prob, alias) = vose(weights)?;
        Ok(AliasTable {
            threshold: prob.iter().map(|&p| threshold_u64(p)).collect(),
            alias,
        })
    }

    pub fn len(&self) -> usize {
        self.alias.len()
    }

    pub fn is_empty(&self) -> bool {
        self.alias.is_empty()
    }

    /// Draws one category (two random numbers).
    #[inline]
    pub fn sample(&self, rng: &mut Rng) -> usize {
        let i = rng.below(self.alias.len() as u64) as usize;
        if rng.next_u64() < self.threshold[i] {
            i
        } else {
            self.alias[i] as usize
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_bad_weights() {
        assert!(AliasTable::new(&[]).is_err());
        assert!(AliasTable::new(&[0.0, 0.0]).is_err());
        assert!(AliasTable::new(&[1.0, -1.0]).is_err());
        assert!(AliasTable::new(&[1.0, f64::NAN]).is_err());
        assert!(AliasTable::new(&[1.0, f64::INFINITY]).is_err());
    }

    #[test]
    fn frequencies_follow_weights() {
        let weights = [5.0, 0.0, 1.0, 3.0, 0.0, 1.0];
        let table = AliasTable::new(&weights).unwrap();
        let mut rng = Rng::for_stream(11, 0);
        let n = 200_000;
        let mut counts = [0u32; 6];
        for _ in 0..n {
            counts[table.sample(&mut rng)] += 1;
        }
        let total: f64 = weights.iter().sum();
        for (i, &w) in weights.iter().enumerate() {
            let expected = n as f64 * w / total;
            if w == 0.0 {
                assert_eq!(counts[i], 0, "zero-weight category {i} was drawn");
            } else {
                let sd = (expected * (1.0 - w / total)).sqrt();
                assert!(
                    (counts[i] as f64 - expected).abs() < 5.0 * sd,
                    "category {i}: {} vs {expected}",
                    counts[i]
                );
            }
        }
    }

    #[test]
    fn zero_weights_are_never_drawn_even_when_many() {
        let mut weights = vec![0.0; 1000];
        weights[17] = 1e-300;
        weights[900] = 1e300;
        let table = AliasTable::new(&weights).unwrap();
        let mut rng = Rng::for_stream(3, 0);
        for _ in 0..100_000 {
            let i = table.sample(&mut rng);
            assert!(i == 17 || i == 900);
        }
    }

    #[test]
    fn single_category() {
        let table = AliasTable::new(&[2.5]).unwrap();
        let mut rng = Rng::for_stream(1, 0);
        assert_eq!(table.sample(&mut rng), 0);
        assert_eq!(table.len(), 1);
    }
}
