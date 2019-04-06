extern crate rand;

use rand::Rng;
use std::error::Error;
use std::hash::{Hash, Hasher};

#[allow(deprecated)]
use std::hash::SipHasher;

pub struct HyperLogLog {
    /// `b` is the number of bit addressing registers, takes the range from 4 to 16.
    b: u8,

    /// `b_mask` is used to take `b` bits from hash value.
    b_mask: usize,

    /// the number of registers, calculated as 2 to the `b`th power.
    m: usize,

    /// ???
    alpha: f64,

    /// Registers with size of `m` bytes.
    registers: Vec<u8>,

    /// Keys are used in initialize SipHasher.
    hasher_key0: u64,
    hasher_key1: u64,
}

#[derive(Debug)]
pub enum Estimator {
    HyperLogLog,
    LinearCounting, // to estimate in small range.
}

impl HyperLogLog {
    pub fn new(b: u8) -> Result<Self, Box<Error>> {
        if b < 4 || b > 16 {
            return Err(From::from(format!("b must be between 4 and 16. b = {}", b)));
        }

        let m = 1 << b;
        let alpha = get_alpha(b)?;

        let mut rng = rand::OsRng::new()
            .map_err(|e| format!("Failed to create and OS RNG: {}", e))?;

        Ok(HyperLogLog {
            b: b,
            b_mask: m - 1,
            m: m,
            alpha: alpha,
            registers: vec![0; m],
            hasher_key0: rng.gen(),
            hasher_key1: rng.gen(),
        })
    }

    pub fn insert<H: Hash>(&mut self, value: &H) {
        let x = self.hash(value);
        let j = x as usize & self.b_mask;
        let w = x >> self.b;

        let p1 = position_of_leftmost_one_bit(w, 64 - self.b);
        let p2 = &mut self.registers[j];
        if *p2 < p1 {
            *p2 = p1
        }
    }

    #[allow(deprecated)]
    fn hash<H: Hash>(&self, value: &H) -> u64 {
        let mut hasher = SipHasher::new_with_keys(self.hasher_key0, self.hasher_key1);
        value.hash(&mut hasher);
        hasher.finish()
    }

    pub fn cardinality(&self) -> f64 {
        estimate_cardinality(self).0
    }
}

fn get_alpha(b: u8) -> Result<f64, Box<Error>> {
    if b < 4 || b > 16 {
        return Err(From::from(format!("b must be between 4 and 16. b = {}", b)));
    } else {
        Ok(match b {
            4 => 0.673,
            5 => 0.697,
            6 => 0.709,
            _ => 0.7213 / (1.0 + 1.079 / (1 << b) as f64),
        })
    }
}

fn position_of_leftmost_one_bit(s: u64, max_width: u8) -> u8 {
    count_leading_zeros(s, max_width) + 1
}

fn count_leading_zeros(mut s: u64, max_width: u8) -> u8 {
    let mut lz = max_width;
    while s != 0 {
        lz -= 1;
        s >>= 1;
    }
    lz
}

fn estimate_cardinality(hll: &HyperLogLog) -> (f64, Estimator) {
    let m_f64 = hll.m as f64;
    let est = raw_hyperloglog_estimate(hll.alpha, m_f64, &hll.registers);

    if est < (5.0 / 2.0 * m_f64) {
        match count_zero_registers(&hll.registers) {
            0 => (est, Estimator::HyperLogLog),
            v => (linear_counting_estimate(m_f64, v as f64), Estimator::LinearCounting),
        }
    } else {
        (est, Estimator::HyperLogLog)
    }
}

fn count_zero_registers(regsiters: &[u8]) -> usize {
    regsiters.iter().filter(|&x| *x == 0).count()
}

fn raw_hyperloglog_estimate(alpha: f64, m: f64, registers: &[u8]) -> f64 {
    let sum = registers.iter()
        .map(|&x| 2.0f64.powi(-(x as i32))).sum::<f64>();
    alpha * m * m / sum
}

fn linear_counting_estimate(m: f64, number_of_zero_registers: f64) -> f64 {
    m * (m / number_of_zero_registers).ln()
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_hll() {
        assert!(HyperLogLog::new(3).is_err());
        assert!(HyperLogLog::new(17).is_err());

        let hll = HyperLogLog::new(4);
        assert!(hll.is_ok());

        let hll = hll.unwrap();
        assert_eq!(hll.b, 4);
        assert_eq!(hll.m, 2_usize.pow(4));
        assert_eq!(hll.alpha, 0.673);
        assert_eq!(hll.registers.len(), 2_usize.pow(4));

        assert!(HyperLogLog::new(16).is_ok());
    }

    #[test]
    fn small_range() {
        let mut hll = HyperLogLog::new(12).unwrap();
        let items = ["test1", "test2", "test3", "test2", "test2", "test2"];

        println!("\n=== Loading {} items.\n", items.len());
        for item in &items {
            hll.insert(item);
        }
        assert_eq!(hll.cardinality().round(), 3.0);
    }
}
