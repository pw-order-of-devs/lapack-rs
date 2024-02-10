/// DLAPY2
///
/// # Documentation
///
/// [Original] Online html documentation available at
/// [http://www.netlib.org/lapack/explore-html/](http://www.netlib.org/lapack/explore-html/)
///
/// Computes sqrt(x² + y²) while trying to avoid unnecessary overflow and underflow.
///
/// # Arguments
///
/// For arguments definitions, please refer to the original documentation.
pub fn dlapy2(x: f64, y: f64) -> f64 {
    if x.is_nan() { return x; }
    if y.is_nan() { return y; }

    let w = x.abs().max(y.abs());
    let z = x.abs().min(y.abs());

    return if z == 0.0 || w > f64::MAX { w }
    else { w * ((1.0 + (z / w).powi(2)).sqrt()) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(0.0, 0.0, 0.0)]
    #[case(3.0, 4.0, 5.0)]
    #[case(-3.0, 4.0, 5.0)]
    #[case(3.0, -4.0, 5.0)]
    #[case(-3.0, -4.0, 5.0)]
    #[case(8.0, 15.0, 17.0)]
    #[case(0.1, 0.0, 0.1)]
    #[case(0.0, 0.1, 0.1)]
    #[case(1e-307, 1e-307, 1.414213562373095e-307)]
    #[case(1e307, 1e307, 1.4142135623730951e307)]
    fn test_dlapy2(#[case] x: f64, #[case] y: f64, #[case] expected: f64) {
        let result = dlapy2(x, y);
        assert!((expected - result).abs() < f64::EPSILON, "Expected {} but got {} for x = {} and y = {}", expected, result, x, y);
    }
}