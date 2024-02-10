/// IEEECK
///
/// [Original] Online HTML documentation available at
/// `http://www.netlib.org/lapack/explore-html/`
///
/// # Arguments
///
/// For arguments definitions, please refer to the original documentation.
pub fn ieeeck(
    ispec: i32,
    zero: f64,
    one: f64,
) -> i32 {
    let mut posinf = one / zero;
    if posinf <= one { return 0 }

    let mut neginf = -one / zero;
    if neginf >= zero { return 0 }

    let negzro = one / (neginf + one);
    if negzro != zero { return 0 }

    neginf = one / negzro;
    if neginf >= zero { return 0 }

    let newzro = negzro + zero;
    if newzro != zero { return 0 }

    posinf = one / newzro;
    if posinf <= one { return 0 }

    neginf = neginf * posinf;
    if neginf >= zero { return 0 }

    posinf = posinf * posinf;
    if posinf <= one { return 0 }

    if ispec == 0 { return 1 }

    let nan1 = posinf + neginf;
    let nan2 = posinf / neginf;
    let nan3 = posinf / posinf;
    let nan4 = posinf * zero;
    let nan5 = neginf * negzro;
    let nan6 = nan5 * zero;

    if nan1.is_nan() || nan2.is_nan() || nan3.is_nan()
        || nan4.is_nan() || nan5.is_nan() || nan6.is_nan() { return 0 }
    1
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(0, 0., 1., 1)] // Basic Test Case
    #[case(1, 0., 1., 0)] // Basic Test Case with NaN checks
    #[case(0, 0., 1e-10, 1)] // Edge Case: Very small positive value
    #[case(0, 0., -1e-10, 0)] // Edge Case: Very small negative value
    #[case(1, 0., 0., 0)] // Zero Division Case: NaN produced for 0. / 0.
    #[case(1, 0., 1e-10, 0)] // Zero Division Case: NaN produced for small positive value / 0.
    #[case(1, 0., -1., 0)] // Negative Infinity Case: NaN produced for sqrt(-1.)
    #[case(1, 0., -1e-10, 0)] // Negative Infinity Case: NaN produced for sqrt(-small positive value)
    #[case(1, 0., 1e20, 0)] // Combination Case: Overflow for very large positive value
    #[case(1, 0., -1e20, 0)] // Combination Case: Underflow for very large negative value
    #[case(1, f64::NAN, 1., 0)] // NaN Arithmetic: NaN input produces NaN output
    #[case(1, 0., f64::NAN, 0)] // NaN Arithmetic: NaN input produces NaN output
    fn test_ieeeck(
        #[case] ispec: i32,
        #[case] zero: f64,
        #[case] one: f64,
        #[case] expected: i32,
    ) {
        assert_eq!(expected, ieeeck(ispec, zero, one))
    }
}
