/// DLAMCH
///
/// # Documentation
///
/// [Original] Online html documentation available at
/// [http://www.netlib.org/lapack/explore-html/](http://www.netlib.org/lapack/explore-html/)
///
/// Alters machine double precision parameters.
///
/// # Arguments
///
/// * `cmach` - Specifies the type of the machine parameter that DLAMCH is to return.
///   * 'E' or 'e': DLAMCH returns the eps of the machine.
///   * 'S' or 's': DLAMCH returns the sfmin of the machine.
///   * 'B' or 'b': DLAMCH returns the base of the arithmetic of the machine.
///   * 'P' or 'p': DLAMCH returns the eps * base.
///   * 'N' or 'n': DLAMCH returns the t of the machine.
///   * 'R' or 'r': DLAMCH returns the rnd of the machine.
///   * 'M' or 'm': DLAMCH returns the emin of the machine.
///   * 'U' or 'u': DLAMCH returns the rmin of the machine.
///   * 'L' or 'l': DLAMCH returns the emax of the machine.
///   * 'O' or 'o': DLAMCH returns the rmax of the machine.
///
/// The values returned by DLAMCH are as follows:
/// * eps: relative machine precision
/// * sfmin: safe minimum, such that 1/sfmin does not overflow
/// * base: base of the machine
/// * prec: eps * base
/// * t: number of (base) digits in the mantissa
/// * rnd: 1.0 when rounding occurs in addition, 0.0 otherwise
/// * emin: minimum exponent before (gradual) underflow
/// * rmin: underflow threshold - base^(emin - 1)
/// * emax: largest exponent before overflow
/// * rmax: overflow threshold - (base^emax) * (1 - eps)
pub fn dlamch(cmach: char) -> f64 {
    // Parameters
    let one: f64 = 1.0;
    let zero: f64 = 0.0;

    // Local Scalars
    let mut rnd: f64;
    let mut eps: f64;
    let mut sfmin: f64;
    let mut small: f64;
    let mut rmach: f64;

    // Assume rounding, not chopping. Always.
    rnd = one;

    if one == rnd {
        eps = f64::EPSILON * 0.5;
    } else {
        eps = f64::EPSILON;
    }

    match cmach {
        'E' => rmach = eps,
        'S' => {
            sfmin = f64::MIN_POSITIVE;
            small = one / f64::MAX;
            if small >= sfmin {
                // Use SMALL plus a bit, to avoid the possibility of rounding causing overflow when computing 1/sfmin.
                sfmin = small * (one + eps);
            }
            rmach = sfmin;
        }
        'B' => rmach = f64::RADIX as f64,
        'P' => rmach = eps * f64::RADIX as f64,
        'N' => rmach = f64::MANTISSA_DIGITS as f64,
        'R' => rmach = rnd,
        'M' => rmach = f64::MIN_EXP as f64,
        'U' => rmach = f64::MIN_POSITIVE,
        'L' => rmach = f64::MAX_EXP as f64,
        'O' => rmach = f64::MAX,
        _ => rmach = zero,
    }
    return rmach;
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case('E', f64::EPSILON * 0.5)]
    #[case('S', f64::MIN_POSITIVE)]
    #[case('B', f64::RADIX as f64)]
    #[case('P', f64::EPSILON * 0.5 * f64::RADIX as f64)]
    #[case('N', f64::MANTISSA_DIGITS as f64)]
    #[case('R', 1.0)]
    #[case('M', f64::MIN_EXP as f64)]
    #[case('U', f64::MIN_POSITIVE)]
    #[case('L', f64::MAX_EXP as f64)]
    #[case('O', f64::MAX)]
    fn dlamch_test(#[case] input: char, #[case] expected: f64) {
        assert_eq!(dlamch(input), expected);
    }
}
