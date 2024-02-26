use crate::array::convert::ToFortranArray;

/// DNRM2
///
/// # Documentation
///
/// [Original] Online html documentation available at
/// [http://www.netlib.org/lapack/explore-html/](http://www.netlib.org/lapack/explore-html/)
///
/// Returns the Euclidean norm of a vector via the function
///
/// `dnrm2` returns sqrt(x' * x)
///
/// # Parameters
///
/// - `n`: The number of elements in input vector(s)
/// - `x`: The double precision array, with dimension `1 + (n - 1) * abs(incx)`
/// - `incx`: The storage spacing between elements of `x`
///
/// If `incx` > 0, `x(1+(i-1)*incx)` = `x(i)` for `1 <= i <= n`
///
/// If `incx` < 0, `x(1-(n-i)*incx)` = `x(i)` for `1 <= i <= n`
///
/// If `incx` = 0, `x` isn't a vector so there is no need to call
/// this function. If you call it anyway, it will count `x(1)`
/// in the vector norm `n` times.
pub fn dnrm2<X>(
    n: i32,
    x: &X,
    incx: i32,
) -> f64 where
    X: ToFortranArray,
{
    let x = x.to_fa();

    // Blue's scaling constants.
    let tsml: f64 = (f64::RADIX as f64).powf((f64::MIN_EXP as f64 - 1.) * 0.5);
    let tbig: f64 = (f64::RADIX as f64).powf((f64::MAX_EXP as f64 - f64::DIGITS as f64 + 1.) * 0.5);
    let ssml: f64 = (f64::RADIX as f64).powf(-(f64::MIN_EXP as f64 - f64::DIGITS as f64) * 0.5);
    let sbig: f64 = (f64::RADIX as f64).powf(-(f64::MAX_EXP as f64 + f64::DIGITS as f64 - 1.) * 0.5);

    // Quick return if possible
    if n <= 0 { return 0.; }

    let mut notbig = true;
    let mut asml = 0.;
    let mut amed = 0.;
    let mut abig = 0.;
    let mut ix = 1;
    if incx < 0 { ix = 1 - (n-1) * incx; }

    for _ in 1..=n {
        let ax = x[ix].abs();
        if ax > tbig {
            abig = abig + (ax*sbig).powf(2.);
            notbig = false;
        } else if ax < tsml {
            if notbig { asml = asml + (ax*ssml).powf(2.); }
        } else {
            amed = amed + ax.powf(2.);
        }
        ix += incx;
    }

    let (scl, sumsq) = if abig > 0. {
        // Combine abig and amed if abig > 0.
        if amed > 0. || amed > f64::MAX || amed != amed {
            abig = abig + (amed * sbig) * sbig;
        }

        (1. / sbig, abig)
    } else if asml > 0. {
        // Combine amed and asml if asml > 0.
        if amed > 0. || amed > f64::MAX || amed != amed {
            (amed, asml) = (amed.sqrt(), asml.sqrt() / ssml);
            let (ymin, ymax) = if asml > amed { (amed, asml) } else { (asml, amed) };

            (1., ymax.powf(2.) * (1. + (ymin / ymax).powf(2.)))
        } else {
            (1. / ssml, asml)
        }
    } else {
        // Otherwise all values are in the mid-range.
        (1., amed)
    };

    scl * sumsq.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(2, vec![2., 3.], 1, 3.605551275463989)]
    #[case(3, vec![4., 3., 2.], 1, 5.385164807134504)]
    #[case(1, vec![5.], 1, 5.)]
    #[case(1, vec![0.], 1, 0.)]
    #[case(2, vec![0., 0.], 1, 0.)]
    #[case(2, vec![-1., -1.], 1, std::f64::consts::SQRT_2)]
    #[case(2, vec![-1., 0.], 1, 1.)]
    #[case(2, vec![2.6, 7.8, 3.2], 2, 4.123105625617661)]
    #[case(3, vec![-1., 0., 1.], 1, std::f64::consts::SQRT_2)]
    #[case(3, vec![-1., 0., 1.], 2, std::f64::consts::SQRT_2)]
    fn test_dnrm2(
        #[case] n: i32,
        #[case] vectors: Vec<f64>,
        #[case] incx: i32,
        #[case] expected: f64) {
        assert!((dnrm2(n, &vectors, incx) - expected).abs() < 1e-15);
    }
}

#[test] fn test() {
    dnrm2(2, &vec![2., 3.], 1);
}
