use crate::array::{AsFortranArray, FortranArray};
use crate::blas::dnrm2::dnrm2;
use crate::blas::dscal::dscal;
use crate::dlamch::dlamch;
use crate::dlapy2::dlapy2;

/// DLARFG
///
/// # Documentation
///
/// [Original] Online html documentation available at
/// [http://www.netlib.org/lapack/explore-html/](http://www.netlib.org/lapack/explore-html/)
///
/// Generates a real elementary reflector H of order n
///
/// H * (alpha) = (beta), H**H * H = I.
///      (  x  )  (  0  )
///
/// Here alpha and beta are scalars, and x is an (n-1)-element real vector.
/// H is represented in the form
/// H = I - tau * (1) * (1 v**H),
///              (v)
/// where tau is a real scalar and v is a real (n-1)-element vector.
///
/// If the elements of x are all zero, then tau = 0 and H is taken to be
/// the unit matrix.
///
/// Otherwise, 1 <= tau <= 2.
///
/// # Arguments
///
/// For arguments definitions, please refer to the original documentation.
pub fn dlarfg<X>(
    n: i32,
    alpha: &mut f64,
    x: &mut X,
    incx: i32,
    mut tau: &mut f64,
) where
    X: AsFortranArray + From<FortranArray>,
{
    if n <= 1 {
        *tau = 0.;
        return;
    }

    let mut xnorm = dnrm2(n - 1, x, incx);
    if xnorm == 0. {
        *tau = 0.;
    } else {
        let mut beta = -alpha.signum() * dlapy2(*alpha, xnorm);
        let safmin = dlamch('S') / dlamch('E');
        let mut knt = 0;

        if beta.abs() < safmin {
            let rsafmn = 1. / safmin;
            loop {
                knt += 1;
                dscal(n - 1, rsafmn, x, incx);
                beta *= rsafmn;
                *alpha *= rsafmn;
                if beta.abs() > safmin || knt > 20 {
                    break;
                }
            }

            xnorm = dnrm2(n - 1, x, incx);
            beta = -alpha.signum() * dlapy2(*alpha, xnorm);
        }

        *tau = (beta - *alpha) / beta;
        dscal(n - 1, 1. / (*alpha - beta), x, incx);
        for _ in 0..knt { beta *= safmin; }
        *alpha = beta;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(1, 2., vec![3.], 1, 0.)]
    #[case(2, 2., vec![4., 5.], 2, 1.4472135954999579)]
    #[case(2, 1.6, vec![2.8, 3.6], 1, 1.4961389383568338)]
    #[case(2, 2.9, vec![3.2, 4.8], 2, 1.6715194247388592)]
    #[case(3, 2.1, vec![2.6, 7.8, 3.2], 2, 1.4538485592273604)]
    #[case(4, 3.4, vec![2.9, 1.0, 4.1, 5.2], 1, 1.5531563983709862)]
    #[case(5, 4.5, vec![3.3, 4.4, 5.5, 6.6, 7.7], 1, 1.403607734816732)]
    #[case(6, 6.8, vec![7.1, 8.2, 9.3, 1.4, 2.5, 3.6], 2, 1.4941176766392699)]
    #[case(7, 8.7, vec![4.8, 2.9, 7.1, 8.2, 9.3, 1.4, 2.5], 4, 1.6392567355572378)]
    fn test_dlarfg(
        #[case] n: i32,
        #[case] mut alpha: f64,
        #[case] mut x: Vec<f64>,
        #[case] incx: i32,
        #[case] expected_tau: f64,
    )
    {
        let mut tau = 0.;
        dlarfg(n, &mut alpha, &mut x, incx, &mut tau);
        assert_eq!(tau, expected_tau);
    }
}
