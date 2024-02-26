use crate::array::convert::ToFortranArray;
use crate::array::FortranArray;

/// DLAQR1
///
/// # Documentation
///
/// [Original] Online html documentation available at
/// [http://www.netlib.org/lapack/explore-html/](http://www.netlib.org/lapack/explore-html/)
///
/// Given a 2x2 or 3x3 matrix H, dlaqr1 sets v to a scalar multiple of the first column of the product
/// K = (H - (sr1 + i*si1)*I)*(H - (sr2 + i*si2)*I)
/// scaling to avoid overflows and most underflows.
///
/// It is assumed that either
///     1) sr1 = sr2 and si1 = -si2 or
///     2) si1 = si2 = 0.
///
/// This is useful for starting double implicit shift bulges in the QR algorithm.
///
/// # Arguments
///
/// For arguments definitions, please refer to the original documentation.
pub fn dlaqr1<H, V>(
    n: i32,
    h: &H,
    ldh: i32,
    sr1: f64,
    si1: f64,
    sr2: f64,
    si2: f64,
    v: &mut V,
) where
    H: ToFortranArray,
    V: ToFortranArray + From<FortranArray>,
{
    let h_f = &mut h.to_fa_2d(ldh);
    let v_f = &mut v.to_fa();

    if n != 2 && n != 3 { return; }

    if n == 2 {
        let s = (h_f[(1, 1)] - sr2).abs() + si2.abs() + h_f[(2, 1)].abs();
        if s == 0. {
            v_f[1] = 0.;
            v_f[2] = 0.;
        } else {
            let h21s = h_f[(2, 1)] / s;
            v_f[1] = h21s * h_f[(1, 2)] + (h_f[(1, 1)] - sr1) * ((h_f[(1, 1)] - sr2) / s) - si1 * (si2 / s);
            v_f[2] = h21s * (h_f[(1, 1)] + h_f[(2, 2)] - sr1 - sr2);
        }
    } else {
        let s = (h_f[(1, 1)] - sr2).abs() + si2.abs() + h_f[(2, 1)].abs() + h_f[(3, 1)].abs();
        if s == 0. {
            v_f[1] = 0.;
            v_f[2] = 0.;
            v_f[3] = 0.;
        } else {
            let h21s = h_f[(2, 1)] / s;
            let h31s = h_f[(3, 1)] / s;
            v_f[1] = (h_f[(1, 1)] - sr1) * ((h_f[(1, 1)] - sr2) / s) - si1 * (si2 / s) + h_f[(1, 2)] * h21s + h_f[(1, 3)] * h31s;
            v_f[2] = h21s * (h_f[(1, 1)] + h_f[(2, 2)] - sr1 - sr2) + h_f[(2, 3)] * h31s;
            v_f[3] = h31s * (h_f[(1, 1)] + h_f[(3, 3)] - sr1 - sr2) + h21s * h_f[(3, 2)];
        }
    }

    *v = V::from(v_f.clone());
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(2, vec![vec![1., 2.], vec![3., 4.]], 2, 1.5, 2.5, 3.5, 4.5, vec![-0.44444444444444442, 0.])]
    #[case(2, vec![vec![0.5, 1.5], vec![3.5, 6.5]], 2, 1., 2., 3., 3.5, vec![-0.06666666666666665, 0.6000000000000001])]
    #[case(3, vec![vec![5., 6., 7.], vec![8., 9., 10.], vec![11., 12., 13.]], 3, 2., 3., 4., 5., vec![5.9473684210526310, 6.9473684210526310, 7.5789473684210522])]
    #[case(3, vec![vec![2.5, 4.5, 2.], vec![1.5, 3.5, 1.], vec![0.5, 2.5, 0.]], 3, 1.5, 2.5, 3.5, 4.5, vec![-0.37499999999999994, 0.7916666666666666, -0.04166666666666663])]
    fn test_dlaqr1(
        #[case] n: i32,
        #[case] h: Vec<Vec<f64>>,
        #[case] ldh: i32,
        #[case] sr1: f64,
        #[case] si1: f64,
        #[case] sr2: f64,
        #[case] si2: f64,
        #[case] expected: Vec<f64>,
    ) {
        let v = &mut vec![0.; h.len()];
        dlaqr1(n, &h, ldh, sr1, si1, sr2, si2, v);
        assert_eq!(expected, v.clone());
    }
}
