use crate::array::{convert::ToFortranArray, FortranArray};

/// DSCAL
///
/// # Documentation
///
/// [Original] Online html documentation available at
/// [http://www.netlib.org/lapack/explore-html/](http://www.netlib.org/lapack/explore-html/)
///
/// Scales a vector by a constant, using unrolled loops for increment equal to 1.
///
/// # Parameters
///
/// - `n`: The number of elements in input vector(s)
/// - `da`: On entry, `da` specifies the scalar alpha
/// - `dx`: `dx` is a mutable double precision array, with dimension `1 + (n - 1) * abs(incx)`
/// - `incx`: The storage spacing between elements of `dx`
///
/// The vector `dx` is scaled by `da`.
pub fn dscal<DX>(
    n: i32,
    da: f64,
    dx: &mut DX,
    incx: i32,
) where
    DX: ToFortranArray + From<FortranArray>,
{
    let dx_f = &mut dx.to_fortran_array();

    if n <= 0 || incx <= 0 || da == 1.0 { return; }
    if incx == 1 {
        // Code for increment equal to 1
        // Clean-up loop
        let m = n % 5;
        if m != 0 {
            for i in 1..=m { dx_f[i] *= da; }
            if n < 5 {
                *dx = DX::from(dx_f.clone());
                return;
            }
        }
        for i in (m+1..=n).step_by(5) {
            dx_f[i] *= da;
            dx_f[i + 1] *= da;
            dx_f[i + 2] *= da;
            dx_f[i + 3] *= da;
            dx_f[i + 4] *= da;
        }
    } else {
        // Code for increment not equal to 1
        let nincx = n * incx;
        let mut i = 1;
        while i < nincx {
            dx_f[i] *= da;
            i += incx;
        }
    }

    *dx = DX::from(dx_f.clone());
}

#[cfg(test)]
mod tests {
    use super::dscal;
    use rstest::rstest;

    #[rstest]
    #[case(0, 1.0, vec![], 1, vec![])]
    #[case(1, 2.0, vec![1.0], 1, vec![2.0])]
    #[case(2, -1.0, vec![1.0, 2.0], 1, vec![-1.0, -2.0])]
    #[case(5, 0.0, vec![1.0, 2.0, 3.0, 4.0, 5.0], 1, vec![0.0, 0.0, 0.0, 0.0, 0.0])]
    #[case(3, 5.0, vec![1.0, 2.0, 3.0], 2, vec![5.0, 2.0, 15.0])]
    #[case(4, -2.0, vec![1.0, 2.0, 3.0, 4.0], 2, vec![-2.0, 2.0, -6.0, 4.0])]
    #[case(4, 3.0, vec![1.0, 2.0, 3.0, 4.0], 3, vec![3.0, 2.0, 3.0, 12.0])]
    #[case(0, 1.0, vec![1.0, 2.0, 3.0, 4.0, 5.0], 2, vec![1.0, 2.0, 3.0, 4.0, 5.0])]
    #[case(5, 1.0, vec![1.0, 2.0, 3.0, 4.0, 5.0], 0, vec![1.0, 2.0, 3.0, 4.0, 5.0])]
    #[case(5, 0.0, vec![1.0, 2.0, 3.0, 4.0, 5.0], -1, vec![1.0, 2.0, 3.0, 4.0, 5.0])]
    fn test_dscal(
        #[case] n: i32,
        #[case] da: f64,
        #[case] dx: Vec<f64>,
        #[case] incx: i32,
        #[case] expected_dx: Vec<f64>,
    ) {
        let mut dx = dx.clone();
        dscal(n, da, &mut dx, incx);
        assert_eq!(dx, expected_dx);
    }
}
