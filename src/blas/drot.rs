use crate::array::{convert::ToFortranArray, FortranArray};

/// DROT
///
/// # Documentation
///
/// [Original] Online html documentation available at
/// [http://www.netlib.org/lapack/explore-html/](http://www.netlib.org/lapack/explore-html/)
///
/// Applies a plane rotation.
///
/// # Arguments
///
/// * `n` - Integer. The number of elements in input vector(s).
/// * `dx` - Vector of `f64` values. Dimension: (1 + (n - 1)*abs(incx)).
/// * `incx` - Integer. The storage spacing between elements of dx.
/// * `dy` - Vector of `f64` values. Dimension: (1 + (n - 1)*abs(incy)).
/// * `incy` - Integer. The storage spacing between elements of dy.
/// * `c` - `f64` scalar.
/// * `s` - `f64` scalar.
pub fn drot<DX, DY>(
    n: i32,
    dx: &mut DX,
    incx: i32,
    dy: &mut DY,
    incy: i32,
    c: f64,
    s: f64,
) where
    DX: ToFortranArray + From<FortranArray>,
    DY: ToFortranArray + From<FortranArray>,
{
    let dx_f = &mut dx.to_fortran_array();
    let dy_f = &mut dy.to_fortran_array();

    if n <= 0 { return; }
    if incx == 1 && incy == 1 {
        // code for both increments equal to 1
        for i in 1..=n {
            let dtemp = c * dx_f[i] + s * dy_f[i];
            dy_f[i] = c * dy_f[i] - s * dx_f[i];
            dx_f[i] = dtemp;
        }
    } else {
        // code for unequal increments or equal increments not equal to 1
        let mut ix = 1;
        let mut iy = 1;
        if incx < 0 { ix = (-n + 1) * incx + 1; }
        if incy < 0 { iy = (-n + 1) * incy + 1; }
        for _ in 1..=n {
            let dtemp = c * dx_f[ix] + s * dy_f[iy];
            dy_f[iy] = c * dy_f[iy] - s * dx_f[ix];
            dx_f[ix] = dtemp;
            ix += incx;
            iy += incy;
        }
    }

    *dx = DX::from(dx_f.clone());
    *dy = DY::from(dy_f.clone());
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(4, vec![1., 2., 3., 4.], 1, vec![5., 6., 7., 8.], 1, 2., 3., vec![17., 22., 27., 32.], vec![7., 6., 5., 4.])]
    #[case(4, vec![2., 4., 6., 8.], 1, vec![4., 5., 4., 5.], 1, 4., 4., vec![24., 36., 40., 52.], vec![8., 4., -8., -12.])]
    #[case(4, vec![0.5, 1.5, 2.5, 3.5], 1, vec![4.5, 5.5, 6.5, 7.5], 1, 2., 3., vec![14.5, 19.5, 24.5, 29.5], vec![7.5, 6.5, 5.5, 4.5])]
    #[case(2, vec![1., 3., 5., 7.], 1, vec![2., 4., 6., 8.], 1, 2., 2., vec![6., 14., 5., 7.], vec![2., 2., 6., 8.])]
    #[case(3, vec![0.1, 0.2, 0.3, 0.4], 1, vec![0.5, 0.6, 0.7, 0.8], 1, 3., 4., vec![2.3, 3., 3.6999999999999997, 0.4], vec![1.1, 0.9999999999999998, 0.8999999999999997, 0.8])]
    fn test_drot(
        #[case] n: i32,
        #[case] mut dx: Vec<f64>,
        #[case] incx: i32,
        #[case] mut dy: Vec<f64>,
        #[case] incy: i32,
        #[case] c: f64,
        #[case] s: f64,
        #[case] expected_dx: Vec<f64>,
        #[case] expected_dy: Vec<f64>,
    ) {
        drot(n, &mut dx, incx, &mut dy, incy, c, s);
        assert_eq!(expected_dx, dx);
        assert_eq!(expected_dy, dy);
    }
}
