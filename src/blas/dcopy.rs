use crate::array::{convert::ToFortranArray, FortranArray};

/// DCOPY
///
/// # Documentation
///
/// [Original] Online html documentation available at
/// [http://www.netlib.org/lapack/explore-html/](http://www.netlib.org/lapack/explore-html/)
///
/// Copy a vector, x, to a vector, y.
/// Use unrolled loops for increments equal to 1.
///
/// # Arguments
///
/// * `n` - number of elements in input vector(s)
/// * `dx` - array
/// * `incx` - storage spacing between elements of dx
/// * `dy` - array
/// * `incy` - storage spacing between elements of dy
pub fn dcopy<DX, DY>(
    n: i32,
    dx: &DX,
    incx: i32,
    dy: &mut DY,
    incy: i32,
) where
    DX: ToFortranArray,
    DY: ToFortranArray + From<FortranArray>,
{
    let dx = dx.to_fa();
    let dy_f = &mut dy.to_fa();

    if n <= 0 { return; }
    if incx == 1 && incy == 1 {
        // Clean-up loop
        let m = n % 7;
        if m != 0 {
            for i in 1..=m { dy_f[i] = dx[i]; }
            if n < 7 {
                *dy = DY::from(dy_f.clone());
                return;
            }
        }
        for i in (m+1..=n).step_by(7) {
            dy_f[i] = dx[i];
            dy_f[i + 1] = dx[i + 1];
            dy_f[i + 2] = dx[i + 2];
            dy_f[i + 3] = dx[i + 3];
            dy_f[i + 4] = dx[i + 4];
            dy_f[i + 5] = dx[i + 5];
            dy_f[i + 6] = dx[i + 6];
        }
    } else {
        let mut ix = 1;
        let mut iy = 1;
        if incx < 0 { ix = (-n + 1) * incx + 1; }
        if incy < 0 { iy = (-n + 1) * incy + 1; }
        for _i in 1..=n {
            dy_f[iy] = dx[ix];
            ix += incx;
            iy += incy;
        }
    }

    *dy = DY::from(dy_f.clone());
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(0, vec![0.], 1, vec![0.], 1, vec![0.])]
    #[case(1, vec![1., 2., 3.], 1, vec![0., 0., 0.], 1, vec![1., 0., 0.])]
    #[case(2, vec![1.,2.,3.,4.,5.,6.,7.,8.], 2, vec![0., 0., 0., 0., 0., 0., 0., 0.], 2, vec![1., 0., 3., 0., 0., 0., 0., 0.])]
    #[case(3, vec![1.1, 2.2, 3.3, 4.4], 1, vec![0., 0., 0., 0.], 1, vec![1.1, 2.2, 3.3, 0.])]
    #[case(4, vec![5.5, 6.6, 7.7, 8.8], 2, vec![0., 0., 0., 0.], 3, vec![5.5, 0., 0., 7.7])]
    #[case(5, vec![9.9, 10.1, 11.1, 12.1], 4, vec![0., 0., 0., 0.], 4, vec![9.9, 0., 0., 0.])]
    fn test_dcopy(
        #[case] n: i32,
        #[case] dx: Vec<f64>,
        #[case] incx: i32,
        #[case] mut dy: Vec<f64>,
        #[case] incy: i32,
        #[case] expected: Vec<f64>,
    ) {
        dcopy(n, &dx, incx, &mut dy, incy);
        assert_eq!(expected, dy)
    }
}