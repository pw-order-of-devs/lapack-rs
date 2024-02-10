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
pub fn dcopy(
    n: i32,
    dx: &Vec<f64>,
    incx: i32,
    mut dy: &mut Vec<f64>,
    incy: i32,
) {
    if n <= 0 { return; }
    if incx == 1 && incy == 1 {
        // Clean-up loop
        let m = n % 7;
        if m != 0 {
            for i in 0..m { dy[i as usize] = dx[i as usize]; }
            if n < 7 { return; }
        }
        for i in (m..n).step_by(7) {
            dy[i as usize] = dx[i as usize];
            dy[i as usize + 1] = dx[i as usize + 1];
            dy[i as usize + 2] = dx[i as usize + 2];
            dy[i as usize + 3] = dx[i as usize + 3];
            dy[i as usize + 4] = dx[i as usize + 4];
            dy[i as usize + 5] = dx[i as usize + 5];
            dy[i as usize + 6] = dx[i as usize + 6];
        }
    } else {
        let mut ix = 0;
        let mut iy = 0;
        if incx < 0 { ix = (-n + 1) * incx + 1; }
        if incy < 0 { iy = (-n + 1) * incy + 1; }
        for _i in 0..n {
            dy[iy as usize] = dx[ix as usize];
            ix += incx;
            iy += incy;
        }
    }
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
    #[case(4, vec![5.5, 6.6, 7.7, 8.8], 1, vec![0., 0., 0., 0.], 1, vec![5.5, 6.6, 7.7, 8.8])]
    // #[case(5, vec![9.9, 10.1, 11.1, 12.1], 1, vec![0., 0., 0., 0.], 1, vec![1.1, 2.2, 3.3, 0.])] # failing case
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