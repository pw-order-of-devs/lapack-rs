/// DROT
///
/// # Documentation
///
/// [Original] Online html documentation available at
/// [http://www.netlib.org/lapack/explore-html/](http://www.netlib.org/lapack/explore-html/)
///
/// Applies a plane rotation.
///
/// The `drot` function applies a plane rotation.
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
pub fn drot(n: i32, dx: &mut Vec<f64>, incx: isize, dy: &mut Vec<f64>, incy: isize, c: f64, s: f64) {
    if n <= 0 { return; }

    if incx == 1 && incy == 1 {
        // code for both increments equal to 1
        for i in 0..n as usize {
            let dtemp = c * dx[i] + s * dy[i];
            dy[i] = c * dy[i] - s * dx[i];
            dx[i] = dtemp;
        }
    } else {
        // code for unequal increments or equal increments not equal to 1
        let mut ix = 0;
        let mut iy = 0;
        if incx < 0 { ix = (n as isize - 1) * incx + 1; }
        if incy < 0 { iy = (n as isize - 1) * incy + 1; }
        for _ in 0..n {
            let dtemp = c * dx[ix as usize] + s * dy[iy as usize];
            dy[iy as usize] = c * dy[iy as usize] - s * dx[ix as usize];
            dx[ix as usize] = dtemp;
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
    #[case(4, vec![1., 2., 3., 4.], 1, vec![5., 6., 7., 8.], 1, 2., 3., vec![17., 22., 27., 32.], vec![7., 6., 5., 4.])]
    #[case(4, vec![2., 4., 6., 8.], 1, vec![4., 5., 4., 5.], 1, 4., 4., vec![24., 36., 40., 52.], vec![8., 4., -8., -12.])]
    #[case(4, vec![0.5, 1.5, 2.5, 3.5], 1, vec![4.5, 5.5, 6.5, 7.5], 1, 2., 3., vec![14.5, 19.5, 24.5, 29.5], vec![7.5, 6.5, 5.5, 4.5])]
    #[case(2, vec![1., 3., 5., 7.], 1, vec![2., 4., 6., 8.], 1, 2., 2., vec![6., 14., 5., 7.], vec![2., 2., 6., 8.])]
    #[case(3, vec![0.1, 0.2, 0.3, 0.4], 1, vec![0.5, 0.6, 0.7, 0.8], 1, 3., 4., vec![2.3, 3., 3.6999999999999997, 0.4], vec![1.1, 0.9999999999999998, 0.8999999999999997, 0.8])]
    fn test_drot(
        #[case] n: i32,
        #[case] mut dx: Vec<f64>,
        #[case] incx: isize,
        #[case] mut dy: Vec<f64>,
        #[case] incy: isize,
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
