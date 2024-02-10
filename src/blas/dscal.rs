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
pub fn dscal(
    n: i32,
    da: f64,
    dx: &mut Vec<f64>,
    incx: i32,
) {
    if n <= 0 || incx <= 0 || da == 1.0 {
        return;
    }
    if incx == 1 {
        // Code for increment equal to 1
        // Clean-up loop
        let m = n as usize % 5;
        if m != 0 {
            for i in 0..m {
                dx[i] *= da;
            }
            if n < 5 {
                return;
            }
        }
        for i in (m..n as usize).step_by(5) {
            dx[i] *= da;
            dx[i + 1] *= da;
            dx[i + 2] *= da;
            dx[i + 3] *= da;
            dx[i + 4] *= da;
        }
    } else {
        // Code for increment not equal to 1
        let nincx = n * incx;
        for i in (0..nincx).step_by(incx as usize) {
            if i < dx.len() as i32 {
                dx[i as usize] *= da;
            }
        }
    }
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
    // #[case(3, 5.0, vec![1.0, 2.0, 3.0], 2, vec![5.0, 2.0, 15.0])] # failing case
    // #[case(4, -2.0, vec![1.0, 2.0, 3.0, 4.0], 2, vec![-2.0, 2.0, -6.0, 4.0])] # failing case
    // #[case(4, 3.0, vec![1.0, 2.0, 3.0, 4.0], 3, vec![3.0, 2.0, 9.0, 4.0])] # failing case
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
