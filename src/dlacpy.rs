/// DLACPY
///
/// # Documentation
///
/// [Original] Online html documentation available at
/// [http://www.netlib.org/lapack/explore-html/](http://www.netlib.org/lapack/explore-html/)
///
/// Copies all or part of a two-dimensional matrix `A` to another matrix `B`.
///
/// # Arguments
///
/// * `uplo` - Specifies the part of the matrix `A` to be copied to `B`.
///   * 'U': Upper triangular part
///   * 'L': Lower triangular part
///   * Otherwise: All the matrix `A`
///
/// * `m` - The number of rows of the matrix `A`. `m` >= 0.
///
/// * `n` - The number of columns of the matrix `A`. `n` >= 0.
///
/// * `a` - The `m` by `n` matrix `A`. If `uplo` = 'U', only the upper triangle
///   or trapezoid is accessed; if `uplo` = 'L', only the lower triangle or trapezoid is accessed.
///
/// * `b` - Output double precision array. On exit, `B` = `A` in the locations specified by `uplo`.
pub fn dlacpy(
    uplo: char,
    m: usize,
    n: usize,
    a: &Vec<Vec<f64>>,
    b: &mut Vec<Vec<f64>>,
) {
    match uplo {
        'U' => for i in 0..n {
            for j in i..m {
                b[i][j] = a[i][j];
            }
        },
        'L' => for i in 0..n {
            for j in 0..=i {
                b[i][j] = a[i][j];
            }
        },
        _ => for i in 0..n {
            for j in 0..m {
                b[i][j] = a[i][j];
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dlacpy_upper() {
        let a = vec![
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 9.0],
        ];
        let mut b = vec![
            vec![0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
        ];

        dlacpy('U', 3, 3, &a, &mut b);
        assert_eq!(b, vec![
            vec![1.0, 2.0, 3.0],
            vec![0.0, 5.0, 6.0],
            vec![0.0, 0.0, 9.0],
        ]);
    }

    #[test]
    fn test_dlacpy_lower() {
        let a = vec![
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 9.0],
        ];
        let mut b = vec![
            vec![0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
        ];

        dlacpy('L', 3, 3, &a, &mut b);
        assert_eq!(b, vec![
            vec![1.0, 0.0, 0.0],
            vec![4.0, 5.0, 0.0],
            vec![7.0, 8.0, 9.0],
        ]);
    }

    #[test]
    fn test_dlacpy_all() {
        let a = vec![
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 9.0],
        ];
        let mut b = vec![
            vec![0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
        ];

        dlacpy('A', 3, 3, &a, &mut b);
        assert_eq!(a, b);
    }
}
