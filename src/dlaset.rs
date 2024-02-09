/// DLASET
///
/// # Documentation
///
/// [Original] Online html documentation available at
/// [http://www.netlib.org/lapack/explore-html/](http://www.netlib.org/lapack/explore-html/)
///
/// Initializes the off-diagonal elements and the diagonal elements of a matrix to given values.
///
/// # Arguments
///
/// * `uplo` - Specifies the part of the matrix `A` to be set.
///   * 'U': Upper triangular part is set; the strictly lower
///          triangular part of `A` is not changed.
///   * 'L': Lower triangular part is set; the strictly upper
///          triangular part of `A` is not changed.
///   * Otherwise:  All the matrix `A` is set.
///
/// * `m` - The number of rows of the matrix `A`. `m` >= 0.
///
/// * `n` - The number of columns of the matrix `A`. `n` >= 0.
///
/// * `alpha` - The constant to which the offdiagonal elements are to be set.
///
/// * `beta` - The constant to which the diagonal elements are to be set.
///
/// * `a` - Matrix `A`, dimension (lda,n).
///   On exit, the leading m-by-n submatrix of `A` is set as follows:
///
///   if `uplo` = 'U', `A(i,j) = alpha`, 1 <= i <= j - 1, 1 <= j <= n,
///
///   if `uplo` = 'L', `A(i,j) = alpha`, j + 1 <= i <= m, 1 <= j <= n,
///
///   otherwise, `A(i,j) = alpha`, 1 <= i <= m, 1 <= j <= n, `i.ne.j`,
///
///   and, for all `uplo`, `A(i,i) = beta`, 1 <= i <= min(m,n).
pub fn dlaset(
    uplo: char,
    m: usize,
    n: usize,
    alpha: f64,
    beta: f64,
    a: &mut Vec<Vec<f64>>,
) {
    match uplo {
        'U' => for i in 0..n {
            for j in i..m {
                a[i][j] = alpha;
            }
        }
        'L' => for i in 0..n {
            for j in 0..=i {
                a[i][j] = alpha;
            }
        }
        _ => for j in 0..n {
            for i in 0..m {
                a[i][j] = alpha;
            }
        }
    }

    let min_dimension = n.min(m);
    for i in 0..min_dimension {
        a[i][i] = beta;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dlaset_upper() {
        let mut matrix = vec![
            vec![0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
        ];
        dlaset('U', 3, 3, 1.0, 2.0, &mut matrix);
        assert_eq!(matrix, vec![
            vec![2.0, 1.0, 1.0],
            vec![0.0, 2.0, 1.0],
            vec![0.0, 0.0, 2.0]
        ]);
    }

    #[test]
    fn test_dlaset_lower() {
        let mut matrix = vec![
            vec![0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
        ];
        dlaset('L', 3, 3, 1.0, 2.0, &mut matrix);
        assert_eq!(matrix, vec![
            vec![2.0, 0.0, 0.0],
            vec![1.0, 2.0, 0.0],
            vec![1.0, 1.0, 2.0]
        ]);
    }

    #[test]
    fn test_dlaset_all() {
        let mut matrix = vec![
            vec![0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
        ];
        dlaset('A', 3, 3, 1.0, 2.0, &mut matrix);
        assert_eq!(matrix, vec![
            vec![2.0, 1.0, 1.0],
            vec![1.0, 2.0, 1.0],
            vec![1.0, 1.0, 2.0]
        ]);
    }
}
