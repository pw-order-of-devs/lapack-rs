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
/// For arguments definitions, please refer to the original documentation.
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
