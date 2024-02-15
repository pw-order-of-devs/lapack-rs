use crate::array::{convert::ToFortranArray, FortranArray};

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
pub fn dlaset<A>(
    uplo: char,
    m: i32,
    n: i32,
    alpha: f64,
    beta: f64,
    a: &mut A,
) where
    A: ToFortranArray + From<FortranArray>,
{
    let a_f = &mut a.to_fortran_array();

    match uplo {
        'U' => for j in 1..=n {
            for i in j..=m {
                a_f[(i, j)] = alpha;
            }
        }
        'L' => for j in 1..=n {
            for i in 1..=j {
                a_f[(i, j)] = alpha;
            }
        }
        _ => for j in 1..=n {
            for i in 1..=m {
                a_f[(i, j)] = alpha;
            }
        }
    }

    for i in 1..=n.min(m) {
        a_f[(i, i)] = beta;
    }

    *a = A::from(a_f.clone());
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
