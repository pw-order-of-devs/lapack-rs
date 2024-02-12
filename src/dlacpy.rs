use crate::array::{AsFortranArray, FortranArray};

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
/// For arguments definitions, please refer to the original documentation.
pub fn dlacpy<A, B>(
    uplo: char,
    m: i32,
    n: i32,
    a: &A,
    b: &mut B,
) where
    A: AsFortranArray,
    B: AsFortranArray + From<FortranArray>,
{
    let a_f = &a.to_fortran_array();
    let b_f = &mut b.to_fortran_array();

    match uplo {
        'U' => for j in 1..=n {
            for i in j..=m {
                b_f[(i, j)] = a_f[(i, j)];
            }
        },
        'L' => for j in 1..=n {
            for i in 1..=j {
                b_f[(i, j)] = a_f[(i, j)];
            }
        },
        _ => for j in 1..=n {
            for i in 1..=m {
                b_f[(i, j)] = a_f[(i, j)];
            }
        }
    }

    *b = B::from(b_f.clone());
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
