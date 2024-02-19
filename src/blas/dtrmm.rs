use crate::array::convert::ToFortranArray;
use crate::array::FortranArray;
use crate::blas::lsame::lsame;
use crate::xerbla::xerbla;

/// DTRMM
///
/// # Documentation
///
/// [Original] Online html documentation available at
/// [http://www.netlib.org/lapack/explore-html/](http://www.netlib.org/lapack/explore-html/)
///
/// Performs one of the matrix-matrix operations
///
/// ```text
/// B := alpha*op( A )*B,   or   B := alpha*B*op( A )
/// ```
///
/// where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
/// non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
///
/// ```text
/// op( A ) = A   or   op( A ) = A**T.
/// ```
///
/// On entry,  SIDE specifies whether  op( A ) multiplies B from
/// the left or right as follows:
///
/// ```text
/// SIDE = 'L' or 'l'   B := alpha*op( A )*B.
///
/// SIDE = 'R' or 'r'   B := alpha*B*op( A ).
/// ```
pub fn dtrmm<A, B>(
    side: char,
    uplo: char,
    transa: char,
    diag: char,
    m: i32,
    n: i32,
    alpha: f64,
    a: &mut A,
    lda: i32,
    b: &mut B,
    ldb: i32,
) where
    A: ToFortranArray + From<FortranArray>,
    B: ToFortranArray + From<FortranArray>,
{
    let a_f = &mut a.to_fortran_array();
    let b_f = &mut b.to_fortran_array();

    let lside = lsame(side, 'L');
    let nrowa = if lside { m } else { n };

    let nounit = lsame(diag, 'N');
    let upper = lsame(uplo, 'U');
    let mut info = 0;
    if  !lside && !lsame(side, 'R') {
        info = 1;
    } else if !upper && !lsame(uplo, 'L') {
        info = 2;
    } else if !lsame(transa, 'N') && !lsame(transa, 'T') && !lsame(transa, 'C') {
        info = 3;
    } else if !lsame(diag, 'U') && !lsame(diag, 'N') {
        info = 4;
    } else if m < 0 {
        info = 5;
    } else if n < 0 {
        info = 6;
    } else if lda < 1.max(nrowa) {
        info = 9;
    } else if ldb < 1.max(m) {
        info = 11;
    }

    if info != 0 {
        xerbla("DTRMM ", info);
        return;
    }

    if m == 0 || n == 0 {
        return;
    }

    if alpha == 0. {
        for j in 1..=n {
            for i in 1..=m {
                b_f[(i, j)] = 0.;
            }
        }
        *b = B::from(b_f.clone());
        return;
    }

    if lside {
        if lsame(transa, 'N') {
            if upper {
                for j in 1..=n {
                    for k in 1..=m {
                        if b_f[(k, j)] != 0. {
                            let mut temp = alpha * b_f[(k, j)];
                            for i in 1..k {
                                b_f[(i, j)] += temp * a_f[(i, k)];
                            }
                            if nounit {
                                temp *= a_f[(k, k)];
                            }
                            b_f[(k, j)] = temp;
                        }
                    }
                }
            } else {
                for j in 1..=n {
                    for k in (1..=m).rev() {
                        if b_f[(k, j)] != 0. {
                            let temp = alpha * b_f[(k, j)];
                            b_f[(k, j)] = temp;
                            if nounit {
                                b_f[(k, j)] *= a_f[(k, k)];
                            }
                            for i in k+1..=m {
                                b_f[(i, j)] += temp * a_f[(i, k)];
                            }
                        }
                    }
                }
            }
        } else {
            if upper {
                for j in 1..=n {
                    for i in (1..=m).rev() {
                        let mut temp = b_f[(i, j)];
                        if nounit {
                            temp *= a_f[(i, i)];
                        }
                        for k in 1..i {
                            temp += a_f[(k, i)] * b_f[(k, j)];
                        }
                        b_f[(i, j)] = alpha * temp;
                    }
                }
            } else {
                for j in 1..=n {
                    for i in 1..=m {
                        let mut temp = b_f[(i, j)];
                        if nounit {
                            temp *= a_f[(i, i)];
                        }
                        for k in i + 1..=m {
                            temp += a_f[(k, i)] * b_f[(k, j)];
                        }
                        b_f[(i, j)] = alpha * temp;
                    }
                }
            }
        }
    } else {
        if lsame(transa, 'N') {
            if upper {
                for j in (1..=n).rev() {
                    let mut temp = alpha;
                    if nounit {
                        temp *= a_f[(j, j)];
                    }
                    for i in 1..=m {
                        b_f[(i, j)] = temp * b_f[(i, j)];
                    }
                    for k in 1..j {
                        if a_f[(k, j)] != 0. {
                            temp = alpha * a_f[(k, j)];
                            for i in 1..=m {
                                b_f[(i, j)] += temp * b_f[(i, k)];
                            }
                        }
                    }
                }
            } else {
                for j in 1..=n {
                    let mut temp = alpha;
                    if nounit {
                        temp *= a_f[(j, j)];
                    }
                    for i in 1..=m {
                        b_f[(i, j)] = temp * b_f[(i, j)];
                    }
                    for k in j + 1..=n {
                        if a_f[(k, j)] != 0. {
                            temp = alpha * a_f[(k, j)];
                            for i in 1..=m {
                                b_f[(i, j)] += temp * b_f[(i, k)];
                            }
                        }
                    }
                }
            }
        } else {
            if upper {
                for k in 1..=n {
                    for j in 1..k {
                        if a_f[(j, k)] != 0. {
                            let temp = alpha * a_f[(j, k)];
                            for i in 1..=m {
                                b_f[(i, j)] = b_f[(i, j)] + temp * b_f[(i, k)];
                            }
                        }
                    }
                    let mut temp = alpha;
                    if nounit {
                        temp *= a_f[(k, k)];
                    }
                    if temp != 1. {
                        for i in 1..=m {
                            b_f[(i, k)] = temp * b_f[(i, k)];
                        }
                    }
                }
            } else {
                for k in (1..=n).rev() {
                    for j in k + 1..=n {
                        if a_f[(j, k)] != 0. {
                            let temp = alpha * a_f[(j, k)];
                            for i in 1..=m {
                                b_f[(i, j)] = b_f[(i, j)] + temp * b_f[(i, k)];
                            }
                        }
                    }
                    let mut temp = alpha;
                    if nounit {
                        temp *= a_f[(k, k)];
                    }
                    if temp != 1. {
                        for i in 1..=m {
                            b_f[(i, k)] = temp * b_f[(i, k)];
                        }
                    }
                }
            }
        }
    }

    *a = A::from(a_f.clone());
    *b = B::from(b_f.clone());
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case('L', 'U', 'N', 'N', 2, 2, 1., [[1., 2.], [0., 1.]], 2, [[3., 4.], [5., 6.]], 2, [[1., 2.], [0., 1.]], [[3., 4.], [5., 6.]])]
    #[case('R', 'U', 'N', 'U', 2, 2, 2., [[1., 0.], [2., 1.]], 2, [[1., 2.], [3., 4.]], 2, [[1., 0.], [2., 1.]], [[2., 4.], [10., 16.]])]
    #[case('L', 'L', 'T', 'N', 2, 2, 3., [[1., 0.], [1., 2.]], 2, [[1., 2.], [3., 4.]], 2, [[1., 0.], [1., 2.]], [[3., 12.], [9., 24.]])]
    #[case('R', 'L', 'N', 'N', 2, 2, 2., [[1., 1.], [0., 2.]], 2, [[1., 2.], [3., 4.]], 2, [[1., 1.], [0., 2.]], [[8., 12.], [12., 16.]])]
    #[case('L', 'U', 'N', 'U', 2, 2, 1., [[1., 2.], [0., 1.]], 2, [[0., 0.], [0., 0.]], 2, [[1., 2.], [0., 1.]], [[0., 0.], [0., 0.]])]
    #[case('R', 'U', 'T', 'N', 2, 2, 2., [[1., 0.], [2., 1.]], 2, [[0., 0.], [0., 0.]], 2, [[1., 0.], [2., 1.]], [[0., 0.], [0., 0.]])]
    fn test_dtrmm(
        #[case] side: char,
        #[case] uplo: char,
        #[case] transa: char,
        #[case] diag: char,
        #[case] m: i32,
        #[case] n: i32,
        #[case] alpha: f64,
        #[case] a: [[f64; 2]; 2],
        #[case] lda: i32,
        #[case] b: [[f64; 2]; 2],
        #[case] ldb: i32,
        #[case] expected_a: [[f64; 2]; 2],
        #[case] expected_b: [[f64; 2]; 2],
    ) {
        let a_mut = &mut a.iter().map(|r| r.to_vec()).collect::<Vec<Vec<f64>>>();
        let b_mut = &mut b.iter().map(|r| r.to_vec()).collect::<Vec<Vec<f64>>>();

        let expected_a = expected_a.iter().map(|r| r.to_vec()).collect::<Vec<Vec<f64>>>();
        let expected_b = expected_b.iter().map(|r| r.to_vec()).collect::<Vec<Vec<f64>>>();

        dtrmm(side, uplo, transa, diag, m, n, alpha, a_mut, lda, b_mut, ldb);

        assert_eq!(expected_a, a_mut.clone());
        assert_eq!(expected_b, b_mut.clone());
    }
}
