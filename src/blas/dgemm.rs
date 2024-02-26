use crate::array::convert::ToFortranArray;
use crate::array::FortranArray;
use crate::blas::lsame::lsame;
use crate::xerbla::xerbla;

/// DGEMM
///
/// # Documentation
///
/// [Original] Online html documentation available at
/// [http://www.netlib.org/lapack/explore-html/](http://www.netlib.org/lapack/explore-html/)
///
/// Performs one of the matrix-matrix operations:
///
/// `C := alpha*op( A )*op( B ) + beta*C,`
///
/// where  `op( X )` is one of:
///
/// `op( X ) = X   or   op( X ) = X**T`
///
/// `alpha` and `beta` are scalars, and `A`, `B` and `C` are matrices, with `op( A )`
/// an `m` by `k` matrix,  `op( B )`  a  `k` by `n` matrix and  `C` an `m` by `n` matrix.
///
/// # Arguments
///
/// For arguments definitions, please refer to the original documentation.
pub fn dgemm<A, B, C>(
    transa: char,
    transb: char,
    m: i32,
    n: i32,
    k: i32,
    alpha: f64,
    a: &A,
    lda: i32,
    b: &B,
    ldb: i32,
    beta: f64,
    c: &mut C,
    ldc: i32,
) where
    A: ToFortranArray,
    B: ToFortranArray,
    C: ToFortranArray + From<FortranArray>,
{
    let a_f = &mut a.to_fa_2d(lda);
    let b_f = &mut b.to_fa_2d(ldb);
    let c_f = &mut c.to_fa_2d(ldc);

    let nota = lsame(transa, 'N');
    let notb = lsame(transb, 'N');
    let nrowa = if nota { m } else { k };
    let nrowb = if notb { k } else { n };

    let mut info = 0;

    if !nota && !lsame(transa, 'C') && !lsame(transa, 'T') {
        info = 1;
    } else if !notb && !lsame(transb, 'C') && !lsame(transb, 'T') {
        info = 2;
    } else if m < 0 {
        info = 3;
    } else if n < 0 {
        info = 4;
    } else if k < 0 {
        info = 5;
    } else if lda < 1.max(nrowa) {
        info = 8;
    } else if ldb < 1.max(nrowb) {
        info = 10;
    } else if ldc < 1.max(m) {
        info = 13;
    }

    if info != 0 {
        xerbla("DGEMM ", info);
        return;
    }

    // Quick return if possible.
    if m == 0 || n == 0 || ((alpha == 0. || k == 0) && beta == 1.) {
        return;
    }

    // And if  alpha.eq.zero.
    if alpha == 0. {
        if beta == 0. {
            for j in 1..=n {
                for i in 1..=m {
                    c_f[(i, j)] = 0.;
                }
            }
        } else {
            for j in 1..=n {
                for i in 1..=m {
                    c_f[(i, j)] *= beta;
                }
            }
        }
        *c = C::from(c_f.clone());
        return;
    }

    // Start the operations.
    if notb {
        if nota {
            // Form  C := alpha*A*B + beta*C.
            for j in 1..=n {
                if beta == 0. {
                    for i in 1..=m {
                        c_f[(i, j)] = 0.;
                    }
                } else if beta != 1. {
                    for i in 1..=m {
                        c_f[(i, j)] *= beta;
                    }
                }
                for l in 1..=k {
                    let temp = alpha * b_f[(l, j)];
                    for i in 1..=m {
                        c_f[(i, j)] += temp * a_f[(i, l)];
                    }
                }
            }
        } else {
            // Form  C := alpha*A**T*B + beta*C
            for j in 1..=n {
                for i in 1..=m {
                    let mut temp = 0.;
                    for l in 1..=k {
                        temp += a_f[(l, i)] * b_f[(l, j)];
                    }
                    if beta == 0. {
                        c_f[(i, j)] = alpha * temp;
                    } else {
                        c_f[(i, j)] = alpha * temp + beta * c_f[(i, j)];
                    }
                }
            }
        }
    } else if nota {
        // Form  C := alpha*A*B**T + beta*C
        for j in 1..=n {
            if beta == 0. {
                for i in 1..=m {
                    c_f[(i, j)] = 0.;
                }
            } else if beta != 1. {
                for i in 1..=m {
                    c_f[(i, j)] *= beta;
                }
            }
            for l in 1..=k {
                let temp = alpha * b_f[(j, l)];
                for i in 1..=m {
                    c_f[(i, j)] += temp * a_f[(i, l)];
                }
            }
        }
    } else {
        // Form  C := alpha*A**T*B**T + beta*C
        for j in 1..=n {
            for i in 1..=m {
                let mut temp = 0.;
                for l in 1..=k {
                    temp += a_f[(l, i)] * b_f[(j, l)];
                }
                if beta == 0. {
                    c_f[(i, j)] = alpha * temp;
                } else {
                    c_f[(i, j)] = alpha * temp + beta * c_f[(i, j)];
                }
            }
        }
    }

    *c = C::from(c_f.clone());
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(
        'T', 'T', 0.5f64, 0.5f64,
        vec![vec![3., 5.], vec![7., 11.]],
        vec![vec![13., 17.], vec![19., 23.]],
        vec![vec![67.5, 150.5], vec![83.5, 186.5]],
    )]
    #[case(
        'T', 'N', 0.5f64, 0.5f64,
        vec![vec![3., 5.], vec![7., 11.]],
        vec![vec![13., 17.], vec![19., 23.]],
        vec![vec![62.5, 139.5], vec![86.5, 193.5]],
    )]
    #[case(
        'N', 'T', 0.5f64, 0.5f64,
        vec![vec![109., 113.], vec![199., 401.]],
        vec![vec![503., 601.], vec![701., 809.]],
        vec![vec![97163.5, 168970.5], vec![113250.5, 196161.5]],
    )]
    #[case(
        'N', 'N', 0.5f64, 0.5f64,
        vec![vec![109., 113.], vec![199., 401.]],
        vec![vec![503., 601.], vec![701., 809.]],
        vec![vec![87213.5, 148920.5], vec![118700.5, 201811.5]],
    )]
    fn dgemm_test(
        #[case] transa: char,
        #[case] transb: char,
        #[case] alpha: f64,
        #[case] beta: f64,
        #[case] a: Vec<Vec<f64>>,
        #[case] b: Vec<Vec<f64>>,
        #[case] expected: Vec<Vec<f64>>,
    ) {
        let (m, n, k) = (2, 2, 2);
        let (lda, ldb, ldc) = (2, 2, 2);

        let c = &mut vec![vec![1.; 2]; 2];
        dgemm(transa, transb, m, n, k, alpha, &a, lda, &b, ldb, beta, c, ldc);
        assert_eq!(expected, c.clone());
    }
}
