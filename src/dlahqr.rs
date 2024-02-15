use crate::array::{convert::ToFortranArray, FortranArray};
use crate::blas::dcopy::dcopy;
use crate::blas::drot::drot;
use crate::dlamch::dlamch;
use crate::dlanv2::dlanv2;
use crate::dlarfg::dlarfg;

/// DLAHQR
///
/// # Documentation
///
/// [Original] Online html documentation available at
/// [http://www.netlib.org/lapack/explore-html/](http://www.netlib.org/lapack/explore-html/)
///
/// Computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using the double-shift/single-shift QR algorithm.
///
/// # Arguments
///
/// For arguments definitions, please refer to the original documentation.
pub fn dlahqr<H, WR, WI, Z>(
    wantt: bool,
    wantz: bool,
    n: i32,
    ilo: i32,
    ihi: i32,
    h: &mut H,
    ldh: i32,
    wr: &mut WR,
    wi: &mut WI,
    iloz: i32,
    ihiz: i32,
    z: &mut Z,
    _ldz: i32,
    info: &mut i32,
) where
    H: ToFortranArray + From<FortranArray>,
    WR: ToFortranArray + From<FortranArray>,
    WI: ToFortranArray + From<FortranArray>,
    Z: ToFortranArray + From<FortranArray>,
{
    let mut h_f = &mut h.to_fortran_array();
    let mut wr_f = &mut wr.to_fortran_array();
    let mut wi_f = &mut wi.to_fortran_array();
    let mut z_f = &mut z.to_fortran_array();

    *info = 0;
    let (dat1, dat2) = (3. / 4., -0.4375);
    let (mut cs, mut sn) = (0., 0.);
    let kexsh = 10;

    if n == 0 { return; }

    if ilo == ihi {
        wr_f[ilo] = h_f[(ilo, ilo)];
        wi_f[ilo] = 0.;

        *wr = WR::from(wr_f.clone());
        *wi = WI::from(wi_f.clone());
        return;
    }

    for j in ilo..=ihi-3 {
        h_f[(j + 2, j)] = 0.;
        h_f[(j + 3, j)] = 0.;
    }
    if ilo <= (ihi - 2) {
        h_f[(ihi, ihi-2)] = 0.;
    }

    let nh = ihi - ilo + 1;
    let nz = ihiz - iloz + 1;

    let safmin = dlamch('S');
    let ulp = dlamch('P');
    let smlnum = safmin * ((nh as f64) / ulp);

    let (mut i1, mut i2) = (0, 0);
    if wantt {
        i1 = 1;
        i2 = n;
    }

    let itmax = 30 * nh.max(10);
    let mut kdefl = 0;
    let mut i = ihi;

    let mut v = FortranArray::vector(&vec![0.; 3]);
    let mut m = 0;
    let mut t1 = 0.;
    let mut nr;
    let mut its = 0;
    let mut l;
    let mut k = i;

    loop { // 'loop_20:
        l = ilo;
        if i < ilo { // GO TO 160
            *h = H::from(h_f.clone());
            *wr = WR::from(wr_f.clone());
            *wi = WI::from(wi_f.clone());
            *z = Z::from(z_f.clone());
            return;
        }

        'loop_140: for curr_its in 0..=itmax {
            its = curr_its;
            for curr_k in (l + 1..=i).rev() { // 'loop_30:
                k = curr_k;
                if h_f[(curr_k, curr_k - 1)].abs() <= smlnum { break; }
                let mut tst = h_f[(curr_k - 1, curr_k - 1)].abs() + h_f[(curr_k, curr_k)].abs();
                if tst == 0. {
                    if (curr_k - 2) >= ilo { tst += h_f[(curr_k - 1, curr_k - 2)].abs(); }
                    if (curr_k + 1) <= ihi { tst += h_f[(curr_k + 1, curr_k)].abs(); }
                }
                if h_f[(curr_k, curr_k - 1)].abs() <= ulp * tst {
                    let ab = h_f[(curr_k, curr_k - 1)].abs().max(h_f[(curr_k - 1, curr_k)].abs());
                    let ba = h_f[(curr_k, curr_k - 1)].abs().min(h_f[(curr_k - 1, curr_k)].abs());
                    let aa = h_f[(curr_k, curr_k)].abs().max((h_f[(curr_k - 1, curr_k - 1)] - h_f[(curr_k, curr_k)]).abs());
                    let bb = h_f[(curr_k, curr_k)].abs().min((h_f[(curr_k - 1, curr_k - 1)] - h_f[(curr_k, curr_k)]).abs());
                    let s = aa + ab;
                    if ba * (ab / s) <= smlnum.max(ulp * (bb * (aa / s))) { break; }
                }
                k -= 1;
            }

            l = k;
            if l > ilo { h_f[(l, l - 1)] = 0.; }

            // go to 150
            if l >= i - 1 { break 'loop_140; }
            kdefl += 1;

            let mut h11;
            let mut h12;
            let mut h21;
            let mut h22;
            let mut i1 = 0;
            let mut i2 = 0;

            if wantt {
                i1 = l;
                i2 = i;
            }

            if kdefl % (2 * kexsh) == 0 {
                // Exceptional shift
                let s = h_f[(i, i - 1)].abs() + h_f[(i - 1, i - 2)].abs();
                h11 = dat1 * s + h_f[(i, i)];
                h12 = dat2 * s;
                h21 = s;
                h22 = h11;
            } else if kdefl % kexsh == 0 {
                // Exceptional shift
                let s = h_f[(l + 1, l)].abs() + h_f[(l + 2, l + 1)].abs();
                h11 = dat1 * s + h_f[(l, l)];
                h12 = dat2 * s;
                h21 = s;
                h22 = h11;
            } else {
                // Prepare to use Francis' double shift
                // (i.e. 2nd degree generalized Rayleigh quotient)
                h11 = h_f[(i - 1, i - 1)];
                h21 = h_f[(i, i - 1)];
                h12 = h_f[(i - 1, i)];
                h22 = h_f[(i, i)];
            }

            let s = h11.abs() + h12.abs() + h21.abs() + h22.abs();
            let (mut rt1r, rt1i, mut rt2r, rt2i);

            if s == 0. {
                rt1r = 0.;
                rt1i = 0.;
                rt2r = 0.;
                rt2i = 0.;
            } else {
                h11 /= s;
                h21 /= s;
                h12 /= s;
                h22 /= s;
                let tr = (h11 + h22) / 2.;
                let det = (h11 - tr) * (h22 - tr) - h12 * h21;
                let rtdisc = det.abs().sqrt();

                if det >= 0. {
                    // complex conjugate shifts
                    rt1r = tr * s;
                    rt2r = rt1r;
                    rt1i = rtdisc * s;
                    rt2i = -rt1i;
                } else {
                    // real shifts (use only one of them)
                    rt1r = tr + rtdisc;
                    rt2r = tr - rtdisc;
                    if (rt1r - h22).abs() <= (rt2r - h22).abs() {
                        rt1r *= s;
                        rt2r = rt1r;
                    } else {
                        rt2r *= s;
                        rt1r = rt2r;
                    }
                    rt1i = 0.;
                    rt2i = 0.;
                }
            }

            for curr_m in (l..=i - 2).rev() {
                m = curr_m;
                // Determine the effect of starting the double-shift QR
                // iteration at row M, and see if this would make H(M,M-1)
                // negligible.  (The following uses scaling to avoid
                // overflows and most underflows.)

                let mut h21s = h_f[(curr_m + 1, curr_m)];
                let mut s = (h_f[(curr_m, curr_m)]-rt2r).abs() + rt2i.abs() + h21s.abs();
                h21s = h_f[(curr_m + 1, curr_m)] / s;
                v[1] = h21s * h_f[(curr_m, curr_m + 1)] + (h_f[(curr_m, curr_m)]-rt1r) * ((h_f[(curr_m, curr_m)]-rt2r) / s) - rt1i * (rt2i / s);
                v[2] = h21s * (h_f[(curr_m, curr_m)] + h_f[(curr_m + 1, curr_m + 1)] - rt1r - rt2r);
                v[3] = h21s * h_f[(curr_m + 2, curr_m + 1)];
                s = v[1].abs() + v[2].abs() + v[3].abs();
                v[1] = v[1]/s;
                v[2] = v[2]/s;
                v[3] = v[3]/s;

                if curr_m == l {
                    // go to 60
                    break;
                }

                if (h_f[(curr_m, curr_m - 1)].abs() * (v[2].abs() + v[3].abs())) <=
                    (ulp * v[1].abs() * (h_f[(curr_m - 1, curr_m - 1)].abs() + h_f[(curr_m, curr_m)].abs() + h_f[(curr_m + 1, curr_m + 1)].abs())) {
                    // go to 60
                    break;
                }
            }

            for k in m..i {
                // The first iteration of this loop determines a reflection G
                // from the vector V and applies it from left and right to H,
                // thus creating a nonzero bulge below the subdiagonal.
                //
                // Each subsequent iteration determines a reflection G to
                // restore the Hessenberg form in the (k - 1)th column, and thus
                // chases the bulge one step toward the bottom of the active
                // submatrix. nr is the order of G.

                nr = 3.min(i-k+1);

                if k > m { dcopy(nr, &h_f[(k, k - 1)..].to_vec(), 1, &mut v, 1); }
                let (v_1, v_2) = (
                    &mut v[1].clone(),
                    &mut FortranArray::vector(&v[2..]));
                dlarfg(nr, v_1, v_2, 1, &mut t1);
                v = v_2.clone();
                v.insert(1, v_1.clone());

                if k > m {
                    h_f[(k, k - 1)] = v[1];
                    h_f[(k + 1, k - 1)] = 0.;
                    if k < i - 1 {
                        h_f[(k + 2, k - 1)] = 0.;
                    }
                } else if m > l {
                    // Use the following instead of h_f[(k, k - 1)] = -h_f[(k, k - 1)], to
                    // avoid a bug when v[1] and v[2] underflow. ====
                    h_f[(k, k - 1)] *= 1. - t1;
                }

                let v2 = v[2];
                let t2 = t1 * v2;
                let mut sum;

                if nr == 3 {
                    let v3 = v[3];
                    let t3 = t1 * v3;

                    // Apply G from the left to transform the rows of the matrix in columns K to I2.
                    for j in k..=i2 { // 'loop_70:
                        sum = h_f[(k, j)] + v2 * h_f[(k + 1, j)] + v3 * h_f[(k + 2, j)];
                        h_f[(k, j)] = h_f[(k, j)] - sum * t1;
                        h_f[(k + 1, j)] = h_f[(k + 1, j)] - sum * t2;
                        h_f[(k + 2, j)] = h_f[(k + 2, j)] - sum * t3;
                    }

                    // Apply G from the right to transform the columns of the matrix in rows I1 to min(K+3,I).
                    for j in i1..=i.min(k + 4) { // 'loop_80:
                        sum = h_f[(j, k)] + v2 * h_f[(j, k + 1)] + v3 * h_f[(j, k + 2)];
                        h_f[(j, k)] = h_f[(j, k)] - sum * t1;
                        h_f[(j, k + 1)] = h_f[(j, k + 1)] - sum * t2;
                        h_f[(j, k + 2)] = h_f[(j, k + 2)] - sum * t3;
                    }

                    if wantz {
                        // Accumulate transformations in the matrix Z
                        for j in iloz..=ihiz { // 'loop_90:
                            sum = z_f[(j, k)] + v2 * z_f[(j, k + 1)] + v3 * z_f[(j, k + 2)];
                            z_f[(j, k)] = z_f[(j, k)] - sum * t1;
                            z_f[(j, k + 1)] = z_f[(j, k + 1)] - sum * t2;
                            z_f[(j, k + 2)] = z_f[(j, k + 2)] - sum * t3;
                        }
                    }

                } else if nr == 2 {
                    // Apply G from the left to transform the rows of the matrix in columns K to i2
                    for j in k..=i2 { // 'loop_100:
                        sum = h_f[(k, j)] + v2 * h_f[(k + 1, j)];
                        h_f[(k, j)] = h_f[(k, j)] - sum * t1;
                        h_f[(k + 1, j)] = h_f[(k + 1, j)] - sum * t2;
                    }

                    // Apply G from the right to transform the columns of the matrix in rows I1 to min(K+3,I).
                    for j in i1..=i { // 'loop_110:
                        sum = h_f[(j, k)] + v2 * h_f[(j, k + 1)];
                        h_f[(j, k)] = h_f[(j, k)] - sum * t1;
                        h_f[(j, k + 1)] = h_f[(j, k + 1)] - sum * t2;
                    }

                    if wantz {
                        // Accumulate transformations in the matrix Z
                        for j in iloz..=ihiz { // 'loop_120:
                            sum = z_f[(j, k)] + v2 * z_f[(j, k + 1)];
                            z_f[(j, k)] = z_f[(j, k)] - sum * t1;
                            z_f[(j, k + 1)] = z_f[(j, k + 1)] - sum * t2;
                        }
                    }
                }
            }
        }

        // condition 150
        if its >= itmax {
            *info = i;

            *h = H::from(h_f.clone());
            *wr = WR::from(wr_f.clone());
            *wi = WI::from(wi_f.clone());
            *z = Z::from(z_f.clone());
            return;
        }

        if l == i {
            wr_f[i] = h_f[(i, i)];
            wi_f[i] = 0.;
        } else if l == i - 1 {
            let (h_00, h_01, h_10, h_11) = (
                &mut h_f[(i-1, i-1)].clone(), &mut h_f[(i-1, i)].clone(), &mut h_f[(i, i-1)].clone(), &mut h_f[(i, i)].clone());
            let (wr_0, wi_0, wr_1, wi_1) = (
                &mut wr_f[i-1].clone(), &mut wi_f[i-1].clone(), &mut wr_f[i].clone(), &mut wi_f[i].clone());
            dlanv2(h_00, h_01, h_10, h_11, wr_0, wi_0, wr_1, wi_1, &mut cs, &mut sn);

            h_f[(i-1, i-1)] = h_00.clone(); h_f[(i-1, i)] = h_01.clone(); h_f[(i, i-1)] = h_10.clone(); h_f[(i, i)] = h_11.clone();
            wr_f[i-1] = wr_0.clone(); wr_f[i] = wr_1.clone();
            wi_f[i-1] = wi_0.clone(); wi_f[i] = wi_1.clone();

            if wantt {
                if i2 > i {
                    let (h_02, h_12) = (&mut h_f[(i-1, i+1)..].to_vec(), &mut h_f[(i, i+1)..].to_vec());
                    drot(i2 - i, h_02, ldh, h_12, ldh, cs, sn);
                    let (h_02_idx, h_12_idx) = (h_f.idx((i-1, i+1)), h_f.idx((i, i+1)));
                    for idx in h_02_idx..=h_f.len() { h_f[idx] = h_02[(idx - h_02_idx) as usize]; }
                    for idx in h_12_idx..=h_f.len() { h_f[idx] = h_12[(idx - h_12_idx) as usize]; }
                }

                let (h_10, h_11) = (&mut h_f[(i1, i-1)..].to_vec(), &mut h_f[(i1, i)..].to_vec());
                drot(i - i1 - 1, &mut h_f[(i1, i-1)..].to_vec(), 1, &mut h_f[(i1, i)..].to_vec(), 1, cs, sn);
                let (h_10_idx, h_11_idx) = (h_f.idx((i1, i-1)), h_f.idx((i1, i)));
                for idx in h_10_idx..=h_f.len() { h_f[idx] = h_10[(idx - h_10_idx) as usize]; }
                for idx in h_11_idx..=h_f.len() { h_f[idx] = h_11[(idx - h_11_idx) as usize]; }
            }

            if wantz {
                let (z_00, z_01) = (&mut z_f[(iloz, i-1)..].to_vec(), &mut z_f[(iloz, i)..].to_vec());
                drot(nz, z_00, 1, z_01, 1, cs, sn);
                let (z_00_idx, z_01_idx) = (z_f.idx((iloz, i-1)), z_f.idx((iloz, i)));
                for idx in z_00_idx..=z_f.len() { z_f[idx] = z_00[(idx - z_00_idx) as usize]; }
                for idx in z_01_idx..=z_f.len() { z_f[idx] = z_01[(idx - z_01_idx) as usize]; }
            }
        }

        kdefl = 0;
        i = l - 1;

        *h = H::from(h_f.clone());
        *wr = WR::from(wr_f.clone());
        *wi = WI::from(wi_f.clone());
        *z = Z::from(z_f.clone());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(
        true, true, 3, 1, 3,
        &mut vec![vec![1., 2., 3.], vec![4., 5., 6.], vec![7., 8., 9.]], 3,
        &mut vec![1., 2., 3.], &mut vec![1., 2., 3.], 1, 3,
        &mut vec![vec![1., 2., 3.], vec![4., 5., 6.], vec![7., 8., 9.]], 3, &mut 0,
        &vec![vec![0.08064984900783352, 1.5862959541840427, 0.0], vec![-0.37825080661482813, 0.08064984900783352, 0.0], vec![-6.221356919060636, -3.7993781215208107, 14.838700301984334]],
        &vec![0.08064984900783352, 0.08064984900783352, 14.838700301984334], &vec![0.774608110078866, -0.774608110078866, 0.0],
        &vec![vec![-0.4651395556132531, -1.2543276506150056, -2.043515745616759], vec![0.34592653350082625, -0.15812789956359435, -0.6621823326280135], vec![8.10333141536393, 9.56042141917825, 11.017511422992568]],
    )]
    #[case(
        false, false, 4, 2, 4,
        &mut vec![vec![7., 8., 9., 10.], vec![11., 12., 13., 14.], vec![15., 16., 17., 18.], vec![19., 20., 21., 22.]], 4,
        &mut vec![4., 5., 6., 7.], &mut vec![4., 5., 6., 7.], 2, 4,
        &mut vec![vec![7., 8., 9., 10.], vec![11., 12., 13., 14.], vec![15., 16., 17., 18.], vec![19., 20., 21., 22.]], 4, &mut 0,
        &vec![vec![7., 8., 9., -3.4983360944002646e-21], vec![-2.0658540465189605, -2.016740798182532, 0., -3.4166091623768073e-28], vec![1.5544763283028376, -0.8003032123903298, 20.692597675634232, 0.0], vec![-26.463481449031356, -28.17468094616857, -35.37092977070222, -11.366325417353767]],
        &vec![4., -2.016740798182532, 20.692597675634232, -11.366325417353767], &vec![4., 0., 0., 0.],
        &mut vec![vec![7., 8., 9., 10.], vec![11., 12., 13., 14.], vec![15., 16., 17., 18.], vec![19., 20., 21., 22.]],
    )]
    #[case(
        false, true, 2, 1, 2,
        &mut vec![vec![2., 3.], vec![3., 2.]], 2,
        &mut vec![2., 3.], &mut vec![2., 3.], 1, 2,
        &mut vec![vec![2., 3.], vec![3., 2.]], 2, &mut 0,
        &vec![vec![5., 0.], vec![0., -1.0000000000000009]],
        &vec![5., -1.0000000000000009], &vec![0., 0.],
        &vec![vec![3.5355339059327378, 3.5355339059327378], vec![0.7071067811865472, -0.7071067811865479]],
    )]
    #[case(
        true, false, 3, 3, 3,
        &mut vec![vec![9., 8., 7.], vec![6., 5., 4.], vec![3., 2., 1.]], 3,
        &mut vec![3., 2., 1.], &mut vec![3., 2., 1.], 3, 3,
        &mut vec![vec![9., 8., 7.], vec![6., 5., 4.], vec![3., 2., 1.]], 3, &mut 0,
        &vec![vec![9.0, 8.0, 7.0], vec![6.0, 5.0, 4.0], vec![3.0, 2.0, 1.0]],
        &vec![3., 2., 1.], &vec![3., 2., 0.],
        &vec![vec![9.0, 8.0, 7.0], vec![6.0, 5.0, 4.0], vec![3.0, 2.0, 1.0]],
    )]
    fn test_dlahqr(
        #[case] wantt: bool,
        #[case] wantz: bool,
        #[case] n: i32,
        #[case] ilo: i32,
        #[case] ihi: i32,
        #[case] h: &mut Vec<Vec<f64>>,
        #[case] ldh: i32,
        #[case] wr: &mut Vec<f64>,
        #[case] wi: &mut Vec<f64>,
        #[case] iloz: i32,
        #[case] ihiz: i32,
        #[case] z: &mut Vec<Vec<f64>>,
        #[case] ldz: i32,
        #[case] info: &mut i32,
        #[case] expected_h: &Vec<Vec<f64>>,
        #[case] expected_wr: &Vec<f64>,
        #[case] expected_wi: &Vec<f64>,
        #[case] expected_z: &Vec<Vec<f64>>,
    ) {
        dlahqr(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, info);

        assert_eq!(expected_h, h);
        assert_eq!(expected_wr, wr);
        assert_eq!(expected_wi, wi);
        assert_eq!(expected_z, z);
    }
}
