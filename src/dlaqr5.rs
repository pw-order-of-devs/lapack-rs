use crate::array::convert::ToFortranArray;
use crate::array::FortranArray;
use crate::blas::dgemm::dgemm;
use crate::dlacpy::dlacpy;
use crate::dlamch::dlamch;
use crate::dlaqr1::dlaqr1;
use crate::dlarfg::dlarfg;
use crate::dlaset::dlaset;

/// DLAQR5
///
/// # Documentation
///
/// [Original] Online html documentation available at
/// [http://www.netlib.org/lapack/explore-html/](http://www.netlib.org/lapack/explore-html/)
///
/// Performs a single small-bulge multi-shift QR sweep.
///
/// # Arguments
///
/// For arguments definitions, please refer to the original documentation.
pub fn dlaqr5<SR, SI, H, Z, V, U, WV, WH>(
    wantt: bool,
    wantz: bool,
    kacc22: i32,
    n: i32,
    ktop: i32,
    kbot: i32,
    nshfts: i32,
    sr: &SR,
    si: &SI,
    h: &mut H,
    ldh: i32,
    iloz: i32,
    ihiz: i32,
    z: &mut Z,
    ldz: i32,
    v: &mut V,
    ldv: i32,
    u: &mut U,
    ldu: i32,
    nv: i32,
    wv: &mut WV,
    ldwv: i32,
    nh: i32,
    wh: &mut WH,
    ldwh: i32,
) where
    SR: ToFortranArray,
    SI: ToFortranArray,
    H: ToFortranArray + From<FortranArray>,
    Z: ToFortranArray + From<FortranArray>,
    V: ToFortranArray + From<FortranArray>,
    U: ToFortranArray + From<FortranArray>,
    WV: ToFortranArray + From<FortranArray>,
    WH: ToFortranArray + From<FortranArray>,
{
    let sr_f = &mut sr.to_fa();
    let si_f = &mut si.to_fa();
    let h_f = &mut h.to_fa_2d(ldh);
    let z_f = &mut z.to_fa_2d(ldz);
    let v_f = &mut v.to_fa_2d(ldv);
    let u_f = &mut u.to_fa_2d(ldu);
    let wv_f = &mut wv.to_fa_2d(ldwv);
    let wh_f = &mut wh.to_fa_2d(ldwh);

    let vt = &mut vec![0.; 3].to_fa();
    if nshfts < 2 { return; }
    if ktop >= kbot { return; }

    for i in (0..nshfts - 2).step_by(2) {
        if si_f[i] != si_f[i + 1] {
            let swap = sr_f[i];
            sr_f[i] = sr_f[i + 1];
            sr_f[i + 1] = sr_f[i + 2];
            sr_f[i + 2] = swap;

            let swap = si_f[i];
            si_f[i] = si_f[i + 1];
            si_f[i + 1] = si_f[i + 2];
            si_f[i + 2] = swap;
        }
    }

    let ns = nshfts - (nshfts % 2);
    let safmin = dlamch('S');
    let ulp = dlamch('P');
    let smlnum = safmin * ((n as f64) / ulp);
    let accum = kacc22 == 1 || kacc22 == 2;

    if ktop + 2 <= kbot {
        h_f[(ktop + 2, ktop)] = 0.;
    }

    let nbmps = ns / 2;
    let kdu = 4 * nbmps;

    for incol in ((ktop - 2*nbmps + 1)..=kbot-2).step_by(2*nbmps as usize) {
        let jtop;
        if accum { jtop = ktop.max(incol); }
        else if wantt { jtop = 1; }
        else { jtop = ktop; }

        let ndcol = incol + kdu;
        if accum { dlaset('A', kdu, kdu, 0., 1., u_f); }

        for krcol in incol..=(incol + 2 * nbmps - 1).min(kbot - 2) {
            let mtop = 1.max((ktop - krcol) / 2 + 1);
            let mbot = nbmps.min((kbot - krcol - 1) / 2);
            let m22 = mbot + 1;
            let bmp22 = mbot < nbmps && krcol + 2 * (m22 - 1) == kbot - 2;

            if bmp22 {
                let k = krcol + 2 * (m22 - 1);
                if k == ktop - 1 {
                    let h_tmp = h_f[(k + 1, k + 1)..].to_vec();
                    let v_tmp = &mut v_f[(1, m22)..].to_vec();
                    dlaqr1(2, &h_tmp, ldh, sr_f[2 * m22 - 1], si_f[2 * m22 - 1],
                           sr_f[2 * m22], si_f[2 * m22], v_tmp);
                    *v_f = v_tmp.to_fa_2d(ldv);

                    let beta = &mut v_f[(1, m22)].clone();
                    let (v_1, v_2) = (
                        &mut v_f[(2, m22)..].to_vec(),
                        &mut v_f[(1, m22)].clone());
                    dlarfg(2, beta, v_1, 1, v_2);
                    *v_f = v_1.to_fa();
                    v_f.insert(1, v_2.clone());
                    v_f.as_2d(ldv);
                } else {
                    let beta = &mut h_f[(k + 1, k)].clone();
                    v_f[(2, m22)] = h_f[(k + 2, k)].clone();
                    let (v_1, v_2) = (
                        &mut v_f[(2, m22)..].to_vec(),
                        &mut v_f[(1, m22)].clone());
                    dlarfg(2, beta, v_1, 1, v_2);
                    *v_f = v_1.to_fa();
                    v_f.insert(1, v_2.clone());
                    v_f.as_2d(ldv);

                    h_f[(k + 1, k)] = beta.clone();
                    h_f[(k + 2, k)] = 0.;
                }

                // Perform update from right within computational window.
                let t1 = v_f[(1, m22)];
                let t2 = t1 * v_f[(2, m22)];
                for j in jtop..=kbot.min(k + 3) {
                    let refsum = h_f[(j, k + 1)] + v_f[(2, m22)] * h_f[(j, k + 2)];
                    h_f[(j, k + 1)] -= refsum * t1;
                    h_f[(j, k + 2)] -= refsum * t2;
                }

                // Perform update from left within computational window.
                let mut jbot;

                if accum { jbot = ndcol.min(kbot); }
                else if wantt { jbot = n; }
                else { jbot = kbot; }

                let t1 = v_f[(1, m22)];
                let t2 = t1 * v_f[(2, m22)];

                for j in k + 1..=jbot {
                    let refsum = h_f[(k + 1, j)] + v_f[(2, m22)] * h_f[(k + 2, j)];
                    h_f[(k + 1, j)] -= refsum * t1;
                    h_f[(k + 2, j)] -= refsum * t2;
                }

                if k >= ktop {
                    if h_f[(k + 1, k)] != 0. {
                        let mut tst1 = h_f[(k, k)].abs() + h_f[(k + 1, k + 1)].abs();
                        if tst1 == 0. {
                            if k >= ktop + 1 { tst1 += h_f[(k, k - 1)].abs(); }
                            if k >= ktop + 2 { tst1 += h_f[(k, k - 2)].abs(); }
                            if k >= ktop + 3 { tst1 += h_f[(k, k - 3)].abs(); }
                            if k <= kbot - 2 { tst1 += h_f[(k + 2, k + 1)].abs(); }
                            if k <= kbot - 3 { tst1 += h_f[(k + 3, k + 1)].abs(); }
                            if k <= kbot - 4 { tst1 += h_f[(k + 4, k + 1)].abs(); }
                        }
                        if h_f[(k + 1, k)].abs() <= smlnum.max(ulp * tst1) {
                            let h12 = h_f[(k + 1, k)].abs().max(h_f[(k, k + 1)].abs());
                            let h21 = h_f[(k + 1, k)].abs().min(h_f[(k, k + 1)].abs());
                            let h11 = h_f[(k + 1, k + 1)].abs().max((h_f[(k, k)] - h_f[(k + 1, k + 1)]).abs());
                            let h22 = h_f[(k + 1, k + 1)].abs().min((h_f[(k, k)] - h_f[(k + 1, k + 1)]).abs());
                            let scl = h11 + h12;
                            let tst2 = h22 * (h11 / scl);

                            if tst2 == 0. || h21 * (h12 / scl) <= (smlnum.max(ulp * tst2)) {
                                h_f[(k + 1, k)] = 0.;
                            }
                        }
                    }
                }

                // Accumulate orthogonal transformations

                if accum {
                    let kms = k - incol;
                    let t1 = v_f[(1, m22)];
                    let t2 = t1 * v_f[(2, m22)];
                    for j in 1.max(ktop - incol)..=kdu {
                        let refsum = u_f[(j, kms + 1)] + v_f[(2, m22)] * u_f[(j, kms + 2)];
                        u_f[(j, kms + 1)] -= refsum * t1;
                        u_f[(j, kms + 2)] -= refsum * t2;
                    }
                } else if wantz {
                    let t1 = v_f[(1, m22)];
                    let t2 = t1 * v_f[(2, m22)];
                    for j in iloz..=ihiz {
                        let refsum = z_f[(j, k + 1)] + v_f[(2, m22)] * z_f[(j, k + 2)];
                        z_f[(j, k + 1)] -= refsum * t1;
                        z_f[(j, k + 2)] -= refsum * t2;
                    }
                }
            }

            // Normal case: Chain of 3-by-3 reflections
            for m in (mtop..=mbot).rev() {
                let k = krcol + 2*(m - 1);
                if k == ktop - 1 {
                    let h_tmp = h_f[(ktop, ktop)..].to_vec();
                    let v_tmp = &mut v_f[(1, m)..].to_vec();
                    dlaqr1(3, &h_tmp, ldh, sr_f[2*m-1], si_f[2*m-1], sr_f[2*m], si_f[2*m], v_tmp);
                    *v_f = v_tmp.to_fa_2d(ldv);

                    let alpha = &mut v_f[(1, m)].clone();
                    let (v_1, v_2) = (
                        &mut v_f[(2, m)..].to_vec(),
                        &mut v_f[(1, m)].clone());
                    dlarfg(3, alpha, v_1, 1, v_2);
                    *v_f = v_1.to_fa();
                    v_f.insert(1, v_2.clone());
                    v_f.as_2d(ldv);
                } else {
                    // Perform delayed transformation of row below Mth bulge. Exploit the fact that the first two elements of the row are actually zero.
                    let t1 = v_f[(1, m)];
                    let t2 = t1 * v_f[(2, m)];
                    let t3 = t1 * v_f[(3, m)];
                    let refsum = v_f[(3, m)] * h_f[(k + 3, k + 2)];
                    h_f[(k + 3, k)] = -refsum * t1;
                    h_f[(k + 3, k + 1)] = -refsum * t2;
                    h_f[(k + 3, k + 2)] -= refsum * t3;

                    // Calculate reflection to move Mth bulge one step
                    let beta = &mut h_f[(k + 1, k)].clone();
                    v_f[(2, m)] = h_f[(k + 2, k)].clone();
                    v_f[(3, m)] = h_f[(k + 3, k)].clone();

                    let (v_1, v_2) = (
                        &mut v_f[(2, m)..].to_vec(),
                        &mut v_f[(1, m)].clone());
                    dlarfg(3, beta, v_1, 1, v_2);
                    *v_f = v_1.to_fa();
                    v_f.insert(1, v_2.clone());
                    v_f.as_2d(ldv);

                    if h_f[(k + 3, k)] != 0. || h_f[(k + 3, k + 1)] != 0. || h_f[(k + 3, k + 2)] == 0. {
                        // Typical case: not collapsed (yet).
                        h_f[(k + 1, k)] = beta.clone();
                        h_f[(k + 2, k)] = 0.;
                        h_f[(k + 3, k)] = 0.;
                    } else {
                        // Atypical case: collapsed.  Attempt to reintroduce ignoring H(K+1, K) and H(K+2, K).
                        // If the fill resulting from the new reflector is too large, then abandon it. Otherwise, use the new one.
                        dlaqr1(3, &h_f[(k+1, k+1)..].to_vec(), ldh, sr_f[2*m-1], si_f[2*m-1], sr_f[2*m], si_f[2*m], vt);

                        let alpha = &mut vt[1].clone();
                        let (vt_1, vt_2) = (
                            &mut vec![vt[2].clone()],
                            &mut vt[1].clone());
                        dlarfg(3, alpha, vt_1, 1, vt_2);
                        *vt = vt_1.to_fa();
                        vt.insert(1, vt_2.clone());

                        let t1 = vt[1];
                        let t2 = t1 * vt[2];
                        let t3 = t1 * vt[3];
                        let refsum = h_f[(k + 1, k)] + vt[2] * h_f[(k + 2, k)];
                        if (h_f[(k + 2, k)] - refsum * t2).abs() + refsum * t3.abs() > ulp * (h_f[(k, k)].abs() + h_f[(k + 1, k + 1)].abs() + h_f[(k + 2, k + 2)].abs()) {
                            // Starting a new bulge here would create non-negligible fill. Use the old one with trepidation.
                            h_f[(k + 1, k)] = beta.clone();
                            h_f[(k + 2, k)] = 0.;
                            h_f[(k + 3, k)] = 0.;
                        } else {
                            // Starting a new bulge here would create only negligible fill. Replace the old reflector with the new one.
                            h_f[(k + 1, k)] -= refsum * t1;
                            h_f[(k + 2, k)] = 0.;
                            h_f[(k + 3, k)] = 0.;
                            v_f[(1, m)] = vt[1];
                            v_f[(2, m)] = vt[2];
                            v_f[(3, m)] = vt[3];
                        }
                    }
                }

                // Apply reflection from the right and the first column of update from the left.
                let t1 = v_f[(1, m)];
                let t2 = t1 * v_f[(2, m)];
                let t3 = t1 * v_f[(3, m)];
                for j in jtop..=kbot.min(k + 3) {
                    let refsum = h_f[(j, k + 1)] + v_f[(2, m)] * h_f[(j, k + 2)] + v_f[(3, m)] * h_f[(j, k + 3)];
                    h_f[(j, k + 1)] -= refsum * t1;
                    h_f[(j, k + 2)] -= refsum * t2;
                    h_f[(j, k + 3)] -= refsum * t3;
                }

                let refsum = h_f[(k + 1, k + 1)] + v_f[(2, m)] * h_f[(k + 2, k + 1)] + v_f[(3, m)] * h_f[(k + 3, k + 1)];
                h_f[(k + 1, k + 1)] -= refsum * t1;
                h_f[(k + 2, k + 1)] -= refsum * t2;
                h_f[(k + 3, k + 1)] -= refsum * t3;

                if k < ktop { continue; }
                if h_f[(k + 1, k)] != 0. {
                    let mut tst1 = h_f[(k, k)].abs() + h_f[(k + 1, k + 1)].abs();
                    if tst1 == 0. {
                        if k >= ktop + 1 { tst1 += h_f[(k, k - 1)].abs(); }
                        if k >= ktop + 2 { tst1 += h_f[(k, k - 2)].abs(); }
                        if k >= ktop + 3 { tst1 += h_f[(k, k - 3)].abs(); }
                        if k <= kbot - 2 { tst1 += h_f[(k + 2, k + 1)].abs(); }
                        if k <= kbot - 3 { tst1 += h_f[(k + 3, k + 1)].abs(); }
                        if k <= kbot - 4 { tst1 += h_f[(k + 4, k + 1)].abs(); }
                    }
                    if h_f[(k + 1, k)].abs() <= smlnum.max(ulp * tst1) {
                        let h12 = h_f[(k + 1, k)].abs().max(h_f[(k, k + 1)].abs());
                        let h21 = h_f[(k + 1, k)].abs().min(h_f[(k, k + 1)].abs());
                        let h11 = h_f[(k + 1, k + 1)].abs().max((h_f[(k, k)] - h_f[(k + 1, k + 1)]).abs());
                        let h22 = h_f[(k + 1, k + 1)].abs().min((h_f[(k, k)] - h_f[(k + 1, k + 1)]).abs());
                        let scl = h11 + h12;
                        let tst2 = h22 * (h11 / scl);
                        if tst2 == 0. || h21 * (h12 / scl) <= smlnum.max(ulp * tst2) {
                            h_f[(k + 1, k)] = 0.;
                        }
                    }
                }
            }

            let jbot =
                if accum { ndcol.min(kbot) }
                else if wantt { n }
                else { kbot };

            for m in (mtop..=mbot).rev() {
                let k = krcol + 2 * (m - 1);
                let t1 = v_f[(1, m)];
                let t2 = t1 * v_f[(2, m)];
                let t3 = t1 * v_f[(3, m)];
                for j in ktop.max(krcol + 2 * m)..=jbot {
                    let refsum = h_f[(k + 1, j)] + v_f[(2, m)] * h_f[(k + 2, j)] + v_f[(3, m)] * h_f[(k + 3, j)];
                    h_f[(k + 1, j)] -= refsum * t1;
                    h_f[(k + 2, j)] -= refsum * t2;
                    h_f[(k + 3, j)] -= refsum * t3;
                }
            }

            if accum {
                // Accumulate U. (If needed, update Z later with an efficient matrix-matrix multiply.)
                for m in (mtop..=mbot).rev() {
                    let k = krcol + 2 * (m - 1);
                    let kms = k - incol;
                    let i2 = 1.max((ktop - incol).max(kms - (krcol - incol) + 1));
                    let i4 = kdu.min(krcol + 2 * (mbot - 1) - incol + 5);
                    let t1 = v_f[(1, m)];
                    let t2 = t1 * v_f[(2, m)];
                    let t3 = t1 * v_f[(3, m)];
                    for j in i2..=i4 {
                        let refsum = u_f[(j, kms + 1)] + v_f[(2, m)] * u_f[(j, kms + 2)] + v_f[(3, m)] * u_f[(j, kms + 3)];
                        u_f[(j, kms + 1)] -= refsum * t1;
                        u_f[(j, kms + 2)] -= refsum * t2;
                        u_f[(j, kms + 3)] -= refsum * t3;
                    }
                }
            } else if wantz {
                // U is not accumulated, so update Z now by multiplying by reflections from the right.
                for m in (mbot..=mtop).rev() {
                    let k = krcol + 2 * (m - 1);
                    let t1 = v_f[(1, m)];
                    let t2 = t1 * v_f[(2, m)];
                    let t3 = t1 * v_f[(3, m)];
                    for j in iloz..=ihiz {
                        let refsum = z_f[(j, k + 1)] + v_f[(2, m)] * z_f[(j, k + 2)] + v_f[(3, m)] * z_f[(j, k + 3)];
                        z_f[(j, k + 1)] -= refsum * t1;
                        z_f[(j, k + 2)] -= refsum * t2;
                        z_f[(j, k + 3)] -= refsum * t3;
                    }
                }
            }
        }

        // Multiply H by reflections from the left
        if accum {
            let (jtop, jbot) =
                if wantt { (1, n) }
                else { (ktop, kbot) };

            let k1 = 1.max(ktop - incol);
            let nu = (kdu - 0.max(ndcol - kbot)) - k1 + 1;

            // Horizontal Multiply
            for jcol in (ndcol.min(kbot) + 1..=jbot).step_by(nh as usize) {
                let jlen = nh.min(jbot - jcol + 1);

                dgemm('C', 'N', nu, jlen, nu, 1., &u_f[(k1, k1)..].to_vec(), ldu,
                      &h_f[(incol + k1, jcol)..].to_vec(), ldh, 0., wh_f, ldwh);
                let h_tmp = &mut h_f[(incol + k1, jcol)..].to_vec().to_fa_2d(ldh);
                dlacpy('A', nu, jlen, wh_f, ldwh, h_tmp, ldh);
                let ini_idx = h_f.idx((incol + k1, jcol)) as usize - 1;
                let mut h_ini = h_f.clone().to_vec()[..ini_idx].to_vec();
                h_ini.extend_from_slice(&h_tmp.clone().to_vec());
                *h_f = h_ini.to_fa_2d(ldz);
            }

            // Vertical multiply
            for jrow in (jtop..=ktop.max(incol) - 1).step_by(nv as usize) {
                let jlen = nv.min(ktop.max(incol) - jrow);

                dgemm('N', 'N', jlen, nu, nu, 1., &h_f[(jrow, incol + k1)..].to_vec(), ldh,
                      &u_f[(k1, k1)..].to_vec(), ldu, 0., wv_f, ldwv);
                let h_tmp = &mut h_f[(jrow, incol + k1)..].to_vec().to_fa_2d(ldh);
                dlacpy('A', jlen, nu, wv_f, ldwv, h_tmp, ldh);
                let ini_idx = h_f.idx((jrow, incol + k1)) as usize - 1;
                let mut h_ini = h_f.clone().to_vec()[..ini_idx].to_vec();
                h_ini.extend_from_slice(&h_tmp.clone().to_vec());
                *h_f = h_ini.to_fa_2d(ldz);
            }

            // Z multiply (also vertical)
            if wantz {
                for jrow in (iloz..=ihiz).step_by(nv as usize) {
                    let jlen = nv.min(ihiz - jrow + 1);

                    dgemm('N', 'N', jlen, nu, nu, 1., &z_f[(jrow, incol + k1)..].to_vec(), ldz,
                          &u_f[(k1, k1)..].to_vec(), ldu, 0., wv_f, ldwv);
                    let z_tmp = &mut z_f[(jrow, incol + k1)..].to_vec().to_fa_2d(ldz);
                    dlacpy('A', jlen, nu, wv_f, ldwv, z_tmp, ldz);
                    let ini_idx = z_f.idx((jrow, incol + k1)) as usize - 1;
                    let mut z_ini = z_f.clone().to_vec()[..ini_idx].to_vec();
                    z_ini.extend_from_slice(&z_tmp.clone().to_vec());
                    *z_f = z_ini.to_fa_2d(ldz);
                }
            }
        }
    }

    *h = H::from(h_f.clone());
    *z = Z::from(z_f.clone());
    *v = V::from(v_f.clone());
    *u = U::from(u_f.clone());
    *wv = WV::from(wv_f.clone());
    *wh = WH::from(wh_f.clone());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dlaqr5_1() {
        let sr = &mut vec![1., 2., 3.];
        let si = &mut vec![0., 0., 0.];
        let h = &mut vec![vec![2., 2., 2.], vec![1., 1., 1.], vec![3., 3., 3.]];
        let z = &mut vec![vec![3., 3., 3.], vec![2., 2., 2.], vec![4., 4., 4.]];
        let v = &mut vec![vec![2., 2., 2.], vec![1., 1., 1.], vec![3., 3., 3.]];
        let u = &mut vec![vec![3., 3., 3.], vec![2., 2., 2.], vec![4., 4., 4.]];
        let wv = &mut vec![vec![0., 0., 0.], vec![0., 0., 0.], vec![0., 0., 0.]];
        let wh = &mut vec![vec![0., 0., 0.], vec![0., 0., 0.], vec![0., 0., 0.]];

        dlaqr5(true, true, 1, 4, 1, 3, 2, sr, si, h, 3,
               1, 3, z, 3, v, 3, u, 3, 3, wv, 5, 3, wh, 3);

        let expected_sr = vec![1., 2., 3.];
        assert_eq!(expected_sr, sr.clone());
        let expected_si = vec![0., 0., 0.];
        assert_eq!(expected_si, si.clone());

        let expected_h = vec![
            vec![3.9999999999999973, 3.6742346141747659, 0.],
            vec![0.81649658092772626, 0.81481481481481510, 5.2378280087892359E-002],
            vec![-2.3094010767585029, -0.65472850109865499, 1.1851851851851853],
        ];
        assert_eq!(expected_h, h.clone());

        let expected_z = vec![
            vec![-4.949747468305832, -4.949747468305832, -4.949747468305832],
            vec![-1.7320508075688774, -1.7320508075688774, -1.7320508075688774],
            vec![5.230171827006869, 5.230171827006869, 5.230171827006869],
        ];
        assert_eq!(expected_z, z.clone());

        let expected_v = vec![
            vec![1.9622504486493764, -0.13870070824202946, 0.4142135623730951],
            vec![1., 1., 1.],
            vec![3., 3., 3.],
        ];
        assert_eq!(expected_v, v.clone());

        let expected_u = vec![
            vec![-0.7071067811865472, 0.0, -0.7071067811865475],
            vec![-0.19245008972987515, -0.9622504486493764, 0.19245008972987515],
            vec![0.6547285010986551, 0.27216552697590857, 0.6804138174397717],
        ];
        assert_eq!(expected_u, u.clone());

        let expected_wv = vec![
            vec![-4.949747468305832, -4.949747468305832, -4.949747468305832],
            vec![-1.7320508075688774, -1.7320508075688774, -1.7320508075688774],
            vec![5.230171827006869, 5.230171827006869, 5.230171827006869],
        ];
        assert_eq!(expected_wv, wv.clone());

        let expected_wh = vec![vec![0.; 3]; 3];
        assert_eq!(expected_wh, wh.clone());
    }

    #[test]
    fn test_dlaqr5_2() {
        let sr = &mut vec![1., 2., 3.];
        let si = &mut vec![0., 0., 0.];
        let h = &mut vec![vec![2., 2., 2.], vec![1., 1., 1.], vec![3., 3., 3.]];
        let z = &mut vec![vec![3., 3., 3.], vec![2., 2., 2.], vec![4., 4., 4.]];
        let v = &mut vec![vec![2., 2., 2.], vec![1., 1., 1.], vec![3., 3., 3.]];
        let u = &mut vec![vec![3., 3., 3.], vec![2., 2., 2.], vec![4., 4., 4.]];
        let wv = &mut vec![vec![0., 0., 0.], vec![0., 0., 0.], vec![0., 0., 0.]];
        let wh = &mut vec![vec![0., 0., 0.], vec![0., 0., 0.], vec![0., 0., 0.]];

        dlaqr5(true, true, 1, 3, 2, 3, 2, sr, si, h, 3,
               1, 3, z, 3, v, 3, u, 3, 3, wv, 3, 3, wh, 3);

        let expected_sr = vec![1., 2., 3.];
        assert_eq!(expected_sr, sr.clone());
        let expected_si = vec![0., 0., 0.];
        assert_eq!(expected_si, si.clone());

        let expected_h = vec![
            vec![2., 2., 2.],
            vec![-1.8973665961010278, 2.4000000000000012, -1.2000000000000002],
            vec![3.1460498941515413, -3.2000000000000006, 1.6000000000000001],
        ];
        assert_eq!(expected_h, h.clone());

        let expected_z = vec![
            vec![3., 3., 3.],
            vec![-3.1622776601683800, -3.1622776601683800, -3.1622776601683800],
            vec![4.3947331922020556, 4.3947331922020556, 4.3947331922020556],
        ];
        assert_eq!(expected_z, z.clone());

        let expected_v = vec![
            vec![1.9486832980505140, 0.16227766016837933, 2.],
            vec![1., 1., 1.],
            vec![3., 3., 3.],
        ];
        assert_eq!(expected_v, v.clone());

        let expected_u = vec![
            vec![-0.94868329805051399, -0.31622776601683794, 0.],
            vec![0.30000000000000010, 0.94868329805051377, 0.],
            vec![0.10000000000000001, 0., 1.],
        ];
        assert_eq!(expected_u, u.clone());

        let expected_wv = vec![
            vec![-3.1622776601683800, -3.1622776601683800, -3.1622776601683800],
            vec![ 4.3947331922020556, 4.3947331922020556, 4.3947331922020556],
            vec![0., 0., 0.],
        ];
        assert_eq!(expected_wv, wv.clone());

        let expected_wh = vec![vec![0.; 3]; 3];
        assert_eq!(expected_wh, wh.clone());
    }
}