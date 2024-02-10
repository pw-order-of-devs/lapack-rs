use crate::dlamch::dlamch;
use crate::dlapy2::dlapy2;

/// DLANV2
///
/// # Documentation
///
/// [Original] Online html documentation available at
/// [http://www.netlib.org/lapack/explore-html/](http://www.netlib.org/lapack/explore-html/)
///
/// Performs Schur factorization on a real 2x2 nonsymmetric matrix in standard form.
///
/// # Arguments
///
/// For arguments definitions, please refer to the original documentation.
pub fn dlanv2(
    a: &mut f64,
    b: &mut f64,
    c: &mut f64,
    d: &mut f64,
    mut rt1r: &mut f64,
    mut rt1i: &mut f64,
    mut rt2r: &mut f64,
    mut rt2i: &mut f64,
    mut cs: &mut f64,
    mut sn: &mut f64,
) {
    let multpl = 4.;
    let safmin = dlamch('S');
    let eps = dlamch('P');
    let safmn2 = dlamch('B').powf((safmin / eps).log(dlamch('B')) / 2.0).trunc();
    let safmx2 = 1.0 / safmn2;

    if *c == 0.0 {
        *cs = 1.0;
        *sn = 0.0;
    } else if *b == 0.0 {
        *cs = 0.0;
        *sn = 1.0;
        let temp = *d;
        *d = *a;
        *a = temp;
        *b = -(*c);
        *c = 0.0;
    } else if (*a - *d).abs() < f64::EPSILON
        && b.signum() != c.signum() {
        *cs = 1.0;
        *sn = 0.0;
    } else {
        let mut temp = *a - *d;
        let p = 0.5 * temp;
        let bcmax = b.abs().max(c.abs());
        let bcmis = b.abs().min(c.abs()) * b.signum() * c.signum();
        let scale = p.abs().max(bcmax);
        let mut z = (p / scale) * p + (bcmax / scale) * bcmis;

        if z >= multpl * eps {
            // Real eigenvalues. Compute A and D.
            z = p + (scale.sqrt() * z.sqrt()).copysign(p);
            *a = *d + z;
            *d = *d - (bcmax / z) * bcmis;

            // Compute B and the rotation matrix.
            let tau = dlapy2(*c, z);
            *cs = z / tau;
            *sn = *c / tau;
            *b = *b - *c;
            *c = 0.0;

            let mut count = 0;
            let mut sigma = *b + *c;

            loop {
                count += 1;
                let scale = temp.abs().max(sigma.abs());

                if scale >= safmx2 {
                    sigma *= safmn2;
                    temp *= safmn2;
                    if count <= 20 {
                        continue;
                    }
                }

                if scale <= safmn2 {
                    sigma *= safmx2;
                    temp *= safmx2;
                    if count <= 20 {
                        continue;
                    }
                }

                break;
            }

            let p = 0.5 * temp;
            let tau = dlapy2(sigma, temp);
            *cs = ((1f64 + sigma.abs() / tau) * 0.5).sqrt();
            *sn = -(p / (tau * *cs)).copysign(1f64 * sigma.signum());

            // Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
            //         [ CC  DD ]   [ C  D ] [ SN  CS ]
            let aa = *a * *cs + *b * *sn;
            let bb = -(*a * *sn) + *b * *cs;
            let cc = *c * *cs + *d * *sn;
            let dd = -(*c * *sn) + *d * *cs;

            // Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
            //         [ C  D ]   [-SN  CS ] [ CC  DD ]
            *a = aa * *cs + cc * *sn;
            *b = bb * *cs + dd * *sn;
            *c = -(aa * *sn) - cc * *cs;
            *d = -(bb * *sn) + dd * *cs;

            temp = 0.5 * (*a + *d);
            *a = temp;
            *d = temp;

            if *c != 0.0 {
                if *b != 0.0 {
                    if 1.0_f64.copysign(*b) == 1.0_f64.copysign(*c) {
                        // Real eigenvalues: reduce to upper triangular form
                        let sab = b.abs().sqrt();
                        let sac = c.abs().sqrt();
                        let p = (sab * sac).copysign(*c);
                        let tau = 1.0 / (*b + *c).abs().sqrt();
                        *a = temp + p;
                        *d = temp - p;
                        *b = *b - *c;
                        *c = 0.0;
                        let cs1 = sab * tau;
                        let sn1 = sac * tau;
                        let temp = *cs * cs1 - *sn * sn1;
                        *sn = *cs * sn1 + *sn * cs1;
                        *cs = temp;
                    }
                } else {
                    *b = -*c;
                    *c = 0.0;
                    let temp = *cs;
                    *cs = -(*sn);
                    *sn = temp;
                }
            }
        }
    }

    *rt1r = *a;
    *rt2r = *d;

    if *c == 0.0 {
        *rt1i = 0.0;
        *rt2i = 0.0;
    } else {
        *rt1i = b.abs().sqrt() * c.abs().sqrt();
        *rt2i = -(*rt1i);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(
        &mut 1.0, &mut 2.0, &mut 3.0, &mut 4.0, &mut 5.0, &mut 6.0, &mut 7.0, &mut 8.0, &mut 9.0, &mut 10.0,
        2.5, 2.066771435679323, -2.0298221281347035, 2.5, 2.5, 2.048213463957949, 2.5, -2.048213463957949, 0.8112421851755609, 0.5847102846637648,
    )]
    #[case(
        &mut 7.0, &mut 10.0, &mut 5.0, &mut 6.0, &mut 11.0, &mut 12.0, &mut 8.0, &mut 9.0, &mut 13.0, &mut 14.0,
        9.28858224651085, 5.11545816559625, 0., 3.7114177534891506, 9.28858224651085, 0., 3.7114177534891506, 0., 0.950617193895293, 0.3103658336071141,
    )]
    #[case(
        &mut 6.0, &mut 12.0, &mut 13.0, &mut 7.0, &mut 10.0, &mut 14.0, &mut 11.0, &mut 8.0, &mut 9.0, &mut 15.0,
        6.500000000000001, 7.985281374238571, -4.449747468305834, 6.500000000000001, 6.500000000000001, 5.960913149738705, 6.500000000000001, -5.960913149738705, 0.9238795325112867, 0.3826834323650898,
    )]
    #[case(
        &mut 10.0, &mut 15.0, &mut 16.0, &mut 17.0, &mut 11.0, &mut 14.0, &mut 18.0, &mut 12.0, &mut 13.0, &mut 19.0,
        13.5, 15.152043533532684, -12.935028842544403, 13.5, 13.5, 13.99971857323331, 13.5, -13.99971857323331, 0.7554539549957063, 0.6552017413601289,
    )]
    fn dlanv2_test(
        #[case] a: &mut f64,
        #[case] b: &mut f64,
        #[case] c: &mut f64,
        #[case] d: &mut f64,
        #[case] rt1r: &mut f64,
        #[case] rt1i: &mut f64,
        #[case] rt2r: &mut f64,
        #[case] rt2i: &mut f64,
        #[case] cs: &mut f64,
        #[case] sn: &mut f64,
        #[case] expected_a: f64,
        #[case] expected_b: f64,
        #[case] expected_c: f64,
        #[case] expected_d: f64,
        #[case] expected_rt1r: f64,
        #[case] expected_rt1i: f64,
        #[case] expected_rt2r: f64,
        #[case] expected_rt2i: f64,
        #[case] expected_cs: f64,
        #[case] expected_sn: f64,
    ) {
        let result = dlanv2(a, b, c, d, rt1r, rt1i, rt2r, rt2i, cs, sn);
        assert_eq!(expected_a, *a);
        assert_eq!(expected_b, *b);
        assert_eq!(expected_c, *c);
        assert_eq!(expected_d, *d);
        assert_eq!(expected_rt1r, *rt1r);
        assert_eq!(expected_rt1i, *rt1i);
        assert_eq!(expected_rt2r, *rt2r);
        assert_eq!(expected_rt2i, *rt2i);
        assert_eq!(expected_cs, *cs);
        assert_eq!(expected_sn, *sn);
    }
}
