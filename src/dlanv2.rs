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
    rt1r: &mut f64,
    rt1i: &mut f64,
    rt2r: &mut f64,
    rt2i: &mut f64,
    cs: &mut f64,
    sn: &mut f64,
) {
    let multpl = 4.;
    let safmin = dlamch('S');
    let eps = dlamch('P');
    let safmn2 = dlamch('B').powf((safmin / eps).log(dlamch('B')) / 2.).trunc();
    let safmx2 = 1. / safmn2;

    if *c == 0. {
        *cs = 1.;
        *sn = 0.;
    } else if *b == 0. {
        *cs = 0.;
        *sn = 1.;
        let temp = *d;
        *d = *a;
        *a = temp;
        *b = -(*c);
        *c = 0.;
    } else if (*a - *d).abs() < f64::EPSILON
        && b.signum() != c.signum() {
        *cs = 1.;
        *sn = 0.;
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
            *c = 0.;
        } else {
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
            *sn = -(p / (tau * *cs)) * sigma.signum();

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
            *c = -(aa * *sn) + cc * *cs;
            *d = -(bb * *sn) + dd * *cs;

            temp = 0.5 * (*a + *d);
            *a = temp;
            *d = temp;

            if *c != 0. {
                if *b != 0. {
                    if b.signum() == c.signum() {
                        // Real eigenvalues: reduce to upper triangular form
                        let sab = b.abs().sqrt();
                        let sac = c.abs().sqrt();
                        let p = c.signum() * sab * sac;
                        let tau = 1. / (*b + *c).abs().sqrt();
                        *a = temp + p;
                        *d = temp - p;
                        *b = *b - *c;
                        *c = 0.;
                        let cs1 = sab * tau;
                        let sn1 = sac * tau;
                        let temp = *cs * cs1 - *sn * sn1;
                        *sn = *cs * sn1 + *sn * cs1;
                        *cs = temp;
                    }
                } else {
                    *b = -*c;
                    *c = 0.;
                    let temp = *cs;
                    *cs = -(*sn);
                    *sn = temp;
                }
            }
        }
    }

    *rt1r = *a;
    *rt2r = *d;

    if *c == 0. {
        *rt1i = 0.;
        *rt2i = 0.;
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
        &mut 1., &mut 2., &mut 3., &mut 4., &mut 5., &mut 6., &mut 7., &mut 8.,
        -0.3722813232690143, -1., 0., 5.372281323269014, -0.3722813232690143, 0., 5.372281323269014, 0., -0.824564840132393848, 0.56576746496899233
    )]
    #[case(
        &mut 7., &mut 10., &mut 5., &mut 6., &mut 11., &mut 12., &mut 8., &mut 9.,
        13.588723439378914, 5., 0., -0.58872343937891181, 13.588723439378914, 0., -0.58872343937891181, 0., 0.83504205669979814, 0.55018611718451338,
    )]
    #[case(
        &mut 6., &mut 12., &mut 13., &mut 7., &mut 10., &mut 14., &mut 11., &mut 8.,
        -6., -1., 0., 19., -6., 0., 19., 0., -0.7071067811865475, 0.7071067811865475,
    )]
    #[case(
        &mut -0.35413764185334262, &mut -0.56298456891977766, &mut 1.4015621918790930, &mut 0.51543733986900975, &mut 1., &mut 1., &mut 2., &mut 2.,
        0.08064984900783352, -0.37825080661482813, 1.5862959541840427, 0.08064984900783352,
        0.08064984900783352, 0.774608110078866, 0.08064984900783352, -0.774608110078866,
        0.9203697166177524, 0.39104933797790564,
    )]
    fn test_dlanv2(
        #[case] a: &mut f64,
        #[case] b: &mut f64,
        #[case] c: &mut f64,
        #[case] d: &mut f64,
        #[case] rt1r: &mut f64,
        #[case] rt1i: &mut f64,
        #[case] rt2r: &mut f64,
        #[case] rt2i: &mut f64,
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
        let (cs, sn) = (&mut 0., &mut 0.);
        dlanv2(a, b, c, d, rt1r, rt1i, rt2r, rt2i, cs, sn);
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
