use crate::ieeeck::ieeeck;
use crate::iparmq::iparmq;

/// ILAENV
///
/// [Original] Online HTML documentation available at
/// `http://www.netlib.org/lapack/explore-html/`
///
/// ILAENV is called to choose problem-dependent parameters for the local environment.
/// See ISPEC for a description of the parameters.
///
/// # Arguments
///
/// For arguments definitions, please refer to the original documentation.
pub(crate) fn ilaenv(
    ispec: i32,
    name: &str,
    _opts: &str,
    n1: i32,
    n2: i32,
    n3: i32,
    n4: i32,
) -> i32 {
    match ispec {
        1..=3 => {
            let mut subnam = name.to_string();
            let ic = subnam.chars().nth(0).unwrap() as u32;
            let iz = 'Z' as u32;

            if iz == 90 || iz == 122 {
                // ASCII character set
                if ic >= 97 && ic <= 122 {
                    let new_char = char::from_u32(ic - 32).unwrap();
                    subnam.replace_range(0..1, &new_char.to_string());
                    for i in 1..6 {
                        let ic = subnam.chars().nth(i).unwrap() as u32;

                        if ic >= 97 && ic <= 122 {
                            let new_char = char::from_u32(ic - 32).unwrap();
                            subnam.replace_range(i..=i, &new_char.to_string());
                        }
                    }
                }
            } else if iz == 233 || iz == 169 {
                // EBCDIC character set
                if (129..=137).contains(&ic) ||
                    (145..=153).contains(&ic) ||
                    (162..=169).contains(&ic) {
                    let new_char = char::from_u32(ic + 64).unwrap();
                    subnam.replace_range(0..1, &new_char.to_string());

                    for i in 1..6 {
                        let ic = subnam.chars().nth(i).unwrap() as u32;

                        if (129..=137).contains(&ic) ||
                            (145..=153).contains(&ic) ||
                            (162..=169).contains(&ic) {
                            let new_char = char::from_u32(ic + 64).unwrap();
                            subnam.replace_range(i..=i, &new_char.to_string());
                        }
                    }
                }
            } else if iz == 218 || iz == 250 {
                // Prime machines: ASCII + 128
                if ic >= 225 && ic <= 250 {
                    let new_char = char::from_u32(ic - 32).unwrap();
                    subnam.replace_range(0..1, &new_char.to_string());
                    for i in 1..6 {
                        let ic = subnam.chars().nth(i).unwrap() as u32;
                        if ic >= 225 && ic <= 250 {
                            let new_char = char::from_u32(ic - 32).unwrap();
                            subnam.replace_range(i..=i, &new_char.to_string());
                        }
                    }
                }
            }

            let c1 = subnam.chars().nth(0).unwrap();
            let sname = c1 == 'S' || c1 == 'D';
            let cname = c1 == 'C' || c1 == 'Z';

            if !(cname || sname) {
                return 1;
            }

            let c2: String = subnam.chars().skip(1).take(2).collect();
            let c3: String = subnam.chars().skip(3).take(3).collect();
            let c4: String = c3.chars().skip(1).take(2).collect();
            let twostage = subnam.len() >= 11 && subnam.chars().nth(10).unwrap() == '2';

            match ispec {
                1 => match (c2.as_str(), c3.as_str()) {
                    _ if (c2 == "LA" && subnam.chars().skip(1).take(5).collect::<String>() == "LAORH")
                        || ((c2 == "SY" || c2 == "HE") && c3 == "TRF" && twostage) => 192,
                    _ if (c2 == "SY" || c2 == "HE") && c3 == "TRF" && !twostage => 64,
                    ("GE", "QRF") | ("GE", "RQF") | ("GE", "LQF") | ("GE", "QLF") | ("GE", "HRD") | ("GE", "BRD") => 32,
                    _ if c2 == "GE" && subnam.chars().skip(3).take(4).collect::<String>() == "QP3RK" => 32,
                    _ if (c2 == "OR" && sname) || (c2 == "UN" && cname)
                        && (c3.starts_with('G') || c3.starts_with('M'))
                        && ["QR", "RQ", "LQ", "QL", "HR", "TR", "BR"].contains(&c4.as_str()) => 32,
                    ("GB", "TRF") | ("PB", "TRF") if n4 > 64 || n2 > 64 => 32,
                    ("TR", "TRI") | ("TR", "EVC") => 64,
                    ("TR", "SYL") if sname => 48.max((n1.min(n2) * 16) / 100).min(240),
                    ("TR", "SYL") if !sname => 24.max((n1.min(n2) * 8) / 100).min(80),
                    ("LA", "UUM") | ("LA", "TRS") => 64,
                    ("ST", "EBZ") if sname => 1,
                    ("GG", _) if c3 == "HD3" => 32,
                    _ => 1,
                },
                2 => match (c2.as_str(), c3.as_str()) {
                    ("GE", "QRF") | ("GE", "RQF") | ("GE", "LQF") | ("GE", "QLF") | ("GE", "HRD") | ("GE", "BRD") => 2,
                    _ if c2 == "GE" && subnam.chars().skip(3).take(4).collect::<String>() == "QP3RK" => 2,
                    ("SY", "TRF") if sname || cname => 8,
                    ("SY", "TRD") | ("HE", "TRD") if sname || cname => 2,
                    _ if (c2 == "OR" && sname) || (c2 == "UN" && cname) => {
                        if c3.starts_with('G') || c3.starts_with('M')
                            && ["QR", "RQ", "LQ", "QL", "HR", "TR", "BR"].contains(&c4.as_str()) { 2 }
                        else { 1 }
                    },
                    ("GG", _) if c3 == "HD3" => 2,
                    _ => 2,
                },
                3 => match (c2.as_str(), c3.as_str()) {
                    ("GE", "QRF") | ("GE", "RQF") | ("GE", "LQF") | ("GE", "QLF") | ("GE", "HRD") | ("GE", "BRD") => 128,
                    _ if c2 == "GE" && subnam.chars().skip(3).take(4).collect::<String>() == "QP3RK" => 128,
                    ("SY", "TRD") | ("HE", "TRD") if sname || cname => 32,
                    ("GG", _) if c3 == "HD3" => 128,
                    _ if (c2 == "OR" && sname) || (c2 == "UN" && cname) =>
                        if c3.starts_with('G')
                            && ["QR", "RQ", "LQ", "QL", "HR", "TR", "BR"].contains(&c4.as_str()) { 128 }
                        else { 0 },
                    _ => 0,
                },
                _ => { -1 },
            }
        },
        4 => 6,
        5 => 2,
        6 => (n1.min(n2) as f64 * 1.6).round() as i32,
        7 => 1,
        8 => 50,
        9 => 25,
        10 => ieeeck(1, 0.0, 1.0),
        11 => ieeeck(0, 0.0, 1.0),
        12..=17 => iparmq(ispec, name, n2, n3),
        _ => -1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(-1, "CGGHRD", 0, 0, 0, 0, -1)]
    #[case(0, "CGGHRD", 0, 0, 0, 0, -1)]
    #[case(1, "CGGHRD", 0, 0, 0, 0, 1)]
    #[case(1, "DGGHRD", 0, 0, 0, 0, 1)]
    #[case(1, "SGBTRF", 0, 128, 0, 128, 32)]
    #[case(1, "DGBTRF", 0, 128, 0, 128, 32)]
    #[case(1, "DGEQRF", 0, 0, 0, 0, 32)]
    #[case(1, "SGEQRF", 0, 0, 0, 0, 32)]
    #[case(2, "CGGHRD", 0, 0, 0, 0, 2)]
    #[case(2, "DGGHRD", 0, 0, 0, 0, 2)]
    #[case(2, "SGBTRF", 0, 0, 0, 0, 2)]
    #[case(2, "DGBTRF", 0, 0, 0, 0, 2)]
    #[case(2, "DGEQRF", 0, 0, 0, 0, 2)]
    #[case(2, "SGEQRF", 0, 0, 0, 0, 2)]
    #[case(3, "CGGHRD", 0, 0, 0, 0, 0)]
    #[case(3, "DGGHRD", 0, 0, 0, 0, 0)]
    #[case(3, "SGBTRF", 0, 0, 0, 0, 0)]
    #[case(3, "DGBTRF", 0, 0, 0, 0, 0)]
    #[case(3, "DGEQRF", 0, 0, 0, 0, 128)]
    #[case(3, "SGEQRF", 0, 0, 0, 0, 128)]
    #[case(4, "CGGHRD", 2, 3, 0, 0, 6)]
    #[case(5, "CGGHRD", 2, 3, 0, 0, 2)]
    #[case(6, "CGGHRD", 2, 3, 0, 0, 3)]
    #[case(7, "CGGHRD", 2, 3, 0, 0, 1)]
    #[case(8, "CGGHRD", 2, 3, 0, 0, 50)]
    #[case(9, "CGGHRD", 2, 3, 0, 0, 25)]
    #[case(10, "CGGHRD", 0, 0, 0, 0, 0)]
    #[case(11, "CGGHRD", 0, 0, 0, 0, 1)]
    #[case(12, "CGGHRD", 2, 3, 4, 5, 75)]
    #[case(13, "CGGHRD", 2, 3, 4, 5, 2)]
    #[case(14, "CGGHRD", 2, 3, 4, 5, 14)]
    #[case(15, "CGGHRD", 2, 3, 4, 5, 2)]
    #[case(16, "CGGHRD", 2, 3, 4, 5, 1)]
    #[case(17, "CGGHRD", 2, 3, 4, 5, 10)]
    #[case(18, "CGGHRD", 0, 0, 0, 0, -1)]
    fn test_ispec_out_of_range(
        #[case] ispec: i32,
        #[case] name: &str,
        #[case] n1: i32,
        #[case] n2: i32,
        #[case] n3: i32,
        #[case] n4: i32,
        #[case] expected: i32,
    ) {
        assert_eq!(expected, ilaenv(ispec, name, "", n1, n2, n3, n4));
    }
}