use crate::errors::LapackError;
use crate::ieeeck::ieeeck;
use crate::iparmq::iparmq;

/// ILAENV
///
/// [Original] Online HTML documentation available at
/// `http://www.netlib.org/lapack/explore-html/`
///
/// # Definition
/// `fn ilaenv(ispec: i32, name: &str, opts: &str, n1: i32, n2: i32, n3: i32, n4: i32) -> Result<i32, LapackError>`
///
/// # Arguments
/// * `ispec: i32`
/// Specifies the parameter to be returned as the value of ILAENV.
/// For more detail on meaning of exact values refer to original documentation.
/// * `name: &str`
/// The name of the calling subroutine, in either upper case or lower case.
/// * `opts: &str`
/// The character options to the subroutine NAME, concatenated into a single character string.
/// * `n1: i32`
/// * `n2: i32`
/// * `n3: i32`
/// * `n4: i32`
/// Problem dimensions for the subroutine NAME; these may not all be required.
///
/// # Returns
/// `Result<i32, LapackError>` Result type returning i32 on success and LapackError on failure.
///
/// # Further Details
/// The following conventions have been used when calling ILAENV from the
/// LAPACK routines:
/// 1) OPTS is a concatenation of all the character options to
/// subroutine NAME, in the same order that they appear in the
/// argument list for NAME, even if they are not used in determining
/// the value of the parameter specified by ISPEC.
/// 2) The problem dimensions N1, N2, N3, N4 are specified in the order
/// that they appear in the argument list for NAME.  N1 is used
/// first, N2 second, and so on, and unused problem dimensions are
/// passed a value of -1.
/// 3) The parameter value returned by ILAENV is checked for validity in
/// the calling subroutine.
pub(crate) fn ilaenv(
    ispec: i32,
    name: &str,
    opts: &str,
    n1: i32,
    n2: i32,
    n3: i32,
    n4: i32,
) -> Result<i32, LapackError> {
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
                // Prime machines: ASCII + Ok(128
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
                return Ok(1);
            }

            let c2: String = subnam.chars().skip(1).take(2).collect();
            let c3: String = subnam.chars().skip(3).take(3).collect();
            let c4: String = c3.chars().skip(1).take(2).collect();
            let twostage = subnam.len() >= 11 && subnam.chars().nth(10).unwrap() == '2';

            println!("{c2}");
            println!("{c3}");
            println!("{c4}");
            match ispec {
                1 => match (c2.as_str(), c3.as_str()) {
                    _ if (c2 == "LA" && subnam.chars().skip(1).take(5).collect::<String>() == "LAORH")
                        || ((c2 == "SY" || c2 == "HE") && c3 == "TRF" && twostage) => Ok(192),
                    _ if (c2 == "SY" || c2 == "HE") && c3 == "TRF" && !twostage => Ok(64),
                    ("GE", "QRF") | ("GE", "RQF") | ("GE", "LQF") | ("GE", "QLF") | ("GE", "HRD") | ("GE", "BRD") => Ok(32),
                    _ if c2 == "GE" && subnam.chars().skip(3).take(4).collect::<String>() == "QP3RK" => Ok(32),
                    _ if (c2 == "OR" && sname) || (c2 == "UN" && cname)
                        && (c3.starts_with('G') || c3.starts_with('M'))
                        && ["QR", "RQ", "LQ", "QL", "HR", "TR", "BR"].contains(&c4.as_str()) => Ok(32),
                    ("GB", "TRF") | ("PB", "TRF") if n4 > 64 || n2 > 64 => Ok(32),
                    ("TR", "TRI") | ("TR", "EVC") => Ok(64),
                    ("TR", "SYL") if sname => Ok(std::cmp::min(std::cmp::max(48, (std::cmp::min(n1, n2) * 16) / 100), 240)),
                    ("TR", "SYL") if !sname => Ok(std::cmp::min(std::cmp::max(24, (std::cmp::min(n1, n2) * 8) / 100), 80)),
                    ("LA", "UUM") | ("LA", "TRS") => Ok(64),
                    ("ST", "EBZ") if sname => Ok(1),
                    ("GG", _) if c3 == "HD3" => Ok(32),
                    _ => Ok(1),
                },
                2 => match (c2.as_str(), c3.as_str()) {
                    ("GE", "QRF") | ("GE", "RQF") | ("GE", "LQF") | ("GE", "QLF") | ("GE", "HRD") | ("GE", "BRD") => Ok(2),
                    _ if c2 == "GE" && subnam.chars().skip(3).take(4).collect::<String>() == "QP3RK" => Ok(2),
                    ("SY", "TRF") if sname || cname => Ok(8),
                    ("SY", "TRD") | ("HE", "TRD") if sname || cname => Ok(2),
                    _ if (c2 == "OR" && sname) || (c2 == "UN" && cname) => {
                        if c3.starts_with('G') || c3.starts_with('M')
                            && ["QR", "RQ", "LQ", "QL", "HR", "TR", "BR"].contains(&c4.as_str()) { Ok(2) }
                        else { Ok(1) }
                    },
                    ("GG", _) if c3 == "HD3" => Ok(2),
                    _ => Ok(2),
                },
                3 => match (c2.as_str(), c3.as_str()) {
                    ("GE", "QRF") | ("GE", "RQF") | ("GE", "LQF") | ("GE", "QLF") | ("GE", "HRD") | ("GE", "BRD") => Ok(128),
                    _ if c2 == "GE" && subnam.chars().skip(3).take(4).collect::<String>() == "QP3RK" => Ok(128),
                    ("SY", "TRD") | ("HE", "TRD") if sname || cname => Ok(32),
                    ("GG", _) if c3 == "HD3" => Ok(128),
                    _ if (c2 == "OR" && sname) || (c2 == "UN" && cname) =>
                        if c3.starts_with('G')
                            && ["QR", "RQ", "LQ", "QL", "HR", "TR", "BR"].contains(&c4.as_str()) { Ok(128) }
                        else { Ok(0) },
                    _ => Ok(0),
                },
                _ => { Ok(-1) },
            }
        },
        4 => Ok(6),
        5 => Ok(2),
        6 => Ok((n1.min(n2) as f64 * 1.6).round() as i32),
        7 => Ok(1),
        8 => Ok(50),
        9 => Ok(25),
        10 => Ok(ieeeck(1, 0.0, 1.0)),
        11 => Ok(ieeeck(0, 0.0, 1.0)),
        12..=17 => iparmq(ispec, name, opts, n1, n2, n3, n4),
        _ => Err(LapackError::InvalidIspecValue { value: ispec }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::errors::LapackError;
    use rstest::rstest;

    #[rstest]
    #[case(-1, "CGGHRD", "test", 0, 0, 0, 0, Err(LapackError::InvalidIspecValue { value: -1 }))]
    #[case(0, "CGGHRD", "test", 0, 0, 0, 0, Err(LapackError::InvalidIspecValue { value: 0 }))]
    #[case(1, "CGGHRD", "test", 0, 0, 0, 0, Ok(1))]
    #[case(1, "DGGHRD", "test", 0, 0, 0, 0, Ok(1))]
    #[case(1, "SGBTRF", "test", 0, 128, 0, 128, Ok(32))]
    #[case(1, "DGBTRF", "test", 0, 128, 0, 128, Ok(32))]
    #[case(1, "DGEQRF", "test", 0, 0, 0, 0, Ok(32))]
    #[case(1, "SGEQRF", "test", 0, 0, 0, 0, Ok(32))]
    #[case(2, "CGGHRD", "test", 0, 0, 0, 0, Ok(2))]
    #[case(2, "DGGHRD", "test", 0, 0, 0, 0, Ok(2))]
    #[case(2, "SGBTRF", "test", 0, 0, 0, 0, Ok(2))]
    #[case(2, "DGBTRF", "test", 0, 0, 0, 0, Ok(2))]
    #[case(2, "DGEQRF", "test", 0, 0, 0, 0, Ok(2))]
    #[case(2, "SGEQRF", "test", 0, 0, 0, 0, Ok(2))]
    #[case(3, "CGGHRD", "test", 0, 0, 0, 0, Ok(0))]
    #[case(3, "DGGHRD", "test", 0, 0, 0, 0, Ok(0))]
    #[case(3, "SGBTRF", "test", 0, 0, 0, 0, Ok(0))]
    #[case(3, "DGBTRF", "test", 0, 0, 0, 0, Ok(0))]
    #[case(3, "DGEQRF", "test", 0, 0, 0, 0, Ok(128))]
    #[case(3, "SGEQRF", "test", 0, 0, 0, 0, Ok(128))]
    #[case(4, "CGGHRD", "test", 2, 3, 0, 0, Ok(6))]
    #[case(5, "CGGHRD", "test", 2, 3, 0, 0, Ok(2))]
    #[case(6, "CGGHRD", "test", 2, 3, 0, 0, Ok(3))]
    #[case(7, "CGGHRD", "test", 2, 3, 0, 0, Ok(1))]
    #[case(8, "CGGHRD", "test", 2, 3, 0, 0, Ok(50))]
    #[case(9, "CGGHRD", "test", 2, 3, 0, 0, Ok(25))]
    #[case(10, "CGGHRD", "test", 0, 0, 0, 0, Ok(0))]
    #[case(11, "CGGHRD", "test", 0, 0, 0, 0, Ok(1))]
    #[case(12, "CGGHRD", "test", 2, 3, 4, 5, Ok(75))]
    #[case(13, "CGGHRD", "test", 2, 3, 4, 5, Ok(2))]
    #[case(14, "CGGHRD", "test", 2, 3, 4, 5, Ok(14))]
    #[case(15, "CGGHRD", "test", 2, 3, 4, 5, Ok(2))]
    #[case(16, "CGGHRD", "test", 2, 3, 4, 5, Ok(1))]
    #[case(17, "CGGHRD", "test", 2, 3, 4, 5, Ok(10))]
    #[case(18, "CGGHRD", "test", 0, 0, 0, 0, Err(LapackError::InvalidIspecValue { value: 18 }))]
    fn test_ispec_out_of_range(
        #[case] ispec: i32,
        #[case] name: &str,
        #[case] opts: &str,
        #[case] n1: i32,
        #[case] n2: i32,
        #[case] n3: i32,
        #[case] n4: i32,
        #[case] expected: Result<i32, LapackError>,
    ) {
        assert_eq!(expected, ilaenv(ispec, name, opts, n1, n2, n3, n4));
    }
}