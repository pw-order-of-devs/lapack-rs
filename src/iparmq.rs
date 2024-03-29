/// IPARMQ
///
/// [Original] Online HTML documentation available at
/// `http://www.netlib.org/lapack/explore-html/`
///
/// # Arguments
///
/// For arguments definitions, please refer to the original documentation.
pub(crate) fn iparmq(
    ispec: i32,
    name: &str,
    ilo: i32,
    ihi: i32,
) -> i32 {
    let inmin = 12;
    let inwin = 13;
    let inibl = 14;
    let ishfts = 15;
    let iacc22 = 16;
    let icost = 17;
    let nmin = 75;
    let k22min = 14;
    let kacmin = 14;
    let nibble = 14;
    let knwswp = 500;
    let rcost = 10;

    let nh = ihi - ilo + 1;
    let mut ns = 2;
    if ispec == ishfts || ispec == inwin || ispec == iacc22 {
        match nh {
            nh if nh < 30 => 2,
            nh if nh < 60 => 4,
            nh if nh < 150 => 10,
            nh if nh < 590 => 10.max(nh / (nh as f64).log(2.).round() as i32),
            nh if nh < 3000 => 64,
            nh if nh < 6000 => 128,
            _ => 256,
        };
        ns = 2.max(ns - ns % 2);
    }

    match ispec {
        _ if ispec == inmin => nmin,
        _ if ispec == inibl => nibble,
        _ if ispec == ishfts => ns,
        _ if ispec == inwin => {
            if nh < knwswp { ns }
            else { 3 * ns / 2 }
        },
        _ if ispec == iacc22 => {
            let mut result = 0;
            let mut subnam = name.to_string();
            let ic = subnam.chars().next().unwrap() as u32;
            let iz = 'Z' as u32;
            if iz == 90 || iz == 122 {
                // ASCII character set
                if ic >= 97 && ic <= 122 {
                    let mut chars: Vec<char> = subnam.chars().collect();
                    for i in 0..chars.len().min(6) {
                        let mut ic = chars[i] as u32;
                        if ic >= 97 && ic <= 122 {
                            ic -= 32;
                            chars[i] = std::char::from_u32(ic).unwrap_or_default();
                        }
                    }
                    subnam = chars.into_iter().collect::<String>();
                }
            } else if iz == 233 || iz == 169 {
                // EBCDIC character set
                let mut chars: Vec<char> = subnam.chars().collect();
                for i in 0..chars.len().min(6) {
                    let mut ic = chars[i] as u32;
                    if (ic >= 129 && ic <= 137) ||
                        (ic >= 145 && ic <= 153) ||
                        (ic >= 162 && ic <= 169) {
                        ic += 64;
                        chars[i] = std::char::from_u32(ic).unwrap_or_default();
                    }
                }
                subnam = chars.into_iter().collect::<String>();
            } else if iz == 218 || iz == 250 {
                // Prime machines: ASCII+128
                let mut chars: Vec<char> = subnam.chars().collect();
                for i in 0..chars.len().min(6) {
                    let mut ic = chars[i] as u32;
                    if ic >= 225 && ic <= 250 {
                        ic -= 32;
                        chars[i] = std::char::from_u32(ic).unwrap_or_default();
                    }
                }
                subnam = chars.into_iter().collect::<String>();
            }

            if &subnam[1..6] == "GGHRD" || &subnam[1..6] == "GGHD3" {
                result = 1;
                if nh >= k22min { return 2 }
            } else if &subnam[3..6] == "EXC" {
                if nh >= kacmin { return 1 }
                else if nh >= k22min { return 2 }
            } else if &subnam[1..6] == "HSEQR" || &subnam[1..5] == "LAQR" {
                if ns >= kacmin { return 1 }
                else if ns >= k22min { return 2 }
            }
            return result
        },
        _ if ispec == icost => rcost,
        _ => -1,
    }
}

#[cfg(test)]
mod tests {
    use super::iparmq;
    use rstest::rstest;

    #[rstest]
    #[case(12, "test", 5, 5, 75)]
    #[case(14, "test", 10, 10, 14)]
    #[case(18, "test", 30, 30, -1)]
    #[case(15, "CGGHRD", 15, 15, 2)]
    #[case(16, "DGGHRD", 20, 20, 1)]
    #[case(17, "CHSEQR", 25, 25, 10)]
    #[case(15, "CTGEXC", 35, 35, 2)]
    #[case(16, "CLAQR0", 40, 40, 0)]
    #[case(17, "CHSEQR", 45, 45, 10)]
    fn test_iparmq(
        #[case] ispec: i32,
        #[case] name: &str,
        #[case] ilo: i32,
        #[case] ihi: i32,
        #[case] expected: i32,
    ) {
        assert_eq!(expected, iparmq(ispec, name, ilo, ihi));
    }
}
