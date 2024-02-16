/// LSAME
///
/// # Documentation
///
/// [Original] Online html documentation available at
/// [http://www.netlib.org/lapack/explore-html/](http://www.netlib.org/lapack/explore-html/)
///
/// Returns true if `ca` is the same letter as `cb` regardless of case.
pub fn lsame(
    ca: char,
    cb: char,
) -> bool {
    let ca_u8 = ca as u8;
    let cb_u8 = cb as u8;

    ca_u8 == cb_u8 || (
        ca_u8.is_ascii_alphabetic()
            && cb_u8.is_ascii_alphabetic()
            && ca_u8.eq_ignore_ascii_case(&cb_u8))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case('A', 'A', true)]
    #[case('A', 'a', true)]
    #[case('a', 'A', true)]
    #[case('a', 'a', true)]
    #[case('a', 'b', false)]
    #[case('A', 'B', false)]
    fn test_lsame(#[case] a: char, #[case] b: char, #[case] expected: bool) {
        assert_eq!(lsame(a, b), expected);
    }
}