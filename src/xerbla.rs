use std::process;

/// XERBLA
///
/// # Documentation
///
/// [Original] Online html documentation available at
/// [http://www.netlib.org/lapack/explore-html/](http://www.netlib.org/lapack/explore-html/)
///
/// # Purpose
///
/// XERBLA is an error handler for the LAPACK routines.
/// It is called by an LAPACK routine if an input parameter has an
/// invalid value. A message is printed and execution stops.
///
/// Installers may consider modifying the STOP statement in order to
/// call system-specific exception-handling facilities.
///
/// # Arguments
///
/// For arguments definitions, please refer to the original documentation.
pub fn xerbla(
    srname: &str,
    info: i32,
) {
    println!("** On entry to {} parameter number {} had an illegal value", srname, info);
    process::exit(1);
}
