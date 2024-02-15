use std::fmt;
use crate::array::FortranArray;

impl fmt::Display for FortranArray {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let p = f.precision().unwrap_or(1);
        if self.is_1d() {
            // 1D Vector
            write!(f, "[")?;
            for (i, value) in self.data.iter().enumerate() {
                if i != 0 { write!(f, ", ")?; }
                write!(f, "{:.p$}", value)?;
            }
            write!(f, "]")
        } else {
            // 2D Matrix
            write!(f, "[")?;
            for row in 0..self.rows {
                write!(f, "[")?;
                for col in 0..self.cols {
                    if col != 0 { write!(f, ", ")?; }
                    let index = (row * self.cols + col) as usize;
                    write!(f, "{:.p$}", self.data[index])?;
                }
                write!(f, "]")?;
            }
            write!(f, "]")
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(FortranArray::vector(& [1.0, 2.0, 3.0, 4.0]), "[1.0, 2.0, 3.0, 4.0]")]
    #[case(FortranArray::matrix(& [1.0, 2.0, 3.0, 4.0], 1, 4), "[[1.0, 2.0, 3.0, 4.0]]")]
    #[case(FortranArray::matrix(& [1.0, 2.0, 3.0, 4.0], 4, 1), "[[1.0][2.0][3.0][4.0]]")]
    #[case(FortranArray::matrix(& [1.0, 2.0, 3.0, 4.0], 2, 2), "[[1.0, 3.0][2.0, 4.0]]")]
    fn test_display(
        #[case] input: FortranArray,
        #[case] expected: &str,
    ) {
        assert_eq!(expected, format!("{input}"))
    }
}
