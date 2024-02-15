use crate::array::FortranArray;

impl From<&[f64]> for FortranArray {
    fn from(data: &[f64]) -> Self {
        let len = data.len();
        FortranArray {
            data: data.to_vec(),
            rows: len as i32,
            cols: 0,
            default_value: f64::NAN,
        }
    }
}

impl From<Vec<f64>> for FortranArray {
    fn from(data: Vec<f64>) -> Self {
        let len = data.len();
        FortranArray {
            data,
            rows: len as i32,
            cols: 0,
            default_value: f64::NAN,
        }
    }
}

impl From<Vec<Vec<f64>>> for FortranArray {
    fn from(data: Vec<Vec<f64>>) -> Self {
        let rows = data.len();
        let cols = if rows > 0 { data[0].len() } else { 0 };
        let data = data.into_iter().flatten().collect();
        FortranArray {
            data,
            rows: rows as i32,
            cols: cols as i32,
            default_value: f64::NAN,
        }
    }
}

impl From<FortranArray> for Vec<f64> {
    fn from(array: FortranArray) -> Self {
        array.data
    }
}

impl From<FortranArray> for Vec<Vec<f64>> {
    fn from(array: FortranArray) -> Self {
        let mut data = vec![vec![0.0; array.rows as usize]; array.cols as usize];
        for (i, value) in array.data.iter().enumerate() {
            let row = i % array.rows as usize;
            let col = i / array.rows as usize;
            data[col][row] = *value;
        }
        data
    }
}

pub trait ToFortranArray {
    fn to_fortran_array(&self) -> FortranArray;
}

impl ToFortranArray for &[f64] {
    fn to_fortran_array(&self) -> FortranArray {
        FortranArray::from(self.to_vec())
    }
}

impl ToFortranArray for Vec<f64> {
    fn to_fortran_array(&self) -> FortranArray {
        FortranArray::from(self.clone())
    }
}

impl ToFortranArray for Vec<Vec<f64>> {
    fn to_fortran_array(&self) -> FortranArray {
        FortranArray::from(self.clone())
    }
}

impl ToFortranArray for FortranArray {
    fn to_fortran_array(&self) -> FortranArray {
        self.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(vec![1.0, 2.0, 3.0, 4.0], FortranArray::vector(&[1.0, 2.0, 3.0, 4.0]))]
    fn test_from_vec(
        #[case] input: Vec<f64>,
        #[case] expected: FortranArray,
    ) {
        let input: FortranArray = input.into();
        assert_eq!(expected, input.into());
    }

    #[rstest]
    #[case(FortranArray::vector(&[1.0, 2.0, 3.0, 4.0]), vec![1.0, 2.0, 3.0, 4.0])]
    fn test_into_vec(
        #[case] input: FortranArray,
        #[case] expected: Vec<f64>,
    ) {
        let input: Vec<f64> = input.into();
        assert_eq!(expected, input);
    }

    #[rstest]
    #[case(vec![vec![1.0, 2.0], vec![3.0, 4.0]], FortranArray::matrix(&[1.0, 2.0, 3.0, 4.0], 2, 2))]
    fn test_from_vec_2d(
        #[case] input: Vec<Vec<f64>>,
        #[case] expected: FortranArray,
    ) {
        let input: FortranArray = input.into();
        assert_eq!(expected, input.into());
    }

    #[rstest]
    #[case(FortranArray::matrix(&[1.0, 2.0, 3.0, 4.0], 2, 2), vec![vec![1.0, 2.0], vec![3.0, 4.0]])]
    fn test_into_vec_2d(
        #[case] input: FortranArray,
        #[case] expected: Vec<Vec<f64>>,
    ) {
        let input: Vec<Vec<f64>> = input.into();
        assert_eq!(expected, input);
    }
}
