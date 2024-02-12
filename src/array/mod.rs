#[derive(Clone)]
pub struct FortranArray {
    data: Vec<f64>,
    rows: i32,
    cols: i32,
    default_value: f64,
}

impl Default for FortranArray {
    fn default() -> Self {
        FortranArray {
            data: vec![],
            rows: 0,
            cols: 0,
            default_value: f64::NAN,
        }
    }
}

impl FortranArray {
    pub fn empty() -> Self {
        FortranArray::default()
    }

    pub fn zeros(rows: i32, cols: i32) -> Self {
        FortranArray {
            data: vec![0.0; (rows * cols) as usize],
            rows,
            cols,
            ..Default::default()
        }
    }

    pub fn vector(data: &[f64]) -> Self {
        FortranArray {
            data: data.to_vec(),
            rows: data.len() as i32,
            ..Default::default()
        }
    }

    pub fn matrix(data: &[f64], rows: i32, cols: i32) -> Self {
        assert_eq!(data.len() as i32, rows * cols);
        FortranArray {
            data: data.to_vec(),
            rows,
            cols,
            ..Default::default()
        }
    }

    pub fn len(&self) -> i32 {
        self.rows
    }

    fn is_1d(&self) -> bool {
        self.cols == 0
    }

    fn is_2d(&self) -> bool {
        self.rows > 0 && self.cols > 0
    }
}

impl std::ops::Index<i32> for FortranArray {
    type Output = f64;

    fn index(&self, index: i32) -> &Self::Output {
        assert!(self.is_1d());
        if index == 0 {
            return &f64::MIN_POSITIVE;
        }

        if let Some(value) = self.data.get(index as usize - 1) {
            &value
        } else {
            &f64::MIN_POSITIVE
        }
    }
}

impl std::ops::IndexMut<i32> for FortranArray {
    fn index_mut(&mut self, index: i32) -> &mut Self::Output {
        assert!(self.is_1d());
        if index == 0 {
            return &mut self.default_value;
        }

        if let Some(value) = self.data.get_mut(index as usize - 1) {
            value
        } else {
            &mut self.default_value
        }
    }
}

impl std::ops::Index<(i32, i32)> for FortranArray {
    type Output = f64;

    fn index(&self, index: (i32, i32)) -> &Self::Output {
        assert!(self.is_2d());
        let index = ((index.1 - 1) * self.cols + (index.0 - 1)) as usize;

        if let Some(value) = self.data.get(index) {
            value
        } else {
            &f64::MIN_POSITIVE
        }
    }
}

impl std::ops::IndexMut<(i32, i32)> for FortranArray {
    fn index_mut(&mut self, index: (i32, i32)) -> &mut Self::Output {
        assert!(self.is_2d());
        let (col, row) = index;
        let index = ((row - 1) * self.cols + (col - 1)) as usize;

        if let Some(value) = self.data.get_mut(index) {
            value
        } else {
            &mut self.default_value
        }
    }
}

impl From<f64> for FortranArray {
    fn from(data: f64) -> Self {
        FortranArray {
            data: vec![data],
            rows: 1,
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

impl From<FortranArray> for Vec<f64> {
    fn from(array: FortranArray) -> Self {
        array.data
    }
}

impl From<FortranArray> for f64 {
    fn from(array: FortranArray) -> Self {
        array.data[0]
    }
}

impl PartialEq for FortranArray {
    fn eq(&self, other: &Self) -> bool {
        self.rows == other.rows
            && self.cols == other.cols
            && self.data == other.data
    }
}

impl std::fmt::Debug for FortranArray {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.is_1d() {
            // 1D Vector
            write!(f, "[")?;
            for (i, value) in self.data.iter().enumerate() {
                if i != 0 { write!(f, ", ")?; }
                write!(f, "{:.1}", value)?;
            }
            write!(f, "]\n")
        } else {
            // 2D Matrix
            for row in 0..self.rows {
                write!(f, "[")?;
                for col in 0..self.cols {
                    if col != 0 { write!(f, ", ")?; }
                    let index = (col * self.rows + row) as usize;
                    write!(f, "{:.1}", self.data[index])?;
                }
                writeln!(f, "]")?;
            }
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(5, 5, 100.0, 100.0)]
    #[case(9, 9, 200.0, 200.0)]
    #[case(1, 1, 0.0, 0.0)]
    #[case(-1, -1, f64::MIN_POSITIVE, f64::MIN_POSITIVE)]
    fn test_set_get(
        #[case] x: i32,
        #[case] y: i32,
        #[case] set_value: f64,
        #[case] expected: f64,
    ) {
        let mut array = FortranArray::zeros(10, 10);
        array[(x, y)] = set_value;
        assert_eq!(expected, array[(x, y)]);
    }

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

    #[rstest]
    #[case(FortranArray::vector(&[1.0, 2.0, 3.0, 4.0]), "[1.0, 2.0, 3.0, 4.0]")]
    #[case(FortranArray::matrix(&[1.0, 2.0, 3.0, 4.0], 1, 4), "[1.0, 2.0, 3.0, 4.0]\n")]
    #[case(FortranArray::matrix(&[1.0, 2.0, 3.0, 4.0], 4, 1), "[1.0]\n[2.0]\n[3.0]\n[4.0]\n")]
    #[case(FortranArray::matrix(&[1.0, 2.0, 3.0, 4.0], 2, 2), "[1.0, 3.0]\n[2.0, 4.0]\n")]
    fn test_debug(
        #[case] input: FortranArray,
        #[case] expected: &str,
    ) {
        assert_eq!(expected, format!("{input:?}"))
    }
}

pub trait AsFortranArray {
    fn to_fortran_array(&self) -> FortranArray;
}

impl AsFortranArray for f64 {
    fn to_fortran_array(&self) -> FortranArray {
        FortranArray::from(vec![self.clone()])
    }
}

impl AsFortranArray for Vec<f64> {
    fn to_fortran_array(&self) -> FortranArray {
        FortranArray::from(self.clone())
    }
}

impl AsFortranArray for Vec<Vec<f64>> {
    fn to_fortran_array(&self) -> FortranArray {
        FortranArray::from(self.clone())
    }
}

impl AsFortranArray for FortranArray {
    fn to_fortran_array(&self) -> FortranArray {
        self.clone()
    }
}
