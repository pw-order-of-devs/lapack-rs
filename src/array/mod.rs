pub mod convert;
pub mod display;
pub mod indexing;

#[derive(Clone, Debug)]
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

    pub fn single(data: f64) -> Self {
        FortranArray {
            data: vec![data],
            rows: 1,
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

    pub(crate) fn insert(&mut self, index: i32, value: f64) {
        if index >= 1 && index < self.len() {
            self.data.insert(index as usize - 1, value)
        }
    }

    pub(crate) fn idx(&self, index: (i32, i32)) -> i32 {
        (index.1 - 1) * self.rows + index.0
    }

    pub(crate) fn len(&self) -> i32 {
        if self.is_2d() { self.rows * self.cols }
        else { self.rows }
    }

    pub(crate) fn is_1d(&self) -> bool {
        self.cols == 0
    }

    pub(crate) fn is_2d(&self) -> bool {
        self.rows > 0 && self.cols > 0
    }

    pub(crate) fn as_2d(&mut self, rows: i32) {
        let len = self.data.len() as i32;
        self.rows = rows;
        self.cols = len / rows;
    }

    pub(crate) fn to_vec(self) -> Vec<f64> {
        self.data
    }
}

impl PartialEq for FortranArray {
    fn eq(&self, other: &Self) -> bool {
        self.rows == other.rows
            && self.cols == other.cols
            && self.data == other.data
    }
}
