use std::ops::{Index, IndexMut, RangeFrom};

use crate::array::FortranArray;

impl Index<i32> for FortranArray {
    type Output = f64;

    fn index(&self, index: i32) -> &Self::Output {
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

impl IndexMut<i32> for FortranArray {

    fn index_mut(&mut self, index: i32) -> &mut Self::Output {
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

impl Index<(i32, i32)> for FortranArray {
    type Output = f64;

    fn index(&self, index: (i32, i32)) -> &Self::Output {
        let index = ((index.1 - 1) * self.rows + (index.0 - 1)) as usize;

        if let Some(value) = self.data.get(index) {
            value
        } else {
            &f64::MIN_POSITIVE
        }
    }
}

impl IndexMut<(i32, i32)> for FortranArray {

    fn index_mut(&mut self, index: (i32, i32)) -> &mut Self::Output {
        let (col, row) = index;
        let index = ((row - 1) * self.rows + (col - 1)) as usize;

        if let Some(value) = self.data.get_mut(index) {
            value
        } else {
            &mut self.default_value
        }
    }
}

impl Index<RangeFrom<i32>> for FortranArray {
    type Output = [f64];

    fn index(&self, range: RangeFrom<i32>) -> &Self::Output {
        if range.start >= 1 && range.start <= self.len() {
            &self.data[(range.start - 1) as usize..]
        } else {
            &[]
        }
    }
}

impl IndexMut<RangeFrom<i32>> for FortranArray {

    fn index_mut(&mut self, range: RangeFrom<i32>) -> &mut Self::Output {
        if range.start >= 1 && range.start <= self.len()  {
            &mut self.data[(range.start - 1) as usize..]
        } else {
            &mut []
        }
    }
}

impl Index<RangeFrom<(i32, i32)>> for FortranArray {
    type Output = [f64];

    fn index(&self, index: RangeFrom<(i32, i32)>) -> &Self::Output {
        let (col, row) = index.start;
        let index = (row - 1) * self.rows + col;
        &self[index..]
    }
}

impl IndexMut<RangeFrom<(i32, i32)>> for FortranArray {

    fn index_mut(&mut self, index: RangeFrom<(i32, i32)>) -> &mut Self::Output {
        let (col, row) = index.start;
        let index = (row - 1) * self.rows + col;
        &mut self[index..]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[rstest]
    #[case(10, 100.0, 100.0)]
    #[case(43, 200.0, 200.0)]
    #[case(66, 0.0, 0.0)]
    #[case(-1, f64::MIN_POSITIVE, f64::MIN_POSITIVE)]
    fn test_set_get(
        #[case] idx: i32,
        #[case] set_value: f64,
        #[case] expected: f64,
    ) {
        let mut array = FortranArray::vector(&vec![0.; 100]);
        array[idx] = set_value;
        assert_eq!(expected, array[idx]);
    }

    #[rstest]
    #[case(5, 5, 100.0, 100.0)]
    #[case(9, 9, 200.0, 200.0)]
    #[case(1, 1, 0.0, 0.0)]
    #[case(-1, -1, f64::MIN_POSITIVE, f64::MIN_POSITIVE)]
    fn test_set_get_double(
        #[case] x: i32,
        #[case] y: i32,
        #[case] set_value: f64,
        #[case] expected: f64,
    ) {
        let mut array = FortranArray::zeros(10, 10);
        array[(x, y)] = set_value;
        assert_eq!(expected, array[(x, y)]);
    }
}
