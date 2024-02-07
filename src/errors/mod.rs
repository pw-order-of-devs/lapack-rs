#[derive(Debug, Clone, PartialOrd, PartialEq, Eq)]
pub enum LapackError {
    InvalidIspecValue { value: i32 },
    CannotBeEmpty { param: &'static str },
}

impl std::error::Error for LapackError {}

impl std::fmt::Display for LapackError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidIspecValue { value } => write!(f, "Invalid ISPEC value: {value}"),
            Self::CannotBeEmpty { param } => write!(f, "Parameter cannot be empty: {param}"),
        }
    }
}
