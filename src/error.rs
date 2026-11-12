//! A minimal error type: every failure is reported to the user as one message.

use std::fmt;

/// An error with a message that is ready to be shown to the user.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Error(pub String);

impl Error {
    pub fn new(msg: impl Into<String>) -> Self {
        Error(msg.into())
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for Error {}

pub type Result<T> = std::result::Result<T, Error>;

/// Adds a description of what was being done to a lower-level error.
pub trait Context<T> {
    fn context<F: FnOnce() -> String>(self, f: F) -> Result<T>;
}

impl<T, E: fmt::Display> Context<T> for std::result::Result<T, E> {
    fn context<F: FnOnce() -> String>(self, f: F) -> Result<T> {
        self.map_err(|e| Error(format!("{}: {}", f(), e)))
    }
}
