//! Optional parts of the app, each behind its own Cargo feature. The viewer
//! stands on its own; a plugin adds jobs and controls around it.

#[cfg(feature = "simplefold")]
pub mod simplefold;
