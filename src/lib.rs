//! Bézier curves and surfaces library.
//!
//! Provides linear, quadratic, and cubic Bézier curves, de Casteljau's algorithm,
//! curve splitting, and arc length parameterization.

pub mod arclength;
pub mod casteljau;
pub mod cubic;
pub mod linear;
pub mod quadratic;

pub use arclength::arc_length;
pub use casteljau::de_casteljau;
pub use cubic::CubicBezier;
pub use linear::LinearBezier;
pub use quadratic::QuadraticBezier;
