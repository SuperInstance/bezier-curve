//! Cubic Bézier curve.

use crate::linear::Point2;

/// A cubic Bézier curve defined by four control points.
#[derive(Debug, Clone, Copy)]
pub struct CubicBezier {
    pub p0: Point2,
    pub p1: Point2,
    pub p2: Point2,
    pub p3: Point2,
}

impl CubicBezier {
    pub fn new(p0: Point2, p1: Point2, p2: Point2, p3: Point2) -> Self {
        Self { p0, p1, p2, p3 }
    }

    /// Evaluate the curve at parameter t ∈ [0, 1].
    pub fn evaluate(&self, t: f64) -> Point2 {
        let u = 1.0 - t;
        let uu = u * u;
        let uuu = uu * u;
        let tt = t * t;
        let ttt = tt * t;

        Point2::new(
            uuu * self.p0.x + 3.0 * uu * t * self.p1.x + 3.0 * u * tt * self.p2.x + ttt * self.p3.x,
            uuu * self.p0.y + 3.0 * uu * t * self.p1.y + 3.0 * u * tt * self.p2.y + ttt * self.p3.y,
        )
    }

    /// First derivative at parameter t.
    pub fn derivative(&self, t: f64) -> Point2 {
        let u = 1.0 - t;
        Point2::new(
            3.0 * u * u * (self.p1.x - self.p0.x)
                + 6.0 * u * t * (self.p2.x - self.p1.x)
                + 3.0 * t * t * (self.p3.x - self.p2.x),
            3.0 * u * u * (self.p1.y - self.p0.y)
                + 6.0 * u * t * (self.p2.y - self.p1.y)
                + 3.0 * t * t * (self.p3.y - self.p2.y),
        )
    }

    /// Second derivative at parameter t.
    pub fn second_derivative(&self, t: f64) -> Point2 {
        let u = 1.0 - t;
        Point2::new(
            6.0 * u * (self.p2.x - 2.0 * self.p1.x + self.p0.x)
                + 6.0 * t * (self.p3.x - 2.0 * self.p2.x + self.p1.x),
            6.0 * u * (self.p2.y - 2.0 * self.p1.y + self.p0.y)
                + 6.0 * t * (self.p3.y - 2.0 * self.p2.y + self.p1.y),
        )
    }

    /// Split the curve at parameter t using de Casteljau.
    pub fn split(&self, t: f64) -> (CubicBezier, CubicBezier) {
        let a = self.p0.lerp(&self.p1, t);
        let b = self.p1.lerp(&self.p2, t);
        let c = self.p2.lerp(&self.p3, t);
        let d = a.lerp(&b, t);
        let e = b.lerp(&c, t);
        let f = d.lerp(&e, t);
        (
            CubicBezier::new(self.p0, a, d, f),
            CubicBezier::new(f, e, c, self.p3),
        )
    }

    /// Approximate arc length using numerical integration.
    pub fn approximate_length(&self, steps: usize) -> f64 {
        let steps = steps.max(1);
        let mut length = 0.0;
        let dt = 1.0 / steps as f64;
        let mut prev = self.evaluate(0.0);
        for i in 1..=steps {
            let t = i as f64 * dt;
            let curr = self.evaluate(t);
            length += prev.distance_to(&curr);
            prev = curr;
        }
        length
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_evaluate_endpoints() {
        let c = CubicBezier::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 3.0),
            Point2::new(3.0, 3.0),
            Point2::new(4.0, 0.0),
        );
        assert_eq!(c.evaluate(0.0), Point2::new(0.0, 0.0));
        assert_eq!(c.evaluate(1.0), Point2::new(4.0, 0.0));
    }

    #[test]
    fn test_evaluate_midpoint_symmetric() {
        let c = CubicBezier::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 3.0),
            Point2::new(3.0, 3.0),
            Point2::new(4.0, 0.0),
        );
        let mid = c.evaluate(0.5);
        assert!((mid.x - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_derivative_at_zero() {
        let c = CubicBezier::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 2.0),
            Point2::new(3.0, 2.0),
            Point2::new(4.0, 0.0),
        );
        let d = c.derivative(0.0);
        assert!((d.x - 3.0).abs() < 1e-10);
        assert!((d.y - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_split_continuity() {
        let c = CubicBezier::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 3.0),
            Point2::new(3.0, 3.0),
            Point2::new(4.0, 0.0),
        );
        let (left, right) = c.split(0.4);
        let end_left = left.evaluate(1.0);
        let start_right = right.evaluate(0.0);
        assert!((end_left.x - start_right.x).abs() < 1e-10);
        assert!((end_left.y - start_right.y).abs() < 1e-10);
    }

    #[test]
    fn test_split_matches_original() {
        let c = CubicBezier::new(
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 8.0),
            Point2::new(6.0, 8.0),
            Point2::new(8.0, 0.0),
        );
        let (left, right) = c.split(0.5);
        // Sample and compare
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let orig = c.evaluate(t);
            let split_pt = if t <= 0.5 {
                left.evaluate(t * 2.0)
            } else {
                right.evaluate((t - 0.5) * 2.0)
            };
            assert!((orig.x - split_pt.x).abs() < 1e-8, "x mismatch at t={}", t);
            assert!((orig.y - split_pt.y).abs() < 1e-8, "y mismatch at t={}", t);
        }
    }

    #[test]
    fn test_length_positive() {
        let c = CubicBezier::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 10.0),
            Point2::new(9.0, 10.0),
            Point2::new(10.0, 0.0),
        );
        let len = c.approximate_length(200);
        assert!(len > 0.0);
    }

    #[test]
    fn test_second_derivative() {
        // Linear cubic (all control points collinear)
        let c = CubicBezier::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, 2.0),
            Point2::new(3.0, 3.0),
        );
        let d2 = c.second_derivative(0.5);
        assert!((d2.x).abs() < 1e-10);
        assert!((d2.y).abs() < 1e-10);
    }
}
