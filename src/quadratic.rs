//! Quadratic Bézier curve.

use crate::linear::Point2;

/// A quadratic Bézier curve defined by three control points.
#[derive(Debug, Clone, Copy)]
pub struct QuadraticBezier {
    pub p0: Point2,
    pub p1: Point2,
    pub p2: Point2,
}

impl QuadraticBezier {
    pub fn new(p0: Point2, p1: Point2, p2: Point2) -> Self {
        Self { p0, p1, p2 }
    }

    /// Evaluate the curve at parameter t ∈ [0, 1].
    /// B(t) = (1-t)²P0 + 2(1-t)tP1 + t²P2
    pub fn evaluate(&self, t: f64) -> Point2 {
        let u = 1.0 - t;
        let uu = u * u;
        let tt = t * t;
        let ut2 = 2.0 * u * t;
        Point2::new(
            uu * self.p0.x + ut2 * self.p1.x + tt * self.p2.x,
            uu * self.p0.y + ut2 * self.p1.y + tt * self.p2.y,
        )
    }

    /// First derivative at parameter t.
    pub fn derivative(&self, t: f64) -> Point2 {
        let u = 1.0 - t;
        Point2::new(
            2.0 * u * (self.p1.x - self.p0.x) + 2.0 * t * (self.p2.x - self.p1.x),
            2.0 * u * (self.p1.y - self.p0.y) + 2.0 * t * (self.p2.y - self.p1.y),
        )
    }

    /// Split the curve at parameter t using de Casteljau.
    pub fn split(&self, t: f64) -> (QuadraticBezier, QuadraticBezier) {
        let q0 = self.p0;
        let q1 = self.p0.lerp(&self.p1, t);
        let q2 = self.p1.lerp(&self.p2, t);
        let r0 = q1.lerp(&q2, t);
        let r1 = q2;
        let r2 = self.p2;
        (
            QuadraticBezier::new(q0, q1, r0),
            QuadraticBezier::new(r0, r1, r2),
        )
    }

    /// Approximate arc length using Gaussian quadrature.
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

    /// Degree-elevate to a cubic Bézier curve.
    pub fn elevate_degree(&self) -> crate::cubic::CubicBezier {
        let p0 = self.p0;
        let p1 = Point2::new(
            self.p0.x + (2.0 / 3.0) * (self.p1.x - self.p0.x),
            self.p0.y + (2.0 / 3.0) * (self.p1.y - self.p0.y),
        );
        let p2 = Point2::new(
            self.p2.x + (2.0 / 3.0) * (self.p1.x - self.p2.x),
            self.p2.y + (2.0 / 3.0) * (self.p1.y - self.p2.y),
        );
        let p3 = self.p2;
        crate::cubic::CubicBezier::new(p0, p1, p2, p3)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_evaluate_endpoints() {
        let c = QuadraticBezier::new(
            Point2::new(0.0, 0.0),
            Point2::new(5.0, 10.0),
            Point2::new(10.0, 0.0),
        );
        assert_eq!(c.evaluate(0.0), Point2::new(0.0, 0.0));
        assert_eq!(c.evaluate(1.0), Point2::new(10.0, 0.0));
    }

    #[test]
    fn test_evaluate_midpoint() {
        let c = QuadraticBezier::new(
            Point2::new(0.0, 0.0),
            Point2::new(5.0, 10.0),
            Point2::new(10.0, 0.0),
        );
        let mid = c.evaluate(0.5);
        assert!((mid.x - 5.0).abs() < 1e-10);
        assert!((mid.y - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_derivative_at_endpoints() {
        let c = QuadraticBezier::new(
            Point2::new(0.0, 0.0),
            Point2::new(5.0, 10.0),
            Point2::new(10.0, 0.0),
        );
        let d0 = c.derivative(0.0);
        assert!((d0.x - 10.0).abs() < 1e-10);
        assert!((d0.y - 20.0).abs() < 1e-10);
    }

    #[test]
    fn test_split_continuity() {
        let c = QuadraticBezier::new(
            Point2::new(0.0, 0.0),
            Point2::new(5.0, 10.0),
            Point2::new(10.0, 0.0),
        );
        let (left, right) = c.split(0.5);
        // End of left should equal start of right
        assert!((left.evaluate(1.0).x - right.evaluate(0.0).x).abs() < 1e-10);
        assert!((left.evaluate(1.0).y - right.evaluate(0.0).y).abs() < 1e-10);
    }

    #[test]
    fn test_approximate_length_positive() {
        let c = QuadraticBezier::new(
            Point2::new(0.0, 0.0),
            Point2::new(0.0, 10.0),
            Point2::new(10.0, 10.0),
        );
        let len = c.approximate_length(100);
        assert!(len > 0.0);
        // Should be at least as long as the straight line from start to end
        let straight = c.p0.distance_to(&c.p2);
        assert!(len >= straight - 0.1);
    }
}
