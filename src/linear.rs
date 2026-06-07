//! Linear Bézier curve (simple lerp between two points).

/// A 2D point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point2 {
    pub x: f64,
    pub y: f64,
}

impl Point2 {
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    pub fn lerp(&self, other: &Point2, t: f64) -> Point2 {
        Point2::new(
            self.x + t * (other.x - self.x),
            self.y + t * (other.y - self.y),
        )
    }

    pub fn distance_to(&self, other: &Point2) -> f64 {
        ((self.x - other.x).powi(2) + (self.y - other.y).powi(2)).sqrt()
    }
}

impl std::ops::Add for Point2 {
    type Output = Point2;
    fn add(self, rhs: Point2) -> Point2 {
        Point2::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl std::ops::Sub for Point2 {
    type Output = Point2;
    fn sub(self, rhs: Point2) -> Point2 {
        Point2::new(self.x - rhs.x, self.y - rhs.y)
    }
}

/// A linear Bézier curve (straight line between P0 and P1).
#[derive(Debug, Clone, Copy)]
pub struct LinearBezier {
    pub p0: Point2,
    pub p1: Point2,
}

impl LinearBezier {
    pub fn new(p0: Point2, p1: Point2) -> Self {
        Self { p0, p1 }
    }

    /// Evaluate the curve at parameter t ∈ [0, 1].
    pub fn evaluate(&self, t: f64) -> Point2 {
        self.p0.lerp(&self.p1, t)
    }

    /// Derivative at parameter t.
    pub fn derivative(&self, _t: f64) -> Point2 {
        Point2::new(self.p1.x - self.p0.x, self.p1.y - self.p0.y)
    }

    /// Split the curve at parameter t, returning two sub-curves.
    pub fn split(&self, t: f64) -> (LinearBezier, LinearBezier) {
        let mid = self.evaluate(t);
        (
            LinearBezier::new(self.p0, mid),
            LinearBezier::new(mid, self.p1),
        )
    }

    /// Arc length (exact for linear = distance between endpoints).
    pub fn length(&self) -> f64 {
        self.p0.distance_to(&self.p1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_evaluate_endpoints() {
        let c = LinearBezier::new(Point2::new(0.0, 0.0), Point2::new(10.0, 0.0));
        assert_eq!(c.evaluate(0.0), Point2::new(0.0, 0.0));
        assert_eq!(c.evaluate(1.0), Point2::new(10.0, 0.0));
    }

    #[test]
    fn test_evaluate_midpoint() {
        let c = LinearBezier::new(Point2::new(0.0, 0.0), Point2::new(10.0, 10.0));
        let mid = c.evaluate(0.5);
        assert!((mid.x - 5.0).abs() < 1e-10);
        assert!((mid.y - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_derivative_constant() {
        let c = LinearBezier::new(Point2::new(0.0, 0.0), Point2::new(5.0, 3.0));
        let d = c.derivative(0.5);
        assert!((d.x - 5.0).abs() < 1e-10);
        assert!((d.y - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_length() {
        let c = LinearBezier::new(Point2::new(0.0, 0.0), Point2::new(3.0, 4.0));
        assert!((c.length() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_split() {
        let c = LinearBezier::new(Point2::new(0.0, 0.0), Point2::new(10.0, 0.0));
        let (left, right) = c.split(0.3);
        assert_eq!(left.p0, Point2::new(0.0, 0.0));
        assert!((left.p1.x - 3.0).abs() < 1e-10);
        assert!((right.p0.x - 3.0).abs() < 1e-10);
        assert_eq!(right.p1, Point2::new(10.0, 0.0));
    }
}
