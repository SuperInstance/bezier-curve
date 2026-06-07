//! De Casteljau's algorithm for evaluating and splitting Bézier curves.

use crate::linear::Point2;

/// Evaluate a Bézier curve of arbitrary degree using de Casteljau's algorithm.
/// `control_points` defines the curve, `t` is the parameter in [0, 1].
pub fn de_casteljau(control_points: &[Point2], t: f64) -> Point2 {
    if control_points.is_empty() {
        return Point2::new(0.0, 0.0);
    }
    if control_points.len() == 1 {
        return control_points[0];
    }

    let mut pts = control_points.to_vec();
    while pts.len() > 1 {
        let mut next = Vec::with_capacity(pts.len() - 1);
        for i in 0..pts.len() - 1 {
            next.push(pts[i].lerp(&pts[i + 1], t));
        }
        pts = next;
    }
    pts[0]
}

/// Get all intermediate points from de Casteljau's algorithm.
/// Returns a vector of layers, where each layer is the points at that level.
pub fn de_casteljau_layers(control_points: &[Point2], t: f64) -> Vec<Vec<Point2>> {
    if control_points.is_empty() {
        return vec![];
    }
    let mut layers = vec![control_points.to_vec()];
    let mut pts = control_points.to_vec();
    while pts.len() > 1 {
        let mut next = Vec::with_capacity(pts.len() - 1);
        for i in 0..pts.len() - 1 {
            next.push(pts[i].lerp(&pts[i + 1], t));
        }
        layers.push(next.clone());
        pts = next;
    }
    layers
}

/// Split a Bézier curve at parameter t using de Casteljau.
/// Returns two sets of control points for the left and right sub-curves.
pub fn de_casteljau_split(control_points: &[Point2], t: f64) -> (Vec<Point2>, Vec<Point2>) {
    if control_points.is_empty() {
        return (vec![], vec![]);
    }

    let layers = de_casteljau_layers(control_points, t);
    let n = layers.len();

    // Left curve: first point of each layer
    let left: Vec<Point2> = (0..n).map(|i| layers[i][0]).collect();

    // Right curve: last point of each layer (reversed)
    let right: Vec<Point2> = (0..n).rev().map(|i| layers[i][layers[i].len() - 1]).collect();

    (left, right)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_linear() {
        let pts = [Point2::new(0.0, 0.0), Point2::new(10.0, 0.0)];
        let p = de_casteljau(&pts, 0.5);
        assert!((p.x - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_quadratic() {
        let pts = [
            Point2::new(0.0, 0.0),
            Point2::new(5.0, 10.0),
            Point2::new(10.0, 0.0),
        ];
        let p = de_casteljau(&pts, 0.5);
        assert!((p.x - 5.0).abs() < 1e-10);
        assert!((p.y - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_cubic() {
        let pts = [
            Point2::new(0.0, 0.0),
            Point2::new(0.0, 10.0),
            Point2::new(10.0, 10.0),
            Point2::new(10.0, 0.0),
        ];
        let p = de_casteljau(&pts, 0.5);
        assert!((p.x - 5.0).abs() < 1e-10);
        assert!((p.y - 7.5).abs() < 1e-10);
    }

    #[test]
    fn test_endpoints() {
        let pts = [
            Point2::new(1.0, 2.0),
            Point2::new(3.0, 4.0),
            Point2::new(5.0, 6.0),
        ];
        assert_eq!(de_casteljau(&pts, 0.0), pts[0]);
        assert_eq!(de_casteljau(&pts, 1.0), pts[2]);
    }

    #[test]
    fn test_single_point() {
        let pts = [Point2::new(5.0, 5.0)];
        assert_eq!(de_casteljau(&pts, 0.7), Point2::new(5.0, 5.0));
    }

    #[test]
    fn test_empty() {
        let p = de_casteljau(&[], 0.5);
        assert_eq!(p, Point2::new(0.0, 0.0));
    }

    #[test]
    fn test_split_matches() {
        let pts = [
            Point2::new(0.0, 0.0),
            Point2::new(3.0, 9.0),
            Point2::new(6.0, 3.0),
            Point2::new(10.0, 0.0),
        ];
        let (left, right) = de_casteljau_split(&pts, 0.5);
        // Left end should be first original point
        assert_eq!(left[0], pts[0]);
        // Right end should be last original point
        assert_eq!(*right.last().unwrap(), pts[3]);
        // Split point should match
        let split_l = de_casteljau(&left, 1.0);
        let split_r = de_casteljau(&right, 0.0);
        assert!((split_l.x - split_r.x).abs() < 1e-10);
        assert!((split_l.y - split_r.y).abs() < 1e-10);
    }

    #[test]
    fn test_higher_degree() {
        // Degree 4 Bézier
        let pts = [
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 5.0),
            Point2::new(5.0, 5.0),
            Point2::new(9.0, 5.0),
            Point2::new(10.0, 0.0),
        ];
        let p = de_casteljau(&pts, 0.5);
        assert!((p.x - 5.0).abs() < 1e-10);
        assert!(p.y > 0.0);
    }
}
