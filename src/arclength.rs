//! Arc length computation and parameterization for Bézier curves.

use crate::linear::Point2;
use crate::casteljau::de_casteljau;

/// Compute the approximate arc length of a Bézier curve defined by control points,
/// using adaptive Simpson's rule on the derivative magnitude.
pub fn arc_length(control_points: &[Point2], tolerance: f64) -> f64 {
    if control_points.len() < 2 {
        return 0.0;
    }
    adaptive_simpson(control_points, 0.0, 1.0, tolerance, 20)
}

fn speed(control_points: &[Point2], t: f64) -> f64 {
    let dt = 1e-6;
    let t0 = (t - dt).max(0.0);
    let t1 = (t + dt).min(1.0);
    let p0 = de_casteljau(control_points, t0);
    let p1 = de_casteljau(control_points, t1);
    let actual_dt = t1 - t0;
    if actual_dt < 1e-12 {
        0.0
    } else {
        p0.distance_to(&p1) / actual_dt
    }
}

fn adaptive_simpson(cps: &[Point2], a: f64, b: f64, tol: f64, depth: usize) -> f64 {
    let c = (a + b) / 2.0;
    let fa = speed(cps, a);
    let fb = speed(cps, b);
    let fc = speed(cps, c);
    let h = b - a;

    let whole = (h / 6.0) * (fa + 4.0 * fc + fb);
    let left = (h / 12.0) * (fa + 4.0 * speed(cps, (a + c) / 2.0) + fc);
    let right = (h / 12.0) * (fc + 4.0 * speed(cps, (c + b) / 2.0) + fb);

    if depth == 0 || (left + right - whole).abs() <= 15.0 * tol {
        left + right + (left + right - whole) / 15.0
    } else {
        let l = adaptive_simpson(cps, a, c, tol / 2.0, depth - 1);
        let r = adaptive_simpson(cps, c, b, tol / 2.0, depth - 1);
        l + r
    }
}

/// Build a lookup table for arc-length parameterization.
/// Returns a vector of (t, arc_length) pairs.
pub fn build_arc_length_lut(control_points: &[Point2], steps: usize) -> Vec<(f64, f64)> {
    let steps = steps.max(1);
    let mut lut = Vec::with_capacity(steps + 1);
    let mut accum = 0.0;
    let mut prev = de_casteljau(control_points, 0.0);
    lut.push((0.0, 0.0));

    for i in 1..=steps {
        let t = i as f64 / steps as f64;
        let curr = de_casteljau(control_points, t);
        accum += prev.distance_to(&curr);
        lut.push((t, accum));
        prev = curr;
    }

    lut
}

/// Given an arc-length parameter s ∈ [0, total_length], find the corresponding t.
pub fn arc_length_to_t(lut: &[(f64, f64)], s: f64) -> f64 {
    if lut.is_empty() {
        return 0.0;
    }
    if s <= 0.0 {
        return 0.0;
    }
    let total = lut.last().unwrap().1;
    if s >= total {
        return 1.0;
    }

    // Binary search for the interval
    let mut lo = 0;
    let mut hi = lut.len() - 1;
    while lo < hi - 1 {
        let mid = (lo + hi) / 2;
        if lut[mid].1 < s {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    let (t0, s0) = lut[lo];
    let (t1, s1) = lut[hi];
    let ds = s1 - s0;
    if ds.abs() < 1e-15 {
        t0
    } else {
        t0 + (s - s0) / ds * (t1 - t0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_arc_length_line() {
        let pts = [Point2::new(0.0, 0.0), Point2::new(3.0, 4.0)];
        let len = arc_length(&pts, 1e-6);
        assert!((len - 5.0).abs() < 0.01, "got {}", len);
    }

    #[test]
    fn test_arc_length_zero_for_point() {
        let pts = [Point2::new(5.0, 5.0)];
        assert!((arc_length(&pts, 1e-6)).abs() < 1e-10);
    }

    #[test]
    fn test_arc_length_curve_longer_than_chord() {
        let pts = [
            Point2::new(0.0, 0.0),
            Point2::new(0.0, 10.0),
            Point2::new(10.0, 10.0),
            Point2::new(10.0, 0.0),
        ];
        let len = arc_length(&pts, 1e-6);
        let chord = pts[0].distance_to(&pts[3]);
        assert!(len > chord, "arc {} should be > chord {}", len, chord);
    }

    #[test]
    fn test_lut_monotonic() {
        let pts = [
            Point2::new(0.0, 0.0),
            Point2::new(5.0, 10.0),
            Point2::new(10.0, 0.0),
        ];
        let lut = build_arc_length_lut(&pts, 50);
        for window in lut.windows(2) {
            assert!(window[0].1 <= window[1].1, "LUT not monotonic");
        }
    }

    #[test]
    fn test_arc_length_to_t_roundtrip() {
        let pts = [
            Point2::new(0.0, 0.0),
            Point2::new(5.0, 10.0),
            Point2::new(10.0, 0.0),
        ];
        let lut = build_arc_length_lut(&pts, 100);
        let total = lut.last().unwrap().1;
        for i in 0..=10 {
            let s = i as f64 / 10.0 * total;
            let t = arc_length_to_t(&lut, s);
            assert!(t >= 0.0 && t <= 1.0, "t={} out of range", t);
        }
    }

    #[test]
    fn test_arc_length_to_t_endpoints() {
        let pts = [
            Point2::new(0.0, 0.0),
            Point2::new(5.0, 10.0),
            Point2::new(10.0, 0.0),
        ];
        let lut = build_arc_length_lut(&pts, 50);
        assert!((arc_length_to_t(&lut, 0.0)).abs() < 1e-10);
        assert!((arc_length_to_t(&lut, lut.last().unwrap().1) - 1.0).abs() < 1e-10);
    }
}
