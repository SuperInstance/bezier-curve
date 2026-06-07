# bezier-curve

Bézier curves and surfaces in pure Rust.

## Features

- Linear, quadratic, and cubic Bézier curves
- De Casteljau's algorithm for arbitrary degree
- Curve splitting at any parameter
- Arc length computation and parameterization
- Degree elevation (quadratic → cubic)

## Usage

```rust
use bezier_curve::{CubicBezier, Point2, arc_length, de_casteljau};

let curve = CubicBezier::new(
    Point2::new(0.0, 0.0),
    Point2::new(1.0, 3.0),
    Point2::new(3.0, 3.0),
    Point2::new(4.0, 0.0),
);
let point = curve.evaluate(0.5);
let (left, right) = curve.split(0.5);
let length = curve.approximate_length(100);
```

Zero external dependencies. Pure `std` Rust.

## License

MIT
