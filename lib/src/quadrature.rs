// SPDX-LICENSE-IDENTIFIER: GPL-3.0-or-later

use std::{borrow::Borrow, f64};

//TODO: support generation of arbitrary gauss-kronrod rules

pub trait Quadrature<T> {
    /// Integrate over the region a..b and return the integral and approximate
    /// error
    fn integrate(
        &self,
        f: impl Fn(T) -> T,
        epsilon: impl Borrow<T>,
        bounds: (impl Borrow<T>, impl Borrow<T>),
    ) -> (T, T);
}

/// A struct providing an implementation of the [`trait@Quadrature`] trait for
/// the Tanh-Sinh quadrature method
#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct TanhSinh {
    /// The upper limit on iteration
    ///
    /// According to <https://www.genivia.com/files/qthsh.pdf>, 6 is "optimal"
    /// and 7 is "just as good". 6 is probably a good starting point
    pub iteration_limit: u64,
}

impl Quadrature<f64> for TanhSinh {
    /// According to page 24 of <https://www.genivia.com/files/qthsh.pdf>,
    /// 1e-9 is a good default for epsilon
    fn integrate(
        &self,
        f: impl Fn(f64) -> f64,
        epsilon: impl Borrow<f64>,
        (a, b): (impl Borrow<f64>, impl Borrow<f64>),
    ) -> (f64, f64) {
        tanh_sinh(f, epsilon, (a, b), self.iteration_limit)
    }
}

/// A struct providing an implementation of the [`trait@Quadrature`] trait for
/// the Gauss-Kronrod quadrature method
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct GaussKronrod<'a> {
    /// The upper limit on intervals
    pub interval_limit: u64,

    /// The quadrature rule to use
    pub rule: &'a [(f64, f64, Option<f64>)],
}

impl<'a> Quadrature<f64> for GaussKronrod<'a> {
    fn integrate(
        &self,
        f: impl Fn(f64) -> f64,
        epsilon: impl Borrow<f64>,
        (a, b): (impl Borrow<f64>, impl Borrow<f64>),
    ) -> (f64, f64) {
        gauss_kronrod(f, self.rule, epsilon, (a, b), self.interval_limit)
    }
}

/// Nodes and weights from G7 / K15 as a triplet of node, Kronrod weight,
/// Gaussian weight (if there is one)
#[allow(clippy::excessive_precision)]
pub const G7_K15: [(f64, f64, Option<f64>); 15] = [
    (
        9.914553711208126392068546975263285e-01,
        2.293532201052922496373200805896959e-02,
        None,
    ),
    (
        9.491079123427585245261896840478513e-01,
        6.309209262997855329070066318920429e-02,
        Some(1.294849661688696932706114326790820e-01),
    ),
    (
        8.648644233597690727897127886409262e-01,
        1.047900103222501838398763225415180e-01,
        None,
    ),
    (
        7.415311855993944398638647732807884e-01,
        1.406532597155259187451895905102379e-01,
        Some(2.797053914892766679014677714237796e-01),
    ),
    (
        5.860872354676911302941448382587296e-01,
        1.690047266392679028265834265985503e-01,
        None,
    ),
    (
        4.058451513773971669066064120769615e-01,
        1.903505780647854099132564024210137e-01,
        Some(3.818300505051189449503697754889751e-01),
    ),
    (
        2.077849550078984676006894037732449e-01,
        2.044329400752988924141619992346491e-01,
        None,
    ),
    (
        0.0,
        2.094821410847278280129991748917143e-01,
        Some(4.179591836734693877551020408163265e-01),
    ),
    (
        -2.077849550078984676006894037732449e-01,
        2.044329400752988924141619992346491e-01,
        None,
    ),
    (
        -4.058451513773971669066064120769615e-01,
        1.903505780647854099132564024210137e-01,
        Some(3.818300505051189449503697754889751e-01),
    ),
    (
        -5.860872354676911302941448382587296e-01,
        1.690047266392679028265834265985503e-01,
        None,
    ),
    (
        -7.415311855993944398638647732807884e-01,
        1.406532597155259187451895905102379e-01,
        Some(2.797053914892766679014677714237796e-01),
    ),
    (
        -8.648644233597690727897127886409262e-01,
        1.047900103222501838398763225415180e-01,
        None,
    ),
    (
        -9.491079123427585245261896840478513e-01,
        6.309209262997855329070066318920429e-02,
        Some(1.294849661688696932706114326790820e-01),
    ),
    (
        -9.914553711208126392068546975263285e-01,
        2.293532201052922496373200805896959e-02,
        None,
    ),
];

pub fn gauss_kronrod(
    f: impl Fn(f64) -> f64,
    rule: &[(f64, f64, Option<f64>)],
    epsilon: impl Borrow<f64>,
    (a, b): (impl Borrow<f64>, impl Borrow<f64>),
    interval_limit: u64,
) -> (f64, f64) {
    let mut n_intervals = 1;

    let mut kahan_t;
    let mut region_width;

    let mut gauss_kronrod_integral = 0.00;
    let mut gauss_kronrod_compensation;

    let mut gauss_integral;
    let mut gauss_compensation;

    let mut gauss_kronrod_acc;
    let mut gauss_acc;

    let mut relative_error = 0.00;
    let mut absolute_region_midpoint;
    let mut half_region_width;

    while n_intervals <= interval_limit {
        region_width = (b.borrow() - a.borrow()) / n_intervals as f64;

        half_region_width = region_width;
        half_region_width /= 2.0;

        gauss_kronrod_integral = 0.00;
        gauss_kronrod_compensation = 0.00;

        gauss_integral = 0.00;
        gauss_compensation = 0.00;

        for interval in 0..n_intervals {
            absolute_region_midpoint = (region_width * interval as f64)
                + a.borrow()
                + half_region_width;
            gauss_kronrod_acc = 0.00;
            gauss_acc = 0.00;

            for (node, gauss_kronrod_weight, gauss_weight) in rule {
                let mut function_input = half_region_width;
                function_input *= node;
                function_input += &absolute_region_midpoint;
                let mut y = f(function_input);

                if let Some(gauss_weight) = gauss_weight {
                    gauss_acc += y * gauss_weight;
                }

                y *= gauss_kronrod_weight;
                gauss_kronrod_acc += y;
            }

            gauss_kronrod_acc *= &half_region_width;
            gauss_acc *= &half_region_width;

            // the following is just the kahan summation algorithm

            gauss_kronrod_acc -= &gauss_kronrod_compensation;
            kahan_t = gauss_kronrod_integral + gauss_kronrod_acc;
            gauss_kronrod_compensation = kahan_t;
            gauss_kronrod_compensation -= &gauss_kronrod_integral;
            gauss_kronrod_compensation -= &gauss_kronrod_acc;
            gauss_kronrod_integral = kahan_t;

            gauss_acc -= &gauss_compensation;
            kahan_t = gauss_integral + gauss_acc;
            gauss_compensation = kahan_t;
            gauss_compensation -= &gauss_integral;
            gauss_compensation -= &gauss_acc;
            gauss_integral = kahan_t;
        }

        relative_error = gauss_kronrod_integral;
        relative_error -= &gauss_integral;
        relative_error /= &gauss_kronrod_integral;
        relative_error = relative_error.abs();
        if &relative_error <= epsilon.borrow() {
            break;
        }

        n_intervals <<= 1;
    }

    (gauss_kronrod_integral, relative_error)
}

// the following function is more or less a rust translation (with a number of
// adjustments made to fit better within the reedbed source code) of the
// c routine on page 24 of https://www.genivia.com/files/qthsh.pdf, which is
// licensed under the mit license. as a consequence, the following routine is
// also licensed under the mit license, which is reproduced below:
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the “Software”),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
pub fn tanh_sinh(
    f: impl Fn(f64) -> f64,
    epsilon: impl Borrow<f64>,
    (a, b): (impl Borrow<f64>, impl Borrow<f64>),
    limit: u64,
) -> (f64, f64) {
    let tolerance = epsilon.borrow() * 10.00;

    let mut region_center = a.borrow() + b.borrow();
    region_center /= 2.0;

    let mut half_region_width = b.borrow() - a.borrow();
    half_region_width /= 2.0;

    let mut s = f(region_center);
    let mut h = 2.0f64;
    let mut iteration = 0;

    //TODO: get more descriptive names for these

    let mut v;

    let mut p;
    let mut q;
    let mut fp;
    let mut fm;
    let mut t;
    let mut eh;

    let mut u;
    let mut r;
    let mut w;
    let mut x;

    let mut temporary;

    loop {
        p = 0.00;
        fp = 0.00;
        fm = 0.00;

        h /= 2.00;
        eh = h.exp();
        t = eh;
        if iteration > 0 {
            eh = eh.powi(2);
        }

        loop {
            u = t.recip();
            u -= &t;
            u = u.exp();

            r = u;
            r += 1.00;
            r = r.recip();
            r *= &u;
            r *= 2.00;

            w = u;
            w += 1.00;
            w = w.recip();
            w *= &r;

            x = t;
            x = x.recip();
            x += &t;

            w *= &x;

            x = half_region_width;
            x *= &r;

            temporary = a.borrow() + x;
            if &temporary > a.borrow() {
                let y = f(temporary);
                if y.is_finite() {
                    fp = y;
                }
            }

            temporary = b.borrow() - x;
            if &temporary < b.borrow() {
                let y = f(temporary);
                if y.is_finite() {
                    fm = y;
                }
            }

            q = fp;
            q += &fm;
            q *= &w;
            p += &q;
            t *= &eh;

            temporary = p.abs() * epsilon.borrow();
            if q.abs() <= temporary {
                break;
            }
        }

        v = s;
        v -= &p;
        s += &p;
        v = v.abs();

        temporary = s.abs() * tolerance;
        iteration += 1;
        if v <= temporary || iteration > limit {
            break;
        }
    }

    let mut e = s.abs();
    e += epsilon.borrow();
    e = e.recip();
    e *= v;

    (half_region_width * s * h, e)
}

#[cfg(test)]
mod tests {
    use super::*;

    //TODO: we should probably parameterize precision here somewhere
    //TODO: test tanh_sinh

    const EPSILON: f64 = 1e-10;

    #[test]
    fn integrate_constant() {
        let (val, _) =
            gauss_kronrod(|_| 3.00, &G7_K15, &EPSILON, (0.00, 5.00), 64);
        assert!((val - 15.00).abs() < EPSILON);
    }

    #[test]
    fn integrate_line() {
        let (val, _) =
            gauss_kronrod(|x| 3.00 + x, &G7_K15, &EPSILON, (0.00, 5.00), 64);
        assert!((val - 15.00 - 25.00 / 2.00).abs() < EPSILON);
    }

    #[test]
    fn integrate_parabola() {
        let (val, _) = gauss_kronrod(
            |x| x + 0.50 * x.powi(2),
            &G7_K15,
            &EPSILON,
            (0.00, 5.00),
            64,
        );
        assert!((val - 25.00 / 2.00 - 125.00 / 6.00).abs() < EPSILON);
    }

    #[test]
    fn integrate_sin_squared() {
        let (val, _) = gauss_kronrod(
            |x| x.sin().powi(2),
            &G7_K15,
            &EPSILON,
            (0.00, 2.00 * f64::consts::PI),
            64,
        );
        assert!((val - f64::consts::PI).abs() < EPSILON);

        let (val, _) = gauss_kronrod(
            |x| x.sin().powi(2),
            &G7_K15,
            &EPSILON,
            (0.00, -2.00 * f64::consts::PI),
            64,
        );
        assert!((val + f64::consts::PI).abs() < EPSILON);
    }
}
