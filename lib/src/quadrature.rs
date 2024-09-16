// SPDX-LICENSE-IDENTIFIER: GPL-3.0-or-later

use rug::{Assign, Float};

//TODO: genericize the parameters here. taking arbitrary-precision floats
//      everywhere is excessive

//TODO: allow selection of nodes and weights

/// Nodes and weights from G7 / K15 as a triplet of node, Kronrod weight,
/// Gaussian weight (if there is one)
const G7_K15: [(f64, f64, Option<f64>); 15] = [
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
    f: impl Fn(Float) -> Float,
    (a, b): (&Float, &Float),
    relative_tolerance: &Float,
    interval_limit: u64,
    precision: u64,
) -> Float {
    let mut n_intervals = 1;

    let mut kahan_t = Float::new_64(precision);
    let mut region_width = Float::new_64(precision);

    let mut gauss_kronrod_integral = Float::new_64(precision);
    let mut gauss_kronrod_compensation = Float::new_64(precision);

    let mut gauss_integral = Float::new_64(precision);
    let mut gauss_compensation = Float::new_64(precision);

    let mut gauss_kronrod_acc = Float::new_64(precision);
    let mut gauss_acc = Float::new_64(precision);

    let mut relative_error = Float::new_64(precision);
    let mut absolute_region_midpoint = Float::new_64(precision);
    let mut half_region_width = Float::new_64(precision);

    while n_intervals <= interval_limit {
        region_width.assign(b - a);
        region_width /= n_intervals;

        half_region_width.assign(&region_width);
        half_region_width /= 2.0;

        gauss_kronrod_integral.assign(0);
        gauss_kronrod_compensation.assign(0);

        gauss_integral.assign(0);
        gauss_compensation.assign(0);

        for interval in 0..n_intervals {
            absolute_region_midpoint.assign(&region_width);
            absolute_region_midpoint *= &interval;
            absolute_region_midpoint += a;
            absolute_region_midpoint += &half_region_width;

            gauss_kronrod_acc.assign(0);
            gauss_acc.assign(0);

            for (node, gauss_kronrod_weight, gauss_weight) in G7_K15 {
                let mut function_input = Float::with_val_64(precision, &half_region_width);
                function_input *= node;
                function_input += &absolute_region_midpoint;
                let mut y = f(function_input);

                if let Some(gauss_weight) = gauss_weight {
                    let mut y = y.clone();
                    y *= gauss_weight;
                    gauss_acc += y;
                }

                y *= gauss_kronrod_weight;
                gauss_kronrod_acc += y;
            }

            gauss_kronrod_acc *= &half_region_width;
            gauss_acc *= &half_region_width;

            // the following is just the kahan summation algorithm

            gauss_kronrod_acc -= &gauss_kronrod_compensation;
            kahan_t.assign(&gauss_kronrod_integral + &gauss_kronrod_acc);
            gauss_kronrod_compensation.assign(&kahan_t);
            gauss_kronrod_compensation -= &gauss_kronrod_integral;
            gauss_kronrod_compensation -= &gauss_kronrod_acc;
            gauss_kronrod_integral.assign(&kahan_t);

            gauss_acc -= &gauss_compensation;
            kahan_t.assign(&gauss_integral + &gauss_acc);
            gauss_compensation.assign(&kahan_t);
            gauss_compensation -= &gauss_integral;
            gauss_compensation -= &gauss_acc;
            gauss_integral.assign(&kahan_t);
        }

        println!("gk compensation: {}, g compensation: {}", &gauss_kronrod_compensation, &gauss_compensation);

        relative_error.assign(&gauss_kronrod_integral);
        relative_error -= &gauss_integral;
        relative_error /= &gauss_kronrod_integral;
        relative_error.abs_mut();

        if &relative_error <= relative_tolerance {
            break;
        }

        n_intervals <<= 1;
    }

    gauss_kronrod_integral
}
