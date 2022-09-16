use std::{process::Command, os::unix::process::CommandExt};

use plotters::prelude::*;
use rand::Rng;
const OUT_FILE_NAME: &'static str = "histogram.svg";
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root = SVGBackend::new(OUT_FILE_NAME, (1920, 1080)).into_drawing_area();

    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .x_label_area_size(35)
        .y_label_area_size(40)
        .caption("Histogram Test", ("sans-serif", 50.0))
        .build_cartesian_2d((0u32..100u32).into_segmented(), 0u32..100u32)?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .bold_line_style(&WHITE.mix(0.3))
        .y_desc("Count")
        .x_desc("Bucket")
        .axis_desc_style(("sans-serif", 15))
        .draw()?;

    let data = {
        let mut a: Vec<u32> = vec![];
        let mut rng = rand::thread_rng();
        for i in 0..100 {
            a.push(rng.gen_range(0..100));
        }
        a
    };

    chart.draw_series(
        Histogram::vertical(&chart)
            .style(RED.mix(0.5).filled())
            .data(data.iter().map(|x: &u32| (*x, 1))),
    )?;

    // To avoid the IO failure being ignored silently, we manually call the present function
    root.present().expect("Unable to write result to file, please make sure 'plotters-doc-data' dir exists under current dir");
    println!("Result has been saved to {}", OUT_FILE_NAME);

    Command::new("nomacs").arg(OUT_FILE_NAME).exec();
    Ok(())
}