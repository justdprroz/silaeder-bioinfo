use std::fs;
use itertools::Itertools;
use plotters::prelude::*;

fn get_dna_string(path: String) -> String {
    // read file with genome
    let raw_genome = fs::read_to_string(path).unwrap();
    // remove header line and replace newline by nothing
    raw_genome[(raw_genome.chars().position(|c| c == '\n').unwrap())..].replace("\n", "")
}

fn main() {
    // get DNA string
    let dna = get_dna_string("genome.fna".to_string());
    // let dna = "TTATGCATGCATAGATAA";

    // codons
    let start_codon = "ATG";
    let stop_codons = ["TAG", "TAA", "TGA"];

    // store found frames
    let mut orfs = Vec::<(usize, usize, bool)>::new();

    // check forward
    for offset in [0, 1, 2] {
        let mut last_start = -1;
        let mut found_start = false;
        for pos in (offset..dna.len()).step_by(3) {
            if pos + 2 < dna.len() {
                let current_triplet = &dna[pos..pos + 3];
                if current_triplet == start_codon && !found_start {
                    found_start = true;
                    last_start = pos as isize
                }
                if stop_codons.contains(&current_triplet) && found_start {
                    found_start = false;
                    orfs.push((last_start as usize + 1, pos as usize + 3, true));
                }
            }
        }
    }

    // reverse and reverse compliment
    let reversed_genome: String = dna.chars().map(|c| match c {
        'A' => 'T',
        'T' => 'A',
        'G' => 'C',
        'C' => 'G',
        _ => {panic!()},
    }).rev().collect();

    // check backwards
    for offset in [0, 1, 2] {
        let mut last_start = -1;
        let mut found_start = false;
        for pos in (offset..reversed_genome.len()).step_by(3) {
            if pos + 2 < reversed_genome.len() {
                let current_triplet = &reversed_genome[pos..pos + 3];
                if current_triplet == start_codon && !found_start {
                    found_start = true;
                    last_start = pos as isize
                }
                if stop_codons.contains(&current_triplet) && found_start {
                    found_start = false;
                    orfs.push((reversed_genome.len() - pos - 2, reversed_genome.len() - last_start as usize, false));
                }
            }
        }
    }

    // print found orfs
    for i in &orfs {
        // println!("{} {} {}", i.0, i.1, ["-", "+"][if i.2 { 1 } else { 0 }])
    }
    // println!("{}", orfs.len());
    let data: Vec<i32> = orfs.iter().map(|o| (o.1 - o.0) as i32).sorted().collect();
    // let unique_lengths = data.iter().map(|f| *f).unique().collect::<Vec<i32>>().len() as i32;
    // let max = *data.iter().max().unwrap();
    // // println!("{} {}", unique_lengths, max);

    // let drawing_area = SVGBackend::new("histogram_vertical.svg", (1000, 1000)).into_drawing_area();
    // drawing_area.fill(&WHITE).unwrap();

    // let mut chart_builder = ChartBuilder::on(&drawing_area);
    // chart_builder.margin(5).set_left_and_bottom_label_area_size(20);

    // let mut chart_context = chart_builder.build_cartesian_2d((0..100).into_segmented(), 0..100).unwrap();
    // chart_context.configure_mesh().draw().unwrap();
    // chart_context.draw_series(Histogram::vertical(&chart_context).style(BLUE.filled()).margin(10)
    //     .data(data.iter().map(|x| ((*x as f32).sqrt() as i32, 1)))).unwrap();
    let root = BitMapBackend::new("plotters-doc-data/0.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let mut chart = ChartBuilder::on(&root)
        .caption("y=x^2", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(-1f32..1f32, -0.1f32..1f32).unwrap();

    chart.configure_mesh().draw().unwrap();

    chart
        .draw_series(LineSeries::new(
            data.iter().enumerate().map(|d| ((*d.1) as f32, d.0 as f32)),
            &RED,
        )).unwrap()
        .label("y = x^2")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw().unwrap();

    root.present().unwrap();
}
