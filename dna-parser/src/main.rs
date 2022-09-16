use std::{fs, collections::HashSet, iter::FromIterator};
use itertools::Itertools;
use plotters::prelude::*;

struct OpenReadingFrame {
    start: usize,
    stop: usize,
    reverse: bool
}

struct CodonsTable {
    start: Vec<&'static str>,
    stop: Vec<&'static str>,
}

// read file with genome
fn get_dna_string(path: String) -> String {
    let raw_genome = fs::read_to_string(path).unwrap();
    // remove header line and replace newline by nothing
    raw_genome[(raw_genome.chars().position(|c| c == '\n').unwrap())..].replace("\n", "")
}

fn create_histogram(path: String, data: Vec<isize>) {
    let unique_data = HashSet::<isize>::from_iter(data.clone().into_iter());
    let data_amount: Vec<isize> = unique_data.iter().map(|f| data.iter().filter(|v| **v == *f).count() as isize).collect();

    let drawing_area = SVGBackend::new(&path, (1920 * 2, 1080 * 2)).into_drawing_area();
    drawing_area.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&drawing_area)
        .x_label_area_size(35)
        .y_label_area_size(40)
        .margin(50)
        .caption("Genome Length", ("sans-serif", 100.0))
        .build_cartesian_2d(0..unique_data.len() as isize, 0..*data_amount.iter().max().unwrap()).unwrap();

    chart
        .configure_mesh()
        .disable_x_mesh()
        .bold_line_style(&WHITE.mix(0.3))
        .y_desc("Length")
        .x_desc("Amount")
        .axis_desc_style(("sans-serif", 50))
        .draw()
        .unwrap();

    chart.draw_series(
        Histogram::vertical(&chart)
            .style(RED.mix(0.5).filled())
            .data(data.iter().map(|x: &isize| (*x, 1)))
            .margin(0),
    ).unwrap();

    drawing_area.present().unwrap();
}

// parse open reading frames
fn parse_orfs(dna: &String, reverse: bool, codons: &CodonsTable) -> Vec<OpenReadingFrame> {
    let mut orfs = vec![];
    for offset in [0, 1, 2] {
        let mut last_start = -1;
        let mut found_start = false;
        for pos in (offset..dna.len()).step_by(3) {
            if pos + 2 < dna.len() {
                let current_triplet = &dna[pos..pos + 3];
                if codons.start.contains(&current_triplet) && !found_start {
                    found_start = true;
                    last_start = pos as isize
                }
                if codons.stop.contains(&current_triplet) && found_start {
                    found_start = false;
                    // orfs.push((last_start as usize + 1, pos as usize + 3, true));
                    orfs.push(
                        OpenReadingFrame {
                            start: if !reverse {last_start as usize + 1} else {dna.len() - pos - 2},
                            stop: if !reverse {pos as usize + 3} else {dna.len() - last_start as usize},
                            reverse
                        }
                    );
                }
            }
        }
    }
    orfs
}

fn main() {
    // get DNA string
    let dna = get_dna_string("genome.fna".to_string());
    // let dna = "TTATGCATGCATAGATAA";

    // codons
    let codons = CodonsTable {
        start: vec!["ATG"],
        stop: vec!["TAG", "TAA", "TGA"] 
    };

    // store found frames
    let mut orfs = Vec::<OpenReadingFrame>::new();

    // parse first strip
    orfs.extend(parse_orfs(&dna, false, &codons));

    // reverse compliment
    let reversed_dna: String = dna.chars().map(|c| match c {
        'A' => 'T',
        'T' => 'A',
        'G' => 'C',
        'C' => 'G',
        _ => {panic!()},
    }).rev().collect();
    orfs.extend(parse_orfs(&reversed_dna, false, &codons));

    // print found orfs
    for i in &orfs {
        println!("{} {} {}", i.start, i.stop, ["-", "+"][if !i.reverse { 1 } else { 0 }])
    }
    println!("{}", orfs.len());

    // create histogram
    let data: Vec<isize> = orfs.iter().map(|o| (o.stop - o.start) as isize).sorted().collect();
    create_histogram(String::from("histogram.svg"), data);
}
