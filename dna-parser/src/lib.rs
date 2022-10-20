use plotters::prelude::*;
use serde::Deserialize;
use std::{
    collections::{BTreeMap, HashSet},
    fs::File,
    io::{BufRead, BufReader},
    iter::FromIterator,
};

#[derive(Debug, Deserialize)]
pub struct OpenReadingFrame {
    pub start: usize,
    pub stop: usize,
    pub reverse: bool,
}

pub struct CodonsTable {
    pub start: Vec<&'static str>,
    pub stop: Vec<&'static str>,
}

#[derive(Debug)]
pub struct GenomeInfo {
    pub header: String,
    pub genome: String,
}

// read all genomes from specified FASTA file
pub fn read_fasta_file(path: String) -> Vec<GenomeInfo> {
    // prepare variables
    let mut parsed_data = Vec::<GenomeInfo>::new();
    let file = File::open(path).unwrap();
    let mut header: Option<String> = None;
    let mut sequence: Option<String> = None;

    // read by lines
    for line_result in BufReader::new(file).lines() {
        if let Ok(line) = line_result {
            // if line red properly
            if line.starts_with(">") {
                // if line is header
                if let Some(hdr) = header.as_mut() {
                    // if header already exists push it to vector
                    parsed_data.push(GenomeInfo {
                        header: hdr.clone(),
                        genome: sequence.as_ref().unwrap().clone(),
                    });
                }
                // reset values
                header = Some(line.trim().to_string());
                sequence = Some(String::new());
            } else {
                // if line is sequence
                if let Some(seq) = sequence.as_mut() {
                    // push line to sequence if not None
                    seq.push_str(line.trim());
                }
            }
        }
    }

    // flush genome info if present
    if let Some(hdr) = header.as_mut() {
        parsed_data.push(GenomeInfo {
            header: hdr.clone(),
            genome: sequence.as_ref().unwrap().clone(),
        });
    }

    // return vector
    parsed_data
}

// Create histogram and save to  file
pub fn create_histogram(path: String, data: Vec<isize>) {
    let unique_data = HashSet::<isize>::from_iter(data.clone().into_iter());
    let data_amount: Vec<isize> = unique_data
        .iter()
        .map(|f| data.iter().filter(|v| **v == *f).count() as isize)
        .collect();

    let drawing_area = SVGBackend::new(&path, (1920 * 2, 1080 * 2)).into_drawing_area();
    drawing_area.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&drawing_area)
        .x_label_area_size(35)
        .y_label_area_size(40)
        .margin(50)
        .caption("Genome Length", ("sans-serif", 100.0))
        .build_cartesian_2d(
            (*unique_data.iter().min().unwrap()..*unique_data.iter().max().unwrap())
                .into_segmented(),
            0..*data_amount.iter().max().unwrap(),
        )
        .unwrap();

    chart
        .configure_mesh()
        .disable_x_mesh()
        .bold_line_style(&WHITE.mix(0.3))
        .y_desc("Amount")
        .x_desc("Length")
        .axis_desc_style(("sans-serif", 50))
        .draw()
        .unwrap();

    chart
        .draw_series(
            Histogram::vertical(&chart)
                .style(RED.mix(0.5).filled())
                .data(data.iter().map(|x: &isize| (*x, 1)))
                .margin(0),
        )
        .unwrap();

    drawing_area.present().unwrap();
}

// parse open reading frames
pub fn parse_orfs(dna: &str, reverse: bool, codons: &CodonsTable) -> Vec<OpenReadingFrame> {
    let mut orfs = vec![]; // vector for storing all found frames
    for offset in [0, 1, 2] {
        // check every offset
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
                    orfs.push(OpenReadingFrame {
                        start: if !reverse {
                            last_start as usize + 1
                        } else {
                            dna.len() - pos - 2
                        },
                        stop: if !reverse {
                            pos as usize + 3
                        } else {
                            dna.len() - last_start as usize
                        },
                        reverse,
                    });
                }
            }
        }
    }
    orfs
}

// parse amino acids from dna
pub fn parse_amino_acid(
    dna: &String,
    rev: &String,
    orf: &OpenReadingFrame,
    lut: &BTreeMap<&str, &str>,
) -> String {
    let mut ret: String = String::new();
    if !orf.reverse {
        for i in (orf.start..orf.stop).step_by(3) {
            let cur_codon = &dna[i - 1..i + 2];
            ret.push_str(lut.get_key_value(cur_codon).unwrap().1);
        }
    } else {
        for i in ((dna.len() - orf.stop + 1)..(dna.len() - orf.start + 1)).step_by(3) {
            let cur_codon = &rev[i - 1..i + 2];
            ret.push_str(lut.get_key_value(cur_codon).unwrap().1);
        }
    }
    ret
}

// reverse compliment dna
pub fn reverse_compliment(dna: &String) -> String {
    dna.chars()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            _ => {
                panic!()
            }
        })
        .rev()
        .collect::<String>()
}
