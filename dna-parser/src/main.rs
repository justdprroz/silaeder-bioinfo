use itertools::Itertools;
use plotters::prelude::*;
use std::{
    collections::{BTreeMap, HashMap, HashSet},
    fs,
    iter::FromIterator,
};

struct OpenReadingFrame {
    start: usize,
    stop: usize,
    reverse: bool,
}

struct CodonsTable {
    start: Vec<&'static str>,
    stop: Vec<&'static str>,
}

// read file with genome
fn get_dna_string(path: String) -> String {
    let raw_genome = fs::read_to_string(path).unwrap();
    raw_genome[(raw_genome.chars().position(|c| c == '\n').unwrap())..].replace("\n", "")
}

// Create histogram and save to  file
fn create_histogram(path: String, data: Vec<isize>) {
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
            0..unique_data.len() as isize,
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
fn parse_orfs(dna: &str, reverse: bool, codons: &CodonsTable) -> Vec<OpenReadingFrame> {
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
fn parse_amino_acid(dna: &String, rev: &String, orf: &OpenReadingFrame, lut: &BTreeMap<&str, &str>) -> String {
    let mut ret: String = String::new();
    if !orf.reverse {
        for i in (orf.start..orf.stop - 3).step_by(3) {
            let cur_codon = &dna[i - 1..i + 2];
            ret.push_str(lut.get_key_value(cur_codon).unwrap().1);
        }
    } else {
        for i in ((dna.len() - orf.stop + 1)..(dna.len() - orf.start - 2)).step_by(3) {
            let cur_codon = &rev[i - 1..i + 2];
            ret.push_str(lut.get_key_value(cur_codon).unwrap().1);
        }
    }
    ret
}

// reverse compliment dna
fn reverse_compliment(dna: &String) -> String {
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

fn main() {
    // get DNA string
    // let dna = get_dna_string("genome.fna".to_string());
    let dna = "TTATGCATGCATAGATAA".to_string();

    // codons
    let codons = CodonsTable {
        start: vec!["ATG"],
        stop: vec!["TAG", "TAA", "TGA"],
    };

    // codon <=> amino acid mapping
    let aminoacid_to_codon: HashMap<&str, Vec<&str>> = HashMap::from([
        ("A", vec!["GCA", "GCC", "GCG", "GCT"]),
        ("C", vec!["TGC", "TGT"]),
        ("D", vec!["GAC", "GAT"]),
        ("E", vec!["GAA", "GAG"]),
        ("F", vec!["TTC", "TTT"]),
        ("G", vec!["GGA", "GGC", "GGG", "GGT"]),
        ("H", vec!["CAC", "CAT"]),
        ("I", vec!["ATA", "ATC", "ATT"]),
        ("K", vec!["AAA", "AAG"]),
        ("L", vec!["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"]),
        ("M", vec!["ATG"]),
        ("N", vec!["AAC", "AAT"]),
        ("P", vec!["CCA", "CCC", "CCG", "CCT"]),
        ("Q", vec!["CAA", "CAG"]),
        ("R", vec!["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"]),
        ("S", vec!["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"]),
        ("T", vec!["ACA", "ACC", "ACG", "ACT"]),
        ("V", vec!["GTA", "GTC", "GTG", "GTT"]),
        ("W", vec!["TGG"]),
        ("Y", vec!["TAC", "TAT"]),
    ]);

    let codon_to_aminoacid: BTreeMap<&str, &str> = {
        let mut tmp: BTreeMap<&str, &str> = BTreeMap::new();
        for aminoacid in aminoacid_to_codon {
            for codon in aminoacid.1 {
                tmp.insert(codon, aminoacid.0);
            }
        }
        tmp
    };

    // store found frames
    let mut orfs = Vec::<OpenReadingFrame>::new();

    // parse first strip
    orfs.extend(parse_orfs(&dna, false, &codons));

    // reverse compliment
    let reversed_dna = reverse_compliment(&dna);

    // parse reversed strip
    orfs.extend(parse_orfs(&reversed_dna, true, &codons));

    // print found orfs
    println!("Current dna: {}", dna);
    for orf in &orfs {
        println!(
            "ORF: {}, {}, {}\t\tAcid: {}",
            orf.start,
            orf.stop,
            ["-", "+"][if !orf.reverse { 1 } else { 0 }],
            parse_amino_acid(&dna, &reversed_dna, &orf, &codon_to_aminoacid)
        );
    }
    println!("ORFs amount: {}", orfs.len());

    // create histogram
    let data: Vec<isize> = orfs
        .iter()
        .map(|o| (o.stop - o.start) as isize)
        .sorted()
        .collect();
    create_histogram(String::from("histogram.svg"), data);
}
