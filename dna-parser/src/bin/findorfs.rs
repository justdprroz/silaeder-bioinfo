use colored::Colorize;
use itertools::Itertools;
use std::{
    collections::{BTreeMap, HashMap},
    fs::File,
    io::Write,
    process,
};

use dna_parser::{
    create_histogram, parse_amino_acid, parse_orfs, read_fasta_file, reverse_compliment,
    CodonsTable, GenomeInfo, OpenReadingFrame,
};

fn main() {
    // Parse args and init settings
    let mut print_to_file = false;
    let mut draw_histogram = false;
    let mut path = String::from("genome.fna");
    let args = std::env::args().collect::<Vec<String>>();
    for arg_index in 0..args.len() {
        let arg = &args[arg_index];
        if arg.starts_with("--") {
            if arg == "--help" {
                eprintln!("{}", "Using:".yellow());
                eprintln!("\t{}", "'--help' - print this message and exit".yellow());
                eprintln!(
                    "\t{}",
                    format!(
                        "'--input {filename}' - read from {filename}",
                        filename = "filename".green()
                    )
                    .yellow()
                );
                eprintln!(
                    "\t{}",
                    format!("'--to-file' - save results in file").yellow()
                );
                eprintln!("\t{}", format!("'--hist' - draw histogram").yellow());
                process::exit(0);
            }
            if arg == "--input" {
                path = args[arg_index + 1].clone();
            }
            if arg == "--to-file" {
                print_to_file = true;
            }
            if arg == "--hist" {
                draw_histogram = true;
            }
        }
    }

    // Print run settings
    eprintln!("{}", "Run settings:".cyan());
    eprintln!(
        "\t{}",
        format!(
            "File to use {}{}",
            path.magenta(),
            if path == String::from("genome.fna") {
                " (Using default genome!)".yellow()
            } else {
                "".clear()
            }
        )
        .blue()
    );
    eprintln!(
        "\t{}",
        if print_to_file {
            "Save to file".blue()
        } else {
            "Print to stdout".blue()
        }
    );
    eprintln!(
        "\t{}",
        if draw_histogram {
            "Draw histogram".blue()
        } else {
            format!("Do {} draw histogram", "NOT".red()).blue()
        }
    );

    // parse genomes from specified file
    let genomes: Vec<GenomeInfo> = read_fasta_file(path);

    if genomes.len() > 0 {
        eprintln!("Processing {} genome(s)", genomes.len())
    } else {
        eprintln!("No genomes found")
    }

    // codons
    let codons = CodonsTable {
        start: vec!["ATG"],
        stop: vec!["TAG", "TAA", "TGA"],
    };

    for genome in genomes {
        // get DNA string
        let dna = genome.genome;
        // eprintln!("{}", &dna[2800..3733]);

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
            ("*", vec!["TAG", "TAA", "TGA"]),
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
        if print_to_file {
            let mut file = File::create(format!("output{}.txt", genome.header)).unwrap();
            writeln!(&mut file, "Current dna: {}", format!("{}", genome.header)).unwrap();
            for orf in &orfs {
                writeln!(
                    &mut file,
                    "ORF: {}, {}, {}\t\tAmino acids: {}",
                    orf.start,
                    orf.stop,
                    ["-", "+"][if !orf.reverse { 1 } else { 0 }],
                    parse_amino_acid(&dna, &reversed_dna, &orf, &codon_to_aminoacid),
                )
                .unwrap();
            }
            writeln!(&mut file, "ORFs amount: {}", format!("{}", orfs.len()),).unwrap();
        } else {
            println!(
                "{} {}",
                "Current dna:".blue(),
                format!("{}", genome.header).cyan(),
            );
            for orf in &orfs {
                println!(
                    "{} {}{c} {}{c} {}\t\t{} {}",
                    "ORF:".green(),
                    format!("{}", orf.start).red(),
                    format!("{}", orf.stop).red(),
                    format!("{}", ["-", "+"][if !orf.reverse { 1 } else { 0 }]).red(),
                    "Acid:".green(),
                    parse_amino_acid(&dna, &reversed_dna, &orf, &codon_to_aminoacid).magenta(),
                    c = ",".green(),
                );
            }
            println!(
                "{} {}",
                "ORFs amount:".blue(),
                format!("{}", orfs.len()).cyan(),
            );
        }

        // create histogram
        if draw_histogram {
            let data: Vec<isize> = orfs
                .iter()
                .map(|o| (o.stop - o.start) as isize)
                .sorted()
                .collect();
            create_histogram(
                String::from(format!("histogram{}.svg", genome.header)),
                data,
            );
        }
    }
}
