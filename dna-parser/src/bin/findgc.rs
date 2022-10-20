use dna_parser::read_fasta_file;

fn gc_content() {
    // println!("\"start\",\"stop\",\"strand\",\"gc_content\"");
    let dna = read_fasta_file("genome.fna".to_string())[0].genome.clone();
    let mut rdr = csv::Reader::from_path("genome_annotation.csv").unwrap();
    let mut gc_contents: Vec<f64> = vec![];
    for result in rdr.records() {
        let record = result.unwrap();
        let start = record[2].parse::<usize>().unwrap() - 1;
        let stop = record[3].parse::<usize>().unwrap() - 1;
        let strand = record[4].as_bytes()[0] as char;
        let mut c_cg = 0f64;
        let mut c_all = 0f64;
        for i in start..stop {
            if dna.as_bytes()[i] as char == 'G' || dna.as_bytes()[i] as char == 'C' {
                c_cg += 1f64;
            }
            c_all += 1f64;
        }
        gc_contents.push(c_cg / c_all);
        println!("{},{},\"{}\",{}", start, stop, strand, c_cg / c_all);
    }
}

fn main() {
    gc_content();
}
