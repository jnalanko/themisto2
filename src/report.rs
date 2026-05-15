use std::io::{BufRead, Write};

#[derive(serde::Deserialize, Debug)]
pub struct PseudoalignmentRecord<'a> {
    #[serde(borrow)]
    pub name: &'a str,
    pub colors: Vec<usize>,
    #[serde(default)]
    pub kmer_hits: Option<Vec<usize>>,
    #[serde(default)]
    pub bases_covered: Option<Vec<usize>>,
}

#[derive(serde::Serialize, Debug)]
pub struct ReportData {
    n_reads: usize,
    n_positive_reads: usize,
    positive_by_color: Vec<usize>,
    unique_positive_by_color: Vec<usize>,
}

impl ReportData {
    fn new(n_colors: usize) -> Self {
        Self { 
            n_reads: 0, 
            n_positive_reads: 0, 
            positive_by_color: vec![n_colors], 
            unique_positive_by_color: vec![n_colors]
        }
    }

    fn add_record(&mut self, rec: &PseudoalignmentRecord) {

        for color in rec.colors.iter() {
            self.positive_by_color[*color] += 1;
        }

        if rec.colors.len() == 1 {
            self.unique_positive_by_color[*rec.colors.first().unwrap()] += 1;
        }

        if rec.colors.len() > 0 {
            self.n_positive_reads += 1;
        }

        self.n_reads += 1;
    }
}

pub fn report(mut reader: impl BufRead, mut writer: impl Write, color_names: &[String]) {
    let mut line = String::new();
    let mut line_no = 0_usize;
    let mut data = ReportData::new(color_names.len());
    loop {
        line.clear();
        let n = reader.read_line(&mut line).unwrap();
        if n == 0 { break; } // EOF
        if line.ends_with('\n') { line.pop(); }
        if line.is_empty() {
            panic!("Empty line on line {} (zero-based indexing)", line_no);
        }
        let record: PseudoalignmentRecord = serde_json::from_str(&line)
            .unwrap_or_else(|e| panic!("Failed to parse JSON on line {}: {} (zero-indexing)", line_no, e));
        data.add_record(&record);

        line_no += 1;
    }

    writer.write(serde_json::to_string_pretty(&data).unwrap().as_bytes());
    writer.flush().unwrap();
}
