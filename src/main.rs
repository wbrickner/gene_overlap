mod args;
use std::{fs::File, io::{BufReader}, ops::Range, thread};
use intervaltree::{Element, IntervalTree};
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};
use serde::{Deserialize, Serialize};
use std::iter::FromIterator;
use arraystring::{ArrayString, typenum::{U4, U8, U16}};

#[derive(Deserialize)]
struct GeneRecord {
  #[serde(rename = "gene")]
  name: ArrayString<U16>,

  strand: ArrayString<U4>,

  #[serde(rename = "chr")]
  chromosome: ArrayString<U8>,

  start: usize,
  end: usize,
}

#[derive(Deserialize)]
struct ReadRecord {
  strand: ArrayString<U4>,
  
  #[serde(rename = "chr")]
  chromosome: ArrayString<U8>,
  
  #[serde(rename = "left")]
  start: usize,
  
  #[serde(rename = "right")]
  end: usize
}

#[derive(Serialize)]
struct GeneResult {
  name: ArrayString<U16>,

  strand: ArrayString<U4>,

  #[serde(rename = "chr")]
  chromosome: ArrayString<U8>,

  intersections: u32
}

fn main() -> std::io::Result<()> {
  let args = args::parse();
  let genes_path = args.genes.clone();
  let genes = thread::spawn(|| {
    println!("Loading genes file...");
    let f = File::open(genes_path).expect("Failed to read gene csv file");
    let reader = csv::Reader::from_reader(BufReader::new(f));
    let records = reader.into_deserialize::<GeneRecord>()
      .map(|record| record.unwrap())
      .collect::<Vec<GeneRecord>>();
    println!("Gene loading complete");
    records
  });

  let reads_path = args.reads.clone();
  let reads = thread::spawn(|| {
    println!("Loading reads file...");
    let f = File::open(reads_path).expect("Failed to read gene csv file");
    let reader = csv::Reader::from_reader(BufReader::new(f));

    let records = reader.into_deserialize::<ReadRecord>()
      .map(|record| record.unwrap())
      .collect::<Vec<ReadRecord>>();

    println!("Read loading complete");
    records
  });

  let (genes, reads) = (genes.join().expect("Failed to load genes"), reads.join().expect("Failed to load reads"));

  println!("Building interval tree...");
  let interval_tree: IntervalTree<usize, usize> = IntervalTree::from_iter(
    reads
      .par_iter()
      .enumerate()
      .map(|(idx, &ReadRecord { start, end, .. })| {
        Element {
          range: {
            if start <= end {
              Range { start, end }
            } else {
              Range { start: end, end: start }
            }
          },
          value: idx
        }
      })
      .collect::<Vec<Element<usize, usize>>>()
  );

  println!("Analyzing gene-read intersections over {} genes and {} reads...", genes.len(), reads.len());
  let genes_output = genes
    .par_iter()
    .map(|GeneRecord { name, strand, chromosome, start, end }| {
      let range = {
        if start <= end {
          Range { start: *start, end: *end }
        } else {
          Range { start: *end, end: *start }
        }
      };

      let intersections = interval_tree.query(range);
      let mut intersection_count: u32 = 0;
      for intersection in intersections {
        let read_index = intersection.value;
        let read = &reads[read_index];
        if &read.strand == strand && &read.chromosome == chromosome {
          intersection_count += 1;
        }
      }

      GeneResult {
        name: *name, 
        strand: *strand,
        chromosome: *chromosome, 
        intersections: intersection_count
      }
    })
    .collect::<Vec<GeneResult>>();

  let mut writer = csv::Writer::from_path(args.output).expect("Failed to write to output file");
  let mut total_intersections = 0;

  for record in genes_output.iter() {
    writer.serialize(record).unwrap();
    total_intersections += record.intersections;
  }
  writer.flush().expect("Failed to write output to disk");

  println!("Total intersections identified: {}", total_intersections);

  Ok(())
}
