#[macro_use]
extern crate clap;
extern crate debruijn;
extern crate dna_io;
extern crate bloom;
extern crate fnv;
extern crate hashbrown;
//extern crate cpuprofiler;
//use cpuprofiler::PROFILER;


use clap::{App};

use hashbrown::{HashMap, HashSet};
//use std::collections::{HashMap};
use fnv::FnvHasher;
use std::hash::BuildHasherDefault;
type FnvHashMap<T,V> = HashMap<T,V, BuildHasherDefault<FnvHasher>>;
type FnvHashSet<V> = HashSet<V, BuildHasherDefault<FnvHasher>>;
use bloom::CountingBloomFilter;
use std::collections::hash_map::RandomState;
type FnvCountingBloomFilter = CountingBloomFilter<BuildHasherDefault<FnvHasher>, RandomState>;

use std::io::{BufWriter, Write};
use std::io::stdout;
use std::fs::File;
use debruijn::*;
use debruijn::kmer::*;

use std::cmp::min;

static mut KMER_SIZE: usize = 21;
static mut KMER_SIZE_MINUS_ONE: usize = 20; // how bad is this. i need speed tho
static mut KMER_MASK:u64=0;

fn main() {
    //PROFILER.lock().unwrap().start("./myprof.profile").expect("Couldn't start");
    let parameters = load_params();
    let kmer_counts = count_kmers_fastq(&parameters);
    detect_pairs(kmer_counts, &parameters);
    //PROFILER.lock().unwrap().stop().expect("Couldn't stop");
}

fn detect_pairs(kmer_counts: FnvHashMap<u64,[u16; 4]>, params: &Params) {
    let mut out_writer: Box<Write> = match params.output {
        Some(ref x) => Box::new(File::create(&x).expect("Unable to open file for writing")),
        None => Box::new(stdout()),
    };
    let mut out_writer = BufWriter::new(out_writer);
    let mut hist = BufWriter::new(
        File::create(&params.output_hist)
        .expect("Unable to create file for writing"));

    //let visited: FnvHashSet<u64> = FnvHashSet::default(); 
    let a_mask = 0;
    let c_mask = ( 1 << ( KX::K() - 1 ) ); // or a C in the middle
    let g_mask = ( 2 << ( KX::K() - 1 ) ); // etc
    let t_mask = ( 3 << ( KX::K() - 1 ) ); // etc
    let masks: [u64; 4] = [ a_mask, c_mask, g_mask, t_mask ];
    
    let mut total = 0;
    let mut best_count = 0;
    let mut second_best_count = 0;
    let mut best_index = 0;
    let mut second_best_index = 0; 
    let mut full_hist = [0u64; 1000];   
    for (middle_base_invariant_kmer, counts) in &kmer_counts {
        for (index, count) in counts.iter().enumerate() {
            total += *count;
            full_hist[min(full_hist.len()-1,*count as usize)] += 1;
            if *count >= best_count {
                second_best_count = best_count;
                second_best_index = best_index;
                best_count = *count;
                best_index = index;
            } else if *count >= second_best_count {
                second_best_count = *count;
                second_best_index = index;
            }
        }
        if best_count as u32 >= params.min_count && second_best_count as u32 >= params.min_count {
            let kmer = KmerX::from_u64(middle_base_invariant_kmer | masks[best_index]);
            let kmer2 = KmerX::from_u64(middle_base_invariant_kmer | masks[second_best_index]);
            writeln!(&mut out_writer,"{}\t{}\t{}\t{}\t{}",kmer.to_string(), best_count,
                kmer2.to_string(), second_best_count, total-best_count-second_best_count)
                .expect("could not write to file");
        }
        total = 0; best_count = 0; best_index = 0; second_best_count = 0; second_best_index = 0;
    }
    out_writer.flush().expect("cannot flush stream");
    for (index, count) in full_hist.iter().enumerate() {
        writeln!(&mut hist, "{}\t{}",index,count).expect("could not write hist");
    }
}

fn count_kmers_fastq(params: &Params) -> FnvHashMap<u64,[u16;4]> {
    //let middle_base_mask: u64 = !( 3 << ( KX::K() - 1 ) ); // make a mask that is 1's outside the two bits at the center of the kmer 
    let estimated_kmers = params.estimated_kmers/params.modimizer;
    if estimated_kmers > std::u32::MAX as u64 { panic!("can't deal with this many kmers in counting bloom filter"); }
    let estimated_kmers = estimated_kmers as u32;
    let mut counting_bits = 0usize;
    while 2u32.pow(counting_bits as u32) <= params.min_count-1 { counting_bits += 1; }
    let mut filter: CountingBloomFilter = 
        CountingBloomFilter::with_rate(counting_bits, 0.05, estimated_kmers);
    let mut counts: FnvHashMap<u64,[u16;4]> = FnvHashMap::default();
    let mut middle_base_mask: u64 = !( 3 << ( KX::K() - 1) );
    //middle_base_mask = middle_base_mask & KmerX::top_mask( KX::K() );
    let middle_base_only: u64 = 3<<(KX::K()-1);
    let a_mask = 0;
    let c_mask = ( 1 << ( KX::K() - 1 ) );
    let g_mask = ( 2 << ( KX::K() - 1 ) );
    let t_mask = ( 3 << ( KX::K() - 1 ) );
    let masks: [u64; 4] = [ a_mask, c_mask, g_mask, t_mask ];
    
    for kmer_file in &params.input_files {
        let reader = dna_io::DnaReader::from_path(kmer_file);
        for record in reader {
            //'kmerloop: for kmeru64 in from_ascii(&record.seq.as_bytes()) {
            'kmerloop: for k in KmerX::kmers_from_ascii(&record.seq.as_bytes()) {
                //let k = KmerX::from_u64(kmeru64);
                let krc = k.rc();
                let middle_base_invariant = min( k.to_u64() & middle_base_mask, krc.to_u64() & middle_base_mask );
                // below is the line that splits thread work and maintains that ...X... and ...Y... kmers end up in the same thread
                if middle_base_invariant % params.modimizer != params.mod_index { continue; }
                let to_hash = min( k.to_u64(), krc.to_u64() );
                //if counts.contains_key(&middle_base_invariant) {
                match counts.get_mut(&middle_base_invariant) {
                    Some(array) => {
                        array[((middle_base_only & to_hash) >> (KX::K()-1)) as usize] += 1;
                        continue 'kmerloop;
                    },
                    None => (),
                }
                match filter.insert_get_count(&to_hash) {
                    a if a >= params.min_count => {
                        let mut array = [0u16; 4];
                        for (index, base_mask) in masks.iter().enumerate() {
                            let kmer = KmerX::from_u64(middle_base_invariant | base_mask);
                            let to_hash = min(kmer.to_u64(), kmer.rc().to_u64());
                            //println!("double check {} {}",k.to_string(),kmer.to_string());
                            array[index] = filter.estimate_count(&to_hash) as u16;
                        }
                        counts.insert(middle_base_invariant, array);
                    },
                    _ => ()
                }
            }
        }
    }
    counts
}



type KmerX = VarIntKmer<u64, KX>;
#[derive(Debug, Hash, Copy, Clone, Ord, PartialOrd, Eq, PartialEq)]
pub struct KX;

impl KmerSize for KX {
    fn K() -> usize {
        unsafe {
            KMER_SIZE
        }
    }
}

#[derive(Clone)]
struct Params {
    min_count: u32,
    estimated_kmers: u64,
    //counting_bits: usize,
    input_files: Vec<String>,
    output_hist: String,
    output: Option<String>,
    mod_index: u64,
    modimizer: u64,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();
    let mut input_files: Vec<String> = Vec::new();
    for input_file in params.values_of("inputs").unwrap() {
        input_files.push(input_file.to_string());
    }
    let output = match params.value_of("output") {
        Some(x) => Some(x.to_string()),
        None => None,
    };
    let min = params.value_of("min_coverage").unwrap();
    let min: u32 = min.to_string().parse::<u32>().unwrap();
    let kmer_size = params.value_of("kmer_size").unwrap_or("21");
    let kmer_size: usize = kmer_size.to_string().parse::<usize>().unwrap();
    if kmer_size % 2 == 0 {
        panic!("kmer size required to be odd");
    }
    unsafe {
        KMER_SIZE = kmer_size;
        KMER_SIZE_MINUS_ONE = kmer_size-1;
        KMER_MASK = 0u64;
        for i in 0..kmer_size {KMER_MASK = (KMER_MASK <<2) | 3u64; }
    }
    let output_hist = params.value_of("output_full_hist").unwrap_or("none");
    let estimated_kmers = params.value_of("estimated_kmers").unwrap_or("1000000000");
    let estimated_kmers: u64 = estimated_kmers.to_string().parse::<u64>().unwrap();
    //let counting_bits = params.value_of("counting_bits").unwrap_or("7");
    //let counting_bits: usize = counting_bits.to_string().parse::<usize>().unwrap();
    let modimizer = params.value_of("modimizer").unwrap_or("1");
    let modimizer: u64 = modimizer.to_string().parse::<u64>().unwrap();
    let mod_index = params.value_of("mod_remainder").unwrap();
    let mod_index = mod_index.to_string().parse::<u64>().unwrap();
    Params{
        min_count: min,
        //threads: threads,
        estimated_kmers: estimated_kmers,
        //counting_bits: counting_bits,
        input_files: input_files,
        output: output,
        output_hist: output_hist.to_string(),
        modimizer: modimizer,
        mod_index: mod_index,
    }
}

//lets make this ugly
#[inline]
pub fn base_to_u64(c: u8) -> u64 {
    match c {
        b'A' => 0u64,
        b'a' => 0u64,
        b'C' => 1u64,
        b'c' => 1u64,
        b'G' => 2u64,
        b'g' => 2u64,
        b'T' => 3u64,
        b't' => 3u64,
        _ => 0u64,
    }
}

pub fn from_ascii(str: &[u8]) -> Vec<u64> {
    unsafe {
        if str.len() < KMER_SIZE {
            return Vec::new();
        }
        let mut kmers: Vec<u64> = Vec::with_capacity(str.len()-KMER_SIZE);
        let mut k0 = 0;
        for i in 0..KMER_SIZE {
            k0 = (k0 << 2) | base_to_u64(str[KMER_SIZE-i-1]);
        }
        kmers.push(k0.clone());
        for index in KMER_SIZE..str.len() {
            k0 = (k0 << 2) | base_to_u64(str[index]) & KMER_MASK;
            kmers.push(k0.clone());
        }
        return kmers;
    }
}
