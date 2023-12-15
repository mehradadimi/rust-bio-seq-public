fn main() {
    // // /* Pattern Matching */

    // BNDM

    use bio::pattern_matching::bndm;
    let pattern = b"GAAAA";
    let text = b"ACGGCTAGAAAAGGCTAGAAAA";
    let bndm = bndm::BNDM::new(pattern);
    let occ: Vec<usize> = bndm.find_all(text).collect();
    assert_eq!(occ, [7, 17]);



    println!("Occurrences BNDM: {:?}", occ);


    // BOM


    // BOM

    use bio::pattern_matching::bom::BOM;
    let text = b"ACGGCTAGGAAAAAGACTGAGGACTGAAAA";
    let pattern = b"GAAAA";
    let bom = BOM::new(pattern);
    let occ: Vec<usize> = bom.find_all(text).collect();
    assert_eq!(occ, [8, 25]);


    // println!("Occurrences BOM: {:?}", occ);


    // MYERS
    use bio::pattern_matching::myers::Myers;

    let text = b"CGGTCCTGAGGGATTAGCAC";
    let pattern = b"TCCTAGGGC";

    let myers = Myers::<u64>::new(pattern);
    let occ: Vec<_> = myers.find_all_end(text, 2).collect();

    assert_eq!(occ, [(11, 2), (12, 2)]);

    println!("Occurrences Myers: {:?}", occ);
    // Ukkonen

    use bio::pattern_matching::ukkonen::{unit_cost, Ukkonen};

    let mut ukkonen = Ukkonen::with_capacity(10, unit_cost);
    let text = b"ACCGTGGATGAGCGCCATAG";
    let pattern = b"TAGCGC";
    let occ: Vec<(usize, usize)> = ukkonen.find_all_end(pattern, text, 1).collect();
    assert_eq!(occ, [(14, 1)]);

    println!("Occurrences UKKONEN: {:?}", occ);



    extern crate bv;
    use bio::data_structures::rank_select::RankSelect;
    use bv::BitVec;
    use bv::BitsMut;

    let mut bits: BitVec<u8> = BitVec::new_fill(false, 64);
    bits.set_bit(5, true);
    bits.set_bit(32, true);
    let rs = RankSelect::new(bits, 1);

    for i in 0..64 {
        println!("RS Rank_0 at position {}: {:?}", i, rs.rank_0(i).unwrap());
    }

    for i in 0..64 {
        println!("RS Rank_1 at position {}: {:?}", i, rs.rank_1(i).unwrap());
    }

    for i in 0..64 {
        println!("RS select_0 at position {}: {:?}", i, rs.select_0(i));
    }

    for i in 0..64 {
        println!("RS select_1 at position {}: {:?}", i, rs.select_1(i));
    }


    // Small Ints
    use bio::data_structures::smallints::SmallInts;
    let mut smallints: SmallInts<u8, usize> = SmallInts::new();
    smallints.push(3);
    smallints.push(4);
    smallints.push(255);
    smallints.push(305093);
    assert_eq!(smallints.get(0).unwrap(), 3);
    smallints.set(0, 50000);
    let values: Vec<usize> = smallints.iter().collect();
    assert_eq!(values, [50000, 4, 255, 305093]);
    println!("{:?}", values);


    // HMM
    use approx::assert_relative_eq;
    use bio::stats::hmm::discrete_emission::Model as DiscreteEmissionHMM;
    use bio::stats::hmm::viterbi;
    use bio::stats::Prob;
    use ndarray::array;

    let transition = array![[0.5, 0.5], [0.4, 0.6]];
    let observation = array![[0.2, 0.3, 0.3, 0.2], [0.3, 0.2, 0.2, 0.3]];
    let initial = array![0.5, 0.5];

    let hmm = DiscreteEmissionHMM::with_float(&transition, &observation, &initial)
        .expect("Dimensions should be consistent");
    let (_path, log_prob) = viterbi(&hmm, &vec![2, 2, 1, 0, 1, 3, 2, 0, 0]);
    let prob = Prob::from(log_prob);
    assert_relative_eq!(4.25e-8_f64, *prob, epsilon = 1e-9_f64);


        
    // Sparse -> only lcskpp

    use bio::alignment::sparse::*;

    let s1 =   b"ACGTACGATAGGTA";
    let s2 = b"TTACGTACGATAGGTATT";
    let k = 8;
    let matches = find_kmer_matches(s1, s2, k);
    let sparse_al = lcskpp(&matches, k);
    let match_path: Vec<(u32,u32)> = sparse_al.path.iter().map(|i| matches[*i]).collect();
    assert_eq!(match_path, vec![(0,2), (1,3), (2,4), (3,5), (4,6), (5,7), (6,8)]);
    assert_eq!(sparse_al.score, 14);

    println!("{:?}", match_path);
    println!("{}",sparse_al.score );


    // Interval Tree

    use bio::data_structures::interval_tree::IntervalTree;
    use bio::utils::Interval;

    let mut tree = IntervalTree::new();
    tree.insert(11..20, "Range_1");
    tree.insert(25..30, "Range_2");
    for r in tree.find(15..25) {
        assert_eq!(r.interval().start, 11);
        println!("{}", r.interval().start);
        assert_eq!(r.interval().end, 20);
        println!("{}", r.interval().end);
        assert_eq!(r.interval(), &(Interval::from(11..20)));
        println!("{:?}", r.interval());
        assert_eq!(r.data(), &"Range_1");
        println!("{}", r.data());
    }   



    // Wavelet matrix




    // Wavelet Matrix data structure for DNA alphabet.



    use bio::data_structures::wavelet_matrix::WaveletMatrix;


    let text = b"AANGGT$ACCNTT$";
    let wm = WaveletMatrix::new(text);  

    println!("Rank of 'A' at position 0: {}", wm.rank(b'A', 0)); // Expected: 1
    println!("Rank of 'G' at position 9: {}", wm.rank(b'G', 9)); // Expected: 2
    println!("Rank of 'T' at position 13: {}", wm.rank(b'T', 13)); // Expected: 3

    for i in 0..64 {
        println!("RS select_1 at position {}: {:?}", i, rs.select_1(i));

    }

   
}
