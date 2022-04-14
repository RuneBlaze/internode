mod internode;
use internode::*;
use ndarray::prelude::*;

fn main() {
    // let tree = parse_newick(&mut TaxonSet::new(),"((A,B)0.58:0.5,(C,D)0.45:0.5);");
    // println!("{:?}", tree);
    // run_fastme()
    let mut ts = TaxonSet::new();
    ts.request("A".to_string());
    ts.request("B".to_string());
    ts.request("C".to_string());
    ts.request("D".to_string());
    // let dm = array![
    //     [0.0, 2.0, 3.0, 3.0],
    //     [2.0, 0.0, 3.0, 3.0],
    //     [3.0, 3.0, 0.0, 2.0],
    //     [3.0, 3.0, 2.0, 0.0],
    // ];
    let dm = array![
        [0.0, 3.0, 2.0, 3.0],
        [3.0, 0.0, 3.0, 2.0],
        [2.0, 3.0, 0.0, 3.0],
        [3.0, 2.0, 3.0, 0.0],
    ];
    for i in 0..4 {
        for j in 0..4 {
            if (dm[[i, j]] != dm[[j, i]]) {
                println!("non-symmetric matrix");
            }
        }
    }
    let res = run_fastme(&ts, &dm);
    println!("{}", res);
}