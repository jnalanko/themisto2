// The input is a weighted compatibility matrix where cell (i,j) has the number of
// compatible genomes in cluster j for read i.
// have size 1, the hit counts are 0 or 1.
// The latent variables (I_n in the paper) are color ids, one for each read.
fn evalulate_model(hit_counts: &Vec<Vec<usize>>, latent_variables: &Vec<usize>, mixing_proportions: &Vec<f64>) -> f64 {
    // Assuming a flat prior on the mixing proportions theta
    // p(r_{n, k} | I_n = k) * p(I_n = k | theta)
    // Is it the case that p(I_n = k | theta) is just theta_n?

    // The likelihood term p(r_{n, k} | I_n = k)...
    // Ordinarily, without clustered references, this would be... I guess...:
    // 0.01 if r_{n,k} = 0, else 0.99
    //
    // But now, with a clustered reference, we have a refinement: Let M_k be the number of genomes in cluster k.
    // Then, we have 0.01 if r_{n,k} as above, but otherwise we have 0.99 * f(r_{n,k}, M_k). The thesis says that
    // a beta-binomial distribution for f(...) would be a reasonable first approximation. It has two hyperparameters
    // alpha and beta. But, then f would not be increasing as a function of hit count, so they define a modified
    // version of the beta-binomial distribution to fix this. There is a lot of math involved in finding the right
    // normalization constant.

    let term1 = |_,k| mixing_proportions[k];
    let term2 = |n: usize, k: usize| {
        if hit_counts[n][k] == 0 {
            0.01
        } else if hit_counts[n][k] == 1 {
            0.99
        } else {
            0.99 // todo: here the beta-binomial thing
        }
    };

    let N = latent_variables.len(); // Number of reads
    let k = hit_counts.first().unwrap().len(); // Number of clusters
    let mut ans = 1_f64;
    for n in 0..N {
        for k in 0..k {
            ans *= term1(n, k)* term2(n,k);
        }
    }
    ans
}

fn main() {
    println!("Hello, world!");
}
