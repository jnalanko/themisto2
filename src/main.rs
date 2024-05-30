#![allow(non_snake_case, clippy::needless_range_loop)] // Using upper-case variable names from the source material

use std::thread::LocalKey;

use statrs::function::gamma::ln_gamma;
use statrs::function::beta::ln_beta;

fn ln_binom(n: usize, k: usize) -> f64 {
    assert!(k <= n); // Otherwise would be minus infinity
    let n = n as f64;
    let k = k as f64;
    ln_gamma(n+1.0) - ln_gamma(k+1.0) - ln_gamma(n-k+1.0)
}

fn ln_beta_binom(k: usize, n: usize, alpha: f64, beta: f64) -> f64 {
    ln_binom(n,k) + ln_beta(k as f64 + alpha, (n - k) as f64 + beta) - ln_beta(alpha, beta)
}

// p^* in Eq. 2.7 in Tommi's thesis
fn ln_pstar_unnormalized(k: usize, n: usize, alpha: f64, beta: f64) -> f64 {
    ln_binom(n,k) + ln_beta(k as f64 + alpha, (n - k) as f64 + beta) - ln_beta(n as f64 + alpha, beta)
}

// Returns p(theta | R)
// Likelihood_matrix[i][j] is the likelihood of read j vs cluster i
// Flat prior on theta
fn theta_posterior(likelihood_matrix: &Vec<Vec<f64>>, theta: &Vec<f64>) -> f64 {
    let k = likelihood_matrix.len();
    let n = likelihood_matrix.first().unwrap().len();
    let mut ans = 0_f64;
    for read in 0..n {
        let likelihood = (0..k).map(|k| likelihood_matrix[k][read] * theta[k]).sum::<f64>();
        ans += likelihood.ln(); // Log-likelihood
    }
    ans
}

// Returns p(theta | R)
// Likelihood_matrix[i][j] is the likelihood of read j vs cluster i
// Flat prior on theta
fn theta_posterior_exponential_formula(likelihood_matrix: &Vec<Vec<f64>>, theta: &Vec<f64>) -> f64 {
    let k = likelihood_matrix.len();
    let n = likelihood_matrix.first().unwrap().len();
    let mut ans = 0_f64;
    
    // Iterate all tuples of length n with values in [0..k).
    for x in 0..k.pow(n as u32) { // Interpreting x as a base-k number
        let mut I = Vec::<usize>::new();
        let mut x_copy = x;
        for _ in 0..n {
            I.push(x_copy % k);
            x_copy /= k;
        }
        let mut I_likelihood = 1_f64;
        for (read, &k) in I.iter().enumerate() {
            I_likelihood *= likelihood_matrix[k][read] * theta[k];
        }
        ans += I_likelihood;
    }

    ans.ln()
}

// Takes in R[0..K)[0..N)
// R[k][n] = number of pseudoalignments from read n to color cluster k
// cluster_sizes[k] = number of colors in cluster k
// Returns a matrix L[0..K)[0..N) such that L[k][n] is the likelihood that read r_n
// is from cluster k.
fn build_likelihood_matrix(R: &Vec<Vec<usize>>, cluster_sizes: &Vec<usize>) -> Vec<Vec<f64>> {
    let K = R.len();
    let N = R.first().unwrap().len();
    let mut L = vec![vec![0_f64; N]; K];

    let pi: f64 = 0.65; // From thesis

    for n in 0..N {
        for k in 0..K {
            L[k][n] = if R[k][n] == 0 {
                0.01
            } else if R[k][n] == 1 && cluster_sizes[k] == 1 {
                0.99
            } else {
                let phi = 1.0 - pi + 0.01 / cluster_sizes[k] as f64;
                let alpha = pi / phi;
                let beta = (1.0 - pi) / phi;
                dbg!(pi, phi, alpha, beta, k, n);
                0.99 * ln_pstar_unnormalized(R[k][n], cluster_sizes[k], alpha, beta).exp()
            }
        }
    }

    L
}

fn main() {

    let cluster_sizes = vec![3,1,4]; // Max number of hits to each cluster
    let R = vec![vec![3,2,1,0], vec![1,0,1,0], vec![0,2,4,4]]; // 4 reads, 3 clusters
    let likelihood_matrix = build_likelihood_matrix(&R, &cluster_sizes);
    let theta = vec![0.3, 0.3, 0.4];
    eprintln!("{}", theta_posterior(&likelihood_matrix, &theta));
    //eprintln!("{}", 

    /*
    let n = 10;
    let pi = 0.65;
    let phi = 1.0 - pi + 0.01 / 10.0;
    let alpha = pi / phi;
    let beta = (1.0 - pi) / phi;
    dbg!(alpha, beta);

    //let alpha = 3_f64;
    //let beta = 2_f64;
    assert!(alpha / (alpha + beta) >= 0.5);
    let mut sum = 0_f64;
    let mut sum2 = 0_f64;
    for k in 0..=n {
        eprintln!("{} {} {}", k, ln_beta_binom(k, n, alpha, beta).exp(), ln_pstar_normalized_thesis_formula(k, n, alpha, beta).exp());
        sum += ln_beta_binom(k, n, alpha, beta).exp();
        sum2 += ln_pstar_normalized_thesis_formula(k, n, alpha, beta).exp();
    }
    eprintln!("{} {}", sum, sum2);
    */
}

#[cfg(test)]
mod tests {
    use statrs::assert_almost_eq;

    use crate::ln_binom;
    #[test]
    fn test_ln_binom(){
        let n = 10;
        let k = 4;
        let correct = 210; // n choose k
        let computed = ln_binom(n,k).exp();
        assert!((correct as f64 - computed).abs() < 1e-6);
    }
}
