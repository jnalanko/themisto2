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

fn pstar_normalization_factor(n: usize, alpha: f64, beta: f64) -> f64 {
    let sum: f64 = (0..=n).map(|k| ln_pstar_unnormalized(k,n,alpha,beta)).map(|x| x.exp()).sum();
    return 1.0 / sum;
}

// p^* in Eq. 2.7 in Tommi's thesis
fn ln_pstar_unnormalized(k: usize, n: usize, alpha: f64, beta: f64) -> f64 {
    ln_binom(n,k) + ln_beta(k as f64 + alpha, (n - k) as f64 + beta) - ln_beta(n as f64 + alpha, beta)
}

fn ln_pstar(k: usize, n: usize, alpha: f64, beta: f64) -> f64 {
    ln_pstar_unnormalized(k, n, alpha, beta) + pstar_normalization_factor(n, alpha, beta).ln()
}

fn ln_pstar_normalized_thesis_formula(k: usize, n: usize, alpha: f64, beta: f64) -> f64 {
    let mut ans = ln_pstar_unnormalized(k, n, alpha, beta);
    for j in 1..=n {
        let n = n as f64;
        let k = k as f64;
        ans += (alpha + n + k - j as f64).ln() - (alpha + beta + 2.0*n - j as f64).ln();
    } 
    ans
}

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

// Likelihood_matrix[i][j] is the likelihood of read i vs cluster j
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

// Likelihood_matrix[i][j] is the likelihood of read i vs cluster j
// Flat prior on theta
fn theta_posterior_with_indicator_variables(likelihood_matrix: &Vec<Vec<f64>>, theta: &Vec<f64>) -> f64 {
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

fn main() {

    let likelihood_matrix = vec![vec![0.1, 0.1], vec![0.6, 0.7], vec![0.5, 0.5]];
    let theta = vec![0.0, 0.5, 0.5];
    eprintln!("{}", theta_posterior(&likelihood_matrix, &theta));
    eprintln!("{}", theta_posterior_with_indicator_variables(&likelihood_matrix, &theta));

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
