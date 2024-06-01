trait Likelihood {

    type Observation;

    // Returns p(x_i | z_i = k, theta)
    // x_i is the i-th observation.
    // z_i is the i-th latent variable. z_i \in {0,...,k-1}
    fn likelihood(&self, x_i: &Self::Observation, k: usize) -> f64;

}

// initial_theta is the initial guess for the mixing fractions.
// \theta_1 + ... + \theta_k = 1.
fn fit_model<L: Likelihood>(likelihood: &L, observations: &Vec<L::Observation>, initital_theta: &Vec<f64>) -> Vec<f64>{
    let n = observations.len();
    let K = initital_theta.len();

    let mut prev_theta = initital_theta.clone();
    loop {
        // Compute latent variable posteriors given previous theta estimate
        let mut Z_posteriors: Vec<Vec<f64>> = vec![vec![0.0; K]; n];
        for i in 0..n {
            let mut denominator: f64 = 0.0;
            for w in 0..K {
                denominator += prev_theta[w] * likelihood.likelihood(&observations[i], w);
            }
            for k in 0..K {
                Z_posteriors[i][k] = prev_theta[k] * likelihood.likelihood(&observations[i], k) / denominator;
            }
        }

        // Estimate a new theta that maximizes the expected likelihood assuming the latent variables are distributed
        // according to their posteriors that were just computed.

        let mut next_theta: Vec<f64> = vec![0.0; K];
        for k in 0..K {
            next_theta[k] = (0..n).fold(0.0, |acc, i| acc + Z_posteriors[i][k]) / n as f64;
        }

        eprintln!("{:?}", next_theta);

        // change = |prev_theta - next_theta|
        let change = (0..K).fold(0.0, |acc, k| acc + (prev_theta[k] - next_theta[k]) * (prev_theta[k] - next_theta[k])).sqrt();

        if change < 1e-9 {
            return next_theta;
        }

        prev_theta = next_theta;

    }
}


#[cfg(test)]
mod tests {
    use super::*;

    struct MyLikelihood {
        true_origins: Vec<usize> 
    }

    impl Likelihood for MyLikelihood {
        type Observation = usize;

        fn likelihood(&self, x_i: &Self::Observation, k: usize) -> f64 {
            match *x_i == k {
                true => 0.99,
                false => 0.01,
            }
        }
    }

    #[test]
    fn test_EM_algo(){
        let observations = vec![0,0,1,2,3,1,2,3,3,2,1,1,2,3,3,1,3,3,3,3,3,1,1,1,1,1,1,1];
        let likelihood = MyLikelihood{true_origins: observations.clone()};
        let initial_theta: Vec<f64> = vec![0.25, 0.25, 0.25, 0.25];
        fit_model(&likelihood, &observations, &initial_theta);
    }
}