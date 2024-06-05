pub trait Likelihood {

    type Observation;

    // Returns p(x_i | z_i = k, theta)
    // x_i is the i-th observation.
    // z_i is the i-th latent variable. z_i \in {0,...,k-1}
    fn likelihood(&self, x_i: &Self::Observation, k: usize) -> f64;

}

// initial_theta is the initial guess for the mixing fractions.
// \theta_1 + ... + \theta_k = 1.
pub fn fit_model<L: Likelihood>(likelihood: &L, observations: &Vec<L::Observation>, initital_theta: &Vec<f64>) -> Vec<f64>{
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

        // Compute the current log-likelihood (just FYI for the user, does not affect the algorithm)
        let mut total_log_likelihood: f64 = 0.0;
        for i in 0..n {
            let prob = (0..K).fold(0.0, |acc, k| acc + next_theta[k] * likelihood.likelihood(&observations[i], k));
            total_log_likelihood += prob.ln();
        }

        log::info!("{}", total_log_likelihood);

        // change = |prev_theta - next_theta|
        let change = (0..K).fold(0.0, |acc, k| acc + (prev_theta[k] - next_theta[k]) * (prev_theta[k] - next_theta[k])).sqrt();
        if change < 1e-9 {
            return next_theta; // Converged
        }

        prev_theta = next_theta;

    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use rand::{rngs::ThreadRng, Rng};
    use rand_distr::{Dirichlet, Distribution};

    struct MyLikelihood {}

    impl Likelihood for MyLikelihood {
        type Observation = usize;

        fn likelihood(&self, x_i: &Self::Observation, k: usize) -> f64 {
            match *x_i == k {
                true => 0.99,
                false => 0.01,
            }
        }
    }

    struct LikelihoodMatrix{
        likelihoods: Vec<Vec<f64>>
    }

    impl Likelihood for LikelihoodMatrix {

        type Observation = usize; // Index of read

        fn likelihood(&self, x_i: &Self::Observation, k: usize) -> f64 {
            self.likelihoods[*x_i][k]
        }
    }

    #[test]
    fn test_EM_algo(){
        let observations = vec![0,0,1,2,3,1,2,3,3,2,1,1,2,3,3,1,3,3,3,3,3,1,1,1,1,1,1,1];
        let likelihood = MyLikelihood{};
        let initial_theta: Vec<f64> = vec![0.25, 0.25, 0.25, 0.25];

        eprintln!("{:?}", (0..=3).map(|k| observations.iter().filter(|x| **x == k).count()));

        fit_model(&likelihood, &observations, &initial_theta);
    }

    #[test]
    fn simulation_test(){
        // Start with some ground truth mixing ratios. Generate observations and a likelihood
        // matrix for them so that the first mixing component tends to have a high likelihood
        // and the others a low likelihood.

        let mut rng = rand::thread_rng();
        let n = 10000;
        let K = 5;
        let mut L: Vec<Vec<f64>> = vec![vec![0.0; K]; n]; // Likelihood matrix
        let initial_theta: Vec<f64> = vec![0.2, 0.2, 0.2, 0.2, 0.2];
        let true_theta: Vec<f64> = vec![0.1, 0.2, 0.3, 0.4, 0.0];
        let cat_theta = rand_distr::WeightedIndex::new(&true_theta).unwrap(); // Categorical distribution with weights theta
        let dirichlets = vec![ // The log likelihoods will be samples from one of these
            Dirichlet::new(&[10.0, 1.0, 1.0, 1.0, 1.0]).unwrap(),
            Dirichlet::new(&[1.0, 10.0, 1.0, 1.0, 1.0]).unwrap(),
            Dirichlet::new(&[1.0, 1.0, 10.0, 1.0, 1.0]).unwrap(),
            Dirichlet::new(&[1.0, 1.0, 1.0, 10.0, 1.0]).unwrap(),
            Dirichlet::new(&[1.0, 1.0, 1.0, 1.0, 10.0]).unwrap(),
        ];
        for i in 0..n {
            let k = cat_theta.sample(&mut rng);
            let likelihoods = dirichlets[k].sample(&mut rng);
            L[i] = likelihoods;
        }

        let M = LikelihoodMatrix{likelihoods: L};
        let observations = (0..n).collect::<Vec::<usize>>();
        fit_model(&M, &observations, &initial_theta);
    }
}