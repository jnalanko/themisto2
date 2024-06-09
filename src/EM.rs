pub trait Likelihood {

    type Observation;

    // Returns p(x_i | z_i = k, theta)
    // x_i is the i-th observation.
    // z_i is the i-th latent variable. z_i \in {0,...,k-1}
    fn likelihood(&self, x_i: &Self::Observation, k: usize) -> f64;

}

fn compute_theta_contributions<O: Sync, L: Likelihood<Observation = O> + Sync + Send>(likelihood: &L, observations: &[L::Observation], observation_counts: &[usize], prev_theta: &[f64]) -> Vec<f64> {
    assert_eq!(observations.len(), observation_counts.len());
    let K = prev_theta.len();

    // Compute latent variable posteriors for each distinct observation given previous theta estimate
    let mut next_theta: Vec<f64> = vec![0.0; K];
    for i in 0..observations.len() {
        let mut denominator: f64 = 0.0; // Normalization factor
        for w in 0..K {
            denominator += prev_theta[w] * likelihood.likelihood(&observations[i], w);
        }
        for k in 0..K {
            let Z_i_k_posterior = prev_theta[k] * likelihood.likelihood(&observations[i], k) / denominator;

            // This contribution should be divided by n_total_observations,
            // but we'll do it in the end for numerical reasons.
            next_theta[k] += observation_counts[i] as f64 * Z_i_k_posterior;
        }
    }
    next_theta
}

// likelihood: see the comments in the trait for an explanation for what it is.
// initial_theta is the initial guess for the mixing fractions.
// \theta_1 + ... + \theta_K = 1.
// observation_counts[i] = number of times observation i was observed
// It should hold that observations.len() == observation_counts.len().
pub fn fit_model<O: Sync, L: Likelihood<Observation = O> + Sync + Send>(likelihood: &L, observations: &[L::Observation], observation_counts: &[usize], initital_theta: &[f64], n_threads: usize) -> Vec<f64>{

    // There is one latent variable per observation. Let's denote the i-th latent variable with Z_i. Each latent variable
    // is assigned a value from the set {0..K-1}, where K is the number of mixture components. The interpretation is that
    // if Z_i = j, then the i-th observation comes from the j-th mixture component.
    //
    // E-step computes the posterior distributions for each latent variable, given the previous estimate for theta.
    // That is, for each Z_i, we compute p(Z_i = k | theta, observations) for all k = 0..K-1.
    // Assuming a flat prior on theta, The formula is:
    // p(Z_i = k | theta, observations) = theta_k * p(observation i | Z_i = k) / N(i)
    // where the likelihood p(observation i | Z_i = k) is queried from the given likelihood model, and N(i)
    // is the normalization factor for observation i to make p(Z_i = k) sum to 1 over k = 0..K-1.
    //
    // In the M-step, we find the theta that maximizes the likelihood of the data, given that the latent variables
    // are distributed according to the posteriors computed in the E-step. There is a closed-form formula for this,
    // which can be derived using lagrange multipliers. The formula just boils down to this: the k-th component of the
    // optimal theta is the expected fraction of all observations assigned to component k. That is, the sum of
    // p(Z_i = k) over all observations i, divided by the number of observations.
    //
    // Since the model is this simple, we can merge the E-step and M-step and do them at the same time.
    // We also compute the contributions to the next theta in parallel by dividing the work into blocks
    // and deferring a division that applies to all contributions until the very end for numerical reasons.
    // The code becomes less readable (sorry about that), but it's fast and parallelized very well.

    let n_distinct_observations = observations.len();
    let n_total_observations: usize = observation_counts.iter().sum();
    let K = initital_theta.len();

    let mut prev_theta = initital_theta.to_owned();
    let slice_len = (n_distinct_observations + n_threads - 1) / n_threads; // ceil n_distinct_observations / n_threads

    loop {

        // Compute the contributions to the next theta for each distinct observation given previous theta estimate.
        // The work is split evenly to n_threads threads.
        let mut next_theta = std::thread::scope(|s| {
            let mut join_handles = Vec::<_>::new();
            for t in 0..n_threads {
                let start = t*slice_len;
                let end = std::cmp::min((t+1)*slice_len, observations.len());
                let ob_slice = &observations[start..end];
                let ob_count_slice = &observation_counts[start..end];
                join_handles.push(s.spawn(|| {
                    compute_theta_contributions(likelihood, ob_slice, ob_count_slice, &prev_theta)
                }));
            }

            // Add up contributions from all threads
            let mut next_theta: Vec<f64> = vec![0.0; K];
            for h in join_handles {
                for (i, theta_i) in h.join().unwrap().iter().enumerate() {
                    next_theta[i] += theta_i;
                }
            }
            next_theta
        });

        for k in 0..K {
            // This divisions should have been done to all contributions, but we do all those divisions
            // here at once for numerical reasons.
            next_theta[k] /= n_total_observations as f64;
        }

        // Compute the current log-likelihood (just FYI for the user, does not affect the algorithm)
        /*
        let mut total_log_likelihood: f64 = 0.0;
        for i in 0..n_distinct_observations {
            let prob = (0..K).fold(0.0, |acc, k| acc + next_theta[k] * likelihood.likelihood(&observations[i], k));
            total_log_likelihood += prob.ln() * observation_counts[i] as f64;
        }
        */

        let change = (0..K).fold(0.0, |acc, k| acc + (prev_theta[k] - next_theta[k]) * (prev_theta[k] - next_theta[k])).sqrt();

        log::info!("{}", change);

        // change = |prev_theta - next_theta|
        if change < 1e-5 {
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
        let observation_counts = vec![1; observations.len()]; 
        let likelihood = MyLikelihood{};
        let initial_theta: Vec<f64> = vec![0.25, 0.25, 0.25, 0.25];

        eprintln!("{:?}", (0..=3).map(|k| observations.iter().filter(|x| **x == k).count()));

        fit_model(&likelihood, &observations, &observation_counts, &initial_theta, 2);
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
        let observation_counts = vec![1; observations.len()];
        let estimated_theta = fit_model(&M, &observations, &observation_counts, &initial_theta, 2);
        eprintln!("{:?}", estimated_theta);
    }
}