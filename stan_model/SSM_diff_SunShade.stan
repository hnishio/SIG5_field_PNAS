
data {
  int N_time;
  int N1[N_time];
  int N2[N_time];
  int sumN1;
  int sumN2;
  vector[sumN1] Y1;
  vector[sumN2] Y2;
}

parameters {
  vector[N_time] mu1;
  real diff0;
  vector<lower=-pi()/2, upper=pi()/2>[N_time-1] diff_unif;   // Uniform distribution
  real<lower=0> s_mu1;
  real<lower=0> s_diff;
  real<lower=0> s_Y;
}

transformed parameters {
  vector[N_time] diff;
  vector[N_time] mu2;
  diff[1] = diff0;
  for (t in 2:N_time) {
    diff[t] = diff[t-1] + s_diff*tan(diff_unif[t-1]);    // diff follows cauchy distribution
  }
  mu2 = mu1 + diff;
}

model {
  // State equation of mu1
  for (t in 3:N_time) {
    mu1[t] ~ normal(2*mu1[t-1] - mu1[t-2], s_mu1);
  }
  // Observation equation of Y1
  for (n in 1:N1[1]) {                   // N1=c(6,6,6,...)
    Y1[n] ~ normal(mu1[1], s_Y);
  }
  for (t in 2:N_time) {
    for (n in 1:N1[t]) {
      Y1[sum(N1[1:(t-1)])+n] ~ normal(mu1[t], s_Y);
    }
  }
  // Observation equation of Y2
  for (n in 1:N2[1]) {                   // N2=c(6,6,6,...)
    Y2[n] ~ normal(mu2[1], s_Y);
  }
  for (t in 2:N_time) {
    for (n in 1:N2[t]) {
      Y2[sum(N2[1:(t-1)])+n] ~ normal(mu2[t], s_Y);
    }
  }
}
