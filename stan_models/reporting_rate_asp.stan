data {
//Observation numbers
int N; // number of obs
int P; //number of points
int C; // number of canyons
int K; //number of obs level predictors
int J; //number of point level predictors
int R; // number of ranges

//Cluster IDs
int <lower=1, upper=P> PtID[N]; // map obs to point
int <lower=1, upper=C> CanID[N]; // map obs to canyon
int <lower=1, upper=R> RgID[N]; // map obs to range

//Lookup IDs
int <lower=1, upper=C> PtInCan[P]; // map point to can
int <lower=1, upper=R> CanInRg[C]; // map can to range

//outcome
int y[N];

//binomial offset
int nsamp[N]; // number of samples per point

//predictors
matrix[P,J] W; //point level predictors
matrix[N,K] X; //obs level predictors
}
parameters{
real mu_alpha;
vector[P] alpha_pt_tilde;
vector[C] alpha_can_tilde;
vector[R] alpha_rg_tilde;

real<lower=0> sigma_alpha;
vector<lower=0>[C] sigma_can;
vector<lower=0>[R] sigma_rg;

vector[J] beta_p;
vector[K] beta_o;
}
transformed parameters{
vector[R] alpha_rg = mu_alpha + sigma_alpha * alpha_rg_tilde;
vector[C] alpha_can = alpha_rg[CanInRg] + sigma_rg[CanInRg] .* alpha_can_tilde[CanInRg];
vector[P] alpha_point = alpha_can[PtInCan] + W * beta_p + sigma_can[PtInCan] .* alpha_pt_tilde[PtInCan];
}
model{
mu_alpha ~ normal(0, 2);
alpha_pt_tilde ~ normal(0, 1);
alpha_can_tilde ~ normal(0, 1);
alpha_rg_tilde ~ normal(0, 1);

sigma_alpha ~ normal(0, 1);
sigma_can ~ normal(0, 1);
sigma_rg ~ normal(0, 1);


beta_p ~ normal(0, 1);
beta_o ~ normal(0, 1);

y ~ binomial_logit(nsamp, X * beta_o + alpha_point[PtID]);
}
generated quantities{
vector[N] y_tilde;
vector[N] obs_dev;
vector[N] obs_dev2;
vector[N] sim_dev;
vector[N] sim_dev2;
vector[N] log_lik;

for (n in 1:N){
y_tilde[n] = binomial_rng(nsamp[n],inv_logit(X[n, 1:K] * beta_o + alpha_point[PtID[n]])); // generate a simulated outcome based on the number of samples and posterior values of p
obs_dev[n] = y[n] - (nsamp[n] * inv_logit(X[n,1:K] * beta_o + alpha_point[PtID[n]])); //estimate the deviation from observed
obs_dev2[n] = obs_dev[n]*obs_dev[n];
sim_dev[n] = y_tilde[n] - (nsamp[n] * inv_logit(X[n,1:K] * beta_o + alpha_point[PtID[n]])); //estimate the deviation for the simulated outcome
sim_dev2[n] = sim_dev[n]*sim_dev[n];
log_lik[n] = binomial_logit_lpmf(y[n] | nsamp[n], X[n, 1:K] * beta_o + alpha_point[PtID[n]]);
  }
}
