data {
  //Observation numbers
  int N; // number of obs
  int P; //number of points
  int C; // number of canyons
  int K; //number of observation predictors
  int J; //number of point predictors
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


  //predictors
  matrix[P,J] W;
  vector[N] X;
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
real beta_o;
}
transformed parameters{
  vector[R] alpha_rg = mu_alpha + sigma_alpha * alpha_rg_tilde;
  vector[C] alpha_can = alpha_rg[CanInRg] + sigma_rg[CanInRg] .* alpha_can_tilde[CanInRg];
  vector[P] alpha_point = alpha_can[PtInCan] + W * beta_p + sigma_can[PtInCan] .* alpha_pt_tilde[PtInCan];
  vector[N] theta = X * beta_o + alpha_point[PtID];
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

y ~ bernoulli_logit(theta);
}
generated quantities{
  vector[N] y_tilde;
for (n in 1:N){
  y_tilde[n] = bernoulli_logit_rng(X[n] * beta_o + alpha_point[PtID[n]]); // generate a simulated outcome based on the number of samples and posterior values of p
  }
}
