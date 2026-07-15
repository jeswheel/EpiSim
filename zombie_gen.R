# ==============================================================================
# Simulation-Based Inference Workshop: The Unobserved Origin
# Instructor Code: S-E-Z-R with reulermultinom and standard rate nomenclature
# ==============================================================================

library(pomp)

# ------------------------------------------------------------------------------
# 1. DEFINE C-SNIPPETS
# ------------------------------------------------------------------------------

rinit_snippet <- Csnippet("
  S = S_0;
  E = E_0;
  Z = 0;
  R_imm = 0; 
  R = 0;     
  // Z_inc = 0;
")

rprocess_snippet <- Csnippet("
  // Active population denominator
  double N = S + E + Z + R_imm; 
  
  // Calculate noise multipliers (dW / dt)
  // If sigma > 0, we draw Gamma noise. Otherwise, multiplier is 1.0 (deterministic)
  double noise_trans = (sigma_sq > 0.0) ? (rgamma(dt / sigma_sq, sigma_sq) / dt) : 1.0;
  double noise_dem   = (sigma_dem > 0.0) ? (rgamma(dt / sigma_dem, sigma_dem) / dt) : 1.0;

  // ---------------------------------------------------------------------------
  // S COMPARTMENT (Competing risks: Infection vs Natural Death)
  // ---------------------------------------------------------------------------
  double rate_S[2], trans_S[2];
  
  // Base rates
  double mu_SE = (beta_E * E / N) + (beta_Z * Z / N);
  double mu_SD = mu_d;
  
  // Apply environmental noise to rates
  rate_S[0] = mu_SE * noise_trans;
  rate_S[1] = mu_SD * noise_dem;
  
  // reulermultinom(number of exit routes, compartment size, rates, dt, output transitions)
  reulermultinom(2, S, &rate_S[0], dt, &trans_S[0]);
  
  double trans_SE = trans_S[0];
  double trans_S_death = trans_S[1];

  // ---------------------------------------------------------------------------
  // E COMPARTMENT (Competing risks: Zombie vs Immune vs Natural Death)
  // ---------------------------------------------------------------------------
  double rate_E[3], trans_E[3];
  
  // The mutation unlocks at day 60
  double mu_E_out = (t < -7 - gamma) ? 0.0 : gamma;
  
  double mu_EZ = mu_E_out * (1.0 - omega);
  double mu_ERimm = mu_E_out * omega;
  double mu_ED = mu_d;

  rate_E[0] = mu_EZ;                 // Progression to Zombie
  rate_E[1] = mu_ERimm;              // Progression to Immune
  rate_E[2] = mu_ED * noise_dem;     // Natural Death
  
  reulermultinom(3, E, &rate_E[0], dt, &trans_E[0]);
  
  double trans_EZ = trans_E[0];
  double trans_E_immune = trans_E[1];
  double trans_E_death = trans_E[2];

  // ---------------------------------------------------------------------------
  // SINGLE EXIT COMPARTMENTS (Z and R_imm)
  // ---------------------------------------------------------------------------
  
  // Zombie removal via combat
  double mu_ZR = alpha * (S + R_imm);
  double trans_ZR = rbinom(Z, 1 - exp(-mu_ZR * dt));

  // Immune natural death
  double mu_RimmD = mu_d * noise_dem;
  double trans_Rimm_death = rbinom(R_imm, 1 - exp(-mu_RimmD * dt));

  // ---------------------------------------------------------------------------
  // BIRTHS & STATE UPDATES
  // ---------------------------------------------------------------------------
  double expected_birth_rate = mu_b * (S + E + R_imm) * noise_dem;
  double births = rpois(expected_birth_rate * dt);

  // Updates
  S += births - trans_SE - trans_S_death;
  E += trans_SE - trans_EZ - trans_E_immune - trans_E_death;
  Z += trans_EZ - trans_ZR;
  R_imm += trans_E_immune - trans_Rimm_death;
  
  // R holds dead bodies and destroyed zombies
  R += trans_ZR + trans_S_death + trans_E_death + trans_Rimm_death;

  // Z_inc += trans_EZ;
")

rmeasure_snippet <- Csnippet("
  double rho_z = rho_0 * exp(-k * Z);
  double mean_cases = rho_z * Z; // Now taking a fraction of total Z

  if (mean_cases > 0.0) {
    cases = rnbinom_mu(1.0 / tau, mean_cases);
  } else {
    cases = 0.0;
  }
")

dmeasure_snippet <- Csnippet("
  double rho_z = rho_0 * exp(-k * Z);
  double mean_cases = rho_z * Z; // Now evaluating against total Z

  if (mean_cases > 0.0) {
    lik = dnbinom_mu(cases, 1.0 / tau, mean_cases, give_log);
  } else {
    lik = (cases == 0.0) ? (give_log ? 0.0 : 1.0) : (give_log ? R_NegInf : 0.0);
  }
")

# rmeasure_snippet <- Csnippet("
#   double rho_z = rho_0 * exp(-k * Z);
#   double mean_cases = rho_z * Z_inc;
# 
#   if (mean_cases > 0.0) {
#     cases = rnbinom_mu(1.0 / tau, mean_cases);
#   } else {
#     cases = 0.0;
#   }
# ")

# dmeasure_snippet <- Csnippet("
#   double rho_z = rho_0 * exp(-k * Z);
#   double mean_cases = rho_z * Z_inc;
# 
#   if (mean_cases > 0.0) {
#     lik = dnbinom_mu(cases, 1.0 / tau, mean_cases, give_log);
#   } else {
#     lik = (cases == 0.0) ? (give_log ? 0.0 : 1.0) : (give_log ? R_NegInf : 0.0);
#   }
# ")

pt_snippet <- parameter_trans(
  log = c("beta_E", "beta_Z", "alpha", "sigma_sq", "sigma_dem", "tau", "mu_b", "mu_d"),
  logit = c("omega") 
)

# ------------------------------------------------------------------------------
# 2. CALIBRATION & POMP OBJECT SETUP
# ------------------------------------------------------------------------------

true_params <- c(
  beta_E = 0.03,       
  beta_Z = 0.025,      
  gamma = 0.5,
  omega = 0.05,        
  alpha = 0.00000185,  
  sigma_sq = 0.05,     
  sigma_dem = 0.1,     
  rho_0 = 0.95,         
  k = 0.0001,          
  tau = 0.1,           
  S_0 = 10000,        
  E_0 = 5,
  mu_b = 0.0000400,    
  mu_d = 0.0000342     
)

# Weekly reports starting at Day 67
weekly_times <- seq(0, 4*364-60, by = 7)
dummy_data <- data.frame(time = weekly_times, cases = NA)

zombie_pomp <- pomp(
  data = dummy_data,
  times = "time",
  t0 = -60, 
  rinit = rinit_snippet,
  rprocess = euler(rprocess_snippet, delta.t = 0.1), 
  rmeasure = rmeasure_snippet,
  dmeasure = dmeasure_snippet,
  partrans = pt_snippet,
  # accumvars = c("Z_inc"),
  statenames = c("S", "E", "Z", "R_imm", "R"),
  paramnames = names(true_params)
)

# ------------------------------------------------------------------------------
# 3. SIMULATE THE OUTBREAK
# ------------------------------------------------------------------------------

set.seed(42)
sim_data <- simulate(
  zombie_pomp, params = true_params, include.data = FALSE, format = 'data.frame'
)

write.csv(sim_data[1:78, c("time", "cases")], file = '')

# plot(sim_data)


# Check the data visually
plot(sim_data$time[1:78], sim_data$cases[1:78], type = "b", col = "darkred", pch = 16, lwd = 2,
     main = "Reported Zombie Spottings", xlab = "Days", ylab = "Reported Zombie Cases")

# write.csv(student_sim[, c("time", "cases")], "outbreak_data.csv", row.names = FALSE)