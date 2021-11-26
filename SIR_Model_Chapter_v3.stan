//------------------------------------------------------------------------
//   SEIQR Model
//------------------------------------------------------------------------
//
// DESCRIPTION
//------------------------------------------------------------------------
// The following code adjusts a SEIQR model to COVID-19 data from Mexico
// and conducts simulations for different scenarios for the paper.
//
// INPUTS
//------------------------------------------------------------------------
// The parameters involved in the model are as follows:
// real t .- Time 
// vector y      
//      y[1]  = Susceptible
//      y[2]  = Exposed
//      y[3]  = Asymptomatic
//      y[4]  = Infected
//      y[5]  = Hospitalized
//      y[6]  = ICU
//      y[7]  = Death
//      y[8]  = Cummulative number of infected
//      y[9]  = Cummulative number of hospitalized
//      y[10] = Cummulative number of icu patients
//      y[11] = Quarantine of suceptibles (1st kind)
//      y[12] = Quarantine of exposed (2nd kind)
//      y[13] = Quarantine of asyntomatic (3rd kind)
//      y[14] = Quarantine of infected (4th kind)
// real beta_0  .- Infection rate for Asymptomatic
// real beta_1  .- Infection rate for Infected
// real beta_2  .- Infection rate for ICU
// real beta_3  .- Infection rate for Hospitalized
// real beta_12 .- Rate of infectious to hospitalized
// real beta_23 .- Rate of hospitalized -> UCI 
// real eta_E,  .- Rate from incubating to symptomatic
// real gamma_0R .- Rate of recovery for asymptomatics
// real gamma_1R .- Rate of recovery for infected 
// real gamma_2R .- Rate of recovery for hospitalized
// real gamma_3R .- Rate of recovery for ICU
// real theta_3M .- Mortality rate for COVID-19
// vector coef   .- Vector of coefficients for alpha(t)
// real p_A,     .- Proportion of asymptomatics
// real p_I2,    .- Proportion of infected -> hospitalized
// real p_I3,    .- Proportion of hospitalized -> UCI
// real p_M      .- Proportion of death from ICU
// real q2       .- Proportion of infected quarantined                    

functions {
  
  //Polynomial of 4th order
  real alpha(real t, vector coef){
    real poly;
    poly = 1  + 
        coef[1]*(t/200) + coef[2]*(t/200)^2 + coef[3]*(t/200)^3 + 
        coef[4]*(t/200)^4 + coef[5]*(t/200)^5;
    return poly;
  }
  
  //dC1 Cummulative number of I1
  real dC1(real p_A, real eta_E, vector y){
    return (1 - p_A)*eta_E*y[2];
  }
  
  //dC2 Cummulative number of I2
  real dC2(real p_I2, real beta_12, vector y){
    return p_I2*beta_12*(y[4] + y[12]);
  }
  
  //dC3 Cummulative number of I3
  real dC3(real p_I3, real beta_23, vector y){
    return p_I3*beta_23*y[5];
  }
  
  //dS = -alpha(t)*(beta_0*Asymptomatic + beta_1*Infected + beta_2*UCI + 
  //         beta_3*Hospitalized)*S 
  real dS(real t, vector coef, real beta_0, real beta_1, real beta_2, 
            real beta_3, vector y){
    return -alpha(t, coef)*(beta_0*y[3] + beta_1*y[4] + beta_2*y[5] 
                + beta_3*y[6])*y[1];
  }
  
  //dE = alpha(t)*(beta_0*Asymptomatic + beta_1*Infected + beta_2*UCI + 
    //         beta_3*Hospitalized)*S - eta_E*E
  real dE(real t, vector coef, real beta_0, real beta_1, real beta_2, 
            real beta_3, real eta_E, vector y){
    return alpha(t, coef)*(beta_0*y[3] + beta_1*y[4] + beta_2*y[5] 
                + beta_3*y[6])*y[1] - eta_E*y[2];
  }
  
  //dA = p_A*eta_E*E - gamma_0R*A
  real dA(real eta_E, real p_A, real gamma_0R, vector y){
    return p_A*eta_E*y[2] - gamma_0R*y[3];
  }
  
  //dI1 = (1 - p_A)*eta_E*E - (p_I2*beta_12 + (1 - p_I2)*gamma_1R)*I1
  real dI1(real eta_E, real p_A, real p_I2, real beta_12, real gamma_1R, real q2,
            vector y){
    return dC1(p_A, eta_E, y) - ((1 - p_I2)*gamma_1R + q2)*y[4] - dC2(p_I2, beta_12, y);
  }
  
  //dI2 = p_I2*beta_12*I_1 - (p_I3*beta_23 + (1 - p_I3)*gamma_2R)*I2
  real dI2(real p_I2, real beta_12, real p_I3, real beta_23, real gamma_2R,
            vector y){
    return dC2(p_I2, beta_12, y) - (1 - p_I3)*gamma_2R*y[5] - dC3(p_I3, beta_23, y);
  }
  
  //dI3 = p_I3*beta_23*I2 - ((1 - pM)*gamma_3R + pM*theta_3M)*I3
  real dI3(real p_I3, real beta_23, real p_M, real gamma_3R, real theta_3M,
            vector y){
    return dC3(p_I3, beta_23, y) - ((1 - p_M)*gamma_3R + p_M*theta_3M)*y[6];
  }
  
  //dM  = pM*theta_3M*I3
  real dM(real p_M, real theta_3M, vector y){
    return p_M*theta_3M*y[6];
  }
  
  //dQ1 = 0 
  real dQ1(){
    return 0;
  }
  
  //dQ2 = -eta_E*Q2
  real dQ2(real eta_E, vector y){
    return -eta_E*y[12];
  }
  
  //dQ3 = p_A*eta_E*Q2 - gamma_0R*Q3;
  real dQ3(real p_A, real eta_E, real gamma_0R, vector y){
    return p_A*eta_E*y[12] - gamma_0R*y[13];
  }
  
  //dQ4 = (1 - p_A)*eta_E*Q2 + q2*I1 - (p_I2*beta_12 + (1 - p_I2)*gamma_1R)*y[14];
  real dQ4(real p_A, real eta_E, real p_I2, real beta_12, real gamma_1R, 
            real q2, vector y){
    return (1 - p_A)*eta_E*y[12] +  q2*y[4] - 
              (p_I2*beta_12 + (1 - p_I2)*gamma_1R)*y[14];
  }
  
  //MODEL 1: Model for fitting SIR to data
  //---------------------------------------------------------
  //Description: In this model no quarantine is assumed. 
  vector SIR_fitted(real t, vector y, real beta_0, real beta_1, real beta_2, 
              real beta_3, real beta_12, real beta_23,  real eta_E, 
              real gamma_0R, real gamma_1R, real gamma_2R, real gamma_3R, 
              real theta_3M, vector coef, real p_A, real p_I2, real p_I3,
              real p_M) {
                    
    vector[14] dydt;
  
    dydt[1]  = dS(t, coef, beta_0, beta_1, beta_2, beta_3, y);
    dydt[2]  = dE(t, coef, beta_0, beta_1, beta_2, beta_3, eta_E, y);
    dydt[3]  = dA(eta_E, p_A, gamma_0R, y);
    dydt[4]  = dI1(eta_E, p_A, p_I2, beta_12, gamma_1R, 0.0, y);
    dydt[5]  = dI2(p_I2, beta_12, p_I3, beta_23, gamma_2R, y);
    dydt[6]  = dI3(p_I3, beta_23, p_M, gamma_3R, theta_3M, y);
    dydt[7]  = dM(p_M, theta_3M, y);
    dydt[8]  = dC1(p_A, eta_E, y);
    dydt[9]  = dC2(p_I2, beta_12, y);
    dydt[10] = dC3(p_I3, beta_23, y);
    dydt[11] = dQ1();
    dydt[12] = dQ2(eta_E,y);
    dydt[13] = dQ3(p_A, eta_E, gamma_0R, y);
    dydt[14] = dQ4(p_A, eta_E, p_I2, beta_12, gamma_1R, 0.0, y);
      
    return dydt;
  }
  
  //MODEL 2
  //---------------------------------------------------------
  //Description: In this model quarantine for a proportion q2 of infected
  vector SIR_quarantine_infected(real t, vector y, real beta_0, real beta_1, 
      real beta_2, real beta_3, real beta_12, real beta_23,  real eta_E, 
      real gamma_0R, real gamma_1R, real gamma_2R, real gamma_3R, real theta_3M, 
      vector coef, real p_A, real p_I2, real p_I3, real p_M, real q2) {
                    
    vector[14] dydt;
  
    dydt[1]  = dS(t, coef, beta_0, beta_1, beta_2, beta_3, y);
    dydt[2]  = dE(t, coef, beta_0, beta_1, beta_2, beta_3, eta_E, y);
    dydt[3]  = dA(eta_E, p_A, gamma_0R, y);
    dydt[4]  = dI1(eta_E, p_A, p_I2, beta_12, gamma_1R, q2, y);
    dydt[5]  = dI2(p_I2, beta_12, p_I3, beta_23, gamma_2R, y);
    dydt[6]  = dI3(p_I3, beta_23, p_M, gamma_3R, theta_3M, y);
    dydt[7]  = dM(p_M, theta_3M, y);
    dydt[8]  = dC1(p_A, eta_E, y);
    dydt[9]  = dC2(p_I2, beta_12, y);
    dydt[10] = dC3(p_I3, beta_23, y);
    dydt[11] = dQ1();
    dydt[12] = dQ2(eta_E,y);
    dydt[13] = dQ3(p_A, eta_E, gamma_0R, y);
    dydt[14] = dQ4(p_A, eta_E, p_I2, beta_12, gamma_1R, q2, y);
      
    return dydt;
  }
  
  //MODEL 3
  //---------------------------------------------------------
  //Description: In this model when saturation_level_hosp or 
  //saturation_level_icu are reached there is a change in 
  //beta_23, p_M and p_I3 
  vector SIR_saturation(real t, vector y, real beta_0, real beta_1, 
      real beta_2, real beta_3, real beta_12, real beta_23_unsaturated,  
      real eta_E,  real gamma_0R, real gamma_1R, real gamma_2R, 
      real gamma_3R, real theta_3M_unsaturated, vector coef, real p_A, 
      real p_I2, real p_I3_unsaturated, real p_M_unsaturated, real q2,
      real beta_23_saturated, real p_I3_saturated, real p_M_saturated,
      real theta_3M_saturated, real saturation_level_hosp, 
      real saturation_level_icu) {
         
    real beta_23;
    real p_M;
    real p_I3;
    real theta_3M;
    vector[14] dydt;
    
    //Check for saturation of hospitals
    if (y[5] < saturation_level_hosp){
      beta_23 = beta_23_unsaturated;
      p_I3    = p_I3_unsaturated; 
    } else {
      beta_23 = beta_23_saturated;
      p_I3    = p_I3_saturated;
    }
    
    //Check for saturation of icu
    if (y[6] < saturation_level_icu){
      p_M      = p_M_unsaturated;
      theta_3M = theta_3M_unsaturated;
    } else {
      p_M      = p_M_saturated;
      theta_3M = theta_3M_saturated;
    }
  
    dydt[1]  = dS(t, coef, beta_0, beta_1, beta_2, beta_3, y);
    dydt[2]  = dE(t, coef, beta_0, beta_1, beta_2, beta_3, eta_E, y);
    dydt[3]  = dA(eta_E, p_A, gamma_0R, y);
    dydt[4]  = dI1(eta_E, p_A, p_I2, beta_12, gamma_1R, q2, y);
    dydt[5]  = dI2(p_I2, beta_12, p_I3, beta_23, gamma_2R, y);
    dydt[6]  = dI3(p_I3, beta_23, p_M, gamma_3R, theta_3M, y);
    dydt[7]  = dM(p_M, theta_3M, y);
    dydt[8]  = dC1(p_A, eta_E, y);
    dydt[9]  = dC2(p_I2, beta_12, y);
    dydt[10] = dC3(p_I3, beta_23, y);
    dydt[11] = dQ1();
    dydt[12] = dQ2(eta_E,y);
    dydt[13] = dQ3(p_A, eta_E, gamma_0R, y);
    dydt[14] = dQ4(p_A, eta_E, p_I2, beta_12, gamma_1R, q2, y);
      
    return dydt;
  }
}

data {
  //Data
  int<lower=0> N;         // Number of times observed
  real t[N];              // Times observed
  vector[14] y0;          // Vector of initial state
  vector[N] Infected;     // Number of infected individuals
  vector[N] Hospitalized; // Number of hospitalized individuals
  vector[N] ICU;          // Number of ICU individuals
  vector[N] Death;        // Number of death individuals
  
  //For simulations:
  int<lower=0> Pop;               // Population size
  int<lower=0> hospital_capacity; // Maximum number of beds before saturation
  int<lower=0> icu_capacity;      // Maximum number of icus before saturation
  real<lower=0,upper=1> q2_low;   // Lower bound for q2 simulation
  real<lower=0,upper=1> q2_up;    // Upper bound for q2 simulation
  real<lower=0,upper=1> exposed_quarantine_low; // Lower bound for %exposed quarantined
  real<lower=0,upper=1> exposed_quarantine_up;  // Upper bound for %exposed quarantined
  
  //Lengths of start and quarantines
  int<lower=0>  pre_first_quarantine_length;
  int<lower=0>  time_first_quarantine_start;
  int<lower=0>  time_first_quarantine_end;
  int<lower=0>  first_quarantine_length;
  int<lower=0>  pre_second_quarantine_length;
  int<lower=0>  time_second_quarantine_start;
  int<lower=0>  time_second_quarantine_end;
  int<lower=0>  second_quarantine_length;
  int<lower=0>  post_second_quarantine_length;
  real t_pre_first[pre_first_quarantine_length];
  real t_first[first_quarantine_length];
  real t_pre_second[pre_second_quarantine_length];
  real t_second[second_quarantine_length];
  real t_post_second[post_second_quarantine_length];
  
  //For periodic quarantine:
  int<lower=0> days_on;  //Days without quarantine
  int<lower=0> days_off; //Days with quarantine
  int<lower=0> number_of_periods; //Number of periods in days_on / days_off
  int<lower=0> periodic_quarantine_start; //When the periodic quarantine is implemented
  int<lower=0> extra_days;
  int<lower=0> extra_days_on;
  int<lower=0> extra_days_off;
  real t_pre_periodic[periodic_quarantine_start];
  real t_periodic[days_on + days_off];
  real t_periodic_on[days_on];
  real t_periodic_off[days_off];
  real t_post_periodic_on[extra_days_on];
  real t_post_periodic_off[extra_days_off];
}

transformed data {
  real t0   = 0;     //Initial time
  int steps = 10000; //Steps for ODE solver
  real tol1 = 1e-10; //Relative tolerance
  real tol2 = 1e-10; //Absolute tolerance
}

parameters {
  
  real<lower=0, upper = 1> beta_0;         //Asymptomatic
  real<lower=0, upper = 1> beta_1;         //Infected
  real<lower=0, upper = 1> beta_2;         //ICU
  real<lower=0, upper = 1> beta_3;         //Hospitalized
  real<lower=0, upper = 1> beta_12;        //Rate of infectious to hospitalized
  real<lower=0, upper = 1> beta_23;        //Rate of hospitalized -> UCI 
  real<lower=0, upper = 1> eta_E;          //Rate from incubating to symptomatic
  real<lower=0, upper = 1> gamma_0R;       //Rate of recovery for asymptomatics
  real<lower=0, upper = 1> gamma_1R;       //Rate of recovery for infected 
  real<lower=0, upper = 1> gamma_2R;       //Rate of recovery for hospitalized
  real<lower=0, upper = 1> gamma_3R;       //Rate of recovery for ICU
  real<lower=0, upper = 1> theta_3M;       //Mortality rate for COVID-19
  vector[5] coef;                          //Vector of coefficients for alpha(t)
  real<lower=0, upper=1> p_A;              //Proportion of asymptomatics
  real<lower=0, upper=1> p_I2;             //Proportion of infected -> hospitalized
  real<lower=0, upper=1> p_I3;             //Proportion of hospitalized -> UCI
  real<lower=0, upper=1> p_M;              //Proportion of death from ICU
  real<lower=0> sigma_death;               //Variance for dead
  real<lower=0> sigma_hospitalized;        //Variance for hospitalized
  real<lower=0> sigma_infected;            //Variance for infected
  real<lower=0> sigma_icu;                 //Variance for ICU
}

transformed parameters {
  vector<lower=0>[14] y[N] = ode_rk45_tol(SIR_fitted, y0, t0, t, tol1, tol2, 
      steps, beta_0, beta_1, beta_2, beta_3, beta_12, beta_23, eta_E, gamma_0R,
      gamma_1R, gamma_2R, gamma_3R, theta_3M, coef, p_A, p_I2, p_I3, p_M);
}
  
model {
  beta_0   ~ normal(0.8276086119485562, 0.001);
  beta_1   ~ normal(0.5515727329250524, 0.001);
  beta_2   ~ normal(0.9205242957230675, 0.001);
  beta_3   ~ normal(0.311219365300846, 0.001);
  beta_12  ~ normal(0.32561971100893067, 0.001);
  beta_23  ~ normal(0.2015508991016397, 0.001);
  eta_E    ~ normal(0.9960149806661017, 0.001);
  gamma_0R ~ normal(0.7530968545973915, 0.001);
  gamma_1R ~ normal(0.8817345287734505, 0.001);
  gamma_2R ~ normal(0.30647950047778616, 0.001);
  gamma_3R ~ normal(0.16014664285960167, 0.001);
  theta_3M ~ normal(0.46510693286688437, 0.001);
  coef[1]  ~ normal(0.0, 0.001);
  coef[2]  ~ normal(0.009881995293258896, 0.001);
  coef[3]  ~ normal(0.18586993515558667, 0.001);
  coef[4]  ~ normal(0.05386425435891228, 0.001);
  coef[5]  ~ normal(0.001711630633523908, 0.001);
  p_A      ~ normal(0.9752388420669302, 0.001);
  p_I2     ~ normal(0.42646170956461865, 0.001);
  p_I3     ~ normal(0.39268337611472104, 0.001);
  p_M      ~ normal(0.8738793655806338, 0.001);
  
  sigma_infected     ~ cauchy(0.001, 1.0);
  sigma_death        ~ cauchy(0.001, 1.0);
  sigma_icu          ~ cauchy(0.001, 1.0);
  sigma_hospitalized ~ cauchy(0.001, 1.0);
  
  Infected     ~ normal(y[, 8],  sigma_infected);
  Hospitalized ~ normal(y[, 9],  sigma_hospitalized);
  ICU          ~ normal(y[, 10], sigma_icu);
  Death        ~ normal(y[, 7],  sigma_death);
}

generated quantities{
  
  //---------------------------------------------------------
  // MODEL 1: Data as fitted in model 
  //---------------------------------------------------------
  
  //Observed variables
  real Infected_model_1[N];
  real Hospitalized_model_1[N];
  real ICU_model_1[N];
  real Death_model_1[N];
  
  vector[14] y_model_1[N] = ode_rk45_tol(SIR_fitted, y0, t0, t, 1e-6, 1e-6, 
      steps, beta_0, beta_1, beta_2, beta_3, beta_12, beta_23, eta_E, gamma_0R,
      gamma_1R, gamma_2R, gamma_3R, theta_3M, coef, p_A, p_I2, p_I3, p_M);
  
  for (j in 1:N) {
    
    Death_model_1[j]        = y_model_1[j,7]  + normal_rng(0,sigma_death);
    Infected_model_1[j]     = y_model_1[j,8]  + normal_rng(0,sigma_infected);
    Hospitalized_model_1[j] = y_model_1[j,9]  + normal_rng(0,sigma_hospitalized);
    ICU_model_1[j]          = y_model_1[j,10] + normal_rng(0,sigma_icu);
    
    if (Death_model_1[j] <= 0){
      Death_model_1[j] = 0;
    }
    
    if (Infected_model_1[j] <= 0){
      Infected_model_1[j] = 0;
    }
    
    if (Hospitalized_model_1[j] <= 0){
      Hospitalized_model_1[j] = 0;
    }
    
    if (ICU_model_1[j] <= 0){
      ICU_model_1[j] = 0;
    }
  }
  
  //---------------------------------------------------------
  // MODEL 2: Assuming ~50% of infected quarantined
  //---------------------------------------------------------
  
  //Simulate quarantine
  real q2_model_2 = uniform_rng(q2_low, q2_up);
  
  vector[14] y_model_2[N] = ode_rk45_tol(SIR_quarantine_infected, y0, t0, t, 
      1e-6, 1e-6, steps, beta_0, beta_1, beta_2, beta_3, beta_12, beta_23, 
      eta_E, gamma_0R, gamma_1R, gamma_2R, gamma_3R, theta_3M, coef, p_A, 
      p_I2, p_I3, p_M, q2_model_2);
  
  //---------------------------------------------------------
  // MODEL 3: Quarantine all
  //---------------------------------------------------------    

  //First we start without quarantine
  vector[14] y_model_3_1[pre_first_quarantine_length] = 
    ode_rk45_tol(SIR_quarantine_infected, y0, t0, t_pre_first, 
      1e-6, 1e-6, steps, beta_0, beta_1, beta_2, beta_3, beta_12, beta_23, 
      eta_E, gamma_0R, gamma_1R, gamma_2R, gamma_3R, theta_3M, coef, p_A, 
      p_I2, p_I3, p_M, 0.0);
        
  //Then we quarantine the population of susceptibles
  vector[14] y_init_first = [
    0.3*y_model_3_1[time_first_quarantine_start, 1], 
    0.3*y_model_3_1[time_first_quarantine_start, 2],
    0.3*y_model_3_1[time_first_quarantine_start, 3],
    y_model_3_1[time_first_quarantine_start, 4],
    y_model_3_1[time_first_quarantine_start, 5],
    y_model_3_1[time_first_quarantine_start, 6],
    y_model_3_1[time_first_quarantine_start, 7],
    y_model_3_1[time_first_quarantine_start, 8],
    y_model_3_1[time_first_quarantine_start, 9],
    y_model_3_1[time_first_quarantine_start, 10],
    0.7*y_model_3_1[time_first_quarantine_start, 1],
    0.7*y_model_3_1[time_first_quarantine_start, 2],
    0.7*y_model_3_1[time_first_quarantine_start, 3],
    y_model_3_1[time_first_quarantine_start, 14]]';
  vector[14] y_model_3_2[first_quarantine_length] = 
    ode_rk45_tol(SIR_quarantine_infected, y_init_first, 
      time_first_quarantine_start, t_first, 1e-6, 1e-6, steps, beta_0, 
      beta_1, beta_2, beta_3, beta_12, beta_23, 
      eta_E, gamma_0R, gamma_1R, gamma_2R, gamma_3R, theta_3M, coef, p_A, 
      p_I2, p_I3, p_M, 0.7);
        
  //We unquarantine the individuals
  vector[14] y_init_pre_second = [
    y_model_3_2[first_quarantine_length, 1] + y_model_3_2[first_quarantine_length, 11], 
    y_model_3_2[first_quarantine_length, 2] + y_model_3_2[first_quarantine_length, 12],
    y_model_3_2[first_quarantine_length, 3] + y_model_3_2[first_quarantine_length, 13],
    y_model_3_2[first_quarantine_length, 4],
    y_model_3_2[first_quarantine_length, 5],
    y_model_3_2[first_quarantine_length, 6],
    y_model_3_2[first_quarantine_length, 7],
    y_model_3_2[first_quarantine_length, 8],
    y_model_3_2[first_quarantine_length, 9],
    y_model_3_2[first_quarantine_length, 10],
    0.0,
    0.0,
    0.0,
    y_model_3_2[first_quarantine_length, 14]]';
  vector[14] y_model_3_3[pre_second_quarantine_length] = 
    ode_rk45_tol(SIR_quarantine_infected, y_init_pre_second, 
      time_first_quarantine_end, t_pre_second, 1e-6, 1e-6, steps, beta_0, 
      beta_1, beta_2, beta_3, beta_12, beta_23, 
      eta_E, gamma_0R, gamma_1R, gamma_2R, gamma_3R, theta_3M, coef, p_A, 
      p_I2, p_I3, p_M, 0.0);   
        
  //Quarantine for second wave
  vector[14] y_init_second = [
    0.4*y_model_3_3[pre_second_quarantine_length, 1], 
    0.4*y_model_3_3[pre_second_quarantine_length, 2],
    0.4*y_model_3_3[pre_second_quarantine_length, 3],
    y_model_3_3[pre_second_quarantine_length, 4],
    y_model_3_3[pre_second_quarantine_length, 5],
    y_model_3_3[pre_second_quarantine_length, 6],
    y_model_3_3[pre_second_quarantine_length, 7],
    y_model_3_3[pre_second_quarantine_length, 8],
    y_model_3_3[pre_second_quarantine_length, 9],
    y_model_3_3[pre_second_quarantine_length, 10],
    0.6*y_model_3_3[pre_second_quarantine_length, 1],
    0.6*y_model_3_3[pre_second_quarantine_length, 12],
    0.6*y_model_3_3[pre_second_quarantine_length, 13],
    y_model_3_3[pre_second_quarantine_length, 14]]';
  vector[14] y_model_3_4[second_quarantine_length] = 
    ode_rk45_tol(SIR_quarantine_infected, y_init_second, 
      time_second_quarantine_start, t_second, 1e-6, 1e-6, steps, beta_0, 
      beta_1, beta_2, beta_3, beta_12, beta_23, 
      eta_E, gamma_0R, gamma_1R, gamma_2R, gamma_3R, theta_3M, coef, p_A, 
      p_I2, p_I3, p_M, 0.6);  
        
  //unquarantine after second wave
  vector[14] y_init_post_second = [
    y_model_3_4[second_quarantine_length, 1] + y_model_3_4[second_quarantine_length, 11], 
    y_model_3_4[second_quarantine_length, 2] + y_model_3_4[second_quarantine_length, 12],
    y_model_3_4[second_quarantine_length, 3] + y_model_3_4[second_quarantine_length, 13],
    y_model_3_4[second_quarantine_length, 4],
    y_model_3_4[second_quarantine_length, 5],
    y_model_3_4[second_quarantine_length, 6],
    y_model_3_4[second_quarantine_length, 7],
    y_model_3_4[second_quarantine_length, 8],
    y_model_3_4[second_quarantine_length, 9],
    y_model_3_4[second_quarantine_length, 10],
    0.0,
    0.0,
    0.0,
    y_model_3_4[second_quarantine_length, 14]]';
  vector[14] y_model_3_5[post_second_quarantine_length] = 
    ode_rk45_tol(SIR_quarantine_infected, y_init_post_second, 
      time_second_quarantine_end, t_post_second, 1e-6, 1e-6, steps, beta_0, 
      beta_1, beta_2, beta_3, beta_12, beta_23, 
      eta_E, gamma_0R, gamma_1R, gamma_2R, gamma_3R, theta_3M, coef, p_A, 
      p_I2, p_I3, p_M, 0.0); 
  
  vector[14] y_model_3[N];
  
  for (i in 1:time_first_quarantine_start){
    y_model_3[i] = y_model_3_1[i];
  }
  
  for (i in (time_first_quarantine_start + 1):time_first_quarantine_end){
    y_model_3[i] = y_model_3_2[i - time_first_quarantine_start];
  }
  
  for (i in (time_first_quarantine_end + 1):time_second_quarantine_start){
    y_model_3[i] = y_model_3_3[i - time_first_quarantine_end];
  }
  
  for (i in (time_second_quarantine_start + 1):time_second_quarantine_end){
    y_model_3[i] = y_model_3_4[i - time_second_quarantine_start];
  }
  
  for (i in (time_second_quarantine_end + 1):N){
    y_model_3[i] = y_model_3_5[i - time_second_quarantine_end];
  }
      
  //---------------------------------------------------------
  // MODEL 4: Quarantine for exposed
  //---------------------------------------------------------    

  //First we start without quarantine
  vector[14] y_model_4_1[pre_first_quarantine_length] = 
    ode_rk45_tol(SIR_quarantine_infected, y0, t0, t_pre_first, 
      1e-6, 1e-6, steps, beta_0, beta_1, beta_2, beta_3, beta_12, beta_23, 
      eta_E, gamma_0R, gamma_1R, gamma_2R, gamma_3R, theta_3M, coef, p_A, 
      p_I2, p_I3, p_M, 0.0);
        
  //Then we quarantine the population of susceptibles
  vector[14] y_init_first_4 = [
    y_model_4_1[time_first_quarantine_start, 1], 
    0.3*y_model_4_1[time_first_quarantine_start, 2],
    y_model_4_1[time_first_quarantine_start, 3],
    y_model_4_1[time_first_quarantine_start, 4],
    y_model_4_1[time_first_quarantine_start, 5],
    y_model_4_1[time_first_quarantine_start, 6],
    y_model_4_1[time_first_quarantine_start, 7],
    y_model_4_1[time_first_quarantine_start, 8],
    y_model_4_1[time_first_quarantine_start, 9],
    y_model_4_1[time_first_quarantine_start, 10],
    y_model_4_1[time_first_quarantine_start, 1],
    0.7*y_model_4_1[time_first_quarantine_start, 2],
    y_model_4_1[time_first_quarantine_start, 13],
    y_model_4_1[time_first_quarantine_start, 14]]';
  vector[14] y_model_4_2[first_quarantine_length] = 
    ode_rk45_tol(SIR_quarantine_infected, y_init_first_4, 
      time_first_quarantine_start, t_first, 1e-6, 1e-6, steps, beta_0, 
      beta_1, beta_2, beta_3, beta_12, beta_23, 
      eta_E, gamma_0R, gamma_1R, gamma_2R, gamma_3R, theta_3M, coef, p_A, 
      p_I2, p_I3, p_M, 0.0);
        
    //We unquarantine the individuals
    vector[14] y_init_pre_second_4 = [
      y_model_4_2[first_quarantine_length, 1], 
      y_model_4_2[first_quarantine_length, 2],
      y_model_4_2[first_quarantine_length, 3],
      y_model_4_2[first_quarantine_length, 4],
      y_model_4_2[first_quarantine_length, 5],
      y_model_4_2[first_quarantine_length, 6],
      y_model_4_2[first_quarantine_length, 7],
      y_model_4_2[first_quarantine_length, 8],
      y_model_4_2[first_quarantine_length, 9],
      y_model_4_2[first_quarantine_length, 10],
      y_model_4_2[first_quarantine_length, 11],
      y_model_4_2[first_quarantine_length, 12],
      y_model_4_2[first_quarantine_length, 13],
      y_model_4_2[first_quarantine_length, 14]]';
    vector[14] y_model_4_3[pre_second_quarantine_length] = 
      ode_rk45_tol(SIR_quarantine_infected, y_init_pre_second_4, 
        time_first_quarantine_end, t_pre_second, 1e-6, 1e-6, steps, beta_0, 
        beta_1, beta_2, beta_3, beta_12, beta_23, 
        eta_E, gamma_0R, gamma_1R, gamma_2R, gamma_3R, theta_3M, coef, p_A, 
        p_I2, p_I3, p_M, 0.0);   
        
    //Quarantine for second wave
    vector[14] y_init_second_4 = [
      y_model_4_3[pre_second_quarantine_length, 1], 
      0.4*y_model_4_3[pre_second_quarantine_length, 2],
      y_model_4_3[pre_second_quarantine_length, 3],
      y_model_4_3[pre_second_quarantine_length, 4],
      y_model_4_3[pre_second_quarantine_length, 5],
      y_model_4_3[pre_second_quarantine_length, 6],
      y_model_4_3[pre_second_quarantine_length, 7],
      y_model_4_3[pre_second_quarantine_length, 8],
      y_model_4_3[pre_second_quarantine_length, 9],
      y_model_4_3[pre_second_quarantine_length, 10],
      0.6*y_model_4_3[pre_second_quarantine_length, 1],
      y_model_4_3[pre_second_quarantine_length, 12],
      y_model_4_3[pre_second_quarantine_length, 13],
      y_model_4_3[pre_second_quarantine_length, 14]]';
    vector[14] y_model_4_4[second_quarantine_length] = 
      ode_rk45_tol(SIR_quarantine_infected, y_init_second_4, 
        time_second_quarantine_start, t_second, 1e-6, 1e-6, steps, beta_0, 
        beta_1, beta_2, beta_3, beta_12, beta_23, 
        eta_E, gamma_0R, gamma_1R, gamma_2R, gamma_3R, theta_3M, coef, p_A, 
        p_I2, p_I3, p_M, 0.0);  
        
    //unquarantine after second wave
    vector[14] y_init_post_second_4 = [
      y_model_4_4[second_quarantine_length, 1], 
      y_model_4_4[second_quarantine_length, 2] + y_model_4_4[second_quarantine_length, 12],
      y_model_4_4[second_quarantine_length, 3],
      y_model_4_4[second_quarantine_length, 4],
      y_model_4_4[second_quarantine_length, 5],
      y_model_4_4[second_quarantine_length, 6],
      y_model_4_4[second_quarantine_length, 7],
      y_model_4_4[second_quarantine_length, 8],
      y_model_4_4[second_quarantine_length, 9],
      y_model_4_4[second_quarantine_length, 10],
      y_model_4_4[second_quarantine_length, 11],
      0.0,
      y_model_4_4[second_quarantine_length, 13],
      y_model_4_4[second_quarantine_length, 14]]';
    vector[14] y_model_4_5[post_second_quarantine_length] = 
      ode_rk45_tol(SIR_quarantine_infected, y_init_post_second_4, 
        time_second_quarantine_end, t_post_second, 1e-6, 1e-6, steps, beta_0, 
        beta_1, beta_2, beta_3, beta_12, beta_23, 
        eta_E, gamma_0R, gamma_1R, gamma_2R, gamma_3R, theta_3M, coef, p_A, 
        p_I2, p_I3, p_M, 0.0); 
  
  vector[14] y_model_4[N];
  
  for (i in 1:time_first_quarantine_start){
    y_model_4[i] = y_model_4_1[i];
  }
  
  for (i in (time_first_quarantine_start + 1):time_first_quarantine_end){
    y_model_4[i] = y_model_4_2[i - time_first_quarantine_start];
  }
  
  for (i in (time_first_quarantine_end + 1):time_second_quarantine_start){
    y_model_4[i] = y_model_4_3[i - time_first_quarantine_end];
  }
  
  for (i in (time_second_quarantine_start + 1):time_second_quarantine_end){
    y_model_4[i] = y_model_4_4[i - time_second_quarantine_start];
  }
  
  for (i in (time_second_quarantine_end + 1):N){
    y_model_4[i] = y_model_4_5[i - time_second_quarantine_end];
  }
  
  //---------------------------------------------------------
  // MODEL 5: Hospital saturation
  //---------------------------------------------------------  
  vector[14] y_model_5[N] = ode_rk45_tol(SIR_saturation, y0, t0, t, 
      1e-6, 1e-6, steps, beta_0, beta_1, beta_2, beta_3, beta_12, beta_23, 
      eta_E, gamma_0R, gamma_1R, gamma_2R, gamma_3R, theta_3M, coef, p_A, 
      p_I2, p_I3, p_M, 0.0, 100.0, 1.0, 1.0, 10.0, hospital_capacity * 1.0 / Pop, 
      icu_capacity * 1.0 / Pop);
      
  //---------------------------------------------------------
  // MODEL 6: Weekly scheme 4 days ON | 3 days off for 70%
  //--------------------------------------------------------- 
  vector[14] y_model_6[N];
  vector[14] y_0_model_6;
  vector[14] y_model_6_aux_on[days_on];
  vector[14] y_model_6_aux_off[days_off];
  
  //First period without periodic quarantine
  vector[14] y_model_6_1[periodic_quarantine_start] = ode_rk45_tol(
      SIR_quarantine_infected, y0, t0, 
      t_pre_periodic, 1e-6, 1e-6, steps, beta_0, beta_1, beta_2, beta_3, 
      beta_12, beta_23, eta_E, gamma_0R, gamma_1R, gamma_2R, gamma_3R, 
      theta_3M, coef, p_A, p_I2, p_I3, p_M, 0.0);
  
  //Assign to model
  y_model_6[1:periodic_quarantine_start] = y_model_6_1;
  
  y_0_model_6 = [
      0.25*y_model_6_1[periodic_quarantine_start, 1], 
      0.25*y_model_6_1[periodic_quarantine_start, 2],
      0.25*y_model_6_1[periodic_quarantine_start, 3],
      y_model_6_1[periodic_quarantine_start, 4],
      y_model_6_1[periodic_quarantine_start, 5],
      y_model_6_1[periodic_quarantine_start, 6],
      y_model_6_1[periodic_quarantine_start, 7],
      y_model_6_1[periodic_quarantine_start, 8],
      y_model_6_1[periodic_quarantine_start, 9],
      y_model_6_1[periodic_quarantine_start, 10],
      0.75*y_model_6_1[periodic_quarantine_start, 1],
      0.75*y_model_6_1[periodic_quarantine_start, 2],
      0.75*y_model_6_1[periodic_quarantine_start, 3],
      y_model_6_1[periodic_quarantine_start, 14]]';
      
  for (period in 1:number_of_periods){
    
    //Days in quarantine
    y_model_6_aux_on = ode_rk45_tol(
      SIR_quarantine_infected, y_0_model_6, 0.0, 
      t_periodic_on, 1e-6, 1e-6, steps, beta_0, beta_1, beta_2, beta_3, 
      beta_12, beta_23, eta_E, gamma_0R, gamma_1R, gamma_2R, gamma_3R, 
      theta_3M, coef, p_A, p_I2, p_I3, p_M, 0.75);
      
    //Assign to model
    y_model_6[(periodic_quarantine_start + (period - 1)*(days_on + days_off) + 1):(periodic_quarantine_start + (period - 1)*(days_on + days_off) + days_on)] = y_model_6_aux_on;
    
    //Unquarantine everyone
    y_0_model_6 = [
      y_model_6_aux_on[days_on, 1] + y_model_6_aux_on[days_on, 11], 
      y_model_6_aux_on[days_on, 2] + y_model_6_aux_on[days_on, 12],
      y_model_6_aux_on[days_on, 3] + y_model_6_aux_on[days_on, 13],
      y_model_6_aux_on[days_on, 4],
      y_model_6_aux_on[days_on, 5],
      y_model_6_aux_on[days_on, 6],
      y_model_6_aux_on[days_on, 7],
      y_model_6_aux_on[days_on, 8],
      y_model_6_aux_on[days_on, 9],
      y_model_6_aux_on[days_on, 10],
      0.0,
      0.0,
      0.0,
      y_model_6_aux_on[days_on, 14]]';
      
    //Days off quarantine
    y_model_6_aux_off = ode_rk45_tol(
      SIR_quarantine_infected, y_0_model_6, 0.0, 
      t_periodic_off, 1e-6, 1e-6, steps, beta_0, beta_1, beta_2, beta_3, 
      beta_12, beta_23, eta_E, gamma_0R, gamma_1R, gamma_2R, gamma_3R, 
      theta_3M, coef, p_A, p_I2, p_I3, p_M, 0.0);  
    
    //Assign to model
    y_model_6[(periodic_quarantine_start + (period - 1)*(days_on + days_off) + days_on + 1):(periodic_quarantine_start + period*(days_on + days_off))] = y_model_6_aux_off;
    
    //Quarantine everyone
    y_0_model_6 = [
      0.25*y_model_6_aux_off[days_off, 1], 
      0.25*y_model_6_aux_off[days_off, 2],
      0.25*y_model_6_aux_off[days_off, 3],
      y_model_6_aux_off[days_off, 4],
      y_model_6_aux_off[days_off, 5],
      y_model_6_aux_off[days_off, 6],
      y_model_6_aux_off[days_off, 7],
      y_model_6_aux_off[days_off, 8],
      y_model_6_aux_off[days_off, 9],
      y_model_6_aux_off[days_off, 10],
      0.75*y_model_6_aux_off[days_off, 1],
      0.75*y_model_6_aux_off[days_off, 2],
      0.75*y_model_6_aux_off[days_off, 3],
      y_model_6_aux_off[days_off, 14]]';  
      
  }
  
  //Last days on:
  if (extra_days_on > 0){
    
    //Days in quarantine
    y_model_6_aux_on = ode_rk45_tol(
      SIR_quarantine_infected, y_0_model_6, 0.0, 
      t_post_periodic_on, 1e-6, 1e-6, steps, beta_0, beta_1, beta_2, beta_3, 
      beta_12, beta_23, eta_E, gamma_0R, gamma_1R, gamma_2R, gamma_3R, 
      theta_3M, coef, p_A, p_I2, p_I3, p_M, 0.75);
      
    //Assign to model
    y_model_6[(periodic_quarantine_start + number_of_periods*(days_on + days_off) + 1):(periodic_quarantine_start + number_of_periods*(days_on + days_off) + extra_days_on)] = y_model_6_aux_on;
  }
  
  //Last days on:
  if (extra_days_off > 0){
    
    //Unquarantine everyone
    y_0_model_6 = [
      y_model_6_aux_on[extra_days_on, 1] + y_model_6_aux_on[extra_days_on, 11], 
      y_model_6_aux_on[extra_days_on, 2] + y_model_6_aux_on[extra_days_on, 12],
      y_model_6_aux_on[extra_days_on, 3] + y_model_6_aux_on[extra_days_on, 13],
      y_model_6_aux_on[extra_days_on, 4],
      y_model_6_aux_on[extra_days_on, 5],
      y_model_6_aux_on[extra_days_on, 6],
      y_model_6_aux_on[extra_days_on, 7],
      y_model_6_aux_on[extra_days_on, 8],
      y_model_6_aux_on[extra_days_on, 9],
      y_model_6_aux_on[extra_days_on, 10],
      0.0,
      0.0,
      0.0,
      y_model_6_aux_on[extra_days_on, 14]]'; 
      
    //Days off quarantine
    y_model_6_aux_off = ode_rk45_tol(
      SIR_quarantine_infected, y_0_model_6, 0.0, 
      t_post_periodic_off, 1e-6, 1e-6, steps, beta_0, beta_1, beta_2, beta_3, 
      beta_12, beta_23, eta_E, gamma_0R, gamma_1R, gamma_2R, gamma_3R, 
      theta_3M, coef, p_A, p_I2, p_I3, p_M, 0.0);  
      
    //Assign to model
    y_model_6[(periodic_quarantine_start + number_of_periods*(days_on + days_off) + 1 + extra_days_on):(periodic_quarantine_start + number_of_periods*(days_on + days_off) + extra_days)] = y_model_6_aux_off;
  }
}