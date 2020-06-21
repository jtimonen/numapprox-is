  // a pharmacokinetic model
  real[] PKPD(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    // STATE VARIABLES
    // - y[1] = Depot
    // - y[2] = Cent
    // - y[3] = Periph_1
    // - y[4] = Periph_2
    // - y[5] = Resp
    //
    // PARAMETERS
    // - theta[1] = k_a
    // - theta[2] = CL
    // - theta[3] = k_m
    // - theta[4] = k_in
    // - theta[5] = k_out
    real dydt[5];
    real V_c = 20;
    real V_p1 = 10;
    real V_p2 = 100;
    real V_max = 0;
    real Q1 = 2;
    real Q2 = 0.5;
    real I_max = 1;
    real IC50 = 2;
    real g = 1;
    
    real k_a = theta[1];
    real CL = theta[2];
    real k_m = theta[3];
    real k_in = theta[4];
    real k_out = theta[5];
    
    real z2 = y[2]/V_c;
    real z3 = y[3]/V_p1;
    real z4 = y[4]/V_p2;
    
    dydt[1] = -k_a*y[1];
    dydt[2] = k_a*y[1] + (CL+V_max/(k_m+z2)+Q1)*z2 + Q1*z3 - Q2*z2 + Q2*z4;
    dydt[3] = Q1*z2 - Q1*z3;
    dydt[4] = Q2*z2 - Q2*z4;
    dydt[5] = k_in*(1-(I_max*z2^g)/(IC50^g+z2^g)) - k_out*y[5];
    return dydt;
  }
  
  // Vector version
  vector odefun(real t, vector y, real[] theta, data real[] x_r, data int[] x_i){
    return to_vector(PKPD(t, to_array_1d(y), theta, x_r, x_i));
  }