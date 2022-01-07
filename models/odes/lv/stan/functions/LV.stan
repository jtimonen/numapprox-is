  // Lotka-Volterra system
  real[] LV(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    real dydt[2];
    dydt[1] = theta[1]*y[1] - y[1]*y[2];
    dydt[2] = y[1]*y[2] - theta[2]*y[2];
    return dydt;
  }
  
  // Vector version of LV
  vector odefun(real t, vector y, real[] theta, data real[] x_r, data int[] x_i){
    return to_vector(LV(t, to_array_1d(y), theta, x_r, x_i));
  }