
// 4th order Runge-Kutta method
real[,] integrate_ode_rk4(real[] y0, real t0, real[] ts, real[] theta,
    data real STEP_SIZE, data int[] INTERP_R, data real[] INTERP_A, 
    data real[] x_r, data int[] x_i){
  
  int d = size(y0);
  int n = size(ts);
  real x[n, d];
  int R_n = INTERP_R[n];
  vector[d] y[R_n+2];
  real t = t0;
  
  vector[d] k1; vector[d] k2; vector[d] k3; vector[d] k4;
  
  y[1] = to_vector(y0);
  for(i in 1:(R_n+1)){
    k1 = STEP_SIZE * odefun(t,                 y[i],          theta, x_r, x_i);
    k2 = STEP_SIZE * odefun(t + 0.5*STEP_SIZE, y[i] + 0.5*k1, theta, x_r, x_i);
    k3 = STEP_SIZE * odefun(t + 0.5*STEP_SIZE, y[i] + 0.5*k2, theta, x_r, x_i);
    k4 = STEP_SIZE * odefun(t + STEP_SIZE,     y[i] + k3,     theta, x_r, x_i);
    y[i+1] = y[i] + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
    t = t + STEP_SIZE;
  }
  
  x = interpolate_cubic(y, INTERP_R, INTERP_A, theta, STEP_SIZE,
      x_r, x_i);
  return(x);
}
