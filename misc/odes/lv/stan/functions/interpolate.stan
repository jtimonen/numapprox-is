
// Linear interpolation from a time grid
real[,] interpolate_linear(vector[] y, data int[] R, data real[] A){
  // INPUT:
  // - size(R) is n
  // - size(A) is n
  // - size(y) is R[n]+2 
  // - y[1] is y0
  // - y[j] is y(t0 + (j-1)*h)
  //
  // OUTPUT:
  // - size(x) is n
  
  int n = size(R);
  int d = num_elements(y[1]);
  real x[n, d];
  for(i in 1:n){
    int R_i = R[i];
    x[i] = to_array_1d(A[i] * y[R_i+1] + (1 - A[i]) * y[R_i+2]);
  }
  return(x);
  
}


// Cubic interpolation from a time grid
real[,] interpolate_cubic(vector[] y, data int[] R, data real[] A,
    real[] theta, real h, data real[] x_r, data int[] x_i){
  // INPUT:
  // - size(R) is n
  // - size(A) is n
  // - size(y) is R[n]+2 
  // - y[1] is y0
  // - y[j] is y(t0 + (j-1)*h)
  //
  // OUTPUT:
  // - size(x) is n
  real a;
  int R_i_prev = -1;
  int n = size(R);
  int d = num_elements(y[1]);
  real x[n, d];
  vector[d] y0;
  vector[d] y1;
  vector[d] f0;
  vector[d] f1;
  vector[d] tmp;
  for(i in 1:n){
    int R_i = R[i];
    if(R_i != R_i_prev){
        y0 = y[R_i+1];
        y1 = y[R_i+2];
        f0 = odefun(0.0, y0, theta, x_r, x_i);
        f1 = odefun(0.0, y1, theta, x_r, x_i);
    }
    a = 1 - A[i];
    tmp = (1.0-2.0*a)*(y1-y0) + (a-1.0)*h*f0+a*h*f1;
    x[i] = to_array_1d((1.0-a)*y0 + a*y1 + a*(a-1.0)*tmp);
    R_i_prev = R_i;
  }
  return(x);
  
}
