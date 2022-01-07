// Linear interpolation from equispaced time grid
real[,] interpolate(vector[] y, data int[] R, data real[] A){
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