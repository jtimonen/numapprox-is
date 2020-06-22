// Air pollution model from [1].
//  
// [1] Ernst Hairer and Gerhard Wanner.
//     Solving Ordinary Differential Equations II - 
//     Stiff and Differential-Algebraic Problems. Springer, 1991.
//

  real[] POLLU(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    real dy[20];
    real p[25] = theta;
    
    dy[1]  = - p[1]*y[1] - p[10]*y[11]*y[1] - p[14]*y[1]*y[6] - p[23]*y[1]*y[4] -
               p[24]*y[19]*y[1] + p[2]*y[2]*y[4] + p[3]*y[5]*y[2] + 
               p[9]*y[11]*y[2] + p[11]*y[13] + p[12]*y[10]*y[2] + 
               p[22]*y[19]+ p[25]*y[20];
    dy[2]  = - p[2]*y[2]*y[4] - p[3]*y[5]*y[2] - p[9]*y[11]*y[2] - 
               p[12]*y[10]*y[2] + p[1]*y[1] + p[21]*y[19];
    dy[3]  = - p[15]*y[3] + p[1]*y[1] + p[17]*y[4] + p[19]*y[16]+ p[22]*y[19];
    dy[4]  = - p[2]*y[2]*y[4] - p[16]*y[4] - p[17]*y[4] - p[23]*y[1]*y[4] +
               p[15]*y[3];
    dy[5]  = - p[3]*y[5]*y[2] + p[4]*y[7] + p[4]*y[7] + p[6]*y[7]*y[6] + 
               p[7]*y[9] + p[13]*y[14] + p[20]*y[17]*y[6];
    dy[6]  = - p[6]*y[7]*y[6] - p[8]*y[9]*y[6] - p[14]*y[1]*y[6] -
               p[20]*y[17]*y[6] + p[3]*y[5]*y[2] + p[18]*y[16]+ p[18]*y[16];
    dy[7]  = - p[4]*y[7] - p[5]*y[7] - p[6]*y[7]*y[6] + p[13]*y[14];
    dy[8]  =   p[4]*y[7] + p[5]*y[7] + p[6]*y[7]*y[6] + p[7]*y[9];
    dy[9]  = - p[7]*y[9] - p[8]*y[9]*y[6];
    dy[10] = - p[12]*y[10]*y[2] + p[7]*y[9] + p[9]*y[11]*y[2];
    dy[11] = - p[9]*y[11]*y[2] - p[10]*y[11]*y[1] + p[8]*y[9]*y[6] + p[11]*y[13];
    dy[12] =   p[9]*y[11]*y[2];
    dy[13] = - p[11]*y[13] + p[10]*y[11]*y[1];
    dy[14] = - p[13]*y[14] + p[12]*y[10]*y[2];
    dy[15] =   p[14]*y[1]*y[6];
    dy[16] = - p[18]*y[16] - p[19]*y[16] + p[16]*y[4];
    dy[17] = - p[20]*y[17]*y[6];
    dy[18] =   p[20]*y[17]*y[6];
    dy[19] = - p[21]*y[19] - p[22]*y[19] - p[24]*y[19]*y[1] + p[23]*y[1]*y[4] +
               p[25]*y[20];
    dy[20] = - p[25]*y[20] + p[24]*y[19]*y[1];
    return dy;
  }
  
  // Vector version
  vector odefun(real t, vector y, real[] theta, data real[] x_r, data int[] x_i){
    return to_vector(POLLU(t, to_array_1d(y), theta, x_r, x_i));
  }
