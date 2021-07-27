/*==========================================================
 * c_GF.c - c-file to accelerate Green's function calculation
 *
 * performs Vector multiplications.
 * 
 *
 *		outMatrix = arrayProduct(multiplier, inMatrix)
 *
 * compile with: mex GCC='/usr/local/bin/gcc-4.7.4' COMPFLAGS='$COMPFLAGS -Wall' c_GF.c
 *
 *========================================================*/

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <time.h>
#include <complex.h>





/* The computational routine */
/* write in such a way, that it works with complex arrays */
void cfunc_3D_product ( double *points, int n_points, double i0p, double *E_p, 
        double *E_n, double *E_q, int p_max, int q_max, int rho_m, int rho_n, double _Complex C1[][p_max], 
        double _Complex C2[][rho_n], double _Complex C3[][q_max], double _Complex C4[][rho_n], 
        double _Complex rho[][rho_n],  double _Complex *G, double _Complex *GK)
{
    
    /* numbers to transmit: #omega points, # x # rho, n of C1, */
    double _Complex test; 
    /*G_test = zeros(numel(points_i), 1);
    GK_test = zeros(numel(points_i),1);*/
    int k_om, n_in, m_in, p_in, q_in;
    int loop_index = 0;
    for (k_om = 0; k_om < n_points; k_om++)
    {
        for (n_in = 0; n_in < rho_m; n_in++)
        {
            for (m_in = 0; m_in < rho_n; m_in++)
            {
                for (p_in = 0; p_in < p_max; p_in++)
                {
                   /* test = rho[n_in][m_in] * C1[m_in][p_in] * C2[p_in][n_in] / (points[k_om] + i0p*I - E_p[p_in] + E_n[m_in]);
                    printf("loop_index: %d,  addends: %f + %f * i\n", creal(test), cimag(test));
                    printf("C1: %f, C2: %f, rho: %f\n", cabs(C1[m_in][p_in] ), cabs(C2[p_in][n_in]), cabs(rho[n_in][m_in]));
                    */
                    G[k_om] = G[k_om] +  rho[n_in][m_in] * C1[m_in][p_in] * C2[p_in][n_in] / (points[k_om] + i0p*I - E_p[p_in] + E_n[m_in]);
                    GK[k_om] = GK[k_om] + -2 *I * i0p * rho[n_in][m_in] * C1[m_in][p_in] * C2[p_in][n_in] / ((points[k_om] - (E_p[p_in]-E_n[m_in]))*(points[k_om] - (E_p[p_in]-E_n[m_in])) + i0p*i0p);
                    loop_index = loop_index + 1;
                }
                
                for (q_in = 0; q_in < q_max; q_in++)
                {
                    G[k_om] = G[k_om] +  rho[n_in][m_in] * C3[m_in][q_in] * C4[q_in][n_in] / (points[k_om] + i0p*I - E_n[n_in] + E_q[q_in]);
                    GK[k_om] = GK[k_om] + 2*I * i0p * rho[n_in][m_in] * C3[m_in][q_in] * C4[q_in][n_in] /((points[k_om] - (E_n[n_in] - E_q[q_in])) * (points[k_om] - (E_n[n_in] - E_q[q_in])) + i0p*i0p);
                    loop_index = loop_index + 1;
                    
                }
            }
        }
    }
        
        
            
           
           
    
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /* input
     points:      omega points
     i0p:         imaginary infinetesimal
     *E_p:        array (vector) with energies 
     *E_n:        array (vector) with energies - equal to m-sector
     *E_q:        array (vector) with energies
     *C1 :        annihilation /creation operator 1
     *C2 :        annihilation /creation operator 2
     *C3 :        annihilation /creation operator 3
     *C4 :        annihilation /creation operator 4
     *rho:        density matrix
     *G:          Greens retarded
     *GK:         Keldysh Greens function
    */
    
    double* points = mxGetPr(prhs[0]);
    double i0p = mxGetScalar(prhs[1]);
    
    
    int p_max, q_max, rho_m, rho_n;
    
    mwSize n_points[1];
    
    *n_points = mxGetNumberOfElements(prhs[0]);
    rho_m = mxGetM(prhs[9]);
    rho_n = mxGetN(prhs[9]);
    p_max = mxGetN(prhs[5]);
    q_max = mxGetN(prhs[7]);
    
    double *E_p, *E_n, *E_q, *C1_real, *C1_imag, *C2_real, *C2_imag;
    double *C3_real, *C3_imag, *C4_real, *C4_imag, *rho_real, *rho_imag;
    
    /* get pointers of input arrays */
    E_p = mxGetPr(prhs[2]);
    E_n = mxGetPr(prhs[3]);
    E_q = mxGetPr(prhs[4]);
    C1_real = mxGetPr(prhs[5]);
    C1_imag = mxGetPi(prhs[5]);
    C2_real = mxGetPr(prhs[6]);
    C2_imag = mxGetPi(prhs[6]);
    C3_real = mxGetPr(prhs[7]);
    C3_imag = mxGetPi(prhs[7]);
    C4_real = mxGetPr(prhs[8]);
    C4_imag = mxGetPi(prhs[8]);
    rho_real = mxGetPr(prhs[9]);
    rho_imag = mxGetPi(prhs[9]);
    
    /*printf("before recasting pointers\n"); */
    
    double (*C1_r)[p_max] = (double (*)[p_max]) C1_real;
    double (*C1_i)[p_max] = (double (*)[p_max]) C1_imag;
    double (*C2_r)[rho_m] = (double (*)[rho_m]) C2_real;
    double (*C2_i)[rho_m] = (double (*)[rho_m]) C2_imag;
    double (*C3_r)[q_max] = (double (*)[q_max]) C3_real;
    double (*C3_i)[q_max] = (double (*)[q_max]) C3_imag;
    double (*C4_r)[rho_m] = (double (*)[rho_m]) C4_real;
    double (*C4_i)[rho_m] = (double (*)[rho_m]) C4_imag;   
    double (*rho_r)[rho_n] = (double (*)[rho_n]) rho_real;
    double (*rho_i)[rho_n] = (double (*)[rho_n]) rho_imag;
    
    
     
    int k,j;
    
    /*printf("before initializing complex arrays\n"); */
    /*initialize Complex arrays*/
    double _Complex C1[rho_n][p_max];
    double _Complex C2[p_max][rho_m];
    double _Complex C3[rho_n][q_max];
    double _Complex C4[q_max][rho_m];
    double _Complex rho[rho_m][rho_n];
    double _Complex G[*n_points];
    memset( G, 0, *n_points*sizeof(double _Complex) );
    
    double _Complex GK[*n_points];
    memset( GK, 0, *n_points*sizeof(double _Complex) );
    
    
  
    
    for( k = 0; k < rho_n; k++)
    {
        for( j = 0; j < p_max; j++)
        {
            C1[k][j] = C1_r[k][j] + I * C1_i[k][j];
        }
        for (j = 0; j < q_max; j++)
        {
            C3[k][j] = C3_r[k][j] + I * C3_i[k][j];
        }
    }
    
    for( k=0; k< rho_m; k++)
    {
        for(j=0;j<rho_n; j++)
        {
            rho[k][j] = rho_r[k][j] + I * rho_i[k][j];
        }
        for(j=0; j< p_max; j++)
        {
            C2[j][k] = C2_r[j][k] + I * C2_i[j][k];
        }
        for(j=0; j < q_max; j++)
        {
            C4[j][k] = C4_r[j][k] + I * C4_i[j][k];
        }
    }
    
    
    /*mxArray* E_p = mxDuplicateArray(prhs[2]);
    mxArray* E_n = mxDuplicateArray(prhs[3]);
    mxArray* E_q = mxDuplicateArray(prhs[4]);
    mxArray* C1 = mxDuplicateArray(prhs[5]);
    mxArray* C2 = mxDuplicateArray(prhs[6]);
    mxArray* C3 = mxDuplicateArray(prhs[7]);
    mxArray* C4 = mxDuplicateArray(prhs[8]);
    mxArray* rho = mxDuplicateArray(prhs[9]);
    */
    
     /* output matrix */
    

    /* check for proper number of arguments */
    
    /* make sure the first input argument is scalar */
    
    /* make sure the second input argument is type double */
    
    /* check that number of rows in second input argument is 1 */
    /*
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    
    if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0]) ||
         mxGetNumberOfElements(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Input multiplier must be a scalar.");
    }
    

    if( !mxIsDouble(prhs[1]) || 
         mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input matrix must be type double.");
    }
    
    
    if(mxGetM(prhs[1])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
    }
    */
    /* get the value of the scalar input  */
    /*multiplier = mxGetScalar(prhs[0]); */

    /* create a pointer to the real data in the input matrix  */
    /*inMatrix = mxGetPr(prhs[1]);*/

    /* get dimensions of the input matrix */
    /*ncols = mxGetN(prhs[1]);*/

    /* create the output matrix */
    /*plhs[0] = mxCreateDoubleMatrix(1,(mwSize)ncols,mxREAL);*/

    /* get a pointer to the real data in the output matrix */
    /*outMatrix = mxGetPr(plhs[0]);*/

    /* call the computational routine */
    /*arrayProduct(multiplier,inMatrix,outMatrix,(mwSize)ncols);*/
    
    /* ( double *points, int n_points, double i0p, double *E_p, 
       double *E_n, double *E_q, int p_max, int q_max, Complex *C1, Complex *C2, Complex *C3, Complex *C4, 
    Complex *rho, int rho_m, int rho_n, Complex *G, Complex *GK) */
     
     
    /* create G and GK */
     
    mxArray *Mex_G_loc_ret = mxCreateNumericArray((mwSize) 1, n_points, mxDOUBLE_CLASS, mxCOMPLEX);
    mxArray *Mex_G_loc_kel = mxCreateNumericArray((mwSize) 1, n_points, mxDOUBLE_CLASS, mxCOMPLEX);
    
    double* Real_G_loc_ret = mxGetPr(Mex_G_loc_ret);
    double* Imag_G_loc_ret = mxGetPi(Mex_G_loc_ret);
    double* Real_G_loc_kel = mxGetPr(Mex_G_loc_kel);
    double* Imag_G_loc_kel = mxGetPi(Mex_G_loc_kel);
    
    
    /* printf("before 3D product function\n"); */
   
    
    
    cfunc_3D_product (points, *n_points, i0p, E_p, E_n, E_q, p_max, q_max, rho_m, rho_n, C1, C2, C3, C4, rho, G, GK);
    /*TODO: Beschreibe die Matlab Arrays mit den Ergebnissen aus G und GK */
    /* printf("after 3D product function\n"); */
    
    for (k=0; k<*n_points; ++k)
    {
        Real_G_loc_ret[k] = creal(G[k]);
        Imag_G_loc_ret[k] = cimag(G[k]);
        Real_G_loc_kel[k] = creal(GK[k]);
        Imag_G_loc_kel[k] = cimag(GK[k]);
    }
    
    
    plhs[0] = Mex_G_loc_ret;
    plhs[1] = Mex_G_loc_kel;
    
}
