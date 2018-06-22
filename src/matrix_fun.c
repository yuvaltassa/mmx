// ========================================================================
// straightforward implementations of matrix multiply Complex times Complex
// ========================================================================
void multCC(double* Cr, double* Ci, double* Ar,double* Ai, double* Br, double* Bi, const int rA, const int cA, const int cB) {
   int i, j, k;
   double *ar,*cr,tmpr;
   double *ai,*ci,tmpi;
   for( i=0; i<cB; i++ ){
      cr   = Cr + i*rA;
        ci   = Ci + i*rA;
      for( k=0; k<cA; k++ ){
         tmpr = Br[i*cA+k];
         tmpi = Bi[i*cA+k];
         ar   = Ar + k*rA;
         ai   = Ai + k*rA;
         for( j=0; j<rA; j++ ){
            cr[j] += tmpr * ar[j] - tmpi * ai[j];
            ci[j] += tmpi * ar[j] + tmpr * ai[j];
         }
      }
   }   
}
void multCtC(double* Cr, double* Ci, double* Ar,double* Ai, double* Br, double* Bi, const int rA, const int cA, const int cB) {
   int i, j, k;
   double *ar, *bi, *cr;  
   double *ai, *br, *ci;  
   for( i=0; i<cB; i++ ){
      br = Br + i*rA;
      bi = Bi + i*rA;
      for( k=0; k<cA; k++ ){
         ar = Ar + k*rA;
         ai = Ai + k*rA;
         cr = Cr + i*cA + k;
         ci = Ci + i*cA + k;
         for( j=0; j<rA; j++ ){
            (*cr) += ar[j]*br[j] - ai[j]*bi[j];
            (*ci) += ai[j]*br[j] + ar[j]*bi[j];
         }
      }
   }
}
void multCCt(double* Cr, double* Ci, double* Ar,double* Ai, double* Br, double* Bi, const int rA, const int cA, const int rB) {
   int i, j, k;
   double *ar, *br;
   double *ai, *bi;
   for( j=0; j<cA; j++ ){
      ar = Ar + j*rA;
      ai = Ai + j*rA;
      br = Br + j*rB;
      bi = Bi + j*rB;      
      for( i=0; i<rB; i++ ){
         for( k=0; k<rA; k++ ){
            Cr[i*rA + k] += ar[k]*br[i] - ai[k]*bi[i];
            Ci[i*rA + k] += ai[k]*br[i] + ar[k]*bi[i];
         }
      }
   }
}
void multCtCt(double* Cr, double* Ci, double* Ar,double* Ai, double* Br, double* Bi, const int rA, const int cA, const int rB) {
   int i, j, k;
   double *br, *cr, tmpr;  
   double *bi, *ci, tmpi;  
   for( i=0; i<cA; i++ ){
      cr   = Cr + i;
      ci   = Ci + i;
      for( k=0; k<rA; k++ ){
         tmpr = Ar[i*rA+k];
         tmpi = Ai[i*rA+k];
         br   = Br + k*rB;
         bi   = Bi + k*rB;
         for( j=0; j<rB; j++ ){
            cr[j*cA] += tmpr * br[j] - tmpi * bi[j];
            ci[j*cA] += tmpi * br[j] + tmpr * bi[j];
         }
      }
   }
}
// ========================================================================
// straightforward implementations of matrix multiply Complex times Real
// ========================================================================
void multCR(double* Cr, double* Ci, double* Ar,double* Ai, double* Br, const int rA, const int cA, const int cB) {
   int i, j, k;
   double *ar,*cr,tmpr;
   double *ai,*ci;
   for( i=0; i<cB; i++ ){
      cr   = Cr + i*rA;
        ci   = Ci + i*rA;
      for( k=0; k<cA; k++ ){
         tmpr = Br[i*cA+k];
         ar   = Ar + k*rA;
         ai   = Ai + k*rA;
         for( j=0; j<rA; j++ ){
            cr[j] += tmpr * ar[j];
            ci[j] += tmpr * ai[j];
         }
      }
   }   
}
void multCtR(double* Cr, double* Ci, double* Ar,double* Ai, double* Br, const int rA, const int cA, const int cB) {
   int i, j, k;
   double *ar, *br, *cr;  
   double *ai,      *ci;  
   for( i=0; i<cB; i++ ){
      br = Br + i*rA;
      for( k=0; k<cA; k++ ){
         ar = Ar + k*rA;
         ai = Ai + k*rA;
         cr = Cr + i*cA + k;
         ci = Ci + i*cA + k;
         for( j=0; j<rA; j++ ){
            (*cr) += ar[j]*br[j];
            (*ci) += ai[j]*br[j];
         }
      }
   }
}
void multCRt(double* Cr, double* Ci, double* Ar,double* Ai, double* Br, const int rA, const int cA, const int rB) {
   int i, j, k;
   double *ar, *br;
   double *ai;
   for( j=0; j<cA; j++ ){
      ar = Ar + j*rA;
      ai = Ai + j*rA;
      br = Br + j*rB;
      for( i=0; i<rB; i++ ){
         for( k=0; k<rA; k++ ){
            Cr[i*rA + k] += ar[k]*br[i];
            Ci[i*rA + k] += ai[k]*br[i];
         }
      }
   }
}
void multCtRt(double* Cr, double* Ci, double* Ar,double* Ai, double* Br, const int rA, const int cA, const int rB) {
   int i, j, k;
   double *br, *cr, tmpr;  
   double      *ci, tmpi;  
   for( i=0; i<cA; i++ ){
      cr   = Cr + i;
      ci   = Ci + i;
      for( k=0; k<rA; k++ ){
         tmpr = Ar[i*rA+k];
         tmpi = Ai[i*rA+k];
         br   = Br + k*rB;
         for( j=0; j<rB; j++ ){
            cr[j*cA] += tmpr * br[j];
            ci[j*cA] += tmpi * br[j];
         }
      }
   }
}
// ========================================================================
// straightforward implementations of matrix multiply Real times Complex
// ========================================================================
void multRC(double* Cr, double* Ci, double* Ar, double* Br, double* Bi, const int rA, const int cA, const int cB) {
   int i, j, k;
   double *ar,*cr,tmpr;
   double     *ci,tmpi;
   for( i=0; i<cB; i++ ){
      cr   = Cr + i*rA;
        ci   = Ci + i*rA;
      for( k=0; k<cA; k++ ){
         tmpr = Br[i*cA+k];
         tmpi = Bi[i*cA+k];
         ar   = Ar + k*rA;
         for( j=0; j<rA; j++ ){
            cr[j] += tmpr * ar[j];
            ci[j] += tmpi * ar[j];
         }
      }
   }   
}
void multRtC(double* Cr, double* Ci, double* Ar, double* Br, double* Bi, const int rA, const int cA, const int cB) {
   int i, j, k;
   double *ar, *bi, *cr;  
   double      *br, *ci;  
   for( i=0; i<cB; i++ ){
      br = Br + i*rA;
      bi = Bi + i*rA;
      for( k=0; k<cA; k++ ){
         ar = Ar + k*rA;
         cr = Cr + i*cA + k;
         ci = Ci + i*cA + k;
         for( j=0; j<rA; j++ ){
            (*cr) += ar[j]*br[j];
            (*ci) += ar[j]*bi[j];
         }
      }
   }
}
void multRCt(double* Cr, double* Ci, double* Ar, double* Br, double* Bi, const int rA, const int cA, const int rB) {
   int i, j, k;
   double *ar, *br;
   double      *bi;
   for( j=0; j<cA; j++ ){
      ar = Ar + j*rA;
      br = Br + j*rB;
      bi = Bi + j*rB;      
      for( i=0; i<rB; i++ ){
         for( k=0; k<rA; k++ ){
            Cr[i*rA + k] += ar[k]*br[i];
            Ci[i*rA + k] += ar[k]*bi[i];
         }
      }
   }
}
void multRtCt(double* Cr, double* Ci, double* Ar, double* Br, double* Bi, const int rA, const int cA, const int rB) {
   int i, j, k;
   double *br, *cr, tmpr;  
   double *bi, *ci;  
   for( i=0; i<cA; i++ ){
      cr   = Cr + i;
      ci   = Ci + i;
      for( k=0; k<rA; k++ ){
         tmpr = Ar[i*rA+k];
         br   = Br + k*rB;
         bi   = Bi + k*rB;
         for( j=0; j<rB; j++ ){
            cr[j*cA] += tmpr * br[j];
            ci[j*cA] += tmpr * bi[j];
         }
      }
   }
}
// ========================================================================
// straightforward implementations of matrix multiply Real times Real
// ========================================================================
void multRR(double* Cr, double* Ar, double* Br, const int rA, const int cA, const int cB) {
   int i, j, k;
   double *ar,*cr,tmpr;
   for( i=0; i<cB; i++ ){
      cr   = Cr + i*rA;
      for( k=0; k<cA; k++ ){
         tmpr = Br[i*cA+k];
         ar   = Ar + k*rA;
         for( j=0; j<rA; j++ ){
            cr[j] += tmpr * ar[j];
         }
      }
   }   
}
void multRtR(double* Cr, double* Ar, double* Br, const int rA, const int cA, const int cB) {
   int i, j, k;
   double *ar, *br, *cr;  
   for( i=0; i<cB; i++ ){
      br = Br + i*rA;
      for( k=0; k<cA; k++ ){
         ar = Ar + k*rA;
         cr = Cr + i*cA + k;
         for( j=0; j<rA; j++ ){
            (*cr) += ar[j]*br[j];
         }
      }
   }
}
void multRRt(double* Cr, double* Ar, double* Br, const int rA, const int cA, const int rB) {
   int i, j, k;
   double *ar, *br;
   for( j=0; j<cA; j++ ){
      ar = Ar + j*rA;
      br = Br + j*rB;
      for( i=0; i<rB; i++ ){
         for( k=0; k<rA; k++ ){
            Cr[i*rA + k] += ar[k]*br[i];
         }
      }
   }
}
void multRtRt(double* Cr, double* Ar, double* Br, const int rA, const int cA, const int rB) {
   int i, j, k;
   double *br, *cr, tmpr;  
   for( i=0; i<cA; i++ ){
      cr   = Cr + i;
      for( k=0; k<rA; k++ ){
         tmpr = Ar[i*rA+k];
         br   = Br + k*rB;
         for( j=0; j<rB; j++ ){
            cr[j*cA] += tmpr * br[j];
         }
      }
   }
}


void scl(double* res, const double* vec, double scl, int n)
{
   int i;
   for( i=0; i<n; i++ )
      res[i] += vec[i] * scl;
}


// same as multAB, but exploits sparsity of B
void multABs(double* C, double* A, double* B, const int a1, const int a2, const int b2) {
   int i, j, k;
   double *a, *c, tmp;
   for( i=0; i<b2; i++ ){
      c  = C + i*a1;
      for( k=0; k<a2; k++ ){
         a = A + k*a1;
         if ((tmp = B[i*a2+k]))
            scl(c, a, tmp, a1);
      }
   }
}


double dot(const double* vec1, const double* vec2, const int n){
   int i;
   double res = 0;
   for( i=0; i<n; i++ ) res += vec1[i] * vec2[i];
   return res;
}

// same as multAB, but exploits sparsity of A
void multAsB(double* C, double* A, double* B, const int a1, const int a2, const int b2) {
   int i, j, k;
   double *a, *b, *c, tmp;
   for( i=0; i<a2; i++ ){
      b  = B+i;
      a  = A + i*a1;
      for( k=0; k<a1; k++ ){
         c  = C+k;
         if ((tmp = a[k])){
            for( j=0; j<b2; j++ ){
               c[j*a1] += tmp * b[j*a2];
            }
         }
      }
   }
}

// =============================
// multiply:   C = op(A) * op(B)
// =============================

void mulCC(   double* Cr, double* Ci,
         double* Ar, double* Ai,
         double* Br, double* Bi,
         const int rA, const int cA, 
         const int rB, const int cB, 
         const char *mod) {
//#ifndef USE_BLAS            
   if ( (mod[0] == 'N') && (mod[1] == 'N') )
      multCC(Cr,Ci,Ar,Ai,Br,Bi,rA, cA, cB);
   else if ( (mod[0] == 'T') && (mod[1] == 'N') )
      multCtC(Cr,Ci,Ar,Ai,Br,Bi,rA, cA, cB);
   else if ( (mod[0] == 'N') && (mod[1] == 'T') )
      multCCt(Cr,Ci,Ar,Ai,Br,Bi,rA, cA, cB);
   else if ( (mod[0] == 'T') && (mod[1] == 'T') )
      multCtCt(Cr,Ci,Ar,Ai,Br,Bi,rA, cA, cB);
//#else
   //TODO: Complex BLAs call
//#endif   
}

void mulCR(   double* Cr, double* Ci,
         double* Ar, double* Ai,
         double* Br,
         const int rA, const int cA, 
         const int rB, const int cB, 
         const char *mod) {
//#ifndef USE_BLAS            
   if ( (mod[0] == 'N') && (mod[1] == 'N') )
      multCR(Cr,Ci,Ar,Ai,Br,rA, cA, cB);
   else if ( (mod[0] == 'T') && (mod[1] == 'N') )
      multCtR(Cr,Ci,Ar,Ai,Br,rA, cA, cB);
   else if ( (mod[0] == 'N') && (mod[1] == 'T') )
      multCRt(Cr,Ci,Ar,Ai,Br,rA, cA, cB);
   else if ( (mod[0] == 'T') && (mod[1] == 'T') )
      multCtRt(Cr,Ci,Ar,Ai,Br,rA, cA, cB);
//#else
   //TODO: Complex BLAs call
//#endif               
}

void mulRC(   double* Cr, double* Ci,
         double* Ar,
         double* Br, double* Bi,
         const int rA, const int cA, 
         const int rB, const int cB, 
         const char *mod) {
//#ifndef USE_BLAS            
   if ( (mod[0] == 'N') && (mod[1] == 'N') )
      multRC(Cr,Ci,Ar,Br,Bi,rA, cA, cB);
   else if ( (mod[0] == 'T') && (mod[1] == 'N') )
      multRtC(Cr,Ci,Ar,Br,Bi,rA, cA, cB);
   else if ( (mod[0] == 'N') && (mod[1] == 'T') )
      multRCt(Cr,Ci,Ar,Br,Bi,rA, cA, cB);
   else if ( (mod[0] == 'T') && (mod[1] == 'T') )
      multRtCt(Cr,Ci,Ar,Br,Bi,rA, cA, cB);
//#else
   //TODO: Complex BLAs call
//#endif      
}

void mulRR(   double* C, 
         double* A, 
         double* B,
         const int a1, const int a2, 
         const int b1, const int b2, 
         const char *mod) {
#ifndef USE_BLAS // naive C implementations, including "half-sparse"

   if ((mod[0] != 'S') && (mod[1] != 'S')){
      if ( (mod[0] == 'N') && (mod[1] == 'N') )
         multRR(C, A, B,a1, a2, b2);
      else if ( (mod[0] == 'T') && (mod[1] == 'N') )
         multRtR(C, A, B, a1, a2, b2);
      else if ( (mod[0] == 'N') && (mod[1] == 'T') )
         multRRt(C, A, B, a1, a2, b1);
      else if ( (mod[0] == 'T') && (mod[1] == 'T') )
         multRtRt(C, A, B, a1, a2, b1);
   } else {  
      if (mod[0] == 'S')
         multAsB(C, A, B, a1, a2, b2);
      else
         multABs(C, A, B, a1, a2, b2);
   }      

#else

   // rows(Op(A)), columns(Op(A)), columns(Op(B)), rows(C)
   ptrdiff_t opA1, opA2, opB2, c1;  
   // can't pass consts to fortran
   ptrdiff_t a10 = a1, b10 = b1;    

   char modA = mod[0], modB = mod[1];
   double one = 1.0, zero = 0.0;

   if ( (mod[0] == 'N') && (mod[1] == 'N') ){
      opA1  = a1;
      opA2  = a2;   
      opB2  = b2;
      c1    = a1;
   } else if ( (mod[0] == 'T') && (mod[1] == 'N') ){
      opA1  = a2;
      opA2  = a1;   
      opB2  = b2;
      c1    = a2;   
   } else if ( (mod[0] == 'N') && (mod[1] == 'T') ){
      opA1  = a1;
      opA2  = a2;   
      opB2  = b1;
      c1    = a1;   
   } else if ( (mod[0] == 'T') && (mod[1] == 'T') ){
      opA1  = a2;
      opA2  = a1;   
      opB2  = b1;
      c1    = a2;   
   }
   dgemm(&modA, &modB, &opA1, &opB2, &opA2, &one, A, &a10, B, &b10, &zero, C, &c1);

#endif
}

// ================================================
// square:   C = A * op(A)  or  C = 0.5*(A*B'+B*A')
// ================================================
void squareCC(   double* Creal, double* Cimag, 
            double* Areal, double* Aimag, 
            double* Breal, double* Bimag, 
            const int a1, const int a2, 
            const int b1, const int b2, 
            const char *mod) {
               
            }
            
void squareCR(   double* Creal, double* Cimag, 
            double* Areal, double* Aimag, 
            double* Breal,
            const int a1, const int a2, 
            const int b1, const int b2, 
            const char *mod) {
               
            }
            
void squareRC(   double* Creal, double* Cimag, 
            double* Areal, 
            double* Breal, double* Bimag, 
            const int a1, const int a2, 
            const int b1, const int b2, 
            const char *mod) {
               
            }


void squareRR(   double* C, 
            double* A, 
            double* B,
            const int a1, const int a2, 
            const int b1, const int b2, 
            const char *mod) {
   // can't pass consts to BLAS
   ptrdiff_t a10 = a1, a20 = a2, b10 = b1; 
   // rows(Op(A)), columns(Op(A)), columns(Op(B)), rows(C)
   ptrdiff_t opA2, c1;  
   int i,j; 

   if ( (mod[0] == 'N') ){
      opA2  = a2;
      c1    = a1;
   } else {
      opA2  = a1;
      c1    = a2;   
   } 

#ifndef USE_BLAS // naive C implementations

   if ((b1 == 0) || (b2 == 0)){  // one input  C = A*A'   
      if ( (mod[0] == 'N') )
         multRRt(C, A, A, a1, a2, a1);
      else
         multRtR(C, A, A, a1, a2, a2);
   }else{
      if ( (mod[0] == 'N') )
         multRRt(C, A, B, a1, a2, b1);
      else
         multRtR(C, A, B, a1, a2, b2);

      double temp;
      // symmetrize
      for( i=0; i<c1; i++ )
         for( j=i; j<c1; j++ ){
            temp = C[i*c1+j] + C[j*c1+i];
            C[i*c1+j] = C[j*c1+i] = 0.5*temp;   
         }
   }

#else
   char  modA = mod[0], modB = mod[1], uplo = 'U';
   double one = 1.0, zero = 0.0, half = 0.5;

   if ((!b1) && (!b2))  // one input  C = A*A'
      dsyrk(&uplo, &modA, &c1, &opA2, &one, A, &a10, &zero, C, &c1);
   else                 // two inputs C = 0.5*(A*B'+B*A')
      dsyr2k(&uplo, &modA, &c1, &opA2, &half, A, &a10, B, &b10, &zero, C, &c1);   

   // symmetrize
   for( i=0; i<c1; i++ )
      for( j=i+1; j<c1; j++ )
          C[i*c1+j] = C[j*c1+i];

#endif
}

// =====================================
// cholesky decomposition:   C = chol(A)
// =====================================
int cholA(double* A, double* scratch, const int n)
{
   int i, j, rank=0;
   double tmp;

   // in-place Cholesky factorization, store 1/L(j,j) in scratch
   for( j=0; j<n; j++ )
   {
      tmp = A[j*n+j];
      if( j )
         tmp -= dot(A+j*n, A+j*n, j);

      if( tmp < 0.000000001 )
         return rank;
      else
      {
         scratch[j] = (double)(1.0/sqrt(tmp));
         rank++;
      }

      // process off-diagonal entries, modify 'A'
      for( i=j+1; i<n; i++ )
      {
         A[i*n+j] -= dot(A+i*n, A+j*n, j);
         A[i*n+j] *= scratch[j];
      }
   }

   // copy 'scratch' to diagonal of A
   for( j=0; j<n; j++ )
      A[j*n+j] = 1./scratch[j];

   return rank;
}



void cholR(double* C, double* A,  const int a1) {
   int i,j;

   // copy upper triangle into C
   for( i=0; i<a1; i++ )
      for( j=0; j<=i; j++ )
          C[i*a1+j] = A[i*a1+j];    
   
#ifndef USE_BLAS // naive C implementations
   double temp = A[0];
   int rank = cholA(C, A, a1);
   // chol used A as scratch, now fix it
   if (rank) A[0] = temp;
   for( i=1; i<rank; i++ )
          A[i] = A[i*a1];     
   //if decomposition failed put -1 in C
   if (rank < a1) C[0] = -1; 
#else
   ptrdiff_t a10 = a1;
   ptrdiff_t info;
   dpotrf("U", &a10, C, &a10, &info );
#endif
}

void cholC(   double* Creal, double* Cimag,
         double* Areal, double* Aimag,
         const int a1) {
   
}


// ================================
// solve linear equations   C = A\B
// ================================
void solveRR(   double* C, 
            double* A, 
            double* B,
            const int a1, const int a2, 
            const int b1, const int b2, 
               const char *mod, double *W, const int LW, ptrdiff_t *S) {
#ifdef USE_BLAS
   int i, j, rank;
   char  uplo = 'U', side = 'L', trans = 'N', unit = 'N';
   double one = 1.0, rcond = 0.000000001;
   ptrdiff_t a10 = a1, a20 = a2,  b20 = b2, Lwork=LW, info;
   ptrdiff_t c10 = (a1>a2) ? a1 : a2;
   //ptrdiff_t ptr_S = S;

   switch (mod[0]){
      case 'L':
      case 'U':
         uplo = mod[0];
         dtrsm(&side, &uplo, &trans, &unit, &c10, &b20, &one, A, &a10, C, &c10);
         break;
      case 'P':
         dposv(&uplo, &a10, &b20, W, &a10, C,  &a10, &info);// A has already been copied into W
         break;
      default:
         if (a1 == a2) {
            //dgesv(&a10, &b20, W, &a10, S, C,  &a10, &info);// A has already been copied into W
            dgesv(&a10, &b20, W, &a10, (ptrdiff_t*)S, C,  &a10, &info);// A has already been copied into W
         }
         else{
            for( i=0; i<b2; i++ )
               for( j=0; j<b1; j++ )
                  C[i*c10+j] = B[i*b1+j];
            dgelsy(&a10, &a20, &b20, A, &a10, C, &c10,
                    (ptrdiff_t*)S, &rcond, (ptrdiff_t*)&rank, W, &Lwork, &info);
         }
   }
#endif
}

void solveCR(   double* Creal, double* Cimag, 
            double* Areal, double* Aimag, 
            double* Breal,
            const int a1, const int a2, 
            const int b1, const int b2, 
            const char *mod, double *W, const int LW, ptrdiff_t *S) {
}

void solveRC(   double* Creal, double* Cimag, 
            double* Areal,
            double* Breal, double* Bimag,
            const int a1, const int a2, 
            const int b1, const int b2, 
            const char *mod, double *W, const int LW, ptrdiff_t *S) {
}

void solveCC(   double* Creal, double* Cimag, 
            double* Areal, double* Aimag, 
            double* Breal, double* Bimag,
            const int a1, const int a2, 
            const int b1, const int b2, 
            const char *mod, double *W, const int LW, ptrdiff_t *S) {
}
