// ========================================================================
// straightforward implementations of matrix multiply Complex times Complex
// ========================================================================
void multCACB(double* Cr, double* Ci, double* Ar,double* Ai, double* Br, double* Bi, const int rA, const int cA, const int cB) {
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
void multCAtCB(double* Cr, double* Ci, double* Ar,double* Ai, double* Br, double* Bi, const int rA, const int cA, const int cB) {
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
void multCACBt(double* Cr, double* Ci, double* Ar,double* Ai, double* Br, double* Bi, const int rA, const int cA, const int rB) {
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
void multCAtCBt(double* Cr, double* Ci, double* Ar,double* Ai, double* Br, double* Bi, const int rA, const int cA, const int rB) {
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
void multCARB(double* Cr, double* Ci, double* Ar,double* Ai, double* Br, const int rA, const int cA, const int cB) {
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
void multCAtRB(double* Cr, double* Ci, double* Ar,double* Ai, double* Br, const int rA, const int cA, const int cB) {
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
void multCARBt(double* Cr, double* Ci, double* Ar,double* Ai, double* Br, const int rA, const int cA, const int rB) {
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
void multCAtRBt(double* Cr, double* Ci, double* Ar,double* Ai, double* Br, const int rA, const int cA, const int rB) {
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
void multRACB(double* Cr, double* Ci, double* Ar, double* Br, double* Bi, const int rA, const int cA, const int cB) {
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
void multRAtCB(double* Cr, double* Ci, double* Ar, double* Br, double* Bi, const int rA, const int cA, const int cB) {
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
void multRACBt(double* Cr, double* Ci, double* Ar, double* Br, double* Bi, const int rA, const int cA, const int rB) {
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
void multRAtCBt(double* Cr, double* Ci, double* Ar, double* Br, double* Bi, const int rA, const int cA, const int rB) {
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
void multRARB(double* Cr, double* Ar, double* Br, const int rA, const int cA, const int cB) {
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
void multRAtRB(double* Cr, double* Ar, double* Br, const int rA, const int cA, const int cB) {
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
void multRARBt(double* Cr, double* Ar, double* Br, const int rA, const int cA, const int rB) {
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
void multRAtRBt(double* Cr, double* Ar, double* Br, const int rA, const int cA, const int rB) {
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


// ===========================================================
// multiply:   C = op(A) * op(B)  with A complex and B complex
// ===========================================================
void mulCMatCMat(double* Cr, double* Ci, double* Ar, double* Ai, double* Br, double* Bi,
				   const int rA, const int cA, const int rB, const int cB, const char *mod) {

	if ( (mod[0] == 'N') && (mod[1] == 'N') )
		multCACB(Cr,Ci,Ar,Ai,Br,Bi,rA, cA, cB);
	else if ( (mod[0] == 'T') && (mod[1] == 'N') )
		multCAtCB(Cr,Ci,Ar,Ai,Br,Bi,rA, cA, cB);
	else if ( (mod[0] == 'N') && (mod[1] == 'T') )
		multCACBt(Cr,Ci,Ar,Ai,Br,Bi,rA, cA, cB);
	else if ( (mod[0] == 'T') && (mod[1] == 'T') )
		multCAtCBt(Cr,Ci,Ar,Ai,Br,Bi,rA, cA, cB);
}
// ===========================================================
// multiply:   C = op(A) * op(B)  with A complex and B real
// ===========================================================
void mulCMatRMat(double* Cr, double* Ci, double* Ar, double* Ai, double* Br,
				   const int rA, const int cA, const int rB, const int cB, const char *mod) {

	if ( (mod[0] == 'N') && (mod[1] == 'N') )
		multCARB(Cr,Ci,Ar,Ai,Br,rA, cA, cB);
	else if ( (mod[0] == 'T') && (mod[1] == 'N') )
		multCAtRB(Cr,Ci,Ar,Ai,Br,rA, cA, cB);
	else if ( (mod[0] == 'N') && (mod[1] == 'T') )
		multCARBt(Cr,Ci,Ar,Ai,Br,rA, cA, cB);
	else if ( (mod[0] == 'T') && (mod[1] == 'T') )
		multCAtRBt(Cr,Ci,Ar,Ai,Br,rA, cA, cB);
}
// ===========================================================
// multiply:   C = op(A) * op(B)  with A real and B complex
// ===========================================================
void mulRMatCMat(double* Cr, double* Ci, double* Ar, double* Br, double* Bi,
				   const int rA, const int cA, const int rB, const int cB, const char *mod) {

	if ( (mod[0] == 'N') && (mod[1] == 'N') )
		multRACB(Cr,Ci,Ar,Br,Bi,rA, cA, cB);
	else if ( (mod[0] == 'T') && (mod[1] == 'N') )
		multRAtCB(Cr,Ci,Ar,Br,Bi,rA, cA, cB);
	else if ( (mod[0] == 'N') && (mod[1] == 'T') )
		multRACBt(Cr,Ci,Ar,Br,Bi,rA, cA, cB);
	else if ( (mod[0] == 'T') && (mod[1] == 'T') )
		multRAtCBt(Cr,Ci,Ar,Br,Bi,rA, cA, cB);
}
// ===========================================================
// multiply:   C = op(A) * op(B)  with A real and B real
// ===========================================================
void mulRMatRMat(double* Cr, double* Ar, double* Br,
				   const int rA, const int cA, const int rB, const int cB, const char *mod) {

	if ( (mod[0] == 'N') && (mod[1] == 'N') )
		multRARB(Cr,Ar,Br,rA,cA,cB);
	else if ( (mod[0] == 'T') && (mod[1] == 'N') )
		multRAtRB(Cr,Ar,Br,rA,cA,cB);
	else if ( (mod[0] == 'N') && (mod[1] == 'T') )
		multRARBt(Cr,Ar,Br,rA,cA,cB);
	else if ( (mod[0] == 'T') && (mod[1] == 'T') )
		multRAtRBt(Cr,Ar,Br,rA,cA,cB);
}

/*
// ================================================
// square:   C = A * op(A)  or  C = 0.5*(A*B'+B*A')
// ================================================
void squareMatMat(double* C, double* A, double* B,
				   const int rA, const int cA, const int rB, const int cB, const char *mod) {
   // can't pass consts to BLAS
   ptrdiff_t rA0 = rA, cA0 = cA, rB0 = rB; 
   // rows(Op(A)), columns(Op(A)), columns(Op(B)), rows(C)
   ptrdiff_t copA, rC;  
   int i,j; 
   double temp; 

   if ( (mod[0] == 'N') ){
      copA  = cA;
      rC    = rA;
   } else {
      copA  = rA;
      rC    = cA;   
   } 

#ifndef USE_BLAS // naive C implementations

   if ((rB == 0) || (cB == 0)){  // one input  C = A*A'   
      if ( (mod[0] == 'N') )
         multABt(C, A, A, rA, cA, rA);
      else
         multAtB(C, A, A, rA, cA, cA);
   }else{
      if ( (mod[0] == 'N') )
         multABt(C, A, B, rA, cA, rB);
      else
         multAtB(C, A, B, rA, cA, cB);

      // symmetrize
      for( i=0; i<rC; i++ )
         for( j=i; j<rC; j++ ){
            temp = C[i*rC+j] + C[j*rC+i];
            C[i*rC+j] = C[j*rC+i] = 0.5*temp;   
         }
   }

#else
   char  modA = mod[0], modB = mod[1], uplo = 'U';
   double one = 1.0, zero = 0.0, half = 0.5;

   if ((!rB) && (!cB))  // one input  C = A*A'
      dsyrk(&uplo, &modA, &rC, &copA, &one, A, &rA0, &zero, C, &rC);
   else                 // two inputs C = 0.5*(A*B'+B*A')
      dsyr2k(&uplo, &modA, &rC, &copA, &half, A, &rA0, B, &rB0, &zero, C, &rC);   

   // symmetrize
   for( i=0; i<rC; i++ )
      for( j=i+1; j<rC; j++ )
          C[i*rC+j] = C[j*rC+i];

#endif
}

// =====================================
// cholesky decomposition:   C = chol(A)
// =====================================
double dot(const double* vec1, const double* vec2, const int n)
{
	int i;
	double res = 0;
	
   for( i=0; i<n; i++ )
      res += vec1[i] * vec2[i];
   
	return res;
}


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



void chol(double* C, double* A,  const int rA) {
   int i,j, rank;
   double temp;

   // copy upper triangle into C
   for( i=0; i<rA; i++ )
      for( j=0; j<=i; j++ )
          C[i*rA+j] = A[i*rA+j];    
   
#ifndef USE_BLAS // naive C implementations
   temp = A[0];
   rank = cholA(C, A, rA);
   // chol used A as scratch, now fix it
   if (rank) A[0] = temp;
   for( i=1; i<rank; i++ )
          A[i] = A[i*rA];     
   //if decomposition failed put -1 in C
   if (rank < rA) C[0] = -1; 
#else
   ptrdiff_t rA0 = rA;
   ptrdiff_t info;
   dpotrf("U", &rA0, C, &rA0, &info );
#endif
}


// ================================
// solve linear equations   C = A\B
// ================================
void solve(double* C, double* A, double* B,
				   const int rA, const int cA, const int rB, const int cB, 
               const char *mod, double *W, const int LW, ptrdiff_t *S) {
#ifdef USE_BLAS
   int i, j, rank;
   char  uplo = 'U', side = 'L', trans = 'N', unit = 'N';
   double one = 1.0, rcond = 0.000000001;
   ptrdiff_t rA0 = rA, cA0 = cA,  cB0 = cB, Lwork=LW, info;
   ptrdiff_t rC0 = (rA>cA) ? rA : cA;
   //ptrdiff_t ptr_S = S;

   switch (mod[0]){
      case 'L':
      case 'U':
         uplo = mod[0];
         dtrsm(&side, &uplo, &trans, &unit, &rC0, &cB0, &one, A, &rA0, C, &rC0);
         break;
      case 'P':
         dposv(&uplo, &rA0, &cB0, W, &rA0, C,  &rA0, &info);// A has already been copied into W
         break;
      default:
         if (rA == cA) {
            //dgesv(&rA0, &cB0, W, &rA0, S, C,  &rA0, &info);// A has already been copied into W
            dgesv(&rA0, &cB0, W, &rA0, (ptrdiff_t*)S, C,  &rA0, &info);// A has already been copied into W
         }
         else{
            for( i=0; i<cB; i++ )
               for( j=0; j<rB; j++ )
                  C[i*rC0+j] = B[i*rB+j];
            //dgelsy(&rA0, &cA0, &cB0, A, &rA0, C, &rC0, S, &rcond, &rank, W, &Lwork, &info);
            dgelsy(&rA0, &cA0, &cB0, A, &rA0, C, &rC0,
                    (ptrdiff_t*)S, &rcond, (ptrdiff_t*)&rank, W, &Lwork, &info);
         }
   }
#endif
}
*/