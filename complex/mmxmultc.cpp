// =======================
// preprocessor directives
// =======================
#include "mex.h"
#include "string.h"
#include <ctype.h>

#ifdef UNIX_SYSTEM
   #include <unistd.h>
   #include <pthread.h>
#endif

#ifdef WIN_SYSTEM
   #include "windows.h"
#endif

#include <math.h>

// #if defined(USE_BLAS)
   // #if defined(MKL_ILP64) or defined(MKL_32)
      // #include "mkl_blas.h"
      // #include "mkl_lapack.h"
      // #define ptrdiff_t MKL_INT
   // #else 
      // #ifdef UNIX_SYSTEM
        // #include "blas.h"
        // #include "lapack.h"
      // #else
        // #include "my_blas.h"
      // #endif
   // #endif
// #endif

#define MAX_THREAD 64

// list of possible values for PARTASK
// positive values are matrix oprations
#define MATMUL 1
#define SQUARE 2
#define CHOL   3
#define BSLASH 4
// negative values are binary-elementwise functions (like bsxfun)

// function declarations
#include "matrix_fun_complex.c"

//#define DEBUG

// ================
// global variables
// ================

// thread related
static bool    INITIALIZED = false;
static int     SCHEDULE[MAX_THREAD][2];   //start and stop index for each thread
static int     NTHREAD = 0;

#ifdef WIN_SYSTEM
static HANDLE  THREAD[MAX_THREAD];
static HANDLE  TSTART[MAX_THREAD];
static HANDLE  TDONE[MAX_THREAD];
#endif


// computation related
int   PARTASK;
double *Areal, *Breal, *Creal; 
double *Aimag, *Bimag, *Cimag;
bool A_is_complex,B_is_complex;
// WORK and C2 are related to LAPACK calls, these are not needed for matrix multiplication
double *WORK, *C2;
int rA, cA, rB ,cB, rC, cC, strideA, strideB, strideC, strideW, strideC2;
int *PAIRS = NULL;
ptrdiff_t *iScratch = NULL;
bool  BSX;
//bool USED_DGELSY = false;
char MODIFY[2];

// ================================
// teval() is called by each thread
// ================================
#ifdef WIN_SYSTEM
DWORD __stdcall teval(void* pn) 
#else
void* teval(void* pn) 
#endif
{
   int i, n = *(int*)pn;
   double *Areali, *Breali;
   double *Aimagi, *Bimagi;
#ifdef WIN_SYSTEM
   while(1){ // thread will be terminated externally
      WaitForSingleObject(TSTART[n], INFINITE); // wait for start signal
#endif
    // loop over data
	if (A_is_complex) {
		if (B_is_complex) {
			for( i=SCHEDULE[n][0]; i<SCHEDULE[n][1]; i++ ){
				// pointers to scheduled data
				if (BSX){      //singleton expansion
					Areali = Areal + strideA*PAIRS[2*i];
					Aimagi = Aimag + strideA*PAIRS[2*i];
					Breali = Breal + strideB*PAIRS[2*i+1];
					Bimagi = Bimag + strideB*PAIRS[2*i+1];
				} else {
					Areali = Areal + strideA*i;
					Aimagi = Aimag + strideA*i;
					Breali = Breal + strideB*i;
					Bimagi = Bimag + strideB*i;
				}
				mulCMatCMat(Creal + strideC*i, Cimag + strideC*i, Areali, Aimagi, Breali, Bimagi , rA, cA, rB, cB, MODIFY);
			}
		} else {
			for( i=SCHEDULE[n][0]; i<SCHEDULE[n][1]; i++ ){
				// pointers to scheduled data
				if (BSX){      //singleton expansion
					Areali = Areal + strideA*PAIRS[2*i];
					Aimagi = Aimag + strideA*PAIRS[2*i];
					Breali = Breal + strideB*PAIRS[2*i+1];
				} else {
					Areali = Areal + strideA*i;
					Aimagi = Aimag + strideA*i;
					Breali = Breal + strideB*i;
				}
				mulCMatRMat(Creal + strideC*i, Cimag + strideC*i, Areali, Aimagi, Breali , rA, cA, rB, cB, MODIFY);
			}
		}
	} else {
		if (B_is_complex) {
			for( i=SCHEDULE[n][0]; i<SCHEDULE[n][1]; i++ ){
				// pointers to scheduled data
				if (BSX){      //singleton expansion
					Areali = Areal + strideA*PAIRS[2*i];
					Breali = Breal + strideB*PAIRS[2*i+1];
					Bimagi = Bimag + strideB*PAIRS[2*i+1];
				} else {
					Areali = Areal + strideA*i;
					Breali = Breal + strideB*i;
					Bimagi = Bimag + strideB*i;
				}
				mulRMatCMat(Creal + strideC*i, Cimag + strideC*i, Areali, Breali, Bimagi , rA, cA, rB, cB, MODIFY);
			}
		} else {
			for( i=SCHEDULE[n][0]; i<SCHEDULE[n][1]; i++ ){
				// pointers to scheduled data
				if (BSX){      //singleton expansion
					Areali = Areal + strideA*PAIRS[2*i];
					Breali = Breal + strideB*PAIRS[2*i+1];
				} else {
					Areali = Areal + strideA*i;
					Breali = Breal + strideB*i;
				}
				mulRMatRMat(Creal + strideC*i, Areali, Breali , rA, cA, rB, cB, MODIFY);
			}
		}
	}
#ifdef WIN_SYSTEM
      //signal that thread is finished
      SetEvent(TDONE[n]);
   }
#endif
   return 0;
}


// =============
// mexFunction()
// =============
void mexFunction(int n_out, mxArray *p_out[], int n_in, const mxArray *p_in[])
{
   mwSize  Andim, Bndim, Cndim;
   mwSize *Adims, *Bdims, *Adims_full, *Bdims_full, *Cdims, *idx;
   mxArray *tArray = NULL;
   char chr;
   char *commandStr;
   int i, j, k, nt, N, iA, iB, iC, iC2, rC2 = 0;
   static int tnum[MAX_THREAD];

   // no input: print documentation
   if( n_in==0 ) {
      mexPrintf(  	"MMX: fast, multithreaded, n-D multiplication \n"
					"Basic usage:\nThe command   C = mmx('mult',A,B);\n"
					"is equivalent to the matlab loop\n"
					" for i=1:N\n    C(:,:,i) = A(:,:,i)*B(:,:,i);\n end\n"
					"Type 'help mmx' for detailed information. \n");   
   }
   // ===================
   // threading machinery
   // ===================

   // single input mmx(nt): set number of threads
   if( n_in==1 ){
      if((!mxIsNumeric(p_in[0]))||(mxGetN(p_in[0])!=1)||(mxGetM(p_in[0])!=1))
         mexErrMsgTxt("A single scalar input specifies the desired thread count. Type 'help mmx' for more info.");
      nt = (int)mxGetScalar(p_in[0]);
   }
   else if (n_in==0 || !INITIALIZED) {
#ifdef WIN_SYSTEM
      SYSTEM_INFO sysinfo;
      GetSystemInfo( &sysinfo );
      nt = sysinfo.dwNumberOfProcessors;
#else
      nt = sysconf(_SC_NPROCESSORS_ONLN);
#endif
   }
   else {
      nt = NTHREAD;
   }

   // if necessary, clear threads
   if ((nt==0) || (INITIALIZED && (nt!=NTHREAD))){
#ifdef WIN_SYSTEM
      mexPrintf("Clearing threads.\n");
      for( i=0; i<NTHREAD; i++ ) {
         TerminateThread(THREAD[i], 0);
         CloseHandle(THREAD[i]);
         CloseHandle(TSTART[i]);
         CloseHandle(TDONE[i]);   
      }
#endif
      NTHREAD     = 0;
      INITIALIZED = false;      
   }

   // start threads
   if( !INITIALIZED && nt  ) {
      NTHREAD = (MAX_THREAD <= nt) ? MAX_THREAD : nt; // set global NTHREAD
      // create events and threads
      for( i=0; i<NTHREAD; i++ ) {
         tnum[i] = i;   // set tnum so CreateThread won't access i
         //mexPrintf("tnum[%d]: %d\n", i, tnum[i]);
      }
#ifdef WIN_SYSTEM
      for( i=0; i<NTHREAD; i++ ) {
         TSTART[i]   = CreateEvent(0, FALSE, FALSE, 0);
         TDONE[i]    = CreateEvent(0, TRUE, FALSE, 0);
         THREAD[i]   = CreateThread(NULL, 0, teval, (void*)(tnum+i), 0, 0);
      }
#endif
      INITIALIZED = true;
      if (n_in == 1) {//print this line only in single-input mode
         mexPrintf("%d threads prepared.\n", NTHREAD);
      }
   }

   // just getting help or setting the thread count, exit now
   if( n_in < 2 ) {
      return;   
   }

   // ==============
   // process inputs
   // ==============   

   // type check
   if ( (!mxIsDouble(p_in[0])) || (!mxIsDouble(p_in[1])) ) {
      mexErrMsgTxt("Only inputs of type 'double' are supported.");
   }
	
	
	// get the parts of A
	A_is_complex = mxIsComplex(p_in[0]);
	Areal     = mxGetPr(p_in[0]);
	if (A_is_complex){
		Aimag = mxGetPi(p_in[0]);
	}
	Andim = mxGetNumberOfDimensions(p_in[0]);
	Adims = (mwSize *) mxGetDimensions(p_in[0]);
	rA    = Adims[0];
	cA    = Adims[1];    

	// get the parts of B
	B_is_complex = mxIsComplex(p_in[1]);
	Breal     = mxGetPr(p_in[1]);
	if (B_is_complex){
		Bimag = mxGetPi(p_in[1]);
	}
	Bndim = mxGetNumberOfDimensions(p_in[1]);
	Bdims = (mwSize *) mxGetDimensions(p_in[1]);
	rB    = Bdims[0];
	cB    = Bdims[1];

	// modifiers
	MODIFY[0] = MODIFY[1] = 'N';
	if( n_in > 2 ){
		if(mxGetClassID(p_in[2]) != mxCHAR_CLASS)
			mexErrMsgTxt("Third argument is a modifier string. Type 'help mmx'.");
		char *modifierStr = mxArrayToString(p_in[2]);
        for( i=0; i<2; i++ ){
            chr = toupper(modifierStr[i]);
            if ((chr == 'N')||(chr == 'T'))
            	MODIFY[i] = chr;
            else if(chr!='\0')
            	mexErrMsgTxt("Unknown modifier.");
        }
		mxFree(modifierStr);
	}

    // catch the bug before it crashes matlab
    if ( (MODIFY[0] == 'T') || (MODIFY[1] == 'T'))
        mexErrMsgTxt("There is a bug in transpose.");
    
	// ================
	// dimension checks
	// ================  
    if ( (MODIFY[0] == 'N') && (MODIFY[1] == 'N') && (cA != rB) )
		mexErrMsgTxt("size(A,2) == size(B,1) should be true.");
    if ( (MODIFY[0] == 'T') && (MODIFY[1] == 'N') && (rA != rB) )
		mexErrMsgTxt("size(A,1) == size(B,1) should be true.");
    if ( (MODIFY[0] == 'N') && (MODIFY[1] == 'T') && (cA != cB) )
		mexErrMsgTxt("size(A,2) == size(B,2) should be true.");
    if ( (MODIFY[0] == 'T') && (MODIFY[1] == 'T') && (rA != cB) )
		mexErrMsgTxt("size(A,1) == size(B,2) should be true.");      

   // ===============
   // process outputs
   // =============== 

   Cndim    = (Andim > Bndim) ? Andim : Bndim;
   Cndim    = (Cndim > 3) ? Cndim : 3;
   Cdims    = (mwSize *) mxMalloc( Cndim * sizeof(mwSize) );
   idx      = (mwSize *) mxMalloc( Cndim * sizeof(mwSize) );

   // set Cdims[0,1]
	rC = (MODIFY[0] == 'N') ? rA : cA;
    cC = (MODIFY[1] == 'N') ? cB : rB;         

	Cdims[0] = rC;
	Cdims[1] = cC; 

	// Adims_full and Bdims_full pad Adims and Bdims with 1s, if necessary
	Adims_full = (mwSize *) mxMalloc( Cndim * sizeof(mwSize) );
	Bdims_full = (mwSize *) mxMalloc( Cndim * sizeof(mwSize) );  

	// get Cdims and check singleton dimensions
	for( i=0; i<Cndim; i++ ) {
		Adims_full[i] = (i < Andim) ? Adims[i] : 1; 
		Bdims_full[i] = (i < Bndim) ? Bdims[i] : 1;
		if (i > 1){//check singleton-expanded dimensions
			Cdims[i] = (Adims_full[i] > Bdims_full[i]) ? Adims_full[i] : Bdims_full[i];
			if ( ( Adims_full[i]!=1 ) && ( Bdims_full[i]!=1 ) && ( Adims_full[i]!=Bdims_full[i] )  ){
				mexErrMsgTxt("Non-singleton dimensions of the two input arrays must match each other.");
			}         
		}
	}

	// stride sizes
	strideA    = rA*cA;
	strideB    = rB*cB;
	strideC    = rC*cC;

	// N is the total number of matrix operations
	N  = 1;
	for( i=2; i<Cndim; i++ ) {
		N *= Cdims[i];
	}

	// if one of the output dimensions is 0 we're done, goodbye
	if ( Cdims[0]*Cdims[1]*N == 0 ) {
		return;
	}

   // =====================================
   // compute pairs for singleton expansion
   // =====================================

   // check if singleton expansion can be avoided
	BSX   = false;
	if ( (rB != 0) && (cB != 0) ) {
		for( j=2; j<Cndim; j++ ) {
			if (Adims_full[j] != Bdims_full[j]) {
				BSX = true;
			}
		}
	}

	if (BSX) {
      // initialze idx
      for( j=2; j<Cndim; j++ ) {
         idx[j] = 0;
      }

      // init PAIRS
      PAIRS = (int *) mxMalloc( 2 * N * sizeof(int) );   
      PAIRS[0] = PAIRS[1] = 0;

      // compute PAIRS
      // (is there a fast way to do this inside the threads ???)
      for( i=1; i<N; i++ ){
         // idx = ind2sub(size(C), i) in C-style indexing
         idx[2]++;
         for( j=2; j<Cndim; j++ ) {
            if (idx[j] > Cdims[j]-1){
               idx[j] = 0;
               idx[j+1]++;
            }
         }
         // {iA,iB} = sub2ind(size({A,B}), idx) while ignoring singletons
         iA = iB = 0;
         for( j=Cndim-1; j>1; j-- ){
            if (Adims_full[j] > 1)  iA = iA*Adims_full[j] + idx[j];
            if (Bdims_full[j] > 1)  iB = iB*Bdims_full[j] + idx[j];         
         }
         PAIRS[2*i]     = iA;
         PAIRS[2*i+1]   = iB;      
      }      
	}

	// allocate C
	if (A_is_complex||B_is_complex){
		p_out[0] = mxCreateNumericArray(Cndim, Cdims, mxDOUBLE_CLASS, mxCOMPLEX);
		Creal  = mxGetPr(p_out[0]);
		Cimag  = mxGetPi(p_out[0]);
	} else {
		p_out[0] = mxCreateNumericArray(Cndim, Cdims, mxDOUBLE_CLASS, mxREAL);
		Creal  = mxGetPr(p_out[0]);
	}
   
	n_out = 1;
   

	// ==================================
	// make schedule, run threads, finish
	// ==================================   

	// set SCHEDULE
	int blksz = N/NTHREAD;
	int extra = N - blksz*NTHREAD;
#ifdef DEBUG   
   mexPrintf("matrix_ops: %d block: %d extra: %d\n", N, blksz, extra);
#endif   
	for( i=0; i<NTHREAD; i++ ) {
		SCHEDULE[i][0] = ((i>0) ? SCHEDULE[i-1][1] : 0);
		SCHEDULE[i][1] = SCHEDULE[i][0] + (blksz + (i<extra));
#ifdef DEBUG 	  
		mexPrintf("SCHEDULE[%d] %d, %d\n", i, SCHEDULE[i][0], SCHEDULE[i][1]);
#endif	  
   }

   // signal threads to start
#ifdef WIN_SYSTEM
   for( i=0; i<NTHREAD; i++ ) {
      SetEvent(TSTART[i]);
   }
   //  wait for all threads to finish
   WaitForMultipleObjects(NTHREAD, TDONE, TRUE, INFINITE);
   // reset TDONE events
   for( i=0; i<NTHREAD; i++ ) {
      ResetEvent(TDONE[i]);
   }
#else
   pthread_t p_threads[NTHREAD];
   for( i=0; i<NTHREAD; i++ ) {
      if (pthread_create(&p_threads[i],
               NULL, teval, (void*)(tnum+i)) != 0) {
         mexPrintf("Could not create thread %d\n", i);
      }
   }
   for( i=0; i<NTHREAD; i++ )
   {
      if (pthread_join(p_threads[i], NULL) != 0) {
         mexPrintf("Could not join thread %d\n", i);    
      }
   }
#endif

   mxFree(Adims_full);
   mxFree(Bdims_full);
   if (BSX) {
      mxFree(PAIRS);
   }
   mxFree(Cdims);
   mxFree(idx);

   if (tArray != NULL) {
      mxDestroyArray(tArray);
   }
}

