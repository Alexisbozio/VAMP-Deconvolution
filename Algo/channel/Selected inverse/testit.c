#include "mex.h"
#include <math.h>
#include <assert.h>
#include "matrix.h"

double * cselinv(int *nnodes, int *nnz, int*ip , int*jp, double*a);

void mexFunction(int nl, mxArray  *pl[], int nr, const mxArray *pr[] ){
  int   nz,nd;
  double *ap,  *res;
  int i,j;
  if( nr != 1 ) mexErrMsgTxt("oneargument\n");
  if( nl != 1 ) mexErrMsgTxt("one result \n");
  mwIndex* Ir=mxGetIr( pr[0] );
  mwIndex* Jc=mxGetJc( pr[0] );
  nd = mxGetN (pr[0]); /*dimension matrix*/
  nz = Jc[nd]; /*number of non-zero entries */
  double *A=mxGetPr( pr[0] ); /* Values in the matrix */

  int *ir= (int*) malloc( nz*sizeof(int));
  int *jc= (int*) malloc( (nd+1)*sizeof(int));

  for(i=0;i< nz;i++){
    ir[i]=  Ir[i]+1;
  }
  for(j=0;j< nd+1;j++){
    jc[j] = Jc[j]+1; 
  }

  res=cselinv(&nd,&nz,jc , ir,A); /* Call Lin Lin code */
  pl[0]=mxCreateDoubleMatrix(nd,1,mxREAL); /* copy result back to output array */
  double*c=mxGetPr(pl[0])	;
  for(j=0;j<nd;j++){
    c[j]=res[j];
  }
  free(ir);
  free(jc);
  free(res);
}

extern int ldlt_preprocess__(int *, int *, int *, int *, int *, int *, int *);
extern int ldlt_fact__(int *, int *, int *, double *);
extern int ldlt_free__(int *);
extern int ldlt_selinv__(int *, double *, int*);

double *cselinv(int *nnodes, int *nnz, int*ia , int*ja, double*a) /*Adapted from the Lin Lin code EXAMPLES/selinv.c*/
{
  int *perm;
  double  *diag;
  int token=0, order=-1; /* order=-1 for default, or 2 for metis*/
  int Lnnz;
  int dumpL=0;
  perm=NULL;

  ldlt_preprocess__(&token, nnodes, ia, ja, &Lnnz, &order, perm);
  diag = (double*)calloc(*nnodes ,sizeof(double));
  ldlt_fact__(&token, ia, ja, a);
  ldlt_selinv__(&token, diag, &dumpL);
  ldlt_free__(&token);
  return diag;
}

