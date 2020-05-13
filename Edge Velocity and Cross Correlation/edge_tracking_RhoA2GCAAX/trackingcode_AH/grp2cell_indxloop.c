#include "mex.h"

/* The matlab code I'm trying to accelerate
 *
 * for i=1:n
 *       cl{i}=M((df(i)+1):df(i+1),:);
 * end
 *
 * Where cl is a cell array, M is a double matrix (must I know its type?) 
 * and df is a list of indexes where in each cell in cl whould be items df(i)+1 to df(i+1)
 * 
 * the matlab call will be grp2cell_indxloop(M,df)
 */

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{     
     mxArray *cl;            
     mwSize ngrp,ncols;
     mxArray *m;
     double *avg;
     int i,j,k;
     unsigned long cnt,strt,fnsh;
     double *df,*df2;
     double *mptr,*M;
     
     M=(double*)mxGetPr(prhs[0]);
     df=(double*)mxGetPr(prhs[1]);
     ngrp = mxGetM(prhs[1])-1; /* Get the size of df-1 */
     ncols = mxGetN(prhs[0]); /* get the number of columns in matrix M */       
     cl = mxCreateCellMatrix(ngrp, (mwSize)1);
     if (nlhs==2) {
         plhs[1]=mxCreateDoubleMatrix(ngrp,1,mxREAL);
         avg=(double*)mxGetPr(plhs[1]);
     }
     cnt=0;
     for (i=0;i<ngrp;i++){
         /*create the array for this cell*/
         strt=(unsigned long)df[i];
         fnsh=(unsigned long)df[i+1];
         cnt=fnsh-strt;
         /*mexPrintf("i: %d cnt: %d strt: %d fnsh: %d \n",i,cnt,strt,fnsh);*/
         m=mxCreateDoubleMatrix((mwSize)cnt,ncols, mxREAL);
         mptr=(double*)mxGetPr(m);
         for (j=0;j<cnt;j++) {
             for (k=0;k<ncols;k++) {
                 mptr[j*ncols+k]=M[(strt+j)*ncols+k];
                 if (nlhs==2) {
                     avg[i]+=mptr[j*ncols+k];
                 }
             }
         }
         /*copy small m into the cell array */
         mxSetCell(cl,i,m);
         if (nlhs==2) {
             avg[i]=avg[i]/cnt/ncols;
         }
     }
     /* retun the created cell arrat*/
     plhs[0]=cl;
}