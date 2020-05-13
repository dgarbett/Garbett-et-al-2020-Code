#include "mex.h"

/*
 * mex implementation of the trajectory building part of ultTracking
 */

/* The matlab code I'm trying to accelerate
 *
 *while any(possibleStart) 
 *   cnt=cnt+1;
 *   nxt=find(possibleStart,1,'first');
 *   while nxt>0
 *       traj(nxt)=cnt;
 *       possibleStart(nxt)=0;
 *       nxt=gr(nxt);
 *   end
 *end
 */

void traverseGraph(unsigned long *gr, bool *possibleStart, double N,unsigned long *traj)
{
  int cnt=0;
  int nxt=1;
  int strtIndx=0;
  int i;
  while (strtIndx < N) {
      cnt++;
      /* next line is the equavalent of find(notvisted,1) */
      while (strtIndx < N && possibleStart[strtIndx]==0) {
          strtIndx++;
      }
      nxt=strtIndx+1;
      while (strtIndx < N && nxt>0) {
          traj[nxt-1]=(double)cnt;
          possibleStart[nxt-1]=0;
          nxt=(int)gr[nxt-1];
      }
  }
}


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  unsigned long *gr; /* pointers for the inputs */
  bool *possibleStart;
  unsigned long *traj; /* pointer for the output */
  mwSize N;
  mxArray *possibleStartCpy;
  
  /*  check for proper number of arguments */
  if(nrhs!=2) 
    mexErrMsgTxt("Two inputs required: traj=traverseGraph(gr,possibleStart);");
  if(nlhs!=1) 
    mexErrMsgTxt("One output required: traj=traverseGraph(gr,possibleStart);");
  
  /*  get the graph gr size */
  N = mxGetN(prhs[0]);
  
  /* create a logical array and copy possibleStart to it*/
  possibleStartCpy = mxDuplicateArray(prhs[1]);
  
  /*  get a pointer for the graph vector input x 
   *since gr isn't changing no need to create a matrix for it (unlike possibleStart)*/
  gr = (unsigned long*) mxGetPr(prhs[0]);
  possibleStart = (bool*) mxGetPr(possibleStartCpy);
 
  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateNumericMatrix(1,N,mxUINT32_CLASS, mxREAL);
  
  /*  create a C pointer to a copy of the output matrix */
  traj = (unsigned long*)mxGetPr(plhs[0]);
  
  /*  call the C subroutine */
  traverseGraph(gr,possibleStart,N,traj);
 
}
