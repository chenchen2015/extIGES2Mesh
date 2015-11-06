/**************************************************************************
 *
 * function u=nrbCrvPlaneIntrsctIGES(nurbs,dnurbs,d2nurbs,normal,d)
 *
 * Computes the parameter value for the closest point on a NURBS curve to a plane
 * using binary search and Newton's method
 *
 * Input:
 * nurbs - NURBS structure
 * dnurbs,d2nurbs - NURBS derivatives (output from nrbDerivativesIGES)
 * normal - 3x1 vector, normal of the plane  (normal'*x=d)
 * d - rhs of the equation of the plane
 *
 * Output:
 * u - NURBS Parameter value at the closest point
 *
 * c-file can be downloaded at
 *
 * http://www.mathworks.com/matlabcentral/fileexchange/13253-iges-toolbox
 *
 * compiling in Matlab is done using the command  "mex nrbCrvPlaneIntrsctIGES.c"
 *
 * See "help mex" for more information
 *
 * written by Per Bergström 2012-02-24
 *
 **************************************************************************/

#include <math.h>
#include "mex.h"

/* Input Arguments */

#define	nurbsstructure	prhs[0]
#define	dnrbsstructure	prhs[1]
#define d2nrbsstructure	prhs[2]
#define normal	prhs[3]
#define planeRHS	prhs[4]

/* Output Arguments */

#define	parametervalue	plhs[0]

/* Misc */

#define MAXITER 20

/* Sub functions (in folder "mexSourceFiles") */

#include "mexSourceFiles/FindSpan.c"
#include "mexSourceFiles/BasisFuns.c"
#include "mexSourceFiles/NURBScurveEval.c"
#include "mexSourceFiles/nrbD1D2eval.c"


/* Main function */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    int i, orderU, degU, ncp;
    double paramtmp[1], minParam[1], paramspacing, umin, umax, funcVal, minFuncVal, bspPnts[4], pnttmp[3], pderu[3], pderuu[3], step;
    double *nrmlPntr, *plconst;
    double *left, *right, *N;
    
    if(nlhs==0){
        mexErrMsgTxt("There must be an output");
    }
    if(nrhs!=5){
        mexErrMsgTxt("Number of inputs must be 5");
    }
    if (mxGetM(mxGetField(nurbsstructure, 0, "coefs"))!=4){
        mexPrintf("Number of rows in nurbs.coefs is %d,\n", mxGetM(mxGetField(nurbsstructure, 0, "coefs")));
        mexErrMsgTxt("nurbs.coefs must have 4 rows.");
    }
    
    nrmlPntr=mxGetPr(normal);
    plconst=mxGetPr(planeRHS);
    
    parametervalue = mxCreateDoubleMatrix(1, 1, mxREAL);
    
    orderU = (int)mxGetScalar(mxGetField(nurbsstructure, 0, "order"));
    degU = orderU-1;
    
    ncp = (int)mxGetScalar(mxGetField(nurbsstructure, 0, "number"));    
    
    left  = malloc(orderU*sizeof(double)+1);
    right = malloc(orderU*sizeof(double)+1);
    N = malloc(orderU*sizeof(double)+1);
    
    umin = (double)mxGetPr(mxGetField(nurbsstructure, 0, "knots"))[0];
    umax = (double)mxGetPr(mxGetField(nurbsstructure, 0, "knots"))[ncp];
    
    paramspacing=(umax-umin)/50.0;
    paramtmp[0]=umin;
    minParam[0]=umin;
    
    NURBScurveEval(degU, mxGetPr(mxGetField(nurbsstructure, 0, "coefs")), ncp, mxGetPr(mxGetField(nurbsstructure, 0, "knots")), &paramtmp[0], 1, &pnttmp[0], left, right, N);
    minFuncVal=fabs(nrmlPntr[0]*pnttmp[0]+nrmlPntr[1]*pnttmp[1]+nrmlPntr[2]*pnttmp[2]-plconst[0]);
    
    if (minFuncVal<1e-12){
        mxGetPr(parametervalue)[0]=paramtmp[0];
        return;
    }
    else {
        
        /* Find the closest point among 50 evaluations */
        
        for (i = 0; i < 50; i++){
            
            paramtmp[0]+=paramspacing;
            NURBScurveEval(degU, mxGetPr(mxGetField(nurbsstructure, 0, "coefs")), ncp, mxGetPr(mxGetField(nurbsstructure, 0, "knots")), &paramtmp[0], 1, &pnttmp[0], left, right, N);
            funcVal=fabs(nrmlPntr[0]*pnttmp[0]+nrmlPntr[1]*pnttmp[1]+nrmlPntr[2]*pnttmp[2]-plconst[0]);
            
            if (funcVal<minFuncVal){
                
                minFuncVal=funcVal;
                minParam[0]=paramtmp[0];
                
                if (minFuncVal<1e-12){
                    mxGetPr(parametervalue)[0]=paramtmp[0];
                    return;
                }
                
            }
            
        }
        
        /* Find a closer point using binary search */
        
        for (i = 0; i < 20; i++){
            
            paramspacing=0.5*paramspacing;
            
            paramtmp[0]=minParam[0]-paramspacing;
            NURBScurveEval(degU, mxGetPr(mxGetField(nurbsstructure, 0, "coefs")), ncp, mxGetPr(mxGetField(nurbsstructure, 0, "knots")), &paramtmp[0], 1, &pnttmp[0], left, right, N);
            funcVal=fabs(nrmlPntr[0]*pnttmp[0]+nrmlPntr[1]*pnttmp[1]+nrmlPntr[2]*pnttmp[2]-plconst[0]);
            
            if (funcVal<minFuncVal){
                
                minFuncVal=funcVal;
                minParam[0]=paramtmp[0];
                
                if (minFuncVal<1e-12){
                    
                    if (paramtmp[0] < umin){
                        paramtmp[0] = umin;
                    }
                    else if (paramtmp[0] > umax){
                        paramtmp[0] = umax;
                    }
                    mxGetPr(parametervalue)[0]=paramtmp[0];
                    return;
                    
                }
                
            }
            else {
                
                paramtmp[0]=minParam[0]+paramspacing;
                NURBScurveEval(degU, mxGetPr(mxGetField(nurbsstructure, 0, "coefs")), ncp, mxGetPr(mxGetField(nurbsstructure, 0, "knots")), &paramtmp[0], 1, &pnttmp[0], left, right, N);
                funcVal=fabs(nrmlPntr[0]*pnttmp[0]+nrmlPntr[1]*pnttmp[1]+nrmlPntr[2]*pnttmp[2]-plconst[0]);
                
                if (funcVal<minFuncVal){
                    
                    minFuncVal=funcVal;
                    minParam[0]=paramtmp[0];
                    
                    if (minFuncVal<1e-12){
                        
                        if (paramtmp[0] < umin){
                            paramtmp[0] = umin;
                        }
                        else if (paramtmp[0] > umax){
                            paramtmp[0] = umax;
                        }
                        mxGetPr(parametervalue)[0]=paramtmp[0];
                        return;
                        
                    }
                    
                }
                
            }
            
        }
        
        /* Find a closer point using Newton's method */
        
        for (i = 0; i < MAXITER; i++){
            
            if (minParam[0] < umin){
                minParam[0] = umin;
            }
            else if (minParam[0] > umax){
                minParam[0] = umax;
            }
            
            nrbD1D2eval(orderU, ncp, nurbsstructure, dnrbsstructure, d2nrbsstructure, &minParam[0], &pnttmp[0], &pderu[0], &pderuu[0], &bspPnts[0], left, right, N);
            funcVal=nrmlPntr[0]*pnttmp[0]+nrmlPntr[1]*pnttmp[1]+nrmlPntr[2]*pnttmp[2]-plconst[0];
            
            if (fabs(funcVal)<1e-12){
                mxGetPr(parametervalue)[0]=minParam[0];
                return;
            }
            else {
                
                bspPnts[0]=nrmlPntr[0]*pderu[0]+nrmlPntr[1]*pderu[1]+nrmlPntr[2]*pderu[2];
                bspPnts[1]=bspPnts[0]*bspPnts[0]+funcVal*(nrmlPntr[0]*pderuu[0]+nrmlPntr[1]*pderuu[1]+nrmlPntr[2]*pderuu[2]);
                
                if (fabs(bspPnts[1])<1e-4){
                    
                    minFuncVal=fabs(funcVal);
                    paramspacing=0.5*paramspacing;
                    
                    paramtmp[0]=minParam[0]-paramspacing;
                    NURBScurveEval(degU, mxGetPr(mxGetField(nurbsstructure, 0, "coefs")), ncp, mxGetPr(mxGetField(nurbsstructure, 0, "knots")), &paramtmp[0], 1, &pnttmp[0], left, right, N);
                    funcVal=fabs(nrmlPntr[0]*pnttmp[0]+nrmlPntr[1]*pnttmp[1]+nrmlPntr[2]*pnttmp[2]-plconst[0]);
                    
                    if (funcVal<minFuncVal){
                        
                        minFuncVal=funcVal;
                        minParam[0]=paramtmp[0];
                        
                        if (minFuncVal<1e-12){
                            
                            if (paramtmp[0] < umin){
                                paramtmp[0] = umin;
                            }
                            else if (paramtmp[0] > umax){
                                paramtmp[0] = umax;
                            }
                            mxGetPr(parametervalue)[0]=paramtmp[0];
                            return;
                            
                        }
                        
                    }
                    else {
                        
                        paramtmp[0]=minParam[0]+paramspacing;
                        NURBScurveEval(degU, mxGetPr(mxGetField(nurbsstructure, 0, "coefs")), ncp, mxGetPr(mxGetField(nurbsstructure, 0, "knots")), &paramtmp[0], 1, &pnttmp[0], left, right, N);
                        funcVal=fabs(nrmlPntr[0]*pnttmp[0]+nrmlPntr[1]*pnttmp[1]+nrmlPntr[2]*pnttmp[2]-plconst[0]);
                        
                        if (funcVal<minFuncVal){
                            
                            minFuncVal=funcVal;
                            minParam[0]=paramtmp[0];
                            
                            if (minFuncVal<1e-12){
                                
                                if (paramtmp[0] < umin){
                                    paramtmp[0] = umin;
                                }
                                else if (paramtmp[0] > umax){
                                    paramtmp[0] = umax;
                                }
                                mxGetPr(parametervalue)[0]=paramtmp[0];
                                return;
                                
                            }
                            
                        }
                        
                    }
                    
                }
                else {
                    
                    paramspacing=0.8*paramspacing;
                    
                    step=(funcVal*bspPnts[0])/(bspPnts[1]);
                    if (step < -paramspacing){
                        step = -paramspacing;
                    }
                    else if (step > paramspacing){
                        step = paramspacing;
                    }
                    minParam[0] -= step;
                }
                
            }
            
        }
        
        if (minParam[0] < umin){
            minParam[0] = umin;
        }
        else if (minParam[0] > umax){
            minParam[0] = umax;
        }
        mxGetPr(parametervalue)[0]=minParam[0];
        
    }
    
    free(left);
    free(right);
    free(N);
    
}
