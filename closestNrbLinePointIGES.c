/**************************************************************************
 *
 * function [P,UV]=closestNrbLinePointIGES(nurbs,dnurbs,d2nurbs,UV0,r0,v)
 *
 * Closest points of NURBS patch and line/point using Newton's method.
 *
 * Line (3D):  r=r0+t*v
 * Point (3D): r0
 *
 * Usage in Matlab:
 *
 * § Closest NURBS-point to line
 * [P,UV]=closestNrbLinePointIGES(nurbs,dnurbs,d2nurbs,UV0,r0,v)
 *
 * § Closest NURBS-point to point
 * [P,UV]=closestNrbLinePointIGES(nurbs,dnurbs,d2nurbs,UV0,r0)
 *
 * Input:
 * nurbs - NURBS structure
 * dnurbs,d2nurbs - NURBS derivatives (output from nrbDerivativesIGES).
 * UV0 - Initial start Parameter values
 * r0,(v) - See Line/Point (3D) above. r0 (and v) must have the dimension (3x{1 or N})
 *          If size (3x1) - same point/line is used, if size (3xN) - different point/lines are used.
 *
 * Curve
 * UV0 - 1xN matrix, N number of u-parameters
 *
 * Surface
 * UV0 - 2xN matrix, N number of (u,v)-parameters
 *
 * Output:
 * P - Closest points on NURBS patch.
 * UV - NURBS Parameter values at closest point. (same dimension as UV0)
 *
 * c-file can be downloaded at
 *
 * http://www.mathworks.com/matlabcentral/fileexchange/13253-iges-toolbox
 *
 * compile in Matlab by using the command  "mex -v closestNrbLinePointIGES.c"
 *
 * See "help mex" for more information
 *
 * written by Per Bergström 2009-12-04
 *
 **************************************************************************/

#include <math.h>
#include "mex.h"

/* Input Arguments */

#define	nurbsstructure	prhs[0]
#define	dnrbsstructure	prhs[1]
#define d2nrbsstructure	prhs[2]
#define initparamvalues	prhs[3]
#define point0	prhs[4]
#define linedirection	prhs[5]

/* Output Arguments */

#define	evaluated_points	plhs[0]
#define	parametervalues	plhs[1]

/* Misc */

#define MAXITER 30

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

/* Sub functions (in folder "mexSourceFiles") */

#include "mexSourceFiles/FindSpan.c"
#include "mexSourceFiles/BasisFuns.c"
#include "mexSourceFiles/BspEval.c"
#include "mexSourceFiles/BspEval2.c"
#include "mexSourceFiles/NURBSsurfaceEval.c"
#include "mexSourceFiles/nrbD1D2eval.c"
#include "mexSourceFiles/nrbD1D2eval2.c"


int updateClosestSrfPnt2Line(int degU, int degV, double *cp, int ncp, int kcp, double *knotU, double *knotV, double *P0, double *matA, double *us, double *ep, double *epu, double *epv, double *epuu, double *epuv, double *epvv, double *leftU, double *rightU, double *NU, double *leftV, double *rightV, double *NV, double *InvHess, double *grad, double *dir, double *ustmp){
    
    int i;
    double f0, k0, s;
    
    InvHess[0]=ep[0]-P0[0];
    InvHess[1]=ep[1]-P0[1];
    InvHess[2]=ep[2]-P0[2];
    
    ep[0]=matA[0]*InvHess[0]+matA[1]*InvHess[1]+matA[2]*InvHess[2];
    ep[1]=matA[1]*InvHess[0]+matA[3]*InvHess[1]+matA[4]*InvHess[2];
    ep[2]=matA[2]*InvHess[0]+matA[4]*InvHess[1]+matA[5]*InvHess[2];
    
    f0=InvHess[0]*ep[0]+InvHess[1]*ep[1]+InvHess[2]*ep[2];
    
    grad[0]=ep[0]*epu[0]+ep[1]*epu[1]+ep[2]*epu[2];
    grad[1]=ep[0]*epv[0]+ep[1]*epv[1]+ep[2]*epv[2];
    
    leftU[0]=matA[0]*epu[0]+matA[1]*epu[1]+matA[2]*epu[2];
    k0=matA[1]*epu[0]+matA[3]*epu[1]+matA[4]*epu[2];
    s=matA[2]*epu[0]+matA[4]*epu[1]+matA[5]*epu[2];
    
    InvHess[2]=ep[0]*epuu[0]+ep[1]*epuu[1]+ep[2]*epuu[2]+leftU[0]*epu[0]+k0*epu[1]+s*epu[2];
    InvHess[1]=-(ep[0]*epuv[0]+ep[1]*epuv[1]+ep[2]*epuv[2]+leftU[0]*epv[0]+k0*epv[1]+s*epv[2]);
    
    leftU[0]=matA[0]*epv[0]+matA[1]*epv[1]+matA[2]*epv[2];
    k0=matA[1]*epv[0]+matA[3]*epv[1]+matA[4]*epv[2];
    s=matA[2]*epv[0]+matA[4]*epv[1]+matA[5]*epv[2];
    
    InvHess[0]=ep[0]*epvv[0]+ep[1]*epvv[1]+ep[2]*epvv[2]+leftU[0]*epv[0]+k0*epv[1]+s*epv[2];
    
    k0=InvHess[0]*InvHess[2]-InvHess[1]*InvHess[1];
    
    if (fabs(k0)<1.0e-20){
        k0=sqrt(grad[0]*grad[0]+grad[1]*grad[1]);
        if (fabs(k0)<1.0e-20){
            return 1;
        }
        else {
            dir[0]=-0.05*grad[0]/k0;
            dir[1]=-0.05*grad[1]/k0;
        }
    }
    else {
        dir[0]=-(InvHess[0]*grad[0]+InvHess[1]*grad[1])/k0;
        dir[1]=-(InvHess[1]*grad[0]+InvHess[2]*grad[1])/k0;
    }
    
    k0=0.7*(dir[0]*grad[0]+dir[1]*grad[1]);
    
    s=1.0;
    ustmp[0]=us[0]+dir[0];
    ustmp[1]=us[1]+dir[1];
    
    for (i = 0; i < 100; i++){
        
        NURBSsurfaceEval(degU, degV, cp, ncp, kcp, knotU, knotV, ustmp, 1, ep, leftU, rightU, NU, leftV, rightV, NV);
        
        InvHess[0]=ep[0]-P0[0];
        InvHess[1]=ep[1]-P0[1];
        InvHess[2]=ep[2]-P0[2];
        
        grad[0]=matA[0]*InvHess[0]+matA[1]*InvHess[1]+matA[2]*InvHess[2];
        grad[1]=matA[1]*InvHess[0]+matA[3]*InvHess[1]+matA[4]*InvHess[2];
        leftU[0]=matA[2]*InvHess[0]+matA[4]*InvHess[1]+matA[5]*InvHess[2];
        
        if ( (InvHess[0]*grad[0]+InvHess[1]*grad[1]+InvHess[2]*leftU[0])<(f0+s*k0) ){
            us[0]=ustmp[0];
            us[1]=ustmp[1];
            return 0;
        }
        else {
            s=0.7*s;
            ustmp[0]=us[0]+s*dir[0];
            ustmp[1]=us[1]+s*dir[1];
        }
        
    }
    
    return 1;
    
}

int updateClosestSrfPnt2Pnt(int degU, int degV, double *cp, int ncp, int kcp, double *knotU, double *knotV, double *P0, double *us, double *ep, double *epu, double *epv, double *epuu, double *epuv, double *epvv, double *leftU, double *rightU, double *NU, double *leftV, double *rightV, double *NV, double *InvHess, double *grad, double *dir, double *ustmp){
    
    int i;
    double f0, k0, s;
    
    ep[0]-=P0[0];
    ep[1]-=P0[1];
    ep[2]-=P0[2];
    
    grad[0]=ep[0]*epu[0]+ep[1]*epu[1]+ep[2]*epu[2];
    grad[1]=ep[0]*epv[0]+ep[1]*epv[1]+ep[2]*epv[2];
    
    InvHess[2]=ep[0]*epuu[0]+ep[1]*epuu[1]+ep[2]*epuu[2]+epu[0]*epu[0]+epu[1]*epu[1]+epu[2]*epu[2];
    InvHess[1]=-(ep[0]*epuv[0]+ep[1]*epuv[1]+ep[2]*epuv[2]+epu[0]*epv[0]+epu[1]*epv[1]+epu[2]*epv[2]);
    InvHess[0]=ep[0]*epvv[0]+ep[1]*epvv[1]+ep[2]*epvv[2]+epv[0]*epv[0]+epv[1]*epv[1]+epv[2]*epv[2];
    
    k0=InvHess[0]*InvHess[2]-InvHess[1]*InvHess[1];
    
    if (fabs(k0)<1.0e-20){
        k0=sqrt(grad[0]*grad[0]+grad[1]*grad[1]);
        if (fabs(k0)<1.0e-20){
            return 1;
        }
        else {
            dir[0]=-0.05*grad[0]/k0;
            dir[1]=-0.05*grad[1]/k0;
        }
    }
    else {
        dir[0]=-(InvHess[0]*grad[0]+InvHess[1]*grad[1])/k0;
        dir[1]=-(InvHess[1]*grad[0]+InvHess[2]*grad[1])/k0;
    }
    
    f0=ep[0]*ep[0]+ep[1]*ep[1]+ep[2]*ep[2];
    k0=0.7*(dir[0]*grad[0]+dir[1]*grad[1]);
    
    s=1.0;
    ustmp[0]=us[0]+dir[0];
    ustmp[1]=us[1]+dir[1];
    
    for (i = 0; i < 100; i++){
        
        NURBSsurfaceEval(degU, degV, cp, ncp, kcp, knotU, knotV, ustmp, 1, ep, leftU, rightU, NU, leftV, rightV, NV);
        
        InvHess[0]=ep[0]-P0[0];
        InvHess[1]=ep[1]-P0[1];
        InvHess[2]=ep[2]-P0[2];
        
        if ( (InvHess[0]*InvHess[0]+InvHess[1]*InvHess[1]+InvHess[2]*InvHess[2])<(f0+s*k0) ){
            us[0]=ustmp[0];
            us[1]=ustmp[1];
            return 0;
        }
        else {
            s=0.7*s;
            ustmp[0]=us[0]+s*dir[0];
            ustmp[1]=us[1]+s*dir[1];
        }
        
    }
    
    return 1;
    
}

/* Main function */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    int i, j;
    double stepl, matHess[3], grad[2], dir[2], ustmp[2], matA[6];
    double t, detH, res[3];
    char TrueFalse=0;
    double paramtmp[2], umin, umax, vmin, vmax, bspPnts[4], pnttmp[3], pderu[3], pderv[3], pderuu[3], pderuv[3], pdervv[3];
    
    if (nlhs!=2){
        mexErrMsgTxt("Number of outputs must be 2.");
    }
    if (nrhs>6 || nrhs<5){
        mexErrMsgTxt("Number of inputs must be 5 or 6.");
    }
    if (mxGetM(mxGetField(nurbsstructure, 0, "coefs"))!=4){
        mexPrintf("Number of rows in nurbs.coefs is %d,\n", mxGetM(mxGetField(nurbsstructure, 0, "coefs")));
        mexErrMsgTxt("nurbs.coefs must have 4 rows.");
    }
    if (mxGetM(point0)!=3){
        mexErrMsgTxt("r0 must have 3 rows.");
    }
    if (mxGetM(initparamvalues)>2 || mxGetM(initparamvalues)==0){
        mexErrMsgTxt("UV0 must be of dim 1xN or 2xN.");
    }
    if (nrhs==6){
        if (mxGetM(linedirection)!=3){
            mexErrMsgTxt("v must have 3 rows.");
        }
        if (mxGetN(linedirection)!=mxGetN(point0)){
            mexErrMsgTxt("r0 and v must be of same size.");
        }
        else if(mxGetN(linedirection)==mxGetN(initparamvalues)){
            if (mxGetN(point0)>1){
                TrueFalse=1;
            }
        }
    }
    else if(nrhs==5){
        if (mxGetN(point0)==mxGetN(initparamvalues)){
            if (mxGetN(point0)>1){
                TrueFalse=1;
            }
        }
    }
    
    evaluated_points = mxCreateDoubleMatrix(3, (int)mxGetN(initparamvalues), mxREAL);
    parametervalues = mxCreateDoubleMatrix((int)mxGetM(initparamvalues), (int)mxGetN(initparamvalues), mxREAL);
    
    if (mxGetM(initparamvalues)==2){
        
        int orderU = (int)mxGetPr(mxGetField(nurbsstructure, 0, "order"))[0];
        int orderV = (int)mxGetPr(mxGetField(nurbsstructure, 0, "order"))[1];
        
        int degU = orderU-1;
        int degV = orderV-1;
        
        int ncp = (int)mxGetPr(mxGetField(nurbsstructure, 0, "number"))[0];
        int kcp = (int)mxGetPr(mxGetField(nurbsstructure, 0, "number"))[1];
        
        double *leftU  = (double*) malloc(orderU*sizeof(double)+1);
        double *rightU = (double*) malloc(orderU*sizeof(double)+1);
        double *NU = (double*) malloc(orderU*sizeof(double)+1);
        double *leftV  = (double*) malloc(orderV*sizeof(double)+1);
        double *rightV = (double*) malloc(orderV*sizeof(double)+1);
        double *NV = (double*) malloc(orderV*sizeof(double)+1);
        
        umin = (double)mxGetPr(mxGetCell(mxGetField(nurbsstructure, 0, "knots"), 0))[0];
        umax = (double)mxGetPr(mxGetCell(mxGetField(nurbsstructure, 0, "knots"), 0))[ncp];
        vmin = (double)mxGetPr(mxGetCell(mxGetField(nurbsstructure, 0, "knots"), 1))[0];
        vmax = (double)mxGetPr(mxGetCell(mxGetField(nurbsstructure, 0, "knots"), 1))[kcp];
        
        if (nrhs==6){
            
            if (!TrueFalse){
                matA[0]=mxGetPr(linedirection)[1]*mxGetPr(linedirection)[1]+mxGetPr(linedirection)[2]*mxGetPr(linedirection)[2];
                matA[1]=-mxGetPr(linedirection)[0]*mxGetPr(linedirection)[1];
                matA[2]=-mxGetPr(linedirection)[0]*mxGetPr(linedirection)[2];
                matA[3]=mxGetPr(linedirection)[0]*mxGetPr(linedirection)[0]+mxGetPr(linedirection)[2]*mxGetPr(linedirection)[2];
                matA[4]=-mxGetPr(linedirection)[1]*mxGetPr(linedirection)[2];
                matA[5]=mxGetPr(linedirection)[0]*mxGetPr(linedirection)[0]+mxGetPr(linedirection)[1]*mxGetPr(linedirection)[1];
            }
            
            for (i = 0; i < mxGetN(initparamvalues); i++){
                
                if (TrueFalse){
                    matA[0]=mxGetPr(linedirection)[3*i+1]*mxGetPr(linedirection)[3*i+1]+mxGetPr(linedirection)[3*i+2]*mxGetPr(linedirection)[3*i+2];
                    matA[1]=-mxGetPr(linedirection)[3*i]*mxGetPr(linedirection)[3*i+1];
                    matA[2]=-mxGetPr(linedirection)[3*i]*mxGetPr(linedirection)[3*i+2];
                    matA[3]=mxGetPr(linedirection)[3*i]*mxGetPr(linedirection)[3*i]+mxGetPr(linedirection)[3*i+2]*mxGetPr(linedirection)[3*i+2];
                    matA[4]=-mxGetPr(linedirection)[3*i+1]*mxGetPr(linedirection)[3*i+2];
                    matA[5]=mxGetPr(linedirection)[3*i]*mxGetPr(linedirection)[3*i]+mxGetPr(linedirection)[3*i+1]*mxGetPr(linedirection)[3*i+1];
                }
                
                paramtmp[0] = (double)mxGetPr(initparamvalues)[2*i];
                paramtmp[1] = (double)mxGetPr(initparamvalues)[2*i+1];
                
                t=0.0;
                for (j = 0; j < MAXITER; j++){
                    
                    if (paramtmp[0] < umin){
                        paramtmp[0] = umin;
                    }
                    else if(paramtmp[0] > umax){
                        paramtmp[0] = umax;
                    }
                    if (paramtmp[1] < vmin){
                        paramtmp[1] = vmin;
                    }
                    else if(paramtmp[1] > vmax){
                        paramtmp[1] = vmax;
                    }
                    
                    nrbD1D2eval2(orderU, orderV, ncp, kcp, nurbsstructure, dnrbsstructure, d2nrbsstructure, paramtmp, pnttmp, pderu, pderv, pderuu, pderuv, pdervv, bspPnts, leftU, rightU, NU, leftV, rightV, NV);
                    
                    if (TrueFalse){
                        if (updateClosestSrfPnt2Line(degU, degV, mxGetPr(mxGetField(nurbsstructure, 0, "coefs")), ncp, kcp, mxGetPr(mxGetCell(mxGetField(nurbsstructure, 0, "knots"), 0)), mxGetPr(mxGetCell(mxGetField(nurbsstructure, 0, "knots"), 1)), &mxGetPr(point0)[3*i], matA, paramtmp, pnttmp, pderu, pderv, pderuu, pderuv, pdervv, leftU, rightU, NU, leftV, rightV, NV, matHess, grad, dir, ustmp)){
                            break;
                        }
                    }
                    else{
                        if (updateClosestSrfPnt2Line(degU, degV, mxGetPr(mxGetField(nurbsstructure, 0, "coefs")), ncp, kcp, mxGetPr(mxGetCell(mxGetField(nurbsstructure, 0, "knots"), 0)), mxGetPr(mxGetCell(mxGetField(nurbsstructure, 0, "knots"), 1)), mxGetPr(point0), matA, paramtmp, pnttmp, pderu, pderv, pderuu, pderuv, pdervv, leftU, rightU, NU, leftV, rightV, NV, matHess, grad, dir, ustmp)){
                            break;
                        }
                    }
                    
                    if (j<4){
                        if (j==0){
                            stepl=dir[0]*dir[0]+dir[1]*dir[1];
                        }
                        else {
                            stepl=MAX(stepl,dir[0]*dir[0]+dir[1]*dir[1]);
                        }
                    }
                    else {
                        if ( 1.0e8*(dir[0]*dir[0]+dir[1]*dir[1])<stepl ){
                            break;
                        }
                    }
                    
                }
                
                if (paramtmp[0] < umin){
                    paramtmp[0] = umin;
                }
                else if(paramtmp[0] > umax){
                    paramtmp[0] = umax;
                }
                if (paramtmp[1] < vmin){
                    paramtmp[1] = vmin;
                }
                else if(paramtmp[1] > vmax){
                    paramtmp[1] = vmax;
                }
                BspEval2(degU, degV, mxGetPr(mxGetField(nurbsstructure, 0, "coefs")), 4, ncp, kcp, mxGetPr(mxGetCell(mxGetField(nurbsstructure, 0, "knots"), 0)), mxGetPr(mxGetCell(mxGetField(nurbsstructure, 0, "knots"), 1)), paramtmp, 1, bspPnts, leftU, rightU, NU, leftV, rightV, NV);
                
                mxGetPr(evaluated_points)[3*i]=(bspPnts[0])/(bspPnts[3]);
                mxGetPr(evaluated_points)[3*i+1]=(bspPnts[1])/(bspPnts[3]);
                mxGetPr(evaluated_points)[3*i+2]=(bspPnts[2])/(bspPnts[3]);
                
                mxGetPr(parametervalues)[2*i]=paramtmp[0];
                mxGetPr(parametervalues)[2*i+1]=paramtmp[1];
                
            }
            
        }
        else if(nrhs==5){
            
            for (i = 0; i < mxGetN(initparamvalues); i++){
                
                paramtmp[0] = (double)mxGetPr(initparamvalues)[2*i];
                paramtmp[1] = (double)mxGetPr(initparamvalues)[2*i+1];
                
                for (j = 0; j < MAXITER; j++){
                    
                    if (paramtmp[0] < umin){
                        paramtmp[0] = umin;
                    }
                    else if(paramtmp[0] > umax){
                        paramtmp[0] = umax;
                    }
                    if (paramtmp[1] < vmin){
                        paramtmp[1] = vmin;
                    }
                    else if(paramtmp[1] > vmax){
                        paramtmp[1] = vmax;
                    }
                    
                    nrbD1D2eval2(orderU, orderV, ncp, kcp, nurbsstructure, dnrbsstructure, d2nrbsstructure, paramtmp, pnttmp, pderu, pderv, pderuu, pderuv, pdervv, bspPnts, leftU, rightU, NU, leftV, rightV, NV);
                    
                    if (TrueFalse){
                        if (updateClosestSrfPnt2Pnt(degU, degV, mxGetPr(mxGetField(nurbsstructure, 0, "coefs")), ncp, kcp, mxGetPr(mxGetCell(mxGetField(nurbsstructure, 0, "knots"), 0)), mxGetPr(mxGetCell(mxGetField(nurbsstructure, 0, "knots"), 1)), &mxGetPr(point0)[3*i], paramtmp, pnttmp, pderu, pderv, pderuu, pderuv, pdervv, leftU, rightU, NU, leftV, rightV, NV, matHess, grad, dir, ustmp)){
                            break;
                        }
                    }
                    else{
                        if (updateClosestSrfPnt2Pnt(degU, degV, mxGetPr(mxGetField(nurbsstructure, 0, "coefs")), ncp, kcp, mxGetPr(mxGetCell(mxGetField(nurbsstructure, 0, "knots"), 0)), mxGetPr(mxGetCell(mxGetField(nurbsstructure, 0, "knots"), 1)), mxGetPr(point0), paramtmp, pnttmp, pderu, pderv, pderuu, pderuv, pdervv, leftU, rightU, NU, leftV, rightV, NV, matHess, grad, dir, ustmp)){
                            break;
                        }
                    }
                    
                    if (j<4){
                        if (j==0){
                            stepl=dir[0]*dir[0]+dir[1]*dir[1];
                        }
                        else {
                            stepl=MAX(stepl,dir[0]*dir[0]+dir[1]*dir[1]);
                        }
                    }
                    else {
                        if ( 1.0e8*(dir[0]*dir[0]+dir[1]*dir[1])<stepl ){
                            break;
                        }
                    }
                    
                }
                
                if (paramtmp[0] < umin){
                    paramtmp[0] = umin;
                }
                else if(paramtmp[0] > umax){
                    paramtmp[0] = umax;
                }
                if (paramtmp[1] < vmin){
                    paramtmp[1] = vmin;
                }
                else if(paramtmp[1] > vmax){
                    paramtmp[1] = vmax;
                }
                BspEval2(degU, degV, mxGetPr(mxGetField(nurbsstructure, 0, "coefs")), 4, ncp, kcp, mxGetPr(mxGetCell(mxGetField(nurbsstructure, 0, "knots"), 0)), mxGetPr(mxGetCell(mxGetField(nurbsstructure, 0, "knots"), 1)), paramtmp, 1, bspPnts, leftU, rightU, NU, leftV, rightV, NV);
                
                mxGetPr(evaluated_points)[3*i]=(bspPnts[0])/(bspPnts[3]);
                mxGetPr(evaluated_points)[3*i+1]=(bspPnts[1])/(bspPnts[3]);
                mxGetPr(evaluated_points)[3*i+2]=(bspPnts[2])/(bspPnts[3]);
                
                mxGetPr(parametervalues)[2*i]=paramtmp[0];
                mxGetPr(parametervalues)[2*i+1]=paramtmp[1];
                
            }
            
        }
        else{
            mexPrintf("Error! Illegeal number of inputs.\n");
        }
        
        free(leftU);
        free(rightU);
        free(NU);
        free(leftV);
        free(rightV);
        free(NV);
        
    }
    
    else if (mxGetM(initparamvalues)==1){
        
        int orderU = (int)mxGetScalar(mxGetField(nurbsstructure, 0, "order"));
        
        int ncp = (int)mxGetScalar(mxGetField(nurbsstructure, 0, "number"));
        
        double *left  = (double*) malloc(orderU*sizeof(double)+1);
        double *right = (double*) malloc(orderU*sizeof(double)+1);
        double *N = (double*) malloc(orderU*sizeof(double)+1);
        
        umin = (double)mxGetPr(mxGetField(nurbsstructure, 0, "knots"))[0];
        umax = (double)mxGetPr(mxGetField(nurbsstructure, 0, "knots"))[ncp];
        
        if (nrhs==6){
            
            for (i = 0; i < mxGetN(initparamvalues); i++){
                
                paramtmp[0] = (double)mxGetPr(initparamvalues)[i];
                
                t=0.0;
                for (j = 0; j < MAXITER; j++){
                    
                    if (paramtmp[0] < umin){
                        paramtmp[0] = umin;
                    }
                    else if(paramtmp[0] > umax){
                        paramtmp[0] = umax;
                    }
                    
                    nrbD1D2eval(orderU, ncp, nurbsstructure, dnrbsstructure, d2nrbsstructure, paramtmp, pnttmp, pderu, pderuu, bspPnts, left, right, N);
                    
                    if (TrueFalse){
                        res[0]=mxGetPr(point0)[3*i]+t*mxGetPr(linedirection)[3*i]-pnttmp[0];
                        res[1]=mxGetPr(point0)[3*i+1]+t*mxGetPr(linedirection)[3*i+1]-pnttmp[1];
                        res[2]=mxGetPr(point0)[3*i+2]+t*mxGetPr(linedirection)[3*i+2]-pnttmp[2];
                        
                        matHess[1]=-(pderu[0]*mxGetPr(linedirection)[3*i]+pderu[1]*mxGetPr(linedirection)[3*i+1]+pderu[2]*mxGetPr(linedirection)[3*i+2]);
                        matHess[2]=  mxGetPr(linedirection)[3*i]*mxGetPr(linedirection)[3*i]+mxGetPr(linedirection)[3*i+1]*mxGetPr(linedirection)[3*i+1]+mxGetPr(linedirection)[3*i+2]*mxGetPr(linedirection)[3*i+2];
                        
                        grad[1]=res[0]*mxGetPr(linedirection)[3*i]+res[1]*mxGetPr(linedirection)[3*i+1]+res[2]*mxGetPr(linedirection)[3*i+2];
                    }
                    else{
                        res[0]=mxGetPr(point0)[0]+t*mxGetPr(linedirection)[0]-pnttmp[0];
                        res[1]=mxGetPr(point0)[1]+t*mxGetPr(linedirection)[1]-pnttmp[1];
                        res[2]=mxGetPr(point0)[2]+t*mxGetPr(linedirection)[2]-pnttmp[2];
                        
                        matHess[1]=-(pderu[0]*mxGetPr(linedirection)[0]+pderu[1]*mxGetPr(linedirection)[1]+pderu[2]*mxGetPr(linedirection)[2]);
                        matHess[2]=  mxGetPr(linedirection)[0]*mxGetPr(linedirection)[0]+mxGetPr(linedirection)[1]*mxGetPr(linedirection)[1]+mxGetPr(linedirection)[2]*mxGetPr(linedirection)[2];
                        
                        grad[1]=res[0]*mxGetPr(linedirection)[0]+res[1]*mxGetPr(linedirection)[1]+res[2]*mxGetPr(linedirection)[2];
                    }
                    matHess[0]=pderu[0]*pderu[0]+pderu[1]*pderu[1]+pderu[2]*pderu[2]-(pderuu[0]*res[0]+pderuu[1]*res[1]+pderuu[2]*res[2]);
                    
                    grad[0]=pderu[0]*res[0]+pderu[1]*res[1]+pderu[2]*res[2];
                    
                    detH = matHess[0]*matHess[2]-matHess[1]*matHess[1];
                    if (fabs(detH)<1e-10){
                        break;
                    }
                    ustmp[0] = (matHess[2]*grad[0]+matHess[1]*grad[1])/detH;
                    ustmp[1] = -(matHess[0]*grad[1]+matHess[1]*grad[0])/detH;
                    
                    paramtmp[0] += 0.7*ustmp[0];
                    t+=0.7*ustmp[1];
                    
                    if ((ustmp[0]*ustmp[0]+ustmp[1]*ustmp[1])<1e-20){
                        break;
                    }
                }
                
                if (paramtmp[0] < umin){
                    paramtmp[0] = umin;
                }
                else if(paramtmp[0] > umax){
                    paramtmp[0] = umax;
                }
                BspEval(orderU-1, mxGetPr(mxGetField(nurbsstructure, 0, "coefs")), 4, ncp, mxGetPr(mxGetField(nurbsstructure, 0, "knots")), paramtmp, 1, bspPnts, left, right, N);
                
                mxGetPr(evaluated_points)[3*i]=(bspPnts[0])/(bspPnts[3]);
                mxGetPr(evaluated_points)[3*i+1]=(bspPnts[1])/(bspPnts[3]);
                mxGetPr(evaluated_points)[3*i+2]=(bspPnts[2])/(bspPnts[3]);
                
                mxGetPr(parametervalues)[i]=paramtmp[0];
                
            }
            
        }
        else if(nrhs==5){
            
            for (i = 0; i < mxGetN(initparamvalues); i++){
                
                paramtmp[0] = (double)mxGetPr(initparamvalues)[i];
                
                for (j = 0; j < MAXITER; j++){
                    
                    if (paramtmp[0] < umin){
                        paramtmp[0] = umin;
                    }
                    else if(paramtmp[0] > umax){
                        paramtmp[0] = umax;
                    }
                    
                    nrbD1D2eval(orderU, ncp, nurbsstructure, dnrbsstructure, d2nrbsstructure, paramtmp, pnttmp, pderu, pderuu, bspPnts, left, right, N);
                    
                    if (TrueFalse){
                        res[0]=mxGetPr(point0)[3*i]-pnttmp[0];
                        res[1]=mxGetPr(point0)[3*i+1]-pnttmp[1];
                        res[2]=mxGetPr(point0)[3*i+2]-pnttmp[2];
                    }
                    else{
                        res[0]=mxGetPr(point0)[0]-pnttmp[0];
                        res[1]=mxGetPr(point0)[1]-pnttmp[1];
                        res[2]=mxGetPr(point0)[2]-pnttmp[2];
                    }
                    matHess[0]=pderu[0]*pderu[0]+pderu[1]*pderu[1]+pderu[2]*pderu[2]-(pderuu[0]*res[0]+pderuu[1]*res[1]+pderuu[2]*res[2]);
                    
                    grad[0]=pderu[0]*res[0]+pderu[1]*res[1]+pderu[2]*res[2];
                    
                    if (fabs(matHess[0])<1e-10){
                        break;
                    }
                    ustmp[0] = grad[0]/matHess[0];
                    
                    paramtmp[0] += 0.7*ustmp[0];
                    
                    if ((ustmp[0]*ustmp[0])<1e-20){
                        break;
                    }
                }
                
                if (paramtmp[0] < umin){
                    paramtmp[0] = umin;
                }
                else if(paramtmp[0] > umax){
                    paramtmp[0] = umax;
                }
                BspEval(orderU-1, mxGetPr(mxGetField(nurbsstructure, 0, "coefs")), 4, ncp, mxGetPr(mxGetField(nurbsstructure, 0, "knots")), paramtmp, 1, bspPnts, left, right, N);
                
                mxGetPr(evaluated_points)[3*i]=(bspPnts[0])/(bspPnts[3]);
                mxGetPr(evaluated_points)[3*i+1]=(bspPnts[1])/(bspPnts[3]);
                mxGetPr(evaluated_points)[3*i+2]=(bspPnts[2])/(bspPnts[3]);
                
                mxGetPr(parametervalues)[i]=paramtmp[0];
                
            }
            
        }
        else{
            mexPrintf("Error! Illegeal number of inputs.\n");
        }
        
        free(left);
        free(right);
        free(N);
        
    }
    
    else{
        mexErrMsgTxt("Wrong dimension of UV");
    }
    
}
