
static void NURBSsurfaceEval(int degU, int degV, double *cp, int ncp, int kcp, double *knotU, double *knotV, double *us, int nus, double *ep, double *leftU, double *rightU, double *NU, double *leftV, double *rightV, double *NV){
    /* Modification of  ALGORITHM A4.3, The NURBS Book, L.Piegl and W. Tiller */
    
    /* Evaluates a NURBS surface at given parameter values */
    
    /* NURBSsurfaceEval( degU - degree of NURBS in u, degV - degree of NURBS in v, cp - pointer to control points, ncp - number of control points in u, kcp - number of control points in v, knotU - pointer to knot sequence in u, knotV - pointer to knot sequence in v, us - pointer to parameter values,  nus - number of parameter values, ep - pointer to evaluated points, leftU - pointer to array for function BasisFuns, rightU - pointer to array for function BasisFuns, NU - pointer to array for function BasisFuns, leftV - pointer to array for function BasisFuns, rightV - pointer to array for function BasisFuns, NV - pointer to array for function BasisFuns) */
    
    int i, j, jj, spanU, spanV, ind, ind2, fourncp, threejj, myint;
    double wgh, NN;
    
    fourncp=4*ncp;
    
    for (jj = 0; jj < nus; jj++){
        
        threejj=3*jj;
        
        if(us[2*jj]<=knotU[degU]){
            
            if(us[2*jj+1]<=knotV[degV]){
                ep[threejj] = cp[0]/cp[3];
                ep[threejj+1] = cp[1]/cp[3];
                ep[threejj+2] = cp[2]/cp[3];
            }
            else if(us[2*jj+1]>=knotV[kcp]){
                myint=fourncp*(kcp-1);
                ep[threejj] = cp[myint]/cp[myint+3];
                ep[threejj+1] = cp[myint+1]/cp[myint+3];
                ep[threejj+2] = cp[myint+2]/cp[myint+3];
            }
            else{
                ep[threejj]=0.0;
                ep[threejj+1]=0.0;
                ep[threejj+2]=0.0;
                wgh=0.0;
                
                spanV = FindSpan(kcp, degV, us[2*jj+1], knotV);
                BasisFuns(spanV, us[2*jj+1], degV, knotV, NV, leftV, rightV);
                
                ind = spanV - degV;
                
                for (i = 0; i <= degV; i++){
                    myint=(ind+i)*fourncp;
                    ep[threejj] += NV[i] * cp[myint];
                    ep[threejj+1] += NV[i] * cp[myint+1];
                    ep[threejj+2] += NV[i] * cp[myint+2];
                    wgh += NV[i] * cp[myint+3];
                }
                ep[threejj]=ep[threejj]/wgh;
                ep[threejj+1]=ep[threejj+1]/wgh;
                ep[threejj+2]=ep[threejj+2]/wgh;
            }
            
        }
        else if(us[2*jj]>=knotU[ncp]){
            
            if(us[2*jj+1]<=knotV[degV]){
                ep[threejj] = cp[fourncp-4]/cp[fourncp-1];
                ep[threejj+1] = cp[fourncp-3]/cp[fourncp-1];
                ep[threejj+2] = cp[fourncp-2]/cp[fourncp-1];
            }
            else if(us[2*jj+1]>=knotV[kcp]){
                myint=fourncp*kcp;
                ep[threejj] = cp[myint-4]/cp[myint-1];
                ep[threejj+1] = cp[myint-3]/cp[myint-1];
                ep[threejj+2] = cp[myint-2]/cp[myint-1];
            }
            else{
                ep[threejj]=0.0;
                ep[threejj+1]=0.0;
                ep[threejj+2]=0.0;
                wgh=0.0;
                
                spanV = FindSpan(kcp, degV, us[2*jj+1], knotV);
                BasisFuns(spanV, us[2*jj+1], degV, knotV, NV, leftV, rightV);
                
                ind = spanV - degV;
                myint=(ind+1)*fourncp;
                
                for (i = 0; i <= degV; i++){
                    ep[threejj] += NV[i] * cp[i*fourncp+myint-4];
                    ep[threejj+1] += NV[i] * cp[i*fourncp+myint-3];
                    ep[threejj+2] += NV[i] * cp[i*fourncp+myint-2];
                    wgh += NV[i] * cp[i*fourncp+myint-1];
                }
                ep[threejj]=ep[threejj]/wgh;
                ep[threejj+1]=ep[threejj+1]/wgh;
                ep[threejj+2]=ep[threejj+2]/wgh;
            }
            
        }
        else{
            
            ep[threejj]=0.0;
            ep[threejj+1]=0.0;
            ep[threejj+2]=0.0;
            wgh=0.0;
            
            if(us[2*jj+1]<=knotV[degV]){
                spanU = FindSpan(ncp, degU, us[2*jj], knotU);
                BasisFuns(spanU, us[2*jj], degU, knotU, NU, leftU, rightU);
                
                ind = spanU - degU;
                
                for (i = 0; i <= degU; i++){
                    ep[threejj] += NU[i] * cp[(ind+i)*4];
                    ep[threejj+1] += NU[i] * cp[(ind+i)*4+1];
                    ep[threejj+2] += NU[i] * cp[(ind+i)*4+2];
                    wgh += NU[i] * cp[(ind+i)*4+3];
                }
                ep[threejj]=ep[threejj]/wgh;
                ep[threejj+1]=ep[threejj+1]/wgh;
                ep[threejj+2]=ep[threejj+2]/wgh;
            }
            else if(us[2*jj+1]>=knotV[kcp]){
                spanU = FindSpan(ncp, degU, us[2*jj], knotU);
                BasisFuns(spanU, us[2*jj], degU, knotU, NU, leftU, rightU);
                
                ind = spanU - degU;
                myint=ind*4+(kcp-1)*fourncp;
                
                for (i = 0; i <= degU; i++){
                    ep[threejj] += NU[i] * cp[i*4+myint];
                    ep[threejj+1] += NU[i] * cp[i*4+myint+1];
                    ep[threejj+2] += NU[i] * cp[i*4+myint+2];
                    wgh += NU[i] * cp[i*4+myint+3];
                }
                ep[threejj]=ep[threejj]/wgh;
                ep[threejj+1]=ep[threejj+1]/wgh;
                ep[threejj+2]=ep[threejj+2]/wgh;
            }
            else{
                spanU = FindSpan(ncp, degU, us[2*jj], knotU);
                BasisFuns(spanU, us[2*jj], degU, knotU, NU, leftU, rightU);
                
                spanV = FindSpan(kcp, degV, us[2*jj+1], knotV);
                BasisFuns(spanV, us[2*jj+1], degV, knotV, NV, leftV, rightV);
                
                ind = spanU - degU;
                ind2 = spanV - degV;
                myint=ind2*fourncp+ind*4;
                
                for (i = 0; i <= degV; i++){
                    for (j = 0; j <= degU; j++){
                        NN=NV[i] * NU[j];
                        ep[threejj] += NN * cp[i*fourncp+j*4+myint];
                        ep[threejj+1] += NN * cp[i*fourncp+j*4+myint+1];
                        ep[threejj+2] += NN * cp[i*fourncp+j*4+myint+2];
                        wgh += NN * cp[i*fourncp+j*4+myint+3];
                    }
                }
                ep[threejj]=ep[threejj]/wgh;
                ep[threejj+1]=ep[threejj+1]/wgh;
                ep[threejj+2]=ep[threejj+2]/wgh;
            }
            
        }
        
    }
    
}
