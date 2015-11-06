
double pwr2(double a){
    return a*a;
}

static void LScurveAppEqPntStartEnd(int deg, int ncp, double *knot, double *us, int nus, double *qmat, unsigned int *ctrlPeq, int numCtrlPeq, double ctrlPeqWgh, double *ctrlpnts, double *startq, double *endq, double *nmat, double *rmat, double *left, double *right, double *N){
    
    /* Originates from section 9.4.1, The NURBS Book, L.Piegl and W. Tiller */
    /* Least Squares Curve Approximation */
    
    int i, j, jj, span, ind;
    double Rk[2], wghtmp;
    
    if (ctrlPeq[numCtrlPeq-1]==(ncp-1)){
        numCtrlPeq--;
    }
    
    for (jj = 0; jj < nus; jj++){
        
        Rk[0]=qmat[2*jj];
        Rk[1]=qmat[2*jj+1];
        
        if(us[jj]<=knot[deg] || deg==0){
            /*
             *          Have no sense
             * Rk[0]-=startq[0];
             * Rk[1]-=startq[1];
             */
        }
        else if(us[jj]>=knot[ncp]){
            /*
             *          Have no sense
             * Rk[0]-=endq[0];
             * Rk[1]-=endq[1];
             */
        }
        else{
            span = FindSpan(ncp, deg, us[jj], knot);
            BasisFuns(span, us[jj], deg, knot, N, left, right);
            
            ind=span-deg-1;
            
            if (span<(ncp-1) || span>deg){
                if (fabs(knot[span+2]-knot[span+1])<1e-12){
                    if (fabs(knot[span]-knot[span-1])>1e-12){
                        wghtmp=(ctrlPeqWgh-1.0)*pwr2((us[jj]-knot[span])/(knot[span+1]-knot[span]))+1.0;
                    }
                    else {
                        wghtmp=1.0;
                    }
                }
                else if (fabs(knot[span]-knot[span-1])<1e-12){
                    wghtmp=(1.0-ctrlPeqWgh)*pwr2((us[jj]-knot[span])/(knot[span+1]-knot[span]))+ctrlPeqWgh;
                }
                else {
                    wghtmp=1.0;
                }
            }
            else {
                wghtmp=1.0;
            }
            
            if(span==deg){
                Rk[0]-=N[0]*startq[0];
                Rk[1]-=N[0]*startq[1];
            }
            if(span==(ncp-1)){
                Rk[0]-=N[deg]*endq[0];
                Rk[1]-=N[deg]*endq[1];
            }
            
            if(span==deg && span<(ncp-1)){
                for (i = 1; i < deg; i++){
                    for (j = 0; j < numCtrlPeq; j++){
                        if ((i+ind+1)==ctrlPeq[j]){
                            N[i+1]+=N[i];
                            N[i]=0.0;
                            break;
                        }
                    }
                }
                for (j = 0; j < numCtrlPeq; j++){
                    if (span==ctrlPeq[j]){
                        nmat[span+span*(ncp-2)] += wghtmp*N[deg]*N[deg];
                        rmat[span] += wghtmp*N[deg]*Rk[0];
                        rmat[span+ncp-2] += wghtmp*N[deg]*Rk[1];
                        for (i = 0; i < deg; i++){
                            nmat[ind+i+span*(ncp-2)] += wghtmp*N[deg]*N[i];
                            nmat[span+(ind+i)*(ncp-2)] += wghtmp*N[i]*N[deg];
                        }
                        N[deg]=0.0;
                        break;
                    }
                }
                for (i = 1; i <= deg; i++){
                    for (j = 1; j <= deg; j++){
                        nmat[j-1+(i-1)*(ncp-2)] += wghtmp*N[i]*N[j];
                    }
                    rmat[i-1] += wghtmp*N[i]*Rk[0];
                    rmat[i-1+ncp-2] += wghtmp*N[i]*Rk[1];
                }
            }
            else if(span==(ncp-1) && span>deg){
                for (i = 0; i < deg; i++){
                    for (j = 0; j < numCtrlPeq; j++){
                        if ((i+ind+1)==ctrlPeq[j]){
                            N[i+1]+=N[i];
                            N[i]=0.0;
                            break;
                        }
                    }
                }
                for (i = 0; i < deg; i++){
                    for (j = 0; j < deg; j++){
                        nmat[j+ind+(i+ind)*(ncp-2)] += wghtmp*N[i]*N[j];
                    }
                    rmat[i+ind]+=wghtmp*N[i]*Rk[0];
                    rmat[i+ind+ncp-2]+=wghtmp*N[i]*Rk[1];
                }
            }
            else if(span==(ncp-1) && span==deg){
                for (i = 1; i < deg; i++){
                    for (j = 0; j < numCtrlPeq; j++){
                        if ((i+ind+1)==ctrlPeq[j]){
                            N[i+1]+=N[i];
                            N[i]=0.0;
                            break;
                        }
                    }
                }
                for (i = 1; i < deg; i++){
                    for (j = 1; j < deg; j++){
                        nmat[j-1+(i-1)*(ncp-2)] += wghtmp*N[i]*N[j];
                    }
                    rmat[i-1]+=wghtmp*N[i]*Rk[0];
                    rmat[i-1+ncp-2]+=wghtmp*N[i]*Rk[1];
                }
            }
            else{
                for (i = 0; i < deg; i++){
                    for (j = 0; j < numCtrlPeq; j++){
                        if ((i+ind+1)==ctrlPeq[j]){
                            N[i+1]+=N[i];
                            N[i]=0.0;
                            break;
                        }
                    }
                }
                for (j = 0; j < numCtrlPeq; j++){
                    if (span==ctrlPeq[j]){
                        nmat[span+span*(ncp-2)] += wghtmp*N[deg]*N[deg];
                        rmat[span] += wghtmp*N[deg]*Rk[0];
                        rmat[span+ncp-2] += wghtmp*N[deg]*Rk[1];
                        for (i = 0; i < deg; i++){
                            nmat[ind+i+span*(ncp-2)] += wghtmp*N[deg]*N[i];
                            nmat[span+(ind+i)*(ncp-2)] += wghtmp*N[i]*N[deg];
                        }
                        N[deg]=0.0;
                        break;
                    }
                }
                for (i = 0; i <= deg; i++){
                    for (j = 0; j <= deg; j++){
                        nmat[j+ind+(i+ind)*(ncp-2)] += wghtmp*N[i]*N[j];
                    }
                    rmat[i+ind]+=wghtmp*N[i]*Rk[0];
                    rmat[i+ind+ncp-2]+=wghtmp*N[i]*Rk[1];
                }
            }
            
        }
        
    }
    
    for (i = 0; i < (ncp-2); i++){
        if (nmat[i+i*(ncp-2)]<1.0){
            nmat[i+i*(ncp-2)] += 1.0;
            rmat[i] += ctrlpnts[4*(i+1)];
            rmat[i+ncp-2] += ctrlpnts[4*(i+1)+1];
        }
    }
    
}
