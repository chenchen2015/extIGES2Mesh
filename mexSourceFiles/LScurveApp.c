
static void LScurveApp(int deg, int ncp, double *knot, double *us, int dim, int nus, double *qmat, double *startq, double *endq, double *nmat, double *rmat, double *left, double *right, double *N){
    
    /* Originates from section 9.4.1, The NURBS Book, L.Piegl and W. Tiller */
    /* Least Squares Curve Approximation */
    
    int i, j, jj, span, sd, ncpm2;
    double Rk[4];
    
    ncpm2=ncp-2;
        
    for (jj = 0; jj < nus; jj++){
        
        for (j = 0; j < dim; j++){
            Rk[j]=qmat[dim*jj+j];
        }
        
        if(us[jj]<=knot[deg] || deg==0){
            /*
             *  Have no sense
             */
        }
        else if(us[jj]>=knot[ncp]){
            /*
             *  Have no sense
             */
        }
        else{
            span = FindSpan(ncp, deg, us[jj], knot);
            BasisFuns(span, us[jj], deg, knot, N, left, right);
            
            sd=span-deg-1;
            
            if(span==deg){
                for (j = 0; j < dim; j++){
                    Rk[j]-=N[0]*startq[j];
                }
            }
            if(span==(ncp-1)){
                for (j = 0; j < dim; j++){
                    Rk[j]-=N[deg]*endq[j];
                }
            }
            
            if(span==deg && span<(ncp-1)){
                for (i = 1; i <= deg; i++){
                    for (j = 1; j <= deg; j++){
                        nmat[j-1+(i-1)*ncpm2] += N[i]*N[j];
                    }
                    for (j = 0; j < dim; j++){
                        rmat[i-1+j*ncpm2]+=N[i]*Rk[j];
                    }
                }
            }
            else if(span==(ncp-1) && span>deg){
                for (i = 0; i < deg; i++){
                    for (j = 0; j < deg; j++){
                        nmat[j+sd+(i+sd)*ncpm2] += N[i]*N[j];
                    }
                    for (j = 0; j < dim; j++){
                        rmat[i+sd+j*ncpm2]+=N[i]*Rk[j];
                    }
                }
            }
            else if(span==(ncp-1) && span==deg){
                for (i = 1; i < deg; i++){
                    for (j = 1; j < deg; j++){
                        nmat[j-1+(i-1)*ncpm2] += N[i]*N[j];
                    }
                    for (j = 0; j < dim; j++){
                        rmat[i-1+j*ncpm2]+=N[i]*Rk[j];
                    }
                }
            }
            else{
                for (i = 0; i <= deg; i++){
                    for (j = 0; j <= deg; j++){
                        nmat[j+sd+(i+sd)*ncpm2] += N[i]*N[j];
                    }
                    for (j = 0; j < dim; j++){
                        rmat[i+sd+j*ncpm2]+=N[i]*Rk[j];
                    }
                }
            }
            
        }
        
    }
    
}
