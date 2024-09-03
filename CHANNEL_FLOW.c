#include<stdio.h>
#include<stdlib.h>
#include<math.h>
int main()
{
    int i,j,m ,n,iteration=0;
    double dx, dy, re, a,b,c, beta, d, f, error_psi, error_w, x, y,L,H;
    double urf = 0.1; //under relaxation factor
    
    //taking user input for grid size and Reynolds number
    printf("enter the value of m in mxn grid: ");
    scanf("%d", &m);
    printf("\nenter the value of n in mxn grid: ");
    scanf("%d", &n); 
    printf("\nenter the value of L: ");
    scanf("%lf", &L);
    printf("\nenter the value of H: ");
    scanf("%lf", &H);
    printf("\nenter the value of Re: ");
    scanf("%lf", &re);
//defining variables
    double psi_new[m+1][n+1], psi_old[m+1][n+1], w_new[m+1][n+1], w_old[m+1][n+1], u[m+1][n+1], v[m+1][n+1];

    dx = L/(m-1.0);
    dy = H/(n-1.0);
    printf("\n%lf , %lf",dx,dy);
    beta = dx/dy;

    d = -2.0/(dx*dx);
    f = -2.0/(dy*dy);
 
    a = 1.0/(2.0*(1+(beta*beta)));
    b = (beta*re)/4.0;
    c = re/(4.0*beta);
   
//B.C. for velocity and streamline
    for(j=1;j<=n;j++){
        for(i=1;i<=m;i++){
            if(i==1){
                u[i][j] = 1.0;
                v[i][j] = 0.0;
            }
            else if(i==m){
                u[i][j] = u[i-1][j];
                v[i][j] = v[i-1][j];
            }
            else{
                u[i][j] = 0.0;
                v[i][j] = 0.0;
            }
        }
    }
    
    for(j=1; j<=n; j++){          
        for(i=1; i<=m; i++){
            if(j==1){                
                psi_new[i][j] = 0.0;                
            }
            else if(i==1){                
                psi_new[i][j+1] = psi_new[i][j] + (dy*u[i][j]);                
            }
            else if(i==m){
                psi_new[i][j] = psi_new[i-1][j];
            }
            else if(j==n){
                psi_new[i][j] = 1.0;
            }            
            else{
                psi_new[i][j] = 0.0;
            }
        }
    }
//BC for vorticity
    for(j=1; j<=n; j++){
        for(i=1; i<=m; i++){
            if(i==1){
                w_new[i][j] = d*(psi_new[i+1][j] - psi_new[i][j]);                
            }
            else if(j==1){
                w_new[i][j] = f*(psi_new[i][j+1] - psi_new[i][j]);                
            }
            else if(i==m){
                w_new[i][j] = w_new[i-1][j];                
            }
            else if(j==n){
                w_new[i][j] = f*(psi_new[i][j-1] - psi_new[i][j]);                
            }
            else {
                w_new[i][j] = 0.0;                
            }
        }
    }

    do{
        for(j=1; j<=n; j++){                 //storing previous values
            for(i=1; i<=m; i++){
                psi_old[i][j] = psi_new[i][j];
                w_old[i][j] = w_new[i][j]; 
            }
        }
        //B.C. for velocity and streamline
        for(j=1;j<=n;j++){
            for(i=1;i<=m;i++){
                if(i==1){
                    u[i][j] = 1.0;
                    v[i][j] = 0.0;
                }
                else if(i==m){
                    u[i][j] = u[i-1][j];
                    v[i][j] = v[i-1][j];
                }
            }
        }
        
        for(j=1; j<=n; j++){          
            for(i=1; i<=m; i++){
                if(j==1){                
                    psi_new[i][j] = 0.0;                
                }
                else if(i==1){                
                    psi_new[i][j+1] = psi_new[i][j] + (dy*u[i][j]);                
                }
                else if(i==m){
                    psi_new[i][j] = psi_new[i-1][j];
                }
                else if(j==n){
                    psi_new[i][j] = 1.0;
                }            
                
            }
        }
//BC for vorticity
        for(j=1; j<=n; j++){
            for(i=1; i<=m; i++){
                if(i==1){
                    w_new[i][j] = d*(psi_new[i+1][j] - psi_new[i][j]);                
                }
                else if(j==1){
                    w_new[i][j] = f*(psi_new[i][j+1] - psi_new[i][j]);                
                }
                else if(i==m){
                    w_new[i][j] = w_new[i-1][j];                
                }
                else if(j==n){
                    w_new[i][j] = f*(psi_new[i][j-1] - psi_new[i][j]);                
                }
            }
        }

        for(j=2; j<n; j++){             //solving for psi values
            for(i=2; i<m; i++){
                psi_new[i][j] = a*((dx*dx*w_new[i][j]) + (beta*beta*(psi_new[i][j+1] + psi_new[i][j-1])) + (psi_new[i+1][j] + psi_new[i-1][j]));
            }
        }

        for(j=2; j<n; j++){             //solving for vorticity values
            for(i=2; i<m; i++){
                w_new[i][j] = (1.0-urf)*w_old[i][j] + (urf*a)*( (w_new[i+1][j]*(1.0-(b*(psi_new[i][j+1] - psi_new[i][j-1])))) + (w_new[i-1][j]*(1.0+(b*(psi_new[i][j+1] - psi_new[i][j-1])))) + (beta*beta*w_new[i][j+1]*(1.0+(c*(psi_new[i+1][j] - psi_new[i-1][j])))) + (beta*beta*w_new[i][j-1]*(1.0-(c*(psi_new[i+1][j] - psi_new[i-1][j])))) );
            }            
        }

        for(j=1; j<=n; j++){           //updating vorticity values at boundary
            for(i=1; i<=m; i++){
                if(i==1){
                    w_new[i][j] = d*(psi_new[i+1][j] - psi_new[i][j]);                
                }
                else if(j==1){
                    w_new[i][j] = f*(psi_new[i][j+1] - psi_new[i][j]);                
                }
                else if(i==m){
                    w_new[i][j] = w_new[i-1][j];                
                }
                else if(j==n){
                    w_new[i][j] = f*(psi_new[i][j-1] - psi_new[i][j]);                
                }                
            }
        }
//finding error
        error_psi = 0.0;
        error_w = 0.0;

        for(j=2; j<n; j++){
            for(i=2; i<m; i++){
                error_psi = error_psi + pow((psi_new[i][j] - psi_old[i][j]),2.0);
                error_w = error_w + pow((w_new[i][j] - w_old[i][j]),2.0);
            }
        }
        error_psi = sqrt(error_psi/ ((m-2)*(n-2)));
        error_w = sqrt(error_w/ ((m-2)*(n-2)));
        printf("\n\niteration = %d\t", iteration);
        printf("error_psi = %.10lf\t error_w = %.10lf",error_psi, error_w);
        iteration++;      

    }while(error_psi>1e-6 || error_w>1e-6);

// updating velocity values
    for(j=2; j<n; j++){
        for(i=2; i<m; i++){
            u[i][j] = (1.0/(2.0*dy))*(psi_new[i][j+1] - psi_new[i][j-1]);
            v[i][j] =  (-1.0/(2.0*dx))*(psi_new[i+1][j] - psi_new[i-1][j]);
        }
    } 

    FILE *f1;
    f1 = fopen("channel_flow_t3.dat", "w");
    if(f1==NULL){
        printf("\n cannot open file");
        exit(0);
    }   

    FILE *f5;
    f5 = fopen("u_5_t3.dat", "w");
    if(f5==NULL){
        printf("\n cannot open file");
        exit(0);
    }   
    FILE *f6;
    f6 = fopen("u_10_t3.dat", "w");
    if(f6==NULL){
        printf("\n cannot open file");
        exit(0);
    }   
    FILE *f7;
    f7 = fopen("u_15_t3.dat", "w");
    if(f7==NULL){
        printf("\n cannot open file");
        exit(0);
    }   
    FILE *f8;
    f8 = fopen("u_1_t3.dat", "w");
    if(f8==NULL){
        printf("\n cannot open file");
        exit(0);
    }   
//extracting data for u at x = 1, 5, 10, 15
    for(j=1;j<=n;j++){        
        fprintf(f5,"%lf", u[5][j]);
        fprintf(f5,"\t%lf",(j-1)*dy);
        fprintf(f5,"\n");
    }
    for(j=1;j<=n;j++){        
        fprintf(f6,"%lf", u[10][j]);
        fprintf(f6,"\t%lf",(j-1)*dy);
        fprintf(f6,"\n");
    }
    for(j=1;j<=n;j++){        
        fprintf(f7,"%lf", u[15][j]);
        fprintf(f7,"\t%lf",(j-1)*dy);
        fprintf(f7,"\n");
    }
    for(j=1;j<=n;j++){        
        fprintf(f8,"%lf", u[1][j]);
        fprintf(f8,"\t%lf",(j-1)*dy);
        fprintf(f8,"\n");
    }


    fprintf(f1, "ZONE I=%d, J=%d\n",m,n);
    for(j=1; j<=n; j++){
        y = (j-1)*dy;
        for(i=1; i<=m; i++){
            x = (i-1)*dx;
            fprintf(f1, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",x,y,u[i][j],v[i][j],psi_new[i][j],w_new[i][j]);
        }
    }

    fclose(f1);
    fclose(f5);
    fclose(f6);
    fclose(f7);
    fclose(f8);
    return 0;
}