#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string.h> //added to compare strings



int ncord(int l, int i, int j, int ix, int iy){
        int nc,i1,i2;

        i1=i+ix;
        i2=j+iy;

        /*reflective boundaries*/
        if (i1>=l) i1=l-1;
        if (i1<0) i1 = 0;



        nc=i1+i2*l;

        return nc;
}



void grow(int ipn1, int ip_togrow, int L, double *Cs, double *Gi, double *S,double *nk, double *V, double *Ik, double *ThT, double *R, double *Rg, double *Ht, double *H, double *distances, double *exp_gammafl_dist, double *exp_gammadif_dist, double gamma_fl, double a_fl, double gamma_dif, double a_dif ){
   
    //grow towards the right
    
    int i;


    Cs[ipn1]=1;
    Gi[ipn1]=Gi[ip_togrow];
    Ik[ipn1]=Ik[ip_togrow];
    S[ipn1]=S[ip_togrow];
    nk[ipn1]=nk[ip_togrow];
    V[ipn1]=V[ip_togrow];
    ThT[ipn1]=ThT[ip_togrow];
    Ht[ipn1]=Ht[ip_togrow];
    H[ipn1]=H[ip_togrow];
    R[ipn1]=R[ip_togrow];
    Rg[ipn1]=Rg[ip_togrow];
    //slowv[ipn1]=slowv[ip_togrow];

    
    if (ipn1>ip_togrow){ //biofilm 1
    for (i=0;i<ipn1;i++){
        distances[i]+=1; //because it is 1D, i==ip
        exp_gammafl_dist[i]=(exp(-gamma_fl*distances[i])+a_fl)/(1+a_fl);
        exp_gammadif_dist[i]=(exp(-gamma_dif*distances[i])+a_dif)/(1+a_dif);     
    }

    }
    else{ //biofilm 2

    for (i=ipn1+1;i<=L;i++){
        distances[i]+=1; //because it is 1D, i==ip
        exp_gammafl_dist[i]=(exp(-gamma_fl*distances[i])+a_fl)/(1+a_fl);
        exp_gammadif_dist[i]=(exp(-gamma_dif*distances[i])+a_dif)/(1+a_dif);
    }
    }

    
}

double laplacian (double *conc, int ip, int *n1, int *n3){
      //laplacian
    return (-2*conc[ip] + conc[n1[ip]]+ conc[n3[ip]]);
}

double Gin(double V, double V0, double gv){
    return 1/(1+exp(gv*(V-V0)));
}

    
double R_S_f(double *Gi, double *S, int ip, double S0, double ug, double Gs0ug, double gamma_s, double *V){
    //S01/(1+pow(A[ip],ua)/As0ua)+
    
    
    return S0/(1+pow(Gi[ip],ug)/Gs0ug)-gamma_s*S[ip]; 
}
    

double R_G_f(double *G, double Ginip, int ip, double alpha_gt, double Kgt){
    return -Ginip*alpha_gt*G[ip]/(Kgt+G[ip]);
}
    
double R_Gi_f_noentry(double *Gi, double *R, int ip, double delta_g, double *V){
    //return Ginip*alpha_gt*G[ip]/(Kgt+G[ip])-(alpha_a*Gi[ip])-delta_g*Gi[ip]*R[ip]-delta_g2*Gi[ip]-delta_g3*Gi[ip]*1/(1+exp(g_dg*(Vth_dg-V[ip])));
    return -delta_g*Gi[ip]*R[ip];

}

double A_production(double *Gi, double H, double alpha_a, int ip){
	return alpha_a * H * Gi[ip];
}
    
double R_R_f(double *A, double *Gi,double *R, int ip, double beta_r, double gamma_r){
    return (beta_r*A[ip]*Gi[ip])-gamma_r*R[ip];
}

double R_rg_f(double *Rg, double *R, int ip, double alpha_Rg, double delta_Rg, double alpha_rRg, double Krgnrg, double nrg){
    double aux;
    aux=pow(R[ip],nrg);
    return alpha_Rg-delta_Rg*Rg[ip]+alpha_rRg*aux/(Krgnrg+aux);

}
double R_Ht_f(double *H, double *Gi, double *Ht, int ip, double alpha_htg,double alpha_hta, double gamma_ht, double K_htA, double K_htG, int nht){
    return alpha_htg*(1/(1+pow(Gi[ip]/K_htG, nht))) -gamma_ht*(Ht[ip]);
}


double R_H_f_act(double *Gi, double *Ht, double *H, int ip, double alpha_h, double gamma_h, double K_hgn, double n_hg){
    return alpha_h*pow(Gi[ip],n_hg)*(Ht[ip])/(K_hgn+pow(Gi[ip],n_hg)); //-gamma_h*H[ip];
}
double R_H_f_deact(double *H, double gamma_h, int ip){
    return gamma_h*H[ip];
}

double R_A_deg(double *A, double *R, int ip, double delta_a){
    return -delta_a*A[ip]*R[ip];
}

double R_V_f(double *nk, double *V, double Vks, double Vls, int ip, double gk, double gl){
    return -gk*pow(nk[ip],4)*(V[ip]-Vks)-gl*(V[ip]-Vls);
}
           
double R_nk_f(double aS,double *nk,int ip, double b){
    return aS*(1-nk[ip])-b*nk[ip];
}
            
double R_Ek_f(double *nk, double *V, double *Ek, double *Ik, double Ikmax, int ip, double Dp, double F_gk, double Vks, double *Gi){
    return F_gk*pow(nk[ip],4)*(V[ip]-Vks)-Dp*Ek[ip]*(Ikmax-Ik[ip])*Gi[ip];
}

double R_T_f(double *V, double *ThT, double alpha_t, double V0, double gamma_t, double V0_tht, int ip, double g_tht){
    
    return alpha_t/(1+exp(g_tht*(V[ip]-(V0_tht))))-gamma_t*ThT[ip];
}


double sum_ar (double *x, int N){
    double s=0;
    int k;
    for (k = 0; k < N; k++){
       s +=  x[k];
    }
    return s;
}
int check_nan (double *x, int N){
    int nanfound=0;
    int k;
    k=0;
    while ((k<N) && (nanfound==0)){
        if ( isnan(x[k]) != 0){
            nanfound=1;
        }
        k++;
    }
    return nanfound;
}

    //#functions

double flow_G(double *G, int ip, double fl_G, double G_media){
    return fl_G*(G_media-G[ip]);
}
    
double flow_Ek(double *Ek,int ip, double fl_Ek, double K_media){
    return fl_Ek*(K_media-Ek[ip]);
}

double flow_A(double *A,int ip, double fl_A, double A_media){
    return fl_A*(A_media-A[ip]);
}

double DG1_f(double D_gdx, double *G, int ip, int *n1, int *n3 ){
    return D_gdx*laplacian(G, ip, n1, n3 );
}

double DK1_f(double D_kedx, double *Ek, int ip, int *n1, int *n3 ){
    return D_kedx*laplacian(Ek, ip, n1, n3 );
}
    
double DA1_f(double D_adx, double *A, int ip, int *n1, int *n3){
return D_adx*laplacian(A, ip, n1, n3);
}



int main(int argc, char *argv[]){
    
    int i, j, ijk, ip, it;
    double t;
    double *G,*G_next,*R_G_k, *G_k1, *G_k2, *G_k3, *G_k4, *D_G1_k; 
    double *Cs,*Cs_next;
    double *Gi,*Gi_next,*R_Gi_k, *Gi_k1, *Gi_k2, *Gi_k3, *Gi_k4; 
    double *S,*S_next,*R_S_k, *S_k1, *S_k2, *S_k3, *S_k4; 
    double *R,*R_next,*R_R_k, *R_k1, *R_k2, *R_k3, *R_k4; 
    double *Rg,*Rg_next,*R_Rg_k, *Rg_k1, *Rg_k2, *Rg_k3, *Rg_k4; 
    double *A, *A_next, *R_A_k, *A_k1, *A_k2, *A_k3, *A_k4, *D_A1_k;
    double *Ht,*Ht_next,*R_Ht_k, *Ht_k1, *Ht_k2, *Ht_k3, *Ht_k4; 
    double *H,*H_next,*R_H_k, *H_k1, *H_k2, *H_k3, *H_k4; 

    

    //double A_k1, A_k2, A_k3, A_k4;
    double *V,*V_next,*R_V_k,*V_k1,*V_k2, *V_k3, *V_k4;
    double *nk,*nk_next,*R_nk_k,*nk_k1,*nk_k2,*nk_k3, *nk_k4;
    double *Ek,*Ek_next,*R_Ek_k,*Ek_k1, *Ek_k2, *Ek_k3, *Ek_k4, *D_Ek1_k;
    double *ThT, *ThT_next, *R_ThT_k, *ThT_k1, *ThT_k2, *ThT_k3, *ThT_k4;
    double *Ik, *Ik_next, *R_Ik_k,*Ik_k1, *Ik_k2, *Ik_k3, *Ik_k4; double *Ginip_print;
    double Vls,Vks,aS;
    double *Vls_pr, *Vks_pr;
    double *SG_pr, *SA_pr;
    double *Ikmax_pr;

    double *distances;
    double *exp_gammafl_dist;
    double *exp_gammadif_dist;

    double *exp_gammaflG_dist;
    double *exp_gammaflEk_dist;
    double *exp_gammadifG_dist;
    double *exp_gammadifEk_dist;
    
    
    /*---- Simulation parameters ---- */

    //be very careful with the ordering!!! if a new parameter is added, it has to be modified here and in the list of variables and pointers below
    
    const char* parnames[] = {"g_tht", "Krg", "nrg", "alpha_Rg", "delta_Rg", "alpha_rRg", "A0", "a_fl", "a_dif", "Gi0", "gamma_dif", "gamma_fl", "std_ic","K_hg", "n_hg", "ht0","h0", "alpha_h", "gamma_h", "alpha_htg","alpha_hta","gamma_ht", "K_htA","K_htG","nht", "Kgt", "dperturbK", "tperturbK", "K_perturb", "tperturbG", "GE_perturb", "V0_tht", "dl", "alpha_t", "gamma_t",  "stopflow1", "stopflow2", "Da", "dx", "tt", "dt", "st", "zero_tol", "gv", "Radius2", "fl_Aext", "Radius", "freeze", "mk", "Ikmax0", "gl", "gk", "fl_Gext", "A_media", "V0", "Dg", "F", "Pgrow", "alpha_gt", "freeze2", "tprint", "Dp", "gamma_r", "Sth", "alpha_a", "K_media", "sameinitial", "a0", "GE", "beta_r", "center1", "read_onlyc", "delta_g", "center2", "delta_a", "seed_", "GS0", "b", "read_input", "S0", "Vk", "Vl", "Dke", "unfreeze", "gamma_s", "ug", "fl_Ekext", "Ls"};

    int npars = 88;

    double *parvals;
    parvals=(double *)calloc(npars,sizeof(double));

    FILE* in_file = fopen(argv[7], "r"); // read only 
        
    if (! in_file ) // equivalent to saying if ( in_file == NULL )
    { 
        printf("Parameter file can't be read. Exiting. \n");
        exit(-1);
    }

    printf("reading\n");
    fflush(stdout);
    char line[256];
    int lnum = 0;
    
    while(fgets(line, 256, in_file) != NULL)
    {
      
       char name[256];
       double value;
       //printf("We just read %s", line);
       sscanf(line, "%s %lg", name, &value);
       //printf("%d, name %s, value %g, theoretical name %s\n", lnum, name, value, parnames[lnum]);
       
       
       
       if (strcmp(name, parnames[lnum])==0)
       {
        parvals[lnum]=(double)value;
       }else{
        printf("Bad parameter file. Exiting\n");
        exit(-1);
          
       }
               
    lnum++;
               
    }
             
   
    //for (i=0;i<npars;i++){
    //    printf("%s,%2.5f\n",parnames[i],parvals[i]);
    //    fflush(stdout);
    //}

    double g_tht, Krg, nrg, alpha_Rg, delta_Rg, alpha_rRg, A0, a_fl, a_dif, Gi0, gamma_dif, gamma_fl, std_ic, K_hg, n_hg, ht0, h0, alpha_h, gamma_h, alpha_htg, alpha_hta, gamma_ht, K_htA, K_htG, nht, Kgt, dperturbK, tperturbK, K_perturb, tperturbG, GE_perturb, V0_tht, dl, alpha_t, gamma_t, stopflow1, stopflow2, Da, dx, tt, dt, st, zero_tol, gv, Radius2, fl_Aext, Radius, freeze, mk, Ikmax0, gl, gk, fl_Gext, A_media, V0, Dg, F, Pgrow, alpha_gt, freeze2, tprint, Dp, gamma_r, Sth, alpha_a, K_media, sameinitial, a0, GE, beta_r, center1, read_onlyc, delta_g, center2, delta_a, seed_, GS0, b, read_input, S0, Vk, Vl, Dke, unfreeze, gamma_s, ug, fl_Ekext, Ls;
    // this is printed in the notebook, in the same order as in the parnames list above

    //list of pointers to the variables. & operator denotes an address in memory
    double* readpars[]={&g_tht, &Krg, &nrg, &alpha_Rg, &delta_Rg, &alpha_rRg, &A0, &a_fl, &a_dif, &Gi0, &gamma_dif, &gamma_fl, &std_ic, &K_hg, &n_hg, &ht0, &h0, &alpha_h, &gamma_h, &alpha_htg, &alpha_hta, &gamma_ht, &K_htA, &K_htG, &nht, &Kgt, &dperturbK, &tperturbK, &K_perturb, &tperturbG, &GE_perturb, &V0_tht, &dl, &alpha_t, &gamma_t, &stopflow1, &stopflow2, &Da, &dx, &tt, &dt, &st, &zero_tol, &gv, &Radius2, &fl_Aext, &Radius, &freeze, &mk, &Ikmax0, &gl, &gk, &fl_Gext, &A_media, &V0, &Dg, &F, &Pgrow, &alpha_gt, &freeze2, &tprint, &Dp, &gamma_r, &Sth, &alpha_a, &K_media, &sameinitial, &a0, &GE, &beta_r, &center1, &read_onlyc, &delta_g, &center2, &delta_a, &seed_, &GS0, &b, &read_input, &S0, &Vk, &Vl, &Dke, &unfreeze, &gamma_s, &ug, &fl_Ekext, &Ls};
             
    printf("assigning values read to variables\n");
    fflush(stdout);
 
    for (i=0;i<npars;i++){
        *readpars[i]=parvals[i]; //assigns the value to the correspondent variable, using the pointer. * returns the value of the variable located at the address specified by its operand.
    }

    printf("Ls is %g", Ls);
    fflush(stdout);

    double Ginip;
    double Ikmax=Ikmax0;
   

    double fl_Ek=fl_Ekext;
    double fl_G=fl_Gext;
    double fl_A=fl_Aext;


    double G_init=GE;
    double G_media=GE;

    double Sthm=pow(Sth, mk);
    double Gs0ug=pow(GS0, ug);

    double K_hgn=pow(K_hg, n_hg);
    double Krgnrg=pow(Krg, nrg);



    // Interval size in xy-direction. changed to x by RMC. h in the original.
    //double R_init=1.0;                             // Initial cell density
    double dth = dt*0.5;
    double dt6=dt/6.0;
    double dx2=1/(pow(dx,2));
    double D_adx=Da*dx2;
    double D_gdx=Dg*dx2;
    double D_kedx=Dke*dx2;
    
    //double dx12=dx2/12.0;
    

    
    int nme =  (int)(tt/dt+0.5);
    int nst = (int)(st/dt+0.5);
    int npat = (int)(tt/st+0.5);
    
    
    int *n1,*n3;
    //int *n1_n1,*n3_n3;

    printf("total time (tt)=%g, dt=%g,dx=%g,dx2=%g \n",tt,dt,dx,dx2);
    fflush(stdout); 
    
    
   
    
    gsl_rng * prng = gsl_rng_alloc(gsl_rng_mt19937); //pointer to instance of random number generator of default type

    
    

    //model parameters:
    
   
    int L=(int)(Ls/dx);//number of cells
    printf("Ls is %g", Ls);
    printf("dx is %g", dx);
    printf("L=%g",L);
    fflush(stdout);

    int N=L*1; //2n
    int seed = (int) seed_;
   
     
    int i_center = (int) (center1/dx);
    int i_center2=(int) (center2/dx);
    int j_center = 0;
    int j_max=1;

    int Rn=(int)(Radius/dx);
    int Rn2=(int)(Radius2/dx);

    //edges of both biofilms
    int left_e = i_center-Rn;//left edge of biofilm.
    int left_e_2= i_center2-Rn2;
    int right_e=i_center+Rn;
    int right_e_2=i_center2+Rn2;
   
    
    int center_chamber=(int)(L/2);
    gsl_rng_set(prng, seed);
    
    double Pgrow_norm = Pgrow * dt;
    
   
    Cs = (double *)calloc(N,sizeof(double));
    Cs_next = (double *)calloc(N,sizeof(double));

    G = (double *)calloc(N,sizeof(double));
    G_next = (double *)calloc(N,sizeof(double));
    
    Gi = (double *)calloc(N,sizeof(double));
    Gi_next = (double *)calloc(N,sizeof(double));

    Ginip_print = (double *)calloc(N,sizeof(double));
    
    R_G_k = (double *)calloc(N,sizeof(double));
    D_G1_k = (double *)calloc(N,sizeof(double));
    G_k1 = (double *)calloc(N,sizeof(double));
    G_k2 = (double *)calloc(N,sizeof(double));
    G_k3 = (double *)calloc(N,sizeof(double));
    G_k4 = (double *)calloc(N,sizeof(double));
    

    R_Gi_k = (double *)calloc(N,sizeof(double));
    Gi_k1 = (double *)calloc(N,sizeof(double));
    Gi_k2 = (double *)calloc(N,sizeof(double));
    Gi_k3 = (double *)calloc(N,sizeof(double));
    Gi_k4 = (double *)calloc(N,sizeof(double));
    
    R=(double *)calloc(N,sizeof(double));
    R_next = (double *)calloc(N,sizeof(double));
    R_R_k = (double *)calloc(N,sizeof(double));
    R_k1 = (double *)calloc(N,sizeof(double));
    R_k2 = (double *)calloc(N,sizeof(double));
    R_k3 = (double *)calloc(N,sizeof(double));
    R_k4 = (double *)calloc(N,sizeof(double));

    Rg=(double *)calloc(N,sizeof(double));
    Rg_next = (double *)calloc(N,sizeof(double));
    R_Rg_k = (double *)calloc(N,sizeof(double));
    Rg_k1 = (double *)calloc(N,sizeof(double));
    Rg_k2 = (double *)calloc(N,sizeof(double));
    Rg_k3 = (double *)calloc(N,sizeof(double));
    Rg_k4 = (double *)calloc(N,sizeof(double));
   
    A= (double *)calloc(N,sizeof(double));
    A_next= (double *)calloc(N,sizeof(double));
    R_A_k=(double *)calloc(N,sizeof(double));
    A_k1= (double *)calloc(N,sizeof(double));
    A_k2= (double *)calloc(N,sizeof(double));
    A_k3=(double *)calloc(N,sizeof(double));
    A_k4=(double *)calloc(N,sizeof(double));
    D_A1_k=(double *)calloc(N,sizeof(double));
    
    Ht= (double *)calloc(N,sizeof(double));
    Ht_next= (double *)calloc(N,sizeof(double));
    R_Ht_k=(double *)calloc(N,sizeof(double));
    Ht_k1= (double *)calloc(N,sizeof(double));
    Ht_k2= (double *)calloc(N,sizeof(double));
    Ht_k3=(double *)calloc(N,sizeof(double));
    Ht_k4=(double *)calloc(N,sizeof(double));

    H= (double *)calloc(N,sizeof(double));
    H_next= (double *)calloc(N,sizeof(double));
    R_H_k=(double *)calloc(N,sizeof(double));
    H_k1= (double *)calloc(N,sizeof(double));
    H_k2= (double *)calloc(N,sizeof(double));
    H_k3=(double *)calloc(N,sizeof(double));
    H_k4=(double *)calloc(N,sizeof(double));
    

    S=(double *)calloc(N,sizeof(double));
    S_next=(double *)calloc(N,sizeof(double));
    
    R_S_k=(double *)calloc(N,sizeof(double));
    S_k1=(double *)calloc(N,sizeof(double));
    S_k2=(double *)calloc(N,sizeof(double));    
    S_k3=(double *)calloc(N,sizeof(double));
    S_k4=(double *)calloc(N,sizeof(double));    

        
    V=(double *)calloc(N,sizeof(double));
    V_next=(double *)calloc(N,sizeof(double));
    R_V_k=(double *)calloc(N,sizeof(double));
    V_k1=(double *)calloc(N,sizeof(double));
    V_k2=(double *)calloc(N,sizeof(double));
    V_k3=(double *)calloc(N,sizeof(double));
    V_k4=(double *)calloc(N,sizeof(double));
    
    nk=(double *)calloc(N,sizeof(double));
    nk_next=(double *)calloc(N,sizeof(double));
    R_nk_k=(double *)calloc(N,sizeof(double));
    nk_k1=(double *)calloc(N,sizeof(double));
    nk_k2=(double *)calloc(N,sizeof(double));
    nk_k3=(double *)calloc(N,sizeof(double));
    nk_k4=(double *)calloc(N,sizeof(double));
    
    Ek=(double *)calloc(N,sizeof(double));
    Ek_next=(double *)calloc(N,sizeof(double));
    R_Ek_k=(double *)calloc(N,sizeof(double));
    D_Ek1_k=(double *)calloc(N,sizeof(double));
    Ek_k1=(double *)calloc(N,sizeof(double));
    Ek_k2=(double *)calloc(N,sizeof(double));
    Ek_k3=(double *)calloc(N,sizeof(double));
    Ek_k4=(double *)calloc(N,sizeof(double));

    ThT=(double *)calloc(N,sizeof(double));
    ThT_next=(double *)calloc(N,sizeof(double));
    R_ThT_k=(double *)calloc(N,sizeof(double));
    ThT_k1=(double *)calloc(N,sizeof(double));
    ThT_k2=(double *)calloc(N,sizeof(double));
    ThT_k3=(double *)calloc(N,sizeof(double));
    ThT_k4=(double *)calloc(N,sizeof(double));


    Ik=(double *)calloc(N,sizeof(double));
    Ik_next=(double *)calloc(N,sizeof(double));
    R_Ik_k=(double *)calloc(N,sizeof(double));
    Ik_k1=(double *)calloc(N,sizeof(double));
    Ik_k2=(double *)calloc(N,sizeof(double));
    Ik_k3=(double *)calloc(N,sizeof(double));
    Ik_k4=(double *)calloc(N,sizeof(double));

    Ikmax_pr=(double *)calloc(N,sizeof(double));
    double *Dp_v, *Dp_v_next, *R_Dp_k, *Dp_k1, *Dp_k2, *Dp_k3, *Dp_k4;

    Dp_v=(double *)calloc(N,sizeof(double));
    Dp_v_next=(double *)calloc(N,sizeof(double));
    R_Dp_k=(double *)calloc(N,sizeof(double));
    Dp_k1=(double *)calloc(N,sizeof(double));
    Dp_k2=(double *)calloc(N,sizeof(double));
    Dp_k3=(double *)calloc(N,sizeof(double));
    Dp_k4=(double *)calloc(N,sizeof(double));
     

    Vls_pr=(double *)calloc(N,sizeof(double));
    Vks_pr=(double *)calloc(N,sizeof(double));

    
    n1=(int *)calloc(N,sizeof(int));
    n3=(int *)calloc(N,sizeof(int));

    distances=(double *)calloc(N,sizeof(double));
    exp_gammafl_dist=(double *)calloc(N,sizeof(double));
    exp_gammadif_dist=(double *)calloc(N,sizeof(double));

    exp_gammaflG_dist=(double *)calloc(N,sizeof(double));
    exp_gammadifG_dist=(double *)calloc(N,sizeof(double));

    exp_gammaflEk_dist=(double *)calloc(N,sizeof(double));
    exp_gammadifEk_dist=(double *)calloc(N,sizeof(double));
    
    

    for (i=0;i<L;i++){
        for (j=0;j<j_max;j++){
            ijk=ncord(L,i,j,0,0);
            n1[ijk]=ncord(L,i,j,1,0);
            //n1_n1[ijk]=ncord(L,i,j,2,0);
            
            n3[ijk]=ncord(L,i,j,-1,0);
            //n3_n3[ijk]=ncord(L,i,j,-2,0);
                
        }
    }    
    
    
    FILE *outpat_c, *outpat_S, *outpat_g, *outpat_a, *outpat_r; 
    FILE *outpat_v, *outpat_nk, *outpat_Ek, *outpat_Ik; 
    FILE *outpat_gi,*outpat_ginip, *outpat_vls, *outpat_vks;
    FILE *outpat_tht;
    FILE *outpat_Dp;
    FILE *outpat_ht, *outpat_h;
    FILE *outpat_dist, *outpat_fl_dist, *outpat_dif_dist;
    FILE *outpat_Rg, *outpat_slowv;
    
    FILE *outpat_c_pr, *outpat_g_pr,*outpat_gi_pr ;
    FILE *outpat_a_pr, *outpat_r_pr;
    FILE *outpat_S_pr;
    FILE *outpat_v_pr, *outpat_nk_pr, *outpat_Ek_pr, *outpat_Ik_pr; 

    FILE *outpat_ikmax;
    FILE *outpat_tht_pr;
    FILE *outpat_Dp_pr;
    FILE *outpat_ht_pr, *outpat_h_pr;

    FILE *inpat_c, *inpat_g, *inpat_gi;
    FILE *inpat_a, *inpat_r;
    FILE *inpat_S;
    FILE *inpat_v, *inpat_nk, *inpat_Ek, *inpat_Ik; 
    FILE *inpat_tht;
    FILE *inpat_Dp;
    FILE *inpat_ht, *inpat_h;
    
    char basepath[200];
    sprintf(basepath,"%s/%s", argv[1], argv[2]); //order changed 18th august 2017

    char basepath_in[200];
    sprintf(basepath_in,"%s/%s", argv[3], argv[4]);

    char filepath[200];
    
    sprintf(filepath,"%s_c.pat",basepath);
    outpat_c = fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_c);
    fwrite(&npat,sizeof(int),1,outpat_c);
        
    sprintf(filepath,"%s_g.pat",basepath);
    outpat_g = fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_g);
    fwrite(&npat,sizeof(int),1,outpat_g);
    
    sprintf(filepath,"%s_gi.pat",basepath);
    outpat_gi = fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_gi);
    fwrite(&npat,sizeof(int),1,outpat_gi);


    sprintf(filepath,"%s_a.pat",basepath);
    outpat_a = fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_a);
    fwrite(&npat,sizeof(int),1,outpat_a);

    sprintf(filepath,"%s_ht.pat",basepath);
    outpat_ht = fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_ht);
    fwrite(&npat,sizeof(int),1,outpat_ht);

    sprintf(filepath,"%s_h.pat",basepath);
    outpat_h = fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_h);
    fwrite(&npat,sizeof(int),1,outpat_h);


    sprintf(filepath,"%s_r.pat",basepath);
    outpat_r = fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_r);
    fwrite(&npat,sizeof(int),1,outpat_r);

    sprintf(filepath,"%s_ginip.pat",basepath);
    outpat_ginip = fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_ginip);
    fwrite(&npat,sizeof(int),1,outpat_ginip);

    sprintf(filepath,"%s_S.pat",basepath);
    outpat_S= fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_S);
    fwrite(&npat,sizeof(int),1,outpat_S);

    sprintf(filepath,"%s_Rg.pat",basepath);
    outpat_Rg= fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_Rg);
    fwrite(&npat,sizeof(int),1,outpat_Rg);

    sprintf(filepath,"%s_Dp.pat",basepath);
    outpat_Dp= fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_Dp);
    fwrite(&npat,sizeof(int),1,outpat_Dp);
    

    sprintf(filepath,"%s_v.pat",basepath);
    outpat_v = fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_v);
    fwrite(&npat,sizeof(int),1,outpat_v);

    sprintf(filepath,"%s_vls.pat",basepath);
    outpat_vls = fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_vls);
    fwrite(&npat,sizeof(int),1,outpat_vls);

    sprintf(filepath,"%s_vks.pat",basepath);
    outpat_vks = fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_vks);
    fwrite(&npat,sizeof(int),1,outpat_vks);
    
    sprintf(filepath,"%s_nk.pat",basepath);
    outpat_nk= fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_nk);
    fwrite(&npat,sizeof(int),1,outpat_nk);
    
    sprintf(filepath,"%s_Ek.pat",basepath);
    outpat_Ek= fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_Ek);
    fwrite(&npat,sizeof(int),1,outpat_Ek);

    sprintf(filepath,"%s_tht.pat",basepath);
    outpat_tht= fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_tht);
    fwrite(&npat,sizeof(int),1,outpat_tht);

    sprintf(filepath,"%s_Ik.pat",basepath);
    outpat_Ik= fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_Ik);
    fwrite(&npat,sizeof(int),1,outpat_Ik);

    sprintf(filepath,"%s_ikmax.pat",basepath);
    outpat_ikmax= fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_ikmax);
    fwrite(&npat,sizeof(int),1,outpat_ikmax);

    sprintf(filepath,"%s_dist.pat",basepath);
    outpat_dist = fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_dist);
    fwrite(&npat,sizeof(int),1,outpat_dist);

    sprintf(filepath,"%s_expfl_dist.pat",basepath);
    outpat_fl_dist = fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_fl_dist);
    fwrite(&npat,sizeof(int),1,outpat_fl_dist);

    sprintf(filepath,"%s_expdif_dist.pat",basepath);
    outpat_dif_dist = fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_dif_dist);
    fwrite(&npat,sizeof(int),1,outpat_dif_dist);
    
    
    
    /*initial conditions*/
    
    printf("going to read ic");
    fflush(stdout);
    

    if (read_input>0.0){    

        sprintf(filepath,"%s_c_pr.inp",basepath_in);
        inpat_c = fopen(filepath,"rb");
        printf("opened input\n");
        fflush(stdout);
        fread(Cs,sizeof(double),N,inpat_c);
        printf("red Cs\n");
        fflush(stdout);
        
        sprintf(filepath,"%s_g_pr.inp",basepath_in);
        inpat_g = fopen(filepath,"rb");
        fread(G,sizeof(double),N,inpat_g);

        printf("g**\n");
        fflush(stdout);

        sprintf(filepath,"%s_S_pr.inp",basepath_in);
        inpat_S = fopen(filepath,"rb");
        fread(S,sizeof(double),N,inpat_S);
        
        sprintf(filepath,"%s_v_pr.inp",basepath_in);
        inpat_v = fopen(filepath,"rb");
        fread(V,sizeof(double),N,inpat_v);
        
        sprintf(filepath,"%s_nk_pr.inp",basepath_in);
        inpat_nk = fopen(filepath,"rb");
        fread(nk,sizeof(double),N,inpat_nk);

        sprintf(filepath,"%s_h_pr.inp",basepath_in);
        inpat_h = fopen(filepath,"rb");
        fread(H,sizeof(double),N,inpat_h);

        sprintf(filepath,"%s_ht_pr.inp",basepath_in);
        inpat_ht = fopen(filepath,"rb");
        fread(Ht,sizeof(double),N,inpat_ht);
        

        sprintf(filepath,"%s_a_pr.inp",basepath_in);
        inpat_a = fopen(filepath,"rb");
        fread(A,sizeof(double),N,inpat_a);

        sprintf(filepath,"%s_r_pr.inp",basepath_in);
        inpat_r = fopen(filepath,"rb");
        fread(R,sizeof(double),N,inpat_r);
        
        sprintf(filepath,"%s_Ik_pr.inp",basepath_in);
        inpat_Ik = fopen(filepath,"rb");
        fread(Ik,sizeof(double),N,inpat_Ik);

        sprintf(filepath,"%s_Ek_pr.inp",basepath_in);
        inpat_Ek = fopen(filepath,"rb");
        fread(Ek,sizeof(double),N,inpat_Ek);

        sprintf(filepath,"%s_gi_pr.inp",basepath_in);
        inpat_gi = fopen(filepath,"rb");
        fread(Gi,sizeof(double),N,inpat_gi);

        sprintf(filepath,"%s_tht_pr.inp",basepath_in);
        inpat_tht = fopen(filepath,"rb");
        fread(ThT, sizeof(double),N,inpat_tht);

        printf("now Sp");
        fflush(stdout);

        sprintf(filepath,"%s_Dp_pr.inp",basepath_in);
        inpat_Dp = fopen(filepath,"rb");
        fread(Dp_v, sizeof(double),N,inpat_Dp);


    
    }else{
    
        if (sameinitial>0.0){
            for (i=0; i<L; i++){
                for (j=0; j<j_max; j++){



                    if ((i>=left_e) && (i<=right_e) && (Rn>0)){

                        
                        ip=ncord(L,i,j,0,0);
                        Cs[ip] = 1; 
                        R[ip]=1;
                        A[ip]=A0;
                        G[ip]=G_init;
                        Gi[ip]=Gi0*exp(gsl_ran_gaussian(prng,std_ic));
                        Ek[ip]=K_media*exp(gsl_ran_gaussian(prng,std_ic));
                        S[ip]=0; //ow(10,-8)
                        V[ip]=-156.0;
                        Ik[ip]=300;
                        ThT[ip]=0;
                        Ht[ip]=ht0;
                        H[ip]=h0;
                        Rg[ip]=1;                     
                        if (i <= i_center){
                            distances[ip]=i-left_e;
                        }else{
                            distances[ip]=right_e-i;
                        }
                        exp_gammafl_dist[ip]=(exp(-gamma_fl*distances[ip])+a_fl)/(1+a_fl);
                        exp_gammadif_dist[ip]=(exp(-gamma_dif*distances[ip])+a_dif)/(1+a_dif);
                    }
                    else if ((i>=left_e_2)&&(i<=right_e_2) && (Rn2>0)) {

                        ip=ncord(L,i,j,0,0);
                        Cs[ip] = 1; 
                        R[ip]=1;
                        A[ip]=A0;
                        G[ip]=G_init;
                        Gi[ip]=Gi0*exp(gsl_ran_gaussian(prng,std_ic));
                        Ek[ip]=K_media*exp(gsl_ran_gaussian(prng,std_ic));
                        S[ip]=pow(10,-8);
                        V[ip]=-156.0;
                        Ik[ip]=300;
                        ThT[ip]=0;
                        Ht[ip]=ht0;
                        H[ip]=h0;
                        Rg[ip]=1;

                        if (i <= i_center2){
                            distances[ip]=i-left_e_2;
                        }else{
                            distances[ip]=right_e_2-i;
                        }
                        exp_gammafl_dist[ip]=(exp(-gamma_fl*distances[ip])+a_fl)/(1+a_fl);
                        exp_gammadif_dist[ip]=(exp(-gamma_dif*distances[ip])+a_dif)/(1+a_dif);
                    }
                    else{
                        ip=ncord(L,i,j,0,0);
                        G[ip]=G_init;
                        Gi[ip]=0;
                        V[ip]=-156.0;
                        Cs[ip]=0;
                        Ek[ip]=K_media*exp(gsl_ran_gaussian(prng,std_ic));
                        Ik[ip]=0;
                        A[ip]=A_media;  
                        ThT[ip]=0;   
                        //Dp_v[ip]=0;
                        Ht[ip]=0;
                        H[ip]=0;
                        Rg[ip]=0;
                        distances[ip]=0;  
                        exp_gammafl_dist[ip]=1;   
                        exp_gammadif_dist[ip]=1; 
                    }
                }
            }
        } else{

            for (i=0; i<L; i++){
                for (j=0; j<j_max; j++){

                    if ((i>=left_e) && (i<=right_e) && (Rn>0)){

                        
                        ip=ncord(L,i,j,0,0);
                        Cs[ip] = 1; 
                        R[ip]=pow(10,-1);
                        A[ip]=A0;
                        G[ip]=G_init;
                        Gi[ip]=Gi0*exp(gsl_ran_gaussian(prng,std_ic));
                        Ek[ip]=K_media*exp(gsl_ran_gaussian(prng,std_ic));
                        S[ip]=0; //pow(10,-8);
                        V[ip]=-156.0;
                        Ik[ip]=300;
                        ThT[ip]=0;
                        Ht[ip]=ht0;
                        H[ip]=h0;
                        Rg[ip]=1;

                        if (i <= i_center){
                            distances[ip]=i-left_e;
                        }else{
                            distances[ip]=right_e-i;
                        }
                        exp_gammafl_dist[ip]=(exp(-gamma_fl*distances[ip])+a_fl)/(1+a_fl);
                        exp_gammadif_dist[ip]=(exp(-gamma_dif*distances[ip])+a_dif)/(1+a_dif);


                    }
                    else if ((i>=left_e_2)&&(i<=right_e_2)&& (Rn2>0)) {

                        ip=ncord(L,i,j,0,0);
                        Cs[ip] = 1; 
                        R[ip]=pow(10,-1);
                        A[ip]=A0;
                        G[ip]=G_init;
                        Gi[ip]=Gi0*exp(gsl_ran_gaussian(prng,std_ic));
                        Ek[ip]=K_media*exp(gsl_ran_gaussian(prng,std_ic));
                        S[ip]=pow(10,-8);
                        V[ip]=-156.0;
                        Ik[ip]=300;
                        ThT[ip]=0;
                        Ht[ip]=ht0;
                        H[ip]=h0;
                        Rg[ip]=1;

                        if (i <= i_center2){
                            distances[ip]=i-left_e_2;
                        }else{
                            distances[ip]=right_e_2-i;
                        }
                        exp_gammafl_dist[ip]=(exp(-gamma_fl*distances[ip])+a_fl)/(1+a_fl);
                        exp_gammadif_dist[ip]=(exp(-gamma_dif*distances[ip])+a_dif)/(1+a_dif);
                    }
                    else{ 
                        ip=ncord(L,i,j,0,0);
                        G[ip]=G_init;
                        Gi[ip]=0;
                        V[ip]=-156.0;
                        Cs[ip]=0;
                        Ek[ip]=K_media*exp(gsl_ran_gaussian(prng,std_ic));
                        Ik[ip]=0;
                        A[ip]=A_media;
                        ThT[ip]=0;   
                        //Dp_v[ip]=0;
                        Ht[ip]=0;
                        H[ip]=0;
                        Rg[ip]=0;
                        distances[ip]=0;  
                        exp_gammafl_dist[ip]=1;   
                        exp_gammadif_dist[ip]=1; 
                    }
                }
            }
        }
    }

    if (tprint>0){
        char basepath_in[200];
        sprintf(basepath_in,"%s/%s", argv[5], argv[6]);
        printf("opening ic to write");
        fflush(stdout);

        
        sprintf(filepath,"%s_c_pr.inp",basepath_in);
        outpat_c_pr = fopen(filepath,"w");
      
        sprintf(filepath,"%s_g_pr.inp",basepath_in);
        outpat_g_pr = fopen(filepath,"w");

        sprintf(filepath,"%s_gi_pr.inp",basepath_in);
        outpat_gi_pr = fopen(filepath,"w");

       
        sprintf(filepath,"%s_a_pr.inp",basepath_in);
        outpat_a_pr = fopen(filepath,"w");

        sprintf(filepath,"%s_ht_pr.inp",basepath_in);
        outpat_ht_pr = fopen(filepath,"w");

        sprintf(filepath,"%s_h_pr.inp",basepath_in);
        outpat_h_pr = fopen(filepath,"w");


        sprintf(filepath,"%s_r_pr.inp",basepath_in);
        outpat_r_pr = fopen(filepath,"w");
      
        sprintf(filepath,"%s_S_pr.inp",basepath_in);
        outpat_S_pr= fopen(filepath,"w");
        
        sprintf(filepath,"%s_v_pr.inp",basepath_in);
        outpat_v_pr = fopen(filepath,"w");
     
        sprintf(filepath,"%s_nk_pr.inp",basepath_in);
        outpat_nk_pr= fopen(filepath,"w");
       
        sprintf(filepath,"%s_Ek_pr.inp",basepath_in);
        outpat_Ek_pr= fopen(filepath,"w");
        
        sprintf(filepath,"%s_Ik_pr.inp",basepath_in);
        outpat_Ik_pr= fopen(filepath,"w");
        
        printf("i_center ==%d, Rn=%d",i_center,Rn);
        fflush(stdout);

        sprintf(filepath,"%s_tht_pr.inp",basepath_in);
        outpat_tht_pr= fopen(filepath,"w");

        sprintf(filepath,"%s_Dp_pr.inp",basepath_in);
        outpat_Dp_pr= fopen(filepath,"w");

    }
   
    G_media=GE;
    
    
    printf("initial glutamate, %g \n",sum_ar(G,N));
    fflush(stdout); 
    
    printf("initial density, %g \n",sum_ar(Cs,N));
    fflush(stdout); 
    printf("N, %d \n",N);
    fflush(stdout);
    
   
    int it_freeze=(int)(freeze/dt);
    //int it_perturb1=(int)(perturb1/dt);
    //int it_perturb2=(int)(perturb2/dt);
    int it_unfreeze=(int)(unfreeze/dt);
    int it_freeze2=(int)(freeze2/dt);
    /*
    int npert=(int)(Tpert/dt+0.5);
    int nwarm=(int)(Twarm/dt+0.5);
    int nstopKpert=(int)(stopKpert/dt+0.5);
    int dist_pK_dx=(int)(dist_pK/dx+0.5);
    
    int delay_onset_it=(int)(delay_onset/dt);
    fflush(stdout);
    int duration_pK_it=(int)(duration_pK/dt+0.5);

    int npertA=(int)(TpertA/dt+0.5);
    int nwarmA=(int)(TwarmA/dt+0.5);
    int nstopApert=(int)(stopApert/dt+0.5);
    int dist_pA_dx=(int)(dist_pA/dx+0.5);
    int duration_pA_it=(int)(duration_pA/dt+0.5);
    
    */
    int it_stopflow1=(int)(stopflow1/dt);
    int it_stopflow2=(int)(stopflow2/dt);
       
    int it_print=(int)(tprint/dt+0.5);

    if (it_print%nst!=0){
        it_print=it_print-(it_print%nst);
    
    }
    if (it_print%nst!=0){
        printf("it print non effective!%d. Exiting. ",it_print);
        fflush(stdout);
        exit(1);
    }
    
    
    double GE_p=GE*GE_perturb;
    double K_media0=K_media;

    int it_perturb_G=(int)(tperturbG/dt);
    int it_perturb_K1=(int)(tperturbK/dt);
    int it_perturb_K2=(int)((tperturbK+dperturbK)/dt);

    /********SIMULATION********************************/

    int ipn;
    int nan_found=0;
    int it_p;
    int ip_togrow;
    int other_ip;
    int ipn1,ipn2;
    

    double W_;
    double F_gk=F*gk;
    double Gentry;
    double Aprod;

    double H_act, H_deact;
   
    for (it=0;it<=nme;it++){

        //printf("%d \n",it);
        //fflush(stdout); 
        if (it > it_perturb_G){

            G_media=GE_p;
        }

        if ((it > it_perturb_K1) && (it < it_perturb_K2)){

            K_media=K_perturb;
        }else{
            K_media=K_media0;
        }


        if ((it==0) || (it%nst==0)){
            
            t=dt*(double)it;
            printf("\nbefore integrating, %g/%g", t, tt);
            fflush(stdout); 

            nan_found=check_nan(A,N);
            if (nan_found>0){
            printf("\n nan found in A");
            fflush(stdout); 
            }


            nan_found=check_nan(G,N);
            if (nan_found>0){
            printf("\n nan found in G");
            fflush(stdout); 
            }

            nan_found=check_nan(V,N);
            if (nan_found>0){
            printf("\n nan found in V");
            fflush(stdout); 
            exit(5);
            }
            fwrite(&t,sizeof(double),1,outpat_g);
            fwrite(G,sizeof(double),N,outpat_g);

            
            fwrite(&t,sizeof(double),1,outpat_gi);
            fwrite(Gi,sizeof(double),N,outpat_gi);

            fwrite(&t,sizeof(double),1,outpat_ginip);
            fwrite(Ginip_print,sizeof(double),N,outpat_ginip);

            fwrite(&t,sizeof(double),1,outpat_c);
            fwrite(Cs,sizeof(double),N,outpat_c);

            
            fwrite(&t,sizeof(double),1,outpat_r);
            fwrite(R,sizeof(double),N,outpat_r);

            fwrite(&t,sizeof(double),1,outpat_Rg);
            fwrite(Rg,sizeof(double),N,outpat_Rg);
            
            fwrite(&t,sizeof(double),1,outpat_a);
            fwrite(A,sizeof(double),N,outpat_a);

            fwrite(&t,sizeof(double),1,outpat_ht);
            fwrite(Ht,sizeof(double),N,outpat_ht);

            fwrite(&t,sizeof(double),1,outpat_h);
            fwrite(H,sizeof(double),N,outpat_h);            

            fwrite(&t,sizeof(double),1,outpat_S);
            fwrite(S,sizeof(double),N,outpat_S);

            fwrite(&t,sizeof(double),1,outpat_v);
            fwrite(V,sizeof(double),N,outpat_v);

            fwrite(&t,sizeof(double),1,outpat_vls);
            fwrite(Vls_pr,sizeof(double),N,outpat_vls);

            fwrite(&t,sizeof(double),1,outpat_vks);
            fwrite(Vks_pr,sizeof(double),N,outpat_vks);
            
            fwrite(&t,sizeof(double),1,outpat_nk);
            fwrite(nk,sizeof(double),N,outpat_nk);
            
            fwrite(&t,sizeof(double),1,outpat_Ek);
            fwrite(Ek,sizeof(double),N,outpat_Ek);
            
            fwrite(&t,sizeof(double),1,outpat_Ik);
            fwrite(Ik,sizeof(double),N,outpat_Ik);

            fwrite(&t,sizeof(double),1,outpat_ikmax);
            fwrite(Ikmax_pr,sizeof(double),N,outpat_ikmax);

            fwrite(&t,sizeof(double),1,outpat_tht);
            fwrite(ThT,sizeof(double),N,outpat_tht);

            fwrite(&t,sizeof(double),1,outpat_dist);
            fwrite(distances,sizeof(double),N,outpat_dist);

            fwrite(&t,sizeof(double),1,outpat_fl_dist);
            fwrite(exp_gammafl_dist,sizeof(double),N,outpat_fl_dist);

            fwrite(&t,sizeof(double),1,outpat_dif_dist);
            fwrite(exp_gammadif_dist,sizeof(double),N,outpat_dif_dist);
           
        
            if (it==it_print){
                //print_ic();
                printf("printing t %d", it);
                fflush(stdout);
                fwrite(G,sizeof(double),N,outpat_g_pr);

                fwrite(Gi,sizeof(double),N,outpat_gi_pr);
                
                fwrite(Cs,sizeof(double),N,outpat_c_pr);

                fwrite(R,sizeof(double),N,outpat_r_pr);

                fwrite(A,sizeof(double),N,outpat_a_pr);

                fwrite(Ht,sizeof(double),N,outpat_ht_pr);

                fwrite(H,sizeof(double),N,outpat_h_pr);
                
                fwrite(S,sizeof(double),N,outpat_S_pr);
                            
                fwrite(V,sizeof(double),N,outpat_v_pr);
                            
                fwrite(nk,sizeof(double),N,outpat_nk_pr);
                            
                fwrite(Ek,sizeof(double),N,outpat_Ek_pr);

                fwrite(Ik,sizeof(double),N,outpat_Ik_pr);

                fwrite(ThT,sizeof(double),N,outpat_tht_pr);

            }
        }

        if (nan_found ==0){

            if ((it>it_stopflow1)&(it<it_stopflow2)){
                fl_Ek=0;
                fl_G=0;
                fl_A=0;
                
            }else{
                    
                fl_Ek=fl_Ekext;
                fl_G=fl_Gext;
                fl_A=fl_Aext;

            }
        
            //calculate *k1

            

            for (ip=0; ip<N; ip++){
                exp_gammadifG_dist[ip]=D_gdx*exp_gammadif_dist[ip];
                exp_gammadifEk_dist[ip]=D_kedx*exp_gammadif_dist[ip];

                exp_gammaflG_dist[ip]=fl_G*exp_gammafl_dist[ip];
                exp_gammaflEk_dist[ip]=fl_Ek*exp_gammafl_dist[ip];
            }


            for (ip=0;ip<N;ip++){
                D_G1_k[ip]=DG1_f(exp_gammadifG_dist[ip],G,ip, n1, n3);
                D_Ek1_k[ip]=DK1_f(exp_gammadifEk_dist[ip],Ek,ip, n1, n3);               
                D_A1_k[ip]=DA1_f(D_adx,A,ip, n1, n3);
              
                if (Cs[ip]>zero_tol){                                     
                    //if ((it>delay_onset_it)| ((it<delay_onset_it)&& (ip<center_chamber))){
                       

                        Ginip=Gin(V[ip], V0, gv)*Rg[ip];

                        R_Rg_k[ip]=R_rg_f(Rg, R, ip, alpha_Rg, delta_Rg, alpha_rRg, Krgnrg, nrg);
                        R_S_k[ip]=R_S_f(Gi,S,ip,S0, ug, Gs0ug, gamma_s, V);                        
                        R_G_k[ip]=R_G_f(G, Ginip,ip, alpha_gt, Kgt);        
                        Aprod=A_production(Gi,H[ip], alpha_a, ip);
                        R_Gi_k[ip]=-R_G_k[ip]-Aprod+R_Gi_f_noentry(Gi, R,ip, delta_g,V);
                        R_G_k[ip]+=flow_G(G,ip,exp_gammaflG_dist[ip],G_media);  //I am adding this here because I use the other value above
                        R_A_k[ip]=Aprod +R_A_deg(A, R,ip, delta_a);

                        
                        H_act=R_H_f_act(Gi, Ht, H, ip, alpha_h, gamma_h, K_hgn, n_hg);
                        H_deact=R_H_f_deact(H, gamma_h, ip);
                        R_H_k[ip]=H_act-H_deact;
                        R_Ht_k[ip]=R_Ht_f(H, Gi, Ht, ip, alpha_htg, alpha_hta, gamma_ht, K_htA, K_htG, nht)-H_act+H_deact;                           
                        R_R_k[ip]=R_R_f(A, Gi,R, ip, beta_r,gamma_r);
 
                        aS=a0*pow(S[ip],mk)/(Sthm+pow(S[ip],mk));
                        Vls=Vl+dl*(K_media-Ek[ip])/(exp((K_media-Ek[ip])/0.1)-1.0);
                        Vks=Vk*log(Ek[ip]/Ik[ip]);
               
                        R_V_k[ip]=R_V_f(nk,V,Vks,Vls,ip, gk, gl);
                        R_nk_k[ip]=R_nk_f(aS,nk,ip, b);
                        W_=R_Ek_f(nk,V,Ek,Ik,Ikmax,ip, Dp,F_gk, Vks,Gi);
                        R_Ek_k[ip]=W_+flow_Ek(Ek,ip, exp_gammaflEk_dist[ip], K_media);
                        R_Ik_k[ip]=-W_;

                        R_ThT_k[ip]=R_T_f(V, ThT, alpha_t, V0, gamma_t,  V0_tht, ip, g_tht);

                        


                        S_k1[ip]=dth*R_S_k[ip];
                                             
                        Gi_k1[ip] =dth*(R_Gi_k[ip]);
                        R_k1[ip] =dth*(R_R_k[ip]);
                        Rg_k1[ip]=dth*R_Rg_k[ip];
                        Dp_k1[ip]=dth*(R_Dp_k[ip]);
                        Ht_k1[ip]=dth*(R_Ht_k[ip]);
                        H_k1[ip]=dth*(R_H_k[ip]);
                        
                        V_k1[ip]=dth*R_V_k[ip];
                        nk_k1[ip]=dth*R_nk_k[ip];                        
                        Ik_k1[ip]=dth*(R_Ik_k[ip]);
                        ThT_k1[ip]=dth*R_ThT_k[ip];

                        Cs_next[ip]=Cs[ip];                 
                        S_next[ip]=S[ip]+S_k1[ip];                       
                        Gi_next[ip] =Gi[ip]+Gi_k1[ip];
                        //A_next[ip] =A[ip]+A_k1[ip];
                        R_next[ip] =R[ip]+R_k1[ip];
                        Rg_next[ip]=Rg[ip]+Rg_k1[ip];
                        //Dp_v_next[ip]=Dp_v[ip]+Dp_k1[ip];
                        Ht_next[ip]=Ht[ip]+Ht_k1[ip];
                        H_next[ip]=H[ip]+H_k1[ip];
                        
                        
                        V_next[ip]=V[ip]+V_k1[ip];
                        nk_next[ip]=nk[ip]+nk_k1[ip];                        
                        Ik_next[ip]=Ik[ip]+Ik_k1[ip];
                        ThT_next[ip]=ThT[ip]+ThT_k1[ip];
                        
                   /* }
                    else{
                    
                        R_S_k[ip]=0;
                        R_G_k[ip]=flow_G(G,ip, fl_G, G_media);
                        R_Gi_k[ip]=0;
                        R_Ik_k[ip]=0;
                        R_V_k[ip]=0;
                        R_Ek_k[ip]=flow_Ek(Ek,ip);
                        R_nk_k[ip]=0;
                        R_R_k[ip]=0;
                        R_A_k[ip]=flow_A(A_unique,ip);
                        //R_Ht_k1[ip]=0;
                        //R_H_k1[ip]=0;
                        
                    }*/
                    
                }else{
                   
                    //R_S_k[ip]=0;
                    R_G_k[ip]=flow_G(G,ip, fl_G, G_media);
                    //R_Gi_k[ip]=0;
                    //R_Ik_k[ip]=0;
                    //R_V_k[ip]=0;
                    R_Ek_k[ip]=flow_Ek(Ek,ip, fl_Ek, K_media);
                    //R_nk_k[ip]=0;
                    //R_R_k[ip]=0;
                    R_A_k[ip]=flow_A(A,ip, fl_A, A_media);

                    //R_ThT_k[ip]=0;
                    //R_Dp_k[ip]=0;

                    //R_Ht_k[ip]=0;
                    //R_H_k[ip]=0;

                    
                }
                A_k1[ip] =dth*(R_A_k[ip]+D_A1_k[ip]); 
                G_k1[ip] =dth*(R_G_k[ip]+D_G1_k[ip]);
                Ek_k1[ip]=dth*(R_Ek_k[ip]+D_Ek1_k[ip]);

                G_next[ip] =G[ip]+G_k1[ip];
                Ek_next[ip]=Ek[ip]+Ek_k1[ip];
                A_next[ip]=A[ip]+A_k1[ip];
            }
            
            //for k2

            

            for (ip=0;ip<N;ip++){
                D_G1_k[ip]=DG1_f(exp_gammadifG_dist[ip],G_next,ip, n1, n3);
                D_Ek1_k[ip]=DK1_f(exp_gammadifEk_dist[ip],Ek_next,ip, n1, n3);               
                D_A1_k[ip]=DA1_f(D_adx,A_next,ip,n1, n3);
                if (Cs_next[ip]>zero_tol){
                    //if ((it>delay_onset_it)| ((it<delay_onset_it)&& (ip<center_chamber))){
                        

                        Ginip=Gin(V_next[ip], V0, gv)*Rg_next[ip];
                        Ginip_print[ip]=Ginip;
                        
                        R_S_k[ip]=R_S_f(Gi_next,S_next,ip,S0, ug, Gs0ug, gamma_s, V_next);
                        R_R_k[ip]=R_R_f(A_next, Gi_next,R_next, ip, beta_r, gamma_r);
                        R_Rg_k[ip]=R_rg_f(Rg_next, R_next, ip, alpha_Rg, delta_Rg, alpha_rRg, Krgnrg,  nrg);
                        R_G_k[ip]=R_G_f(G_next, Ginip,ip, alpha_gt, Kgt);     
                        Aprod=A_production(Gi_next, H_next[ip],alpha_a, ip);
                        R_Gi_k[ip]=-R_G_k[ip]-Aprod+R_Gi_f_noentry(Gi_next, R_next,ip, delta_g, V_next);
                        R_G_k[ip]+=flow_G(G_next,ip,exp_gammaflG_dist[ip],G_media);    
                        R_A_k[ip]=Aprod +R_A_deg(A_next, R_next,ip, delta_a);
                        
                        H_act=R_H_f_act(Gi_next, Ht_next, H_next, ip, alpha_h, gamma_h, K_hgn, n_hg);
                        H_deact=R_H_f_deact(H_next, gamma_h,ip);
                        R_H_k[ip]=H_act-H_deact;
                        R_Ht_k[ip]=R_Ht_f(H_next, Gi_next, Ht_next, ip, alpha_htg, alpha_hta, gamma_ht, K_htA, K_htG, nht)-H_act+H_deact;
                        
                        aS=a0*pow(S_next[ip],mk)/(Sthm+pow(S_next[ip],mk));
                        Vls=Vl+dl*(K_media-Ek_next[ip])/(exp((K_media-Ek_next[ip])/0.1)-1.0);
                        Vks=Vk*log(Ek_next[ip]/Ik_next[ip]);
                        
                        R_V_k[ip]=R_V_f(nk_next,V_next,Vks,Vls,ip, gk, gl);
                        R_nk_k[ip]=R_nk_f(aS,nk_next,ip, b);
                        W_=R_Ek_f(nk_next,V_next,Ek_next,Ik_next,Ikmax,ip, Dp,F_gk, Vks,Gi_next);
                        R_Ek_k[ip]=W_+flow_Ek(Ek_next,ip, exp_gammaflEk_dist[ip], K_media);
                        R_Ik_k[ip]=-W_;

                        R_ThT_k[ip]=R_T_f(V_next, ThT_next, alpha_t, V0, gamma_t, V0_tht, ip,g_tht);
                       

                        
                        S_k2[ip]=dt*R_S_k[ip];                                             
                        Gi_k2[ip] =dt*(R_Gi_k[ip]);
                        R_k2[ip] =dt*(R_R_k[ip]);
                        Rg_k2[ip] =dt*(R_Rg_k[ip]);
                        Dp_k2[ip]=dt*R_Dp_k[ip];
                        Ht_k2[ip]=dt*R_Ht_k[ip];
                        H_k2[ip]=dt*R_H_k[ip];
                        
                        V_k2[ip]=dt*R_V_k[ip];
                        nk_k2[ip]=dt*R_nk_k[ip];                       
                        Ik_k2[ip]=dt*(R_Ik_k[ip]);
                        ThT_k2[ip]=dt*(R_ThT_k[ip]);

                        Cs_next[ip]=Cs[ip];                
                        S_next[ip]=S[ip]+S_k2[ip]/2;                       
                        Gi_next[ip] =Gi[ip]+Gi_k2[ip]/2;
                        R_next[ip] =R[ip]+R_k2[ip]/2;
                        Rg_next[ip]=Rg[ip]+Rg_k2[ip]/2;
                        Ht_next[ip]=Ht[ip]+Ht_k2[ip]/2;
                        H_next[ip]=H[ip]+H_k2[ip]/2;
                                                
                        V_next[ip]=V[ip]+V_k2[ip]/2;
                        nk_next[ip]=nk[ip]+nk_k2[ip]/2;                       
                        Ik_next[ip]=Ik[ip]+Ik_k2[ip]/2;
                        ThT_next[ip]=ThT[ip]+ThT_k2[ip]/2;
    



                   
                }else{
                
 
                    R_G_k[ip]=flow_G(G_next,ip, fl_G, G_media);
                    R_A_k[ip]=flow_A(A_next,ip, fl_A, A_media);
                    R_Ek_k[ip]=flow_Ek(Ek_next,ip, fl_Ek, K_media);
                    

                }
                A_k2[ip]=dt*(R_A_k[ip]+D_A1_k[ip]);
                G_k2[ip] =dt*(R_G_k[ip]+D_G1_k[ip]);
                Ek_k2[ip]=dt*(R_Ek_k[ip]+D_Ek1_k[ip]);

                G_next[ip] =G[ip]+G_k2[ip]/2;
                Ek_next[ip]=Ek[ip]+Ek_k2[ip]/2;
                A_next[ip]=A[ip]+A_k2[ip]/2;
            }
        
           
            
            

            //k3

            for (ip=0;ip<N;ip++){
                D_G1_k[ip]=DG1_f(exp_gammadifG_dist[ip],G_next,ip, n1, n3);
                D_Ek1_k[ip]=DK1_f(exp_gammadifEk_dist[ip],Ek_next,ip, n1, n3); 
                D_A1_k[ip]=DA1_f(D_adx,A_next,ip,n1, n3);
                if (Cs_next[ip]>zero_tol){
                    //if ((it>delay_onset_it)| ((it<delay_onset_it)&& (ip<center_chamber))){
                        

                        Ginip=Gin(V_next[ip], V0, gv)*Rg_next[ip];
                        //Ginip_print[ip]=Ginip;
                        
                        R_S_k[ip]=R_S_f(Gi_next,S_next,ip,S0, ug, Gs0ug, gamma_s, V_next);                        
                        R_R_k[ip]=R_R_f(A_next, Gi_next,R_next, ip, beta_r, gamma_r);
                        R_Rg_k[ip]=R_rg_f(Rg_next, R_next, ip, alpha_Rg, delta_Rg, alpha_rRg, Krgnrg,  nrg);
                        
                        R_G_k[ip]=R_G_f(G_next, Ginip,ip, alpha_gt, Kgt);    
                        Aprod=A_production(Gi_next,H_next[ip], alpha_a, ip);
                        R_Gi_k[ip]=-R_G_k[ip]-Aprod+R_Gi_f_noentry(Gi_next, R_next,ip, delta_g, V_next);
                        R_G_k[ip]+=flow_G(G_next,ip,exp_gammaflG_dist[ip],G_media);  
                        R_A_k[ip]=Aprod +R_A_deg(A_next, R_next,ip, delta_a);
                        
                        H_act=R_H_f_act(Gi_next, Ht_next, H_next, ip, alpha_h, gamma_h, K_hgn, n_hg);
                        H_deact=R_H_f_deact(H_next, gamma_h,ip);
                        R_H_k[ip]=H_act-H_deact;
                        R_Ht_k[ip]=R_Ht_f(H_next, Gi_next, Ht_next, ip, alpha_htg, alpha_hta, gamma_ht, K_htA, K_htG, nht)-H_act+H_deact;
     
                        aS=a0*pow(S_next[ip],mk)/(Sthm+pow(S_next[ip],mk));
                        Vls=Vl+dl*(K_media-Ek_next[ip])/(exp((K_media-Ek_next[ip])/0.1)-1.0);
                        Vks=Vk*log(Ek_next[ip]/Ik_next[ip]);

                        Vls_pr[ip]=Vls;
                        Vks_pr[ip]=Vks;
                        
                        
                        R_V_k[ip]=R_V_f(nk_next,V_next,Vks,Vls,ip, gk, gl);
                        R_nk_k[ip]=R_nk_f(aS,nk_next,ip, b);
                        W_=R_Ek_f(nk_next,V_next,Ek_next,Ik_next,Ikmax,ip, Dp,F_gk, Vks,Gi_next);
                        R_Ek_k[ip]=W_+flow_Ek(Ek_next,ip, exp_gammaflEk_dist[ip], K_media);
                        R_Ik_k[ip]=-W_;

                        R_ThT_k[ip]=R_T_f(V_next, ThT_next, alpha_t, V0, gamma_t, V0_tht, ip,g_tht);

                        

                        S_k3[ip]=dt*R_S_k[ip];
                                                
                        Gi_k3[ip] =dt*(R_Gi_k[ip]);
                        R_k3[ip] =dt*(R_R_k[ip]);
                        Rg_k3[ip]=dt*R_Rg_k[ip];
                        Dp_k3[ip]=dt*R_Dp_k[ip];
                        Ht_k3[ip]=dt*R_Ht_k[ip];
                        H_k3[ip]=dt*R_H_k[ip];
                        
                        
                        V_k3[ip]=dt*R_V_k[ip];
                        nk_k3[ip]=dt*R_nk_k[ip];                       
                        Ik_k3[ip]=dt*(R_Ik_k[ip]);
                        ThT_k3[ip]=dt*R_ThT_k[ip];

                        Cs_next[ip]=Cs[ip];                
                        S_next[ip]=S[ip]+S_k3[ip];                        
                        Gi_next[ip] =Gi[ip]+Gi_k3[ip];
                        
                        R_next[ip] =R[ip]+R_k3[ip];
                        Rg_next[ip]=Rg[ip]+Rg_k3[ip];
                        Ht_next[ip]=Ht[ip]+Ht_k3[ip];
                        H_next[ip]=H[ip]+H_k3[ip];
                        
                        
                        V_next[ip]=V[ip]+V_k3[ip];
                        nk_next[ip]=nk[ip]+nk_k3[ip];                        
                        Ik_next[ip]=Ik[ip]+Ik_k3[ip];
                        ThT_next[ip]=ThT[ip]+ThT_k3[ip];



                   
                }else{

                    R_G_k[ip]=flow_G(G_next,ip, fl_G, G_media);
                    R_A_k[ip]=flow_A(A_next,ip, fl_A, A_media);
                    R_Ek_k[ip]=flow_Ek(Ek_next,ip, fl_Ek, K_media);

                }
                A_k3[ip] =dt*(R_A_k[ip]+D_A1_k[ip]);
                G_k3[ip] =dt*(R_G_k[ip]+D_G1_k[ip]);
                Ek_k3[ip]=dt*(R_Ek_k[ip]+D_Ek1_k[ip]);

                G_next[ip] =G[ip]+G_k3[ip];
                Ek_next[ip]=Ek[ip]+Ek_k3[ip];
                A_next[ip] =A[ip]+A_k3[ip];
            } //k3 
        
            //for k4        
            
            for (ip=0;ip<N;ip++){
                D_G1_k[ip]=DG1_f(exp_gammadifG_dist[ip],G_next,ip, n1, n3);
                D_Ek1_k[ip]=DK1_f(exp_gammadifEk_dist[ip],Ek_next,ip, n1, n3);              
                D_A1_k[ip]=DA1_f(D_adx,A_next,ip,n1, n3);
                if (Cs_next[ip]>zero_tol){
                    //if ((it>delay_onset_it)| ((it<delay_onset_it)&& (ip<center_chamber))){
                       

                        Ginip=Gin(V_next[ip], V0, gv)*Rg_next[ip];
                         
                        R_S_k[ip]=R_S_f(Gi_next,S_next,ip,S0, ug, Gs0ug, gamma_s, V_next);
                        R_R_k[ip]=R_R_f(A_next, Gi_next,R_next, ip, beta_r, gamma_r);
                        R_Rg_k[ip]=R_rg_f(Rg_next, R_next, ip, alpha_Rg, delta_Rg, alpha_rRg, Krgnrg,  nrg);
                        R_G_k[ip]=R_G_f(G_next, Ginip,ip, alpha_gt, Kgt);
                        Aprod=A_production(Gi_next,H_next[ip], alpha_a, ip);
                        R_Gi_k[ip]=-R_G_k[ip]-Aprod+R_Gi_f_noentry(Gi_next, R_next,ip, delta_g, V_next);
                        R_G_k[ip]+=flow_G(G_next,ip,exp_gammaflG_dist[ip],G_media);         
                        R_A_k[ip]=Aprod +R_A_deg(A_next, R_next,ip, delta_a);
                       
                        H_act=R_H_f_act(Gi_next, Ht_next, H_next, ip, alpha_h, gamma_h, K_hgn, n_hg);
                        H_deact=R_H_f_deact(H_next, gamma_h,ip);
                        R_H_k[ip]=H_act-H_deact;
                        R_Ht_k[ip]=R_Ht_f(H_next, Gi_next, Ht_next, ip, alpha_htg, alpha_hta, gamma_ht, K_htA, K_htG, nht)-H_act+H_deact;
    
                        
                        aS=a0*pow(S_next[ip],mk)/(Sthm+pow(S_next[ip],mk));
                        Vls=Vl+dl*(K_media-Ek_next[ip])/(exp((K_media-Ek_next[ip])/0.1)-1.0);
                        Vks=Vk*log(Ek_next[ip]/Ik_next[ip]);
                        
                        R_V_k[ip]=R_V_f(nk_next,V_next,Vks,Vls,ip, gk, gl);
                        R_nk_k[ip]=R_nk_f(aS,nk_next,ip, b);
                        W_=R_Ek_f(nk_next,V_next,Ek_next,Ik_next,Ikmax,ip, Dp,F_gk, Vks,Gi_next);
                        R_Ek_k[ip]=W_+flow_Ek(Ek_next,ip, exp_gammaflEk_dist[ip], K_media);
                        R_Ik_k[ip]=-W_;
                        R_ThT_k[ip]=R_T_f(V_next, ThT_next, alpha_t, V0, gamma_t, V0_tht, ip,g_tht);

                        
                        S_k4[ip]=dt6*R_S_k[ip];
                        
                                                
                        Gi_k4[ip] =dt6*(R_Gi_k[ip]);
                        R_k4[ip] =dt6*(R_R_k[ip]);
                        Rg_k4[ip]=dt6*R_Rg_k[ip];
                        Ht_k4[ip]=dt6*R_Ht_k[ip];
                        H_k4[ip]=dt6*R_H_k[ip];                        
                        V_k4[ip]=dt6*R_V_k[ip];
                        nk_k4[ip]=dt6*R_nk_k[ip];                        
                        Ik_k4[ip]=dt6*(R_Ik_k[ip]);
                        ThT_k4[ip]=dt6*R_ThT_k[ip];


                   
                }else{

                    R_G_k[ip]=flow_G(G_next,ip, fl_G, G_media);
                    R_A_k[ip]=flow_A(A_next,ip, fl_A, A_media);
                    R_Ek_k[ip]=flow_Ek(Ek_next,ip, fl_Ek, K_media);
                }

                A_k4[ip] =dt6*(R_A_k[ip]+D_A1_k[ip]);
                G_k4[ip] =dt6*(R_G_k[ip]+D_G1_k[ip]);
                Ek_k4[ip]=dt6*(R_Ek_k[ip]+D_Ek1_k[ip]);
            }


            
            for (ip=0;ip<N;ip++){

               
                S[ip]=S[ip]+S_k1[ip]/3+S_k2[ip]/3+S_k3[ip]/3+S_k4[ip];
                G[ip]=G[ip]+G_k1[ip]/3+G_k2[ip]/3+G_k3[ip]/3+G_k4[ip];
                A[ip]=A[ip]+A_k1[ip]/3+A_k2[ip]/3+A_k3[ip]/3+A_k4[ip];
                R[ip]=R[ip]+R_k1[ip]/3+R_k2[ip]/3+R_k3[ip]/3+R_k4[ip];
                Rg[ip]=Rg[ip]+Rg_k1[ip]/3+Rg_k2[ip]/3+Rg_k3[ip]/3+Rg_k4[ip];
                Gi[ip]=Gi[ip]+Gi_k1[ip]/3+Gi_k2[ip]/3+Gi_k3[ip]/3+Gi_k4[ip];
                Ht[ip]=Ht[ip]+Ht_k1[ip]/3+Ht_k2[ip]/3+Ht_k3[ip]/3+Ht_k4[ip];
                H[ip]=H[ip]+H_k1[ip]/3+H_k2[ip]/3+H_k3[ip]/3+H_k4[ip];
                V[ip]=V[ip]+V_k1[ip]/3+V_k2[ip]/3+V_k3[ip]/3+V_k4[ip];
                nk[ip]=nk[ip]+nk_k1[ip]/3+nk_k2[ip]/3+nk_k3[ip]/3+nk_k4[ip];
                Ek[ip]=Ek[ip]+Ek_k1[ip]/3+Ek_k2[ip]/3+Ek_k3[ip]/3+Ek_k4[ip]; 
                Ik[ip]=Ik[ip]+Ik_k1[ip]/3+Ik_k2[ip]/3+Ik_k3[ip]/3+Ik_k4[ip]; 

                ThT[ip]=ThT[ip]+ThT_k1[ip]/3+ThT_k2[ip]/3+ThT_k3[ip]/3+ThT_k4[ip];
                //slowv[ip]=slowv[ip]+slowv_k1[ip]/3+slowv_k2[ip]/3+slowv_k3[ip]/3+slowv_k4[ip];               
                

            }




            if ((it<it_freeze) |  ((it> it_unfreeze)&&(it<it_freeze2))){
               //stochastic growth

                if (Rn>0){ 
                    //printf("checking growth R[right_e]= %g\n", R[right_e]);
                    //fflush(stdout);

            

                    if (gsl_rng_uniform(prng) < (Pgrow_norm * R[right_e])){
                        grow(right_e+1, right_e,L, Cs, Gi, S, nk, V, Ik, ThT,R, Rg, Ht, H, distances, exp_gammafl_dist, exp_gammadif_dist, gamma_fl, a_fl, gamma_dif, a_dif);
                        right_e=right_e+1;     

                        }
                }

                if (Rn2>0){ 

                    
                    

                    if (gsl_rng_uniform(prng) < (Pgrow_norm * R[left_e_2])){

                        grow(left_e_2-1, left_e_2, L, Cs, Gi, S, nk, V, Ik, ThT,R, Rg, Ht, H, distances, exp_gammafl_dist, exp_gammadif_dist,gamma_fl, a_fl, gamma_dif, a_dif);
                                                                    
                        left_e_2=left_e_2-1;
                        
                    }
                }
             
            }//if it<
        } //if nan_found   
    }//for it
   
    fclose(outpat_c);
    fclose(outpat_g);
    fclose(outpat_gi);
    fclose(outpat_ginip);
    fclose(outpat_S);
    fclose(outpat_a);
    fclose(outpat_r);
    fclose(outpat_Rg);
    fclose(outpat_Dp);
    fclose(outpat_v);
    fclose(outpat_vls); 
    fclose(outpat_vks);
    fclose(outpat_nk);
    fclose(outpat_Ek);
    fclose(outpat_Ik);
    fclose(outpat_tht);
    fclose(outpat_ht);
    fclose(outpat_h);
    fclose(outpat_slowv);
    //fclose(outpat_Sv);
    printf("Quitting ");
    fflush(stdout); 

}// main
