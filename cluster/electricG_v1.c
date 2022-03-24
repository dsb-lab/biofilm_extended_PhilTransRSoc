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


void grow(int ipn1, int ip_togrow, int L, double *Cs, double *Gi, double *Rg, double *S, double *nk, double *V, double *Ik, double *ThT, double *distances, double *exp_gammafl_dist, double *exp_gammadif_dist, double gamma_fl, double a_fl, double gamma_dif, double a_dif ){
   
    //grow towards the right
    
    int i;

    Cs[ipn1]=1;
    Gi[ipn1]=Gi[ip_togrow];
    Ik[ipn1]=Ik[ip_togrow];
    S[ipn1]=S[ip_togrow];
    nk[ipn1]=nk[ip_togrow];
    V[ipn1]=V[ip_togrow];
    ThT[ipn1]=ThT[ip_togrow];
    
    Rg[ipn1]=Rg[ip_togrow];
    

    
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

    

int main(int argc, char *argv[]){
    
    int i, j, ijk, ip, it;
    double t;
    double *G,*G_next,*R_G_k, *G_k1, *G_k2, *G_k3, *G_k4, *D_G1_k; 
    double *Cs,*Cs_next;
    double *Gi,*Gi_next,*R_Gi_k, *Gi_k1, *Gi_k2, *Gi_k3, *Gi_k4; 
    double *Rg,*Rg_next,*R_Rg_k, *Rg_k1, *Rg_k2, *Rg_k3, *Rg_k4;  
    double *S,*S_next,*R_S_k, *S_k1, *S_k2, *S_k3, *S_k4; 
    double *slowv,*slowv_next,*R_slowv_k, *slowv_k1, *slowv_k2, *slowv_k3, *slowv_k4; 
    double *V,*V_next,*R_V_k,*V_k1,*V_k2, *V_k3, *V_k4, *Vls_pr, *Vks_pr;
    double *nk,*nk_next,*R_nk_k,*nk_k1,*nk_k2,*nk_k3, *nk_k4;
    double *Ek,*Ek_next,*R_Ek_k,*Ek_k1, *Ek_k2, *Ek_k3, *Ek_k4, *D_Ek1_k;
    double *ThT, *ThT_next, *R_ThT_k, *ThT_k1, *ThT_k2, *ThT_k3, *ThT_k4;
    double *Ik, *Ik_next, *R_Ik_k,*Ik_k1, *Ik_k2, *Ik_k3, *Ik_k4; 
    double Vls,Vks,aS, Gentry, gkr, glr; //these are auxilliarly variables to store values during the computation
    double *distances;
    double *exp_gammafl_dist;
    double *exp_gammadif_dist;


    double *exp_gammaflG_dist;
    double *exp_gammaflEk_dist;
    double *exp_gammadifG_dist;
    double *exp_gammadifEk_dist;

     int *n1,*n3;
    
    
    /*---- Simulation parameters ---- */

    //be very careful with the ordering!!! if a new parameter is added, it has to be modified here and in the list of variables and pointers below
    
    const char* parnames[] = {"Krg", "nrg", "alpha_Rg", "delta_Rg", "alpha_rRg", "a_fl", "a_dif", "gamma_dif", 
    "gamma_fl", "std_ic","Gi01", "Gi02", "Kgt", "dperturbK", "tperturbK", "K_perturb", "tperturbG", "GE_perturb", 
    "V0_tht", "g_tht", "dl", "alpha_t", "gamma_t",  "stopflow1", "stopflow2", "dx", "tt", "dt", "st", "zero_tol", "gv", 
    "Radius2", "Radius1", "freeze", "mk", "Ikmax", "gl", "gk", "fl_Gext", "V0", "Dg", "F", "Pgrow", "alpha_gt", "freeze2", 
    "tprint", "Dp", "Sth", "K_media", "a0", "GE", "center1", "delta_g", "center2", "seed_", "GS0", "b", "read_input", "S0", 
    "Vk", "Vl", "Dke", "unfreeze", "gamma_s", "ug", "fl_Ekext", "Ls", "dist_pK", "alpha_sv"};

    int npars = 69;

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

    double Krg, nrg, alpha_Rg, delta_Rg, alpha_rRg, a_fl, a_dif, gamma_dif, gamma_fl, std_ic, Gi01, Gi02, Kgt, dperturbK, tperturbK, K_perturb, tperturbG, GE_perturb, V0_tht, g_tht, dl, alpha_t, gamma_t, stopflow1, stopflow2, dx, tt, dt, st, zero_tol, gv, Radius2, Radius1, freeze, mk, Ikmax, gl, gk, fl_Gext, V0, Dg, F, Pgrow, alpha_gt, freeze2, tprint, Dp, Sth, K_media, a0, GE, center1, delta_g, center2, seed_, GS0, b, read_input, S0, Vk, Vl, Dke, unfreeze, gamma_s, ug, fl_Ekext, Ls, dist_pK, alpha_sv;
    

    //list of pointers to the variables. & operator denotes an address in memory
    // this has to be in the same order as in the parnames list above
    double* readpars[]={&Krg, &nrg, &alpha_Rg, &delta_Rg, &alpha_rRg, &a_fl, &a_dif, &gamma_dif, &gamma_fl, &std_ic, &Gi01, &Gi02, &Kgt, &dperturbK, &tperturbK, &K_perturb, &tperturbG, &GE_perturb, &V0_tht, &g_tht, &dl, &alpha_t, &gamma_t, &stopflow1, &stopflow2, &dx, &tt, &dt, &st, &zero_tol, &gv, &Radius2, &Radius1, &freeze, &mk, &Ikmax, &gl, &gk, &fl_Gext, &V0, &Dg, &F, &Pgrow, &alpha_gt, &freeze2, &tprint, &Dp, &Sth, &K_media, &a0, &GE, &center1,  &delta_g, &center2, &seed_, &GS0, &b, &read_input, &S0, &Vk, &Vl, &Dke, &unfreeze, &gamma_s, &ug, &fl_Ekext, &Ls, &dist_pK, &alpha_sv};
             
    printf("assigning values read to variables\n");
    fflush(stdout);
 
    for (i=0;i<npars;i++){
        *readpars[i]=parvals[i]; //assigns the value to the correspondent variable, using the pointer. * returns the value of the variable located at the address specified by its operand.
    }

    printf("Ls is %g", Ls);
    fflush(stdout);

    double Ginip;

    double fl_Ek=fl_Ekext;
    double fl_G=fl_Gext;

    double G_media=GE;

    double Sthm=pow(Sth, mk);
    double Gs0ug=pow(GS0, ug);
    
    double dth = dt*0.5;
    double dt6=dt/6.0;
    double dx2=1/(pow(dx,2));

    double D_gdx=Dg*dx2;
    double D_kedx=Dke*dx2;
    double Krgnrg=pow(Krg, nrg);
    
    
     int nme =  (int)(tt/dt+0.5);
     int nst = (int)(st/dt+0.5);
     int npat = (int)(tt/st+0.5);
    
   

    printf("total time (tt)=%g, dt=%g,dx=%g,dx2=%g \n",tt,dt,dx,dx2);
    fflush(stdout); 

    
    gsl_rng * prng = gsl_rng_alloc(gsl_rng_mt19937); //pointer to instance of random number generator of default type

    
     int i_center = (int) (center1/dx);
     int i_center2=(int) (center2/dx);
     int j_center = 0;
     int j_max=1;

     int Rn=(int)(Radius1/dx);
     int Rn2=(int)(Radius2/dx);

    //edges of both biofilms
     int left_e = i_center-Rn;//left edge of biofilm.
     int left_e_2= i_center2-Rn2;
     int right_e=i_center+Rn;
     int right_e_2=i_center2+Rn2;
    //int left, centre, right, ipn1, ip_togrow; 
    
     int L=(int)(Ls/dx);//number of cells
    printf("Ls is %g", Ls);
    printf("dx is %g", dx);
    printf("L=%g",L);
    fflush(stdout);

     int N=L*1; //2n
    int seed = (int) seed_;
   
     int center_chamber=(int) L/2;
    gsl_rng_set(prng, seed);
    
    double Pgrow_norm = Pgrow * dt;

    //calloc initializes to 0
   
    Cs = (double *)calloc(N,sizeof(double));
    Cs_next = (double *)calloc(N,sizeof(double));

    G = (double *)calloc(N,sizeof(double));
    G_next = (double *)calloc(N,sizeof(double));
    R_G_k = (double *)calloc(N,sizeof(double));
    D_G1_k = (double *)calloc(N,sizeof(double));
    G_k1 = (double *)calloc(N,sizeof(double));
    G_k2 = (double *)calloc(N,sizeof(double));
    G_k3 = (double *)calloc(N,sizeof(double));
    G_k4 = (double *)calloc(N,sizeof(double));
    
    Gi = (double *)calloc(N,sizeof(double));
    Gi_next = (double *)calloc(N,sizeof(double));
    R_Gi_k = (double *)calloc(N,sizeof(double));
    Gi_k1 = (double *)calloc(N,sizeof(double));
    Gi_k2 = (double *)calloc(N,sizeof(double));
    Gi_k3 = (double *)calloc(N,sizeof(double));
    Gi_k4 = (double *)calloc(N,sizeof(double));

    Rg=(double *)calloc(N,sizeof(double));
    Rg_next = (double *)calloc(N,sizeof(double));
    R_Rg_k = (double *)calloc(N,sizeof(double));
    Rg_k1 = (double *)calloc(N,sizeof(double));
    Rg_k2 = (double *)calloc(N,sizeof(double));
    Rg_k3 = (double *)calloc(N,sizeof(double));
    Rg_k4 = (double *)calloc(N,sizeof(double));

    S=(double *)calloc(N,sizeof(double));
    S_next=(double *)calloc(N,sizeof(double));
    R_S_k=(double *)calloc(N,sizeof(double));
    S_k1=(double *)calloc(N,sizeof(double));
    S_k2=(double *)calloc(N,sizeof(double));    
    S_k3=(double *)calloc(N,sizeof(double));
    S_k4=(double *)calloc(N,sizeof(double));  

    // slowv=(double *)calloc(N,sizeof(double));
    // slowv_next=(double *)calloc(N,sizeof(double));
    // R_slowv_k=(double *)calloc(N,sizeof(double));
    // slowv_k1=(double *)calloc(N,sizeof(double));
    // slowv_k2=(double *)calloc(N,sizeof(double));    
    // slowv_k3=(double *)calloc(N,sizeof(double));
    // slowv_k4=(double *)calloc(N,sizeof(double));     
    
        
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

    Vls_pr=(double *)calloc(N,sizeof(double));
    Vks_pr=(double *)calloc(N,sizeof(double));

    n1=( int *)calloc(N,sizeof( int));
    n3=( int *)calloc(N,sizeof( int));

    distances=(double *)calloc(N,sizeof(double));
    exp_gammafl_dist=(double *)calloc(N,sizeof(double));
    exp_gammadif_dist=(double *)calloc(N,sizeof(double));
    exp_gammaflG_dist=(double *)calloc(N,sizeof(double));
    exp_gammadifG_dist=(double *)calloc(N,sizeof(double));
    exp_gammaflEk_dist=(double *)calloc(N,sizeof(double));
    exp_gammadifEk_dist=(double *)calloc(N,sizeof(double));

    
    double *aux=(double *)calloc(10,sizeof(double));


    
    for (i=0;i<L;i++){
        for (j=0;j<j_max;j++){
            ijk=ncord(L,i,j,0,0);
            n1[ijk]=ncord(L,i,j,1,0);
            n3[ijk]=ncord(L,i,j,-1,0);
        }
    }    
    
    
    FILE *outpat_c, *outpat_S, *outpat_g; 
    FILE *outpat_v, *outpat_nk, *outpat_Ek, *outpat_Ik; 
    FILE *outpat_gi, *outpat_vls, *outpat_vks;
    FILE *outpat_tht;
    FILE *outpat_dist, *outpat_fl_dist, *outpat_dif_dist;
    FILE *outpat_Rg;
    FILE *outpat_slowv;
    
    FILE *outpat_c_pr, *outpat_g_pr,*outpat_gi_pr ;
    FILE *outpat_S_pr;
    FILE *outpat_v_pr, *outpat_nk_pr, *outpat_Ek_pr, *outpat_Ik_pr; 
    FILE *outpat_tht_pr;
    FILE *outpat_dist_pr;
    
    FILE *inpat_c, *inpat_g, *inpat_gi;
    FILE *inpat_S;
    FILE *inpat_v, *inpat_nk, *inpat_Ek, *inpat_Ik; 
    FILE *inpat_tht;
    FILE *inpat_dist;
    
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

    sprintf(filepath,"%s_S.pat",basepath);
    outpat_S= fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_S);
    fwrite(&npat,sizeof(int),1,outpat_S);

    // sprintf(filepath,"%s_slowv.pat",basepath);
    // outpat_slowv= fopen(filepath,"w");
    // fwrite(&L,sizeof(int),1,outpat_slowv);
    // fwrite(&npat,sizeof(int),1,outpat_slowv);

    sprintf(filepath,"%s_Rg.pat",basepath);
    outpat_Rg= fopen(filepath,"w");
    fwrite(&L,sizeof(int),1,outpat_Rg);
    fwrite(&npat,sizeof(int),1,outpat_Rg);
    
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
    
    
    /*initial conditions*/
    
    printf("setting ic");
    fflush(stdout);
   
    if (read_input>0.0){
    
            
        printf("reading input");
        fflush(stdout);

        sprintf(filepath,"%s_c_pr.inp",basepath_in);
        inpat_c = fopen(filepath,"rb");
        printf("opened input\n");
        fflush(stdout);
        fread(Cs,sizeof(double),N,inpat_c);
        
        sprintf(filepath,"%s_g_pr.inp",basepath_in);
        inpat_g = fopen(filepath,"rb");
        fread(G,sizeof(double),N,inpat_g);

        sprintf(filepath,"%s_gi_pr.inp",basepath_in);
        inpat_gi = fopen(filepath,"rb");
        fread(Gi,sizeof(double),N,inpat_gi);

        sprintf(filepath,"%s_S_pr.inp",basepath_in);
        inpat_S = fopen(filepath,"rb");
        fread(S,sizeof(double),N,inpat_S);
        
        sprintf(filepath,"%s_v_pr.inp",basepath_in);
        inpat_v = fopen(filepath,"rb");
        fread(V,sizeof(double),N,inpat_v);
        
        sprintf(filepath,"%s_nk_pr.inp",basepath_in);
        inpat_nk = fopen(filepath,"rb");
        fread(nk,sizeof(double),N,inpat_nk);
        
        sprintf(filepath,"%s_Ik_pr.inp",basepath_in);
        inpat_Ik = fopen(filepath,"rb");
        fread(Ik,sizeof(double),N,inpat_Ik);

        sprintf(filepath,"%s_Ek_pr.inp",basepath_in);
        inpat_Ek = fopen(filepath,"rb");
        fread(Ek,sizeof(double),N,inpat_Ek);

        sprintf(filepath,"%s_tht_pr.inp",basepath_in);
        inpat_tht = fopen(filepath,"rb");
        fread(ThT, sizeof(double),N,inpat_tht);
        printf("read tht\n");
        fflush(stdout);

        sprintf(filepath,"%s_dist_pr.inp",basepath_in);
        inpat_dist = fopen(filepath,"rb");
        fread(distances, sizeof(int),N,inpat_dist);
        printf("read dist\n");
        fflush(stdout);

        for (ip=0;ip<N;ip++){

        	exp_gammafl_dist[ip]=(exp(-gamma_fl*distances[ip])+a_fl)/(1+a_fl);
            exp_gammadif_dist[ip]=(exp(-gamma_dif*distances[ip])+a_dif)/(1+a_dif);
        }

        left_e=-1;
        right_e=-1;
        left_e_2=-1;
        right_e_2=-1;

        for (ip=0;ip<N-1;ip++){
        	if ((distances[ip]<0.5)&&(distances[ip+1]>0.5)){
        		if (left_e<0){
        			left_e=ip; //here ip==i because it is 1 D
        		}else{
        			left_e_2=ip;
        		}

        	}

        	if ((distances[ip]<0.5) && (distances[ip-1]>0.5)){
			    if (right_e<0){
			    	right_e=ip;
			    }else{
			    	right_e_2=ip;
			    }
			}
        }
    
    }else{ //no reading initial conditions
        for (i=0; i<L; i++){
            for (j=0; j<j_max; j++){
            	if ((i>=left_e) && (i<=right_e) && (Rn>0)){
                        
                    ip=ncord(L,i,j,0,0);
                    Cs[ip] = 1; 
                    G[ip]=G_media;
                    Gi[ip]=Gi01*exp(gsl_ran_gaussian(prng,std_ic));
                    Ek[ip]=K_media*exp(gsl_ran_gaussian(prng,std_ic));
                    S[ip]=0; //pow(10,-8);
                    V[ip]=-156.0;
                    Ik[ip]=300;
                    ThT[ip]=0;
                    Rg[ip]=1;
                    //slowv[ip]=GS0;
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
                    G[ip]=G_media;
                    Gi[ip]=Gi02*exp(gsl_ran_gaussian(prng,std_ic));
                    Ek[ip]=K_media*exp(gsl_ran_gaussian(prng,std_ic));
                    S[ip]=pow(10,-8);
                    V[ip]=-156.0;
                    Ik[ip]=300;
                    ThT[ip]=0;
                    Rg[ip]=1;
                    //slowv[ip]=GS0;
                    if (i <= i_center2){
                        distances[ip]=i-left_e_2;
                    }else{
                        distances[ip]=right_e_2-i;
                    }
                    exp_gammafl_dist[ip]=(exp(-gamma_fl*distances[ip])+a_fl)/(1+a_fl);
                    exp_gammadif_dist[ip]=(exp(-gamma_dif*distances[ip])+a_dif)/(1+a_dif);
                }else{
                    ip=ncord(L,i,j,0,0);
                    G[ip]=G_media;
                    Gi[ip]=0;
                    V[ip]=-156; //I'm setting this to basal value to ease the plotting
                    Cs[ip]=0;
                    Ek[ip]=K_media*exp(gsl_ran_gaussian(prng,std_ic));
                    Ik[ip]=0;
                    ThT[ip]=0;  
                    Rg[ip]=0;
                    //slowv[ip]=0;
                    distances[ip]=0;  
                    exp_gammafl_dist[ip]=1;   
                    exp_gammadif_dist[ip]=1;                   
                } 
            }
        }   
    }

    printf("finished initialisation");
    fflush(stdout);

    if (tprint>0){
        char basepath_in[150];
        sprintf(basepath_in,"%s/%s", argv[5], argv[6]);
        printf("opening ic to write");
        fflush(stdout);
        
        sprintf(filepath,"%s_c_pr.inp",basepath_in);
        outpat_c_pr = fopen(filepath,"w");
      
        sprintf(filepath,"%s_g_pr.inp",basepath_in);
        outpat_g_pr = fopen(filepath,"w");

        sprintf(filepath,"%s_gi_pr.inp",basepath_in);
        outpat_gi_pr = fopen(filepath,"w");

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

        sprintf(filepath,"%s_tht_pr.inp",basepath_in);
        outpat_tht_pr= fopen(filepath,"w");

        sprintf(filepath,"%s_dist_pr.inp",basepath_in);
        outpat_dist_pr= fopen(filepath,"w");

    }
   
    
    
    
    printf("initial glutamate, %g \n",sum_ar(G,N));
    fflush(stdout); 
    
    printf("initial number of cells, %g \n",sum_ar(Cs,N));
    fflush(stdout); 
    printf("N, %d \n",N);
    fflush(stdout);
   

    //Perturbations 

    int it_freeze=(int)(freeze/dt);
    int it_unfreeze=(int)(unfreeze/dt);
    int it_freeze2=(int)(freeze2/dt);
    
    //int npert=(int)(tKpert/dt+0.5);
    //int nwarm=(int)(Twarm/dt+0.5);
    //int nstopKpert=(int)(stopKpert/dt+0.5);
    int dist_pK_dx=(int)(dist_pK/dx+0.5);
    /*
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
    
       
    double W_;

    double aux_pow;
    
   
    for (it=0;it<=nme;it++){

        //printf("%d \n",it);
        //fflush(stdout); 
         if (it > it_perturb_G){

             G_media=GE_p;
        }

        K_media = K_media0;

         if ((it > it_perturb_K1) && (it < it_perturb_K2)){
            for (ip=right_e+dist_pK_dx; ip< N; ip++){
                    Ek[ip]=K_perturb;
                }

            }    

        //     K_media=K_perturb;
        // }else{
        //     K_media=K_media0;
        // }


        if ((it==0) || (it%nst==0)){

            
            t=dt*(double)it;
            printf("\nbefore integrating, %g/%g", t, tt);
            fflush(stdout); 
            
            nan_found=check_nan(G,N);
            if (nan_found>0){
            printf("\n nan found in G");
            fflush(stdout); 
            }

            
            fwrite(&t,sizeof(double),1,outpat_g);
            fwrite(G,sizeof(double),N,outpat_g);
            
            fwrite(&t,sizeof(double),1,outpat_gi);
            fwrite(Gi,sizeof(double),N,outpat_gi);

            fwrite(&t,sizeof(double),1,outpat_Rg);
            fwrite(Rg,sizeof(double),N,outpat_Rg);

            fwrite(&t,sizeof(double),1,outpat_c);
            fwrite(Cs,sizeof(double),N,outpat_c);

            fwrite(&t,sizeof(double),1,outpat_dist);
            fwrite(distances,sizeof(double),N,outpat_dist);

            fwrite(&t,sizeof(double),1,outpat_fl_dist);
            fwrite(exp_gammafl_dist,sizeof(double),N,outpat_fl_dist);

            fwrite(&t,sizeof(double),1,outpat_dif_dist);
            fwrite(exp_gammadif_dist,sizeof(double),N,outpat_dif_dist);

            fwrite(&t,sizeof(double),1,outpat_S);
            fwrite(S,sizeof(double),N,outpat_S);

            //fwrite(&t,sizeof(double),1,outpat_slowv);
            //fwrite(slowv,sizeof(double),N,outpat_slowv);
            
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

            fwrite(&t,sizeof(double),1,outpat_tht);
            fwrite(ThT,sizeof(double),N,outpat_tht);
          
            if (it==it_print){
                //print_ic();
                printf("printing t %d", it);
                fflush(stdout);
                fwrite(G,sizeof(double),N,outpat_g_pr);

                fwrite(Gi,sizeof(double),N,outpat_gi_pr);
                
                fwrite(Cs,sizeof(double),N,outpat_c_pr);
                
                fwrite(S,sizeof(double),N,outpat_S_pr);
                            
                fwrite(V,sizeof(double),N,outpat_v_pr);
                            
                fwrite(nk,sizeof(double),N,outpat_nk_pr);
                            
                fwrite(Ek,sizeof(double),N,outpat_Ek_pr);

                fwrite(Ik,sizeof(double),N,outpat_Ik_pr);

                fwrite(ThT,sizeof(double),N,outpat_tht_pr);

                fwrite(distances,sizeof(double),N,outpat_dist_pr);

            }
        }

        if (nan_found ==0){

            
            if ((it>it_stopflow1)&(it<it_stopflow2)){
                fl_Ek=0;
                fl_G=0;
                
            }else{
                    
                fl_Ek=fl_Ekext;
                fl_G=fl_Gext;

            }
                    
            //calculate *k1
            for (ip=0; ip<N; ip++){
                exp_gammadifG_dist[ip]=D_gdx*exp_gammadif_dist[ip];
                exp_gammadifEk_dist[ip]=D_kedx*exp_gammadif_dist[ip];

                exp_gammaflG_dist[ip]=fl_G*exp_gammafl_dist[ip];
                exp_gammaflEk_dist[ip]=fl_Ek*exp_gammafl_dist[ip];
            }

            for (ip=0;ip<N;ip++){

                D_G1_k[ip]=exp_gammadifG_dist[ip]*(-2*G[ip] + G[n1[ip]]+ G[n3[ip]]);
                D_Ek1_k[ip]=exp_gammadifEk_dist[ip]*(-2*Ek[ip] + Ek[n1[ip]]+ Ek[n3[ip]]);
              
                if (Cs[ip]>zero_tol){

                	//aux=[Ginip, Gentry, aS, Vls, Vks, gkr, W_]

                    Ginip=Rg[ip]/(1+exp(gv*(V[ip]-V0))); //Ginip
                    Gentry=-Ginip*alpha_gt*G[ip]/(Kgt+G[ip]);  //G entry
                    R_G_k[ip]=Gentry+ exp_gammaflG_dist[ip]*(G_media-G[ip]);    
                    R_Gi_k[ip]=-Gentry-delta_g*Gi[ip];
                    aux_pow=pow(Gi[ip],nrg);
                    R_Rg_k[ip]=alpha_Rg-delta_Rg*Rg[ip]+alpha_rRg*aux_pow/(Krgnrg+aux_pow); 
                    
                    R_S_k[ip]=S0/(1+pow(Gi[ip],ug)/Gs0ug)-gamma_s*S[ip]; 
                    
                    aux_pow=pow(S[ip],mk);
                    aS=a0*aux_pow/(Sthm+aux_pow);
                    Vls=Vl+dl*(K_media-Ek[ip])/(exp((K_media-Ek[ip])/0.1)-1.0);
                    Vks=Vk*log(Ek[ip]/Ik[ip]);

                    gkr=-gk*pow(nk[ip],4)*(V[ip]-Vks);
                    W_=-F*gkr-Dp*Ek[ip]*(Ikmax-Ik[ip])*Gi[ip];

                    R_V_k[ip]=gkr-gl*(V[ip]-Vls);
                    R_nk_k[ip]=aS*(1-nk[ip])-b*nk[ip];
                    R_Ek_k[ip]=W_+exp_gammaflEk_dist[ip]*(K_media-Ek[ip]); 
                    R_Ik_k[ip]=-W_;
                    R_ThT_k[ip]=alpha_t/(1+exp(g_tht*(V[ip]-(V0_tht))))-gamma_t*ThT[ip];


                    S_k1[ip]=dth*R_S_k[ip];
                    Gi_k1[ip] =dth*(R_Gi_k[ip]);
                    V_k1[ip]=dth*R_V_k[ip];
                    nk_k1[ip]=dth*R_nk_k[ip];
                    Ik_k1[ip]=dth*(R_Ik_k[ip]);
                    ThT_k1[ip]=dth*R_ThT_k[ip];
                    Rg_k1[ip]=dth*R_Rg_k[ip];

                    S_next[ip]=S[ip]+S_k1[ip];
	                Gi_next[ip] =Gi[ip]+Gi_k1[ip];
	                V_next[ip]=V[ip]+V_k1[ip];
	                nk_next[ip]=nk[ip]+nk_k1[ip];	                
	                Ik_next[ip]=Ik[ip]+Ik_k1[ip];
	                ThT_next[ip]=ThT[ip]+ThT_k1[ip];
                    Rg_next[ip]=Rg[ip]+Rg_k1[ip];
                   
                    
                }else{
                    
                    R_G_k[ip]=exp_gammaflG_dist[ip]*(G_media-G[ip]);
                    R_Ek_k[ip]=exp_gammaflEk_dist[ip]*(K_media-Ek[ip]); 
 
                    //_k1 will be 0 for the others so no need to do again the calculation
                }

                G_k1[ip] =dth*(R_G_k[ip]+D_G1_k[ip]);
                Ek_k1[ip]=dth*(R_Ek_k[ip]+D_Ek1_k[ip]);

                G_next[ip] =G[ip]+G_k1[ip];
                Ek_next[ip]= Ek[ip]+Ek_k1[ip];
                Cs_next[ip]=Cs[ip]; 


            }
            
            //calculate k2
            for (ip=0;ip<N;ip++){

                D_G1_k[ip]=exp_gammadifG_dist[ip]*(-2*G_next[ip] + G_next[n1[ip]]+ G_next[n3[ip]]);
                D_Ek1_k[ip]=exp_gammadifEk_dist[ip]*(-2*Ek_next[ip] + Ek_next[n1[ip]]+ Ek_next[n3[ip]]);
              
                if (Cs_next[ip]>zero_tol){

                    Ginip=Rg_next[ip]/(1+exp(gv*(V_next[ip]-V0)));
                    Gentry=-Ginip*alpha_gt*G_next[ip]/(Kgt+G_next[ip]); 
                    R_G_k[ip]=Gentry+ exp_gammaflG_dist[ip]*(G_media-G_next[ip]);
                    R_Gi_k[ip]=-Gentry-delta_g*Gi_next[ip];
                    aux_pow=pow(Gi_next[ip],nrg);
                    R_Rg_k[ip]=alpha_Rg-delta_Rg*Rg_next[ip]+alpha_rRg*aux_pow/(Krgnrg+aux_pow);                     
                    R_S_k[ip]=S0/(1+pow(Gi_next[ip],ug)/Gs0ug)-gamma_s*S_next[ip]; 
                    
                    aux_pow=pow(S_next[ip],mk);
                    aS=a0*aux_pow/(Sthm+aux_pow);
                    Vls=Vl+dl*(K_media-Ek_next[ip])/(exp((K_media-Ek_next[ip])/0.1)-1.0);
                    Vks=Vk*log(Ek_next[ip]/Ik_next[ip]);

                    gkr=-gk*pow(nk_next[ip],4)*(V_next[ip]-Vks);
                    W_=-F*gkr-Dp*Ek_next[ip]*(Ikmax-Ik_next[ip])*Gi_next[ip];

                    R_V_k[ip]=gkr-gl*(V_next[ip]-Vls);
                    R_nk_k[ip]=aS*(1-nk_next[ip])-b*nk_next[ip];
                    R_Ek_k[ip]=W_+exp_gammaflEk_dist[ip]*(K_media-Ek_next[ip]); 
                    R_Ik_k[ip]=-W_;
                    R_ThT_k[ip]=alpha_t/(1+exp(g_tht*(V_next[ip]-(V0_tht))))-gamma_t*ThT_next[ip];

                    S_k2[ip]=dt*R_S_k[ip];	                
	                Gi_k2[ip] =dt*(R_Gi_k[ip]);
	                V_k2[ip]=dt*R_V_k[ip];
	                nk_k2[ip]=dt*R_nk_k[ip];	                
	                Ik_k2[ip]=dt*(R_Ik_k[ip]);
	                ThT_k2[ip]=dt*(R_ThT_k[ip]);
                    Rg_k2[ip]=dt*R_Rg_k[ip];

	                S_next[ip]=S[ip]+S_k2[ip]/2;	                
	                Gi_next[ip] =Gi[ip]+Gi_k2[ip]/2;
	                V_next[ip]=V[ip]+V_k2[ip]/2;
	                nk_next[ip]=nk[ip]+nk_k2[ip]/2;	                
	                Ik_next[ip]=Ik[ip]+Ik_k2[ip]/2;
	                ThT_next[ip]=ThT[ip]+ThT_k2[ip]/2;
                    Rg_next[ip]=Rg[ip]+Rg_k2[ip]/2;
	                   
                    
                }else{
                    
                    R_G_k[ip]=exp_gammaflG_dist[ip]*(G_media-G_next[ip]);
                    R_Ek_k[ip]=exp_gammaflEk_dist[ip]*(K_media-Ek_next[ip]); 
                }
                G_k2[ip] =dt*(R_G_k[ip]+D_G1_k[ip]);
                Ek_k2[ip]=dt*(R_Ek_k[ip]+D_Ek1_k[ip]);

                G_next[ip] =G[ip]+G_k2[ip]/2;
                Ek_next[ip]= Ek[ip]+Ek_k2[ip]/2;
            }
                  

            //k3
            for (ip=0;ip<N;ip++){

                D_G1_k[ip]=exp_gammadifG_dist[ip]*(-2*G_next[ip] + G_next[n1[ip]]+ G_next[n3[ip]]);
                D_Ek1_k[ip]=exp_gammadifEk_dist[ip]*(-2*Ek_next[ip] + Ek_next[n1[ip]]+ Ek_next[n3[ip]]);
              
                if (Cs_next[ip]>zero_tol){

                    Ginip=Rg_next[ip]/(1+exp(gv*(V_next[ip]-V0)));
                    Gentry=-Ginip*alpha_gt*G_next[ip]/(Kgt+G_next[ip]); 
                    R_G_k[ip]=Gentry+ exp_gammaflG_dist[ip]*(G_media-G_next[ip]);
                    R_Gi_k[ip]=-Gentry-delta_g*Gi_next[ip];
                    aux_pow=pow(Gi_next[ip],nrg);
                    R_Rg_k[ip]=alpha_Rg-delta_Rg*Rg_next[ip]+alpha_rRg*aux_pow/(Krgnrg+aux_pow); 
                    
                    R_S_k[ip]=S0/(1+pow(Gi_next[ip],ug)/Gs0ug)-gamma_s*S_next[ip]; 
                    
                    aux_pow=pow(S_next[ip],mk);
                    aS=a0*aux_pow/(Sthm+aux_pow);
                    Vls=Vl+dl*(K_media-Ek_next[ip])/(exp((K_media-Ek_next[ip])/0.1)-1.0);
                    Vks=Vk*log(Ek_next[ip]/Ik_next[ip]);

                    gkr=-gk*pow(nk_next[ip],4)*(V_next[ip]-Vks);
                    W_=-F*gkr-Dp*Ek_next[ip]*(Ikmax-Ik_next[ip])*Gi_next[ip];

                    R_V_k[ip]=gkr-gl*(V_next[ip]-Vls);
                    R_nk_k[ip]=aS*(1-nk_next[ip])-b*nk_next[ip];
                    R_Ek_k[ip]=W_+exp_gammaflEk_dist[ip]*(K_media-Ek_next[ip]); 
                    R_Ik_k[ip]=-W_;
                    R_ThT_k[ip]=alpha_t/(1+exp(g_tht*(V_next[ip]-(V0_tht))))-gamma_t*ThT_next[ip];

                    S_k3[ip]=dt*R_S_k[ip];                
	                Gi_k3[ip] =dt*(R_Gi_k[ip]);
	                V_k3[ip]=dt*R_V_k[ip];
	                nk_k3[ip]=dt*R_nk_k[ip];	                
	                Ik_k3[ip]=dt*(R_Ik_k[ip]);
	                ThT_k3[ip]=dt*R_ThT_k[ip];
                    Rg_k3[ip]=dt*R_Rg_k[ip];

	                S_next[ip]=S[ip]+S_k3[ip];	                
	                Gi_next[ip] =Gi[ip]+Gi_k3[ip];
	                V_next[ip]=V[ip]+V_k3[ip];
	                nk_next[ip]=nk[ip]+nk_k3[ip];	                
	                Ik_next[ip]=Ik[ip]+Ik_k3[ip];
	                ThT_next[ip]=ThT[ip]+ThT_k3[ip];
                    Rg_next[ip]=Rg[ip]+Rg_k3[ip];
                   
                    
                }else{
                    
                    R_G_k[ip]=exp_gammaflG_dist[ip]*(G_media-G_next[ip]);
                    R_Ek_k[ip]=exp_gammaflEk_dist[ip]*(K_media-Ek_next[ip]); 
                    
                }
                G_k3[ip] =dt*(R_G_k[ip]+D_G1_k[ip]);
                Ek_k3[ip]=dt*(R_Ek_k[ip]+D_Ek1_k[ip]);
                
                //Cs_next[ip]=Cs_next[ip]; 
                G_next[ip] =G[ip]+G_k3[ip];
                Ek_next[ip]= Ek[ip]+Ek_k3[ip];
            } //k3

            //for k4
            for (ip=0;ip<N;ip++){

                D_G1_k[ip]=exp_gammadifG_dist[ip]*(-2*G_next[ip] + G_next[n1[ip]]+ G_next[n3[ip]]);
                D_Ek1_k[ip]=exp_gammadifEk_dist[ip]*(-2*Ek_next[ip] + Ek_next[n1[ip]]+ Ek_next[n3[ip]]);
              
                if (Cs_next[ip]>zero_tol){

                    Ginip=Rg_next[ip]/(1+exp(gv*(V_next[ip]-V0)));
                    Gentry=-Ginip*alpha_gt*G_next[ip]/(Kgt+G_next[ip]); 
                    R_G_k[ip]=Gentry+ exp_gammaflG_dist[ip]*(G_media-G_next[ip]);
                    R_Gi_k[ip]=-Gentry-delta_g*Gi_next[ip];
                    aux_pow=pow(Gi_next[ip],nrg);
                    R_Rg_k[ip]=alpha_Rg-delta_Rg*Rg_next[ip]+alpha_rRg*aux_pow/(Krgnrg+aux_pow); 
                    
                    R_S_k[ip]=S0/(1+pow(Gi_next[ip],ug)/Gs0ug)-gamma_s*S_next[ip]; 
                    
                    aux_pow=pow(S_next[ip],mk);
                    aS=a0*aux_pow/(Sthm+aux_pow);                    
                    Vls=Vl+dl*(K_media-Ek_next[ip])/(exp((K_media-Ek_next[ip])/0.1)-1.0);
                    Vks=Vk*log(Ek_next[ip]/Ik_next[ip]);

                    gkr=-gk*pow(nk_next[ip],4)*(V_next[ip]-Vks);
                    W_=-F*gkr-Dp*Ek_next[ip]*(Ikmax-Ik_next[ip])*Gi_next[ip];

                    R_V_k[ip]=gkr-gl*(V_next[ip]-Vls);
                    R_nk_k[ip]=aS*(1-nk_next[ip])-b*nk_next[ip]; 
                    R_Ek_k[ip]=W_+exp_gammaflEk_dist[ip]*(K_media-Ek_next[ip]); 
                    R_Ik_k[ip]=-W_;
                    R_ThT_k[ip]=alpha_t/(1+exp(g_tht*(V_next[ip]-(V0_tht))))-gamma_t*ThT_next[ip];

                    S_k4[ip]=dt6*R_S_k[ip];               
                	Gi_k4[ip] =dt6*(R_Gi_k[ip]);
                	V_k4[ip]=dt6*R_V_k[ip];
                	nk_k4[ip]=dt6*R_nk_k[ip];                
	                Ik_k4[ip]=dt6*(R_Ik_k[ip]);
	                ThT_k4[ip]=dt6*R_ThT_k[ip];
                    Rg_k4[ip]=dt6*R_Rg_k[ip];
                   
                    
                }else{
                    
                    R_G_k[ip]=exp_gammaflG_dist[ip]*(G_media-G_next[ip]);
                    R_Ek_k[ip]=exp_gammaflEk_dist[ip]*(K_media-Ek_next[ip]); 
                    
                }

                G_k4[ip] =dt6*(R_G_k[ip]+D_G1_k[ip]);
                Ek_k4[ip]=dt6*(R_Ek_k[ip]+D_Ek1_k[ip]);
            }

            
            
            for (ip=0;ip<N;ip++){
               
                S[ip]=S[ip]+S_k1[ip]/3+S_k2[ip]/3+S_k3[ip]/3+S_k4[ip];
                G[ip]=G[ip]+G_k1[ip]/3+G_k2[ip]/3+G_k3[ip]/3+G_k4[ip];
                Gi[ip]=Gi[ip]+Gi_k1[ip]/3+Gi_k2[ip]/3+Gi_k3[ip]/3+Gi_k4[ip];
                V[ip]=V[ip]+V_k1[ip]/3+V_k2[ip]/3+V_k3[ip]/3+V_k4[ip];
                nk[ip]=nk[ip]+nk_k1[ip]/3+nk_k2[ip]/3+nk_k3[ip]/3+nk_k4[ip];
                Ek[ip]= Ek[ip]+Ek_k1[ip]/3+Ek_k2[ip]/3+Ek_k3[ip]/3+Ek_k4[ip]; 
                Ik[ip]=Ik[ip]+Ik_k1[ip]/3+Ik_k2[ip]/3+Ik_k3[ip]/3+Ik_k4[ip]; 
                ThT[ip]=ThT[ip]+ThT_k1[ip]/3+ThT_k2[ip]/3+ThT_k3[ip]/3+ThT_k4[ip];
                Rg[ip]=Rg[ip]+Rg_k1[ip]/3+Rg_k2[ip]/3+Rg_k3[ip]/3+Rg_k4[ip];

            }

            if ((it<it_freeze) |  ((it> it_unfreeze)&&(it<it_freeze2))){
               
               // I have calculated the left and right edges at the beginning and I am going to keep track of it all the time
                
                if (Rn>0){

                    if (gsl_rng_uniform(prng) < (Pgrow_norm * Gi[right_e])){                        
                       
                        grow(right_e+1, right_e,L, Cs, Gi, Rg, S, nk, V, Ik, ThT, distances, exp_gammafl_dist, exp_gammadif_dist, gamma_fl, a_fl, gamma_dif, a_dif);
                        //left_e=left_e-1;
                        right_e=right_e+1;                  
                    }
                }

                if (Rn2>0){                 

                    if (gsl_rng_uniform(prng) < (Pgrow_norm * Gi[right_e_2])){   

                        grow(left_e_2-1, left_e_2,L, Cs, Gi, Rg, S,nk, V, Ik, ThT, distances, exp_gammafl_dist, exp_gammadif_dist,gamma_fl, a_fl, gamma_dif, a_dif);
                                                                    
                        left_e_2=left_e_2-1;
                        //right_e_2=right_e_2+1;                  
                    }
                }
            }//if it<
        } //if nan_found   
    }//for it
   
    fclose(outpat_c);
    fclose(outpat_g);
    fclose(outpat_gi);
    fclose(outpat_Rg);
    fclose(outpat_S);
    //fclose(outpat_slowv);
    fclose(outpat_v);
    fclose(outpat_vls); 
    fclose(outpat_vks);
    fclose(outpat_nk);
    fclose(outpat_Ek);
    fclose(outpat_Ik);
    fclose(outpat_tht);
    fclose(outpat_dist);
    fclose(outpat_dif_dist);
    fclose(outpat_fl_dist);
    printf("Quitting ");
    fflush(stdout); 
}// main
