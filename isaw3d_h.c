/*
The idea is to run in parallel different Boltzmann distributions (chains) 
each characterized by a given temperature (so there is a 1-1 
correspondence between chain and temperature) and during
the run to make some swapping between any chain pairs.
The swap between two chain, let's say ind1 and ind2,
will put the configuration previosuly located in the chain ind1
(labeled by the index bb[ind1]) into the chain ind2 and
viceversa.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "marsagliazo.h"
#include "array_alloc.h"
#include "hash_chain.h"
#include "pivot_chain_ete.h"
#include "local_chain_ete.h"
#include "ising_saw3d_h.h"

#define TRUE 1
#define FALSE 0
#define MINI -2147483648


FILE *input;
FILE *output[30];
FILE *outconf[30];
FILE *out_save;
FILE *out_dia;         /* output file with the diagnostic   */
FILE *out_stat;
FILE *out_dat;
FILE *out_tab;
FILE *out_ave;
FILE *out_mag_series;
FILE *out_ene_series;
FILE *out_cont_series;
FILE *out_RG_series;
FILE *out_x_series;
int n,mm; //Qui è definito n

static int ***iw,***htab,**sigma,**iv;

int npi;         /* check on lenght of  pivot (non used in this program)*/
int nlocr,nloca; /* chech on the very local moves (not used here )*/
int mloc;        /* how many local moves we want to perform */
int mpibi;       /* every how many piv we want local moves */
int lx,ly,lz;             /* span of the walk in a box */
unsigned int iseed;
char str_conf[30][30];

void FindAutoTau(double *var,int no,double *ave,double *err,double *tau,double *ertau);
void write_chain (int ***iw,int n,int **sigma,int chain,int kc);

int main( int argn, char **argc) //Main
{
	int Delta = 5; 
	int K = 3;
	char *str_in;
	char str_out[30][30]; 
	char str_dia[30];
	char str_stat[30];
        char str_dat[30];
        char str_tab[30];
        char str_ave[30];
	int ivv[3];
        int int_dia = 250000;        
	int *bb;
	double *bet;               /* Array with the temperatures */
	double *bet_old;               /* Array with the temperatures */
	double betta,betty,h;
        int *c;                             
        int *cen;                     /* number of hydrophobic contacts */
	double ***aver,***erro;
	double *amed,*vamed;
	double base,expo;
	double blow,bhigh,db,deno,dno;
        double cden,cdensw,ca;
	double dnorm,amedr,amedr2;
	int conf,conf1,conf2;  /* indices of the configurations
	                          respectively in the chain
				  k,k-1 and k+1 */
        int cc,icc,ccen;
	int ind,ind1,ind2;
        int i,k,j,nc,ii,m,ik,im,l,kl,kr;
	int nobs,no;
        int quantum = 200 ;   
        int cycles = 0;
        int cycle_limit = 10000000 ;
        int relax;              /* number of iterations to reach equilibrium */


        double xcm,ycm,zcm;   	          /* coord. of c.m. for the walk  */

	double dx,dy,dz,dn,dn2;
	double rx,ry,rz,r; 
	double xx,yy,zz,xyz;
	double xy,xz,yz;
	double x2,y2,z2;
        double r2,tr2;                  
        double rx2,ry2,rz2,rxz,ryz,rxy;       
	double dr;
	int Mn,lp1,lp2;

        double e2,te2;                    /* squared end-to end */	
        double **vr2,*ve2,*vc,**vce,**vma,*vma2,*vma22,*vama, *vx, *vx2, **vxt;
        double *vr2pho,*vr2phi;
		int *vfreq, *vfreqt;
        double *ve22,*vc2,*vc22,*vce2,*vce22;
        double *vr2pho2,*vr2phi2;

	double ave,err,tau,ertau;
	double **av,**er,**ta,**erta,*var;
        double dperc,den;
	int icont,swapint;
	int nchain, is,mi,mc,ibl;
        int iprint = 0;
	int npho,nphi;

	double cspe,ercspe,cspc,ercspc,erchi,chi, chix, erchix;
        double qtmp;
	double xxpho,yypho,zzpho;
	double xypho,xzpho,yzpho;
	double x2pho,y2pho,z2pho;
        double r2pho;                   /* square of the radius of gyration */
        double rx2pho,ry2pho,rxypho,rz2pho,rxzpho,ryzpho; /*Comp of Inert tens*/

        double mag,psig,amag, sqa, x;
	double xxphi,yyphi,zzphi,xyzphi;
	double xyphi,xzphi,yzphi; 
	double x2phi,y2phi,z2phi;
        double r2phi;                   /* square of the radius of gyration */
        double rx2phi,ry2phi,rz2phi,rxyphi,rxzphi,ryzphi;/* Component of Inertia tensor*/

        str_in = (char *)cvector(0,29);

        printf("File with the walk please ..\n");
        scanf("%s",str_in);
        if((input = fopen(str_in,"r"))==NULL)
        {
          printf("non riesco ad aprire il file di input!!\n");
          exit(1);
        }
	
        printf("Number of chains please.... ?\n");
	scanf("%d",&nchain);
	bet = (double *)dvector(0,nchain-1);    
	bet_old = (double *)dvector(0,nchain-1);    
	printf("blow = ?, bhigh = ?\n");
	scanf("%lf %lf",&blow, &bhigh);


	for(i=0;i<nchain;i++)
	{
	  printf(" beta for the %d chain ?\n",i);
	  scanf("%lf",&betta);
	  printf(" beta for the %d chain = %lf\n",i,betta);
	  bet[i] = exp(betta);
          bet_old[i]=betta;
        }


	for (i = 0;i<nchain;i++)
	{
        printf("File in output for the %d th chain ?\n",i);
          scanf("%s",str_conf[i]);
	  
          printf("beta[%d] (contacts) = ? %lf\n",i,bet[i]);
	}

        printf("file name for the diagnostic?\n");
        scanf("%s",str_dia);

        printf("file name for the statistic?\n");
        scanf("%s",str_stat);

        printf("file name for the data?\n");
        scanf("%s",str_dat);

        printf("file name for the plot?\n");
        scanf("%s",str_tab);

        printf("file name for the ave?\n");
        scanf("%s",str_ave);

        printf("quantum=?\n");
        scanf("%d",&quantum);
        printf("cycle_limit=?\n");
        scanf("%d",&cycle_limit);
        nobs = cycle_limit/quantum; 
	printf("relaxation = ? \n");
	scanf("%d",&relax);
        printf("seed=?\n");
        scanf("%d",&iseed);
        printf("how many local moves we want to perform = ?\n");
        scanf("%d",&mloc);
        printf("every how many piv  we want local moves = ?\n");
        scanf("%d",&mpibi);
	printf("length of the walk = ?\n");
	scanf("%d",&n);  
	printf("length of the hash-table = ?\n");
	scanf("%d",&mm);  
        printf("every how many moves we want a swap = ?\n");
        scanf("%d",&swapint);
        printf("after how many cycles we wanna diagnostic=?\n");
        scanf("%d",&int_dia);
	printf("psig\n");
	scanf("%lf",&psig);

	Mn = n*0.50;
	dr = 0.5;
        dn = 1./((double)n);
	dn2 = dn*dn;

        fprintf(stderr,"sto allocando... \n");
        fprintf(stderr,"nobs = %d  \n",nobs);

        iw = (int ***)imatrix3(n,3,nchain);
        iv = (int **)imatrix(0,n,0,3);
        htab = (int ***)imatrix3(mm,4,nchain);
        bb = (int *)ivector(0,nchain-1);
	c = (int *)ivector(0,nchain-1);    
	cen = (int *)ivector(0,nchain-1);    
	sigma = (int **)imatrix(0,n,0,nchain);    

        vr2 = (double **)dmatrix(0,nobs,0,nchain-1);
        vce = (double **)dmatrix(0,nobs,0,nchain-1);
        vma = (double **)dmatrix(0,nobs,0,nchain-1);
        vxt = (double **)dmatrix(0,nobs,0,nchain-1);
        vx = (double *)dvector(0,nchain-1);
        vx2 = (double *)dmatrix(0,nobs,0,nchain-1);
        vama = (double *)dvector(0,nchain-1);
        vma2 = (double *)dvector(0,nchain-1);
        vma22 = (double *)dvector(0,nchain-1);

        ve2 = (double *)dvector(0,nchain-1);
        vc = (double *)dvector(0,nchain-1);
        vr2pho = (double *)dvector(0,nchain-1);
        vr2phi = (double *)dvector(0,nchain-1);

        ve22 = (double *)dvector(0,nchain-1);
        vc2 = (double *)dvector(0,nchain-1);
        vc22 = (double *)dvector(0,nchain-1);
        vce2 = (double *)dvector(0,nchain-1);
        vce22 = (double *)dvector(0,nchain-1);
        vr2pho2 = (double *)dvector(0,nchain-1);
        vr2phi2 = (double *)dvector(0,nchain-1);
		vfreq = (int *)dvector(0, nchain-1);
		vfreqt = (int *)dvector(0, nchain-1);

        av = (double **)dmatrix(0,14,0,nchain-1);
        er = (double **)dmatrix(0,14,0,nchain-1);
        ta= (double **)dmatrix(0,14,0,nchain-1);
        erta = (double **)dmatrix(0,14,0,nchain-1);
	var = (double *)dvector(0,nobs); 


       fprintf(stderr,"Ho finito l'allocazione\n");
	       for(k=0;k<nchain;k++)
	       {
                vr2pho[k] = vr2phi[k] = 0;
	        ve2[k] = vc[k]=  0;
                vr2pho2[k] = vr2phi2[k] = 0;
	        ve22[k] = vc2[k] = vc22[k] = 0;
	        vce2[k] = vce22[k] =  0;
                vma2[k] = vama[k]= vma22[k]=0;
	       }
       


       for(k=0;k<n;k++) for(j=0;j<3;j++)fscanf(input,"%d",&iv[k][j]);

       initran_(&iseed); 
       dinuovo: no = 0;
       cycles = 0;
       icc = 0;

		double randomNumber;
       for(k=0;k<nchain;k++)
       {
        sigma[0][k] = 1;
        ccen = 0 - Delta*sigma[0][k]*sigma[0][k];
        for(i=1;i<n;i++){
        randomNumber = ran1double_();
         if(randomNumber < 0.333){sigma[i][k] = -1;}
         else if (randomNumber < 0.666 && randomNumber >= 0.333 ) {sigma[i][k] = 0;}
         else {sigma[i][k] = 1;}
         ccen = ccen + sigma[i][k]*sigma[i-1][k] + K*pow(sigma[i][k]*sigma[i-1][k], 2) - Delta*sigma[i][k]*sigma[i][k]; //Perché inizialmente la catena è dritta
        
        }
         cen[k] = ccen;
       }

       for(i = 0; i<nchain;i++) bb[i] = i;
       sethash_c(htab,mm,nchain); /* initialization of the hash tables */

       icont = 0;
	mc=0;

/* I'm going to read the vertice's coordinates from the input file */
fprintf(stderr," I'm going to read the vertice's coordinates from the input file\n"); 

       for(i=0;i<nchain;i++)for(k=0;k<n;k++)for(j=0;j<3;j++)iw[k][j][i] = iv[k][j];


/* filling the hash tables with the vertices of the first configuration  */
        for(i=0;i<n;i++)
		{
          for(j=0;j<3;j++) ivv[j] = iw[i][j][0];
          for(k=0;k<nchain;k++)addhash_c(ivv,htab,mm,i,k);
		}		                       

       for (k=0;k<nchain;k++)c[k] = 0;/*true if the initial config. is a line */
       for (k=0;k<nchain;k++)fprintf(stderr,"cen[%d] = %d\n",k,cen[k]);

        fprintf(stderr,"Ho  finito l'assegnazione\n");

	setpivot(n);
        setlocal(n);
        setising(n);
       
	fprintf(stderr,"length of the walk = %d\n",n);
	fprintf(stderr,"length of the hash-table = %d\n",mm);
        fprintf(stderr,"%d %d %d %d %d\n", n,quantum,iseed,cycle_limit,swapint);



/* perform  relaxation iterations for the nchain chains: */


	for(i=0;i<relax;i++) 
	{
	  for(k=0;k<nchain;k++)
	  {
	    cc = c[k];
	    ccen = cen[k];
	    pivot_chain_ete(iw,n,htab,mm,&cc,&ccen,sigma,bet[k],k, K); 
      if(mloc)im =local_chain_ete(iw,n,htab,mm,&cc,&ccen,sigma,bet[k],k, K); 
           ising_saw3d_h(iw,n,htab,mm,&ccen,sigma,bet[k],h,k, Delta, K);
           cen[k] = ccen;
	   c[k] = cc;
          }
	}


/* Beginning of the big cycle */

           fprintf(stderr,"prima di big cyckle\n");
           
            printf("prima di big cyckle\n");
           if ((out_mag_series= fopen("datiM.csv" ,"w"))==NULL) {
    	perror("Error opening file");
    	return 1; 
		}
		
		if ((out_ene_series= fopen("datiE.csv" ,"w"))==NULL) {
    	perror("Error opening file");
    	return 1; 
		}
		
		if ((out_cont_series= fopen("datiC.csv" ,"w"))==NULL) {
    	perror("Error opening file");
    	return 1; 
		}
		
		if ((out_RG_series= fopen("datiRG.csv" ,"w"))==NULL) {
    	perror("Error opening file");
    	return 1; 
		}
		
		if ((out_x_series= fopen("datiX.csv" ,"w"))==NULL) {
    	perror("Error opening file");
    	return 1; 
		}
                
	while (TRUE) 
	{
	   for (i=0;i<quantum;i++) 
	   {
	     icont += 1;
	     
	 
	     if(icont%swapint == 0)
	     {                                  /* Try the swapping */
	      icc = icc +1;						
              do
	      {
	      m=nchain;
	      ind1 = ranint_(&m);
	      m=nchain;
	      ind2 = ranint_(&m);
	      }while((ind1-ind2) != 1 && (ind1-ind2) != -1 ); 
 		  vfreqt[ind1] =  vfreqt[ind1] + 1;
 		  vfreqt[ind2] =  vfreqt[ind2] + 1;
/*
You need the part above if you want to swap between nearest
neighbours chains; if not, the statements: 
              do
              {
	      m=nchain;
	      ind1 = ranint_(&m);
	      m=nchain;
	      ind2 = ranint_(&m);
              }while(ind1 == ind2);

are enough.
*/

	       expo = (double)(cen[bb[ind2]] - cen[bb[ind1]]);
	       base = bet[ind1]/bet[ind2];
	       if(pow(base,expo) > ran1double_())
	       {                                 /* make the swapping */
               conf  = bb[ind1];			
	       bb[ind1] = bb[ind2];
	       bb[ind2] = conf;
	       vfreq[ind1] = vfreq[ind1] + 1;
	       vfreq[ind2] = vfreq[ind2] + 1;
	       }

	     }    

        for(k=0;k<nchain;k++) 
	{       
	  cc = c[bb[k]];
	  ccen = cen[bb[k]];
          pivot_chain_ete(iw,n,htab,mm,&cc,&ccen,sigma,bet[k],bb[k], K); 
          c[bb[k]] = cc;
          cen[bb[k]] = ccen;

          /* local moves every mpibi attempetd pivots    */ 

          if((i+1)%mpibi == 0)
          {
	   for (ii=0;ii<mloc;ii++) 
	   {
	    cc = c[bb[k]];
	    ccen = cen[bb[k]];
            im = local_chain_ete(iw,n,htab,mm,&cc,&ccen,sigma,bet[k],bb[k], K);
            c[bb[k]] = cc;
            cen[bb[k]] = ccen;
           }               
          }

	 ccen = cen[bb[k]];
         ising_saw3d_h(iw,n,htab,mm,&ccen,sigma,bet[k],h,bb[k], Delta, K);
         cen[bb[k]] = ccen;

         } 
       }   

           cycles = cycles + quantum;
	   no++;
 
/* Printing on the output files after quantum trials */

          for (k=0;k<nchain;k++)
	  {
/*   Compute the radius of gyration for the polygon */ 
            xx = yy = zz = xyz= 0;
            xy = yz = xz = 0;
	    x2 = y2 = z2 = 0;

            mag = amag= 0;
            sqa = 0.0;
            xxpho = yypho = zzpho  = 0;
            xypho = yzpho = xzpho =0;
	    x2pho = y2pho = z2pho = 0;

            xxphi = yyphi = zzphi = xyzphi =0;
            xyphi = yzphi = xzphi =0;
	    x2phi = y2phi = z2phi = 0;
            npho = 0;
            nphi = 0;

	    for(i=0;i<n;i++)
	    { 
               mag = mag + sigma[i][bb[k]];
               sqa = sqa + sigma[i][bb[k]]*sigma[i][bb[k]];
               }
            mag = mag /((double) n);
            x = 1. - sqa/((double) n);
        
            if(mag < 0)amag = -mag;
            else amag = mag;
            
              fprintf(out_mag_series,"%.4lf",mag);
		fprintf(out_mag_series,"\t");
		
		fprintf(out_x_series,"%.4lf",x);
		fprintf(out_x_series,"\t");
            
		
        
             
	    for(i=0;i<n;i++)
	    {
	       dx = (double)(iw[i][0][bb[k]]);
	       dy = (double)(iw[i][1][bb[k]]);
	       dz = (double)(iw[i][2][bb[k]]);
	       xx += dx;
	       yy += dy;
               zz += dz;
	       x2 += dx*dx;
	       y2 += dy*dy;
               z2 += dz*dz; 

	       if(sigma[i][bb[k]]< 0)
	       {
               nphi ++;
	       xxphi += dx;
	       yyphi += dy;
               zzphi += dz;
               z2phi += dz*dz;
	       x2phi += dx*dx;
	       y2phi += dy*dy;
	       }
	       else
	       {
               npho ++;
	       xxpho += dx;
	       yypho += dy;
	       zzpho += dz;
	       x2pho += dx*dx;
	       y2pho += dy*dy;
               z2pho += dz*dz; 
	       }
	     }

	     xcm = xx*dn;
	     ycm = yy*dn;
             zcm = zz*dn; 

             rx2 = x2*dn-(xx*xx)*dn2; 
             ry2 = y2*dn-(yy*yy)*dn2; 
             rz2 = z2*dn-(zz*zz)*dn2; 

             if(npho != 0){
             rx2pho = x2pho/(double)(npho)-(xxpho*xxpho)/(double)(npho*npho); 
             ry2pho = y2pho/(double)(npho)-(yypho*yypho)/(double)(npho*npho); 
             rz2pho = z2pho/(double)(npho)-(zzpho*zzpho)/(double)(npho*npho); 
             }
             else {rx2pho = ry2pho = rz2pho = 0;}

             if(nphi != 0){
             rx2phi = x2phi/(double)(nphi)-(xxphi*xxphi)/(double)(nphi*nphi); 
             ry2phi = y2phi/(double)(nphi)-(yyphi*yyphi)/(double)(nphi*nphi); 
             rz2phi = z2phi/(double)(nphi)-(zzphi*zzphi)/(double)(nphi*nphi); 
             }
             else {rx2phi = ry2phi = rz2phi = 0;}

	     r2 = rx2 + ry2 + rz2  ;
	     r2pho = rx2pho + ry2pho +rz2pho  ;
	     r2phi = rx2phi + ry2phi + rz2phi ;

             /* span_c(iw,n,bb[k]);*/

              e2 = (iw[n-1][0][bb[k]]-iw[0][0][bb[k]])*
                      (iw[n-1][0][bb[k]]-iw[0][0][bb[k]])+
                      (iw[n-1][1][bb[k]]-iw[0][1][bb[k]])*
                      (iw[n-1][1][bb[k]]-iw[0][1][bb[k]])+
                     (iw[n-1][2][bb[k]]-iw[0][2][bb[k]])*
                      (iw[n-1][2][bb[k]]-iw[0][2][bb[k]]);

               vr2[no][k] = r2;
	       vce[no][k] = cen[bb[k]];
	       vma[no][k] = mag;
	       vxt[no][k] = x;
	    

               vr2pho[k] = vr2pho[k] + r2pho;
               vr2phi[k] = vr2phi[k] + r2phi;
	       ve2[k] = ve2[k] + e2;
	       vc[k] = vc[k] + c[bb[k]];
               vama[k] = vama[k] + amag;

              vr2pho2[k] = vr2pho2[k] + r2pho*r2pho;
              vr2phi2[k] = vr2phi2[k] + r2phi*r2phi;
	      ve22[k] = ve22[k] + e2*e2;
	      vc2[k] = vc2[k] + c[bb[k]]*c[bb[k]];
	      vce2[k] = vce2[k] + cen[bb[k]]*cen[bb[k]];
	      vce22[k] = vce22[k] + cen[bb[k]]*cen[bb[k]]*cen[bb[k]]*cen[bb[k]];
	      vc22[k] = vc22[k] + c[bb[k]]*c[bb[k]]*c[bb[k]]*c[bb[k]];
	      vma2[k] = vma2[k] + mag*mag;
	      vma22[k] = vma22[k] + mag*mag*mag*mag;
	      vx[k] = vx[k] + x;
	      vx2[k] = vx2[k] + pow(x,2);



  		fprintf(out_cont_series,"%d",c[bb[k]]);
		fprintf(out_cont_series,"\t");
		
		fprintf(out_ene_series,"%d",cen[bb[k]]);
		fprintf(out_ene_series,"\t");
	      
	    fprintf(out_RG_series,"%.4lf",r2);
		fprintf(out_RG_series,"\t");
                

	  }  
	  
	     fprintf(out_cont_series,"\n");
     fprintf(out_ene_series,"\n");
      fprintf(out_RG_series,"\n");
      fprintf(out_mag_series,"\n");
      fprintf(out_x_series,"\n");
	  


          if(!(cycles%int_dia))
          {

	    dno = 1/(double)no;
            for(k=0;k<nchain;k++)
	    {
	    for(ik=0;ik<no;ik++)var[ik] = vce[ik][k]; 
            FindAutoTau(var,no,&ave,&err,&tau,&ertau); 
	    av[1][k] = ave;er[1][k]=err;ta[1][k]=tau;erta[1][k]=ertau;

	    for(ik=0;ik<no;ik++)var[ik] = vr2[ik][k];
            FindAutoTau(var,no,&ave,&err,&tau,&ertau);
	    av[2][k] = ave;er[2][k]=err;ta[2][k]=tau;erta[2][k]=ertau;

	    av[0][k] = vc[k]*dno;
	    er[0][k]= sqrt(2*ta[1][k]*dno*(vc2[k]*dno-av[0][k]*av[0][k])) ;

	    av[3][k] = ve2[k]*dno;
	    er[3][k]= sqrt(2*ta[2][k]*dno*(ve22[k]*dno-av[3][k]*av[3][k]));

	    av[4][k] = vr2pho[k]*dno;
	    er[4][k]= sqrt(2*ta[2][k]*dno*(vr2pho2[k]*dno-av[4][k]*av[4][k]));

	    av[5][k] = vr2phi[k]*dno;
	    er[5][k]= sqrt(2*ta[2][k]*dno*(vr2phi2[k]*dno-av[5][k]*av[5][k]));

	    for(ik=0;ik<no;ik++)var[ik] = vma[ik][k]; 
            FindAutoTau(var,no,&ave,&err,&tau,&ertau); 
	    av[6][k] = ave;er[6][k]=err;ta[6][k]=tau;erta[6][k]=ertau;

	    av[7][k] = vc2[k]*dno;
	    er[7][k]= sqrt(2*ta[1][k]*dno*(vc22[k]*dno-av[7][k]*av[7][k]));

	    av[8][k] = vce2[k]*dno;
	    er[8][k]=sqrt(2*ta[1][k]*dno*(vce22[k]*dno-av[8][k]*av[8][k]));

	    av[9][k] = vce22[k]*dno;
	    er[9][k]=sqrt(2*ta[1][k]*dno*(av[9][k]*av[9][k]));

	    av[10][k] = vma2[k]*dno;
	    er[10][k]=sqrt(2*ta[6][k]*dno*(vma2[k]*dno-av[10][k]*av[10][k]));

        av[11][k] = vama[k]*dno;
        er[11][k]=sqrt(2*ta[6][k]*dno*(vma2[k]*dno-av[11][k]*av[11][k]));
        
	    av[12][k] = vma22[k]*dno;
	    er[12][k]=sqrt(2*ta[6][k]*dno*(vma22[k]*dno-av[12][k]*av[12][k]));
	    
	     for(ik=0;ik<no;ik++)var[ik] = vxt[ik][k]; 
            FindAutoTau(var,no,&ave,&err,&tau,&ertau); 
	    av[13][k] = ave;er[13][k]=err;ta[13][k]=tau;erta[13][k]=ertau;
	    av[14][k] = vx2[k]*dno;
	    er[14][k] = sqrt(2*ta[14][k]*dno*(vx2[k]*dno-av[14][k]*av[14][k]));
	    

	    }
            if((out_dia= fopen(str_dia,"w"))==NULL)
            {
              printf("non riesco ad aprire il file di output !!\n");
              exit(1);
            }

            if((out_ave= fopen(str_ave,"w"))==NULL)
            {
              printf("non riesco ad aprire il file di output !!\n");
              exit(1);
            }
            if((out_tab= fopen(str_tab,"w"))==NULL)
            {
              printf("non riesco ad aprire il file di output !!\n");
              exit(1);
            }

	    fprintf(out_ave,"#nobs = %d\n",no);
            fprintf(out_ave,"#b  ");
	    fprintf(out_ave,"\n");
	    fprintf(out_dia,"#nobs = %d\n",no);
            fprintf(out_dia,"#b   en   er  tau   er   r2   er   tau    er M er tau er x er tau er \n");

	    fprintf(out_tab,"#nobs = %d\n",no);
            fprintf(out_tab,"#b      c    er     en    er   csp   er    cspen   er    r2   er    e2    er     r2pho   er     r2phi   er    M   er AM er gr  chi  erchi  x er   Chix    ErChiX \n");
	  for(k=0;k<nchain;k++)
	  {
	   cspc = (av[7][k]-av[0][k]*av[0][k])*dn;
	   ercspc = (er[7][k]+2.*er[0][k]*av[0][k])*dn;
	   cspe = (av[8][k]-av[1][k]*av[1][k])*dn;
	   ercspe = (er[8][k]+2.*er[1][k]*av[1][k])*dn;
	   chi = (av[10][k]-av[11][k]*av[11][k])*n;
	   erchi = (er[10][k]+2.*er[11][k]*av[11][k])*n;
	    chix = (av[14][k]-av[13][k]*av[13][k])*dn;
	   erchix = (er[14][k]+2.*er[13][k]*av[13][k])*dn;
	

	   fprintf(out_tab,"%.3lf ",log(bet[k]));
	   fprintf(out_ave,"%.3lf ",log(bet[k]));
	   fprintf(out_dia,"%.3lf ",log(bet[k]));

	   for(i=0;i<2;i++)fprintf(out_tab,"%.2lf %.3lf ",av[i][k],er[i][k]);
	   fprintf(out_tab, "%.4lf %.4lf %.4lf %.4lf ",cspc,ercspc,cspe,ercspe);
	   for(i=2;i<7;i++)fprintf(out_tab, "%.2lf %.3lf ",av[i][k],er[i][k]);
           fprintf(out_tab, "%.2lf %.3lf ",av[11][k],er[11][k]);
           fprintf(out_tab, "%.4lf ",av[12][k]/(av[10][k]*av[10][k])-3);
	   fprintf(out_tab, "%.6lf %.6lf ",chi,erchi);
	    fprintf(out_tab, "%.4lf %.4lf ",av[13][k], er[13][k]);
	     fprintf(out_tab, "%.4lf %.4lf ",av[14][k], er[14][k]);
	
	  
	   for(i=1;i<3;i++)fprintf(out_dia,"%.3lf %.4lf %.3lf %.4lf ",av[i][k],er[i][k],ta[i][k],erta[i][k]);
	   fprintf(out_dia,"%.3lf %.4lf %.3lf %.4lf ",av[6][k],er[6][k],ta[6][k],erta[6][k]);
	    fprintf(out_dia,"%.3lf %.4lf %.3lf %.4lf ",av[13][k],er[13][k],ta[13][k],erta[13][k]);
	   fprintf(out_dia,"\n");
           fprintf(out_ave,"\n");
	   fprintf(out_tab,"\n");
	  }
          fclose(out_tab);
	  fclose(out_ave);
	  fclose(out_dia);

          }

/* End part of the saving */

       if (cycles >= cycle_limit) 
       {

            if((out_ave= fopen(str_ave,"w"))==NULL)
            {
              printf("non riesco ad aprire il file di output !!\n");
              exit(1);
            }
            if((out_tab= fopen(str_tab,"w"))==NULL)
            {
              printf("non riesco ad aprire il file di output !!\n");
              exit(1);
            }

	    fprintf(out_ave,"#nobs = %d\n",no);
        fprintf(out_ave,"#b  ");
	    fprintf(out_ave,"\n");
	    fprintf(out_dia,"#nobs = %d\n",no);
        fprintf(out_dia,"#b     en    er    tau   er    r2   er    tau    er x er tau er \n");

	    fprintf(out_tab,"#nobs = %d\n",no);
        fprintf(out_tab,"#b      c    er     en    er   csp   er    cspen   er    r2   er    e2    er     r2pho   er     r2phi   er    M   er AM er gr  chi  erchi  x er   Chix    ErChiX \n");
	  for(k=0;k<nchain;k++)
	  {
	   cspc = (av[7][k]-av[0][k]*av[0][k])*dn;
	   ercspc = (er[7][k]+2.*er[0][k]*av[0][k])*dn;
	   cspe = (av[8][k]-av[1][k]*av[1][k])*dn;
	   ercspe = (er[8][k]+2.*er[1][k]*av[1][k])*dn;
	   chi = (av[10][k]-av[11][k]*av[11][k])*n;
	   erchi = (er[10][k]+2.*er[11][k]*av[11][k])*n;
		chix = (av[14][k]-av[13][k]*av[13][k])*dn;
	   erchix = (er[14][k]+2.*er[13][k]*av[13][k])*dn;

	   fprintf(out_tab,"%.3lf ",log(bet[k]));
	   fprintf(out_ave,"%.3lf ",log(bet[k]));
	   fprintf(out_dia,"%.3lf ",log(bet[k]));

	   for(i=0;i<2;i++)fprintf(out_tab,"%.2lf %.3lf ",av[i][k],er[i][k]);
	   fprintf(out_tab, "%.4lf %.4lf %.4lf %.4lf ",cspc,ercspc,cspe,ercspe);
	   for(i=2;i<7;i++)fprintf(out_tab, "%.2lf %.3lf ",av[i][k],er[i][k]);
           fprintf(out_tab, "%.2lf %.3lf ",av[11][k],er[11][k]);
          fprintf(out_tab, "%.4lf ",av[12][k]/(av[10][k]*av[10][k])-3);
	   fprintf(out_tab, "%.4lf %.4lf ",chi,erchi);
	   fprintf(out_tab, "%.4lf %.4lf ",av[13][k], er[13][k]);
	   	     fprintf(out_tab, "%.4lf %.4lf ",av[14][k], er[14][k]);

	    for(i=1;i<3;i++)fprintf(out_dia,"%.3lf %.4lf %.3lf %.4lf ",av[i][k],er[i][k],ta[i][k],erta[i][k]);
	   fprintf(out_dia,"%.3lf %.4lf %.3lf %.4lf ",av[6][k],er[6][k],ta[6][k],erta[6][k]);
	  fprintf(out_dia,"\n");
	   fprintf(out_ave,"\n");
	   fprintf(out_tab,"\n");
	   
	   printf("k: ");
	   printf("%d %.4lf %4lf\n", k, log(bet[k]), vfreq[k]/(double)vfreqt[k]);
	   
	   
	  }
          fclose(out_tab);
	  fclose(out_ave);
    //      fclose(out_dia);
          fclose(out_mag_series);


          for(k=0;k<nchain;k++)write_chain(iw, n,sigma,bb[k],k);
          exit(0); 
         }

 	}   
 	
 	
 	  fclose(out_mag_series);
          fclose(out_ene_series);
          fclose(out_cont_series);
          fclose(out_RG_series);
          fclose(out_x_series);
 	
 	
 	
 return(1);
} 



void FindAutoTau(double *var, int no,double *ave,double *err,double *tau,double *ertau)
{
static double *texp,x,y;
static double e,w,z;
static double s,sx,sy,sxy,sxx;

 int i,j,k,l,c,i0;
 double a,t,sa,sb;
 float xx;
i = (int) ((double) no/ 10.) + 1;
texp = (double *)malloc(i * sizeof(double));
if(!texp) { fprintf(stderr,"\n ...Memory allocation fails...\n");
   exit(1); }
/* Calculate the mean of the stochastic variable... */
x = 0; 
for(i=0;i<no;i++) x += var[i]; x = x / (double) no;
/* Calculate the mean of the square stochastic variable... */
y = 0; for(i=0;i<no;i++) y += var[i] * var[i]; y = y / (double) no;
/* Calculate the autocorrelation function now... */
y -= x * x;
e = 0.;s = 0.; sx = 0.; sy = 0.; sxy = 0.; sxx = 0.; i0 = 1;
for(k=0;k<no;k++){
 w = 0.; z = 0.;
 for(l=0;l<no-1-k;l++){ w += 1.0;
   z += 2.0 * (var[l] - x) * (var[l+k+1] - x);}
 z = z / w;
   if(z > 0.) e += z;
   else break;
 a = z / (2.0 * y); /* This is the value of the autocorrelation function.  If we in its
               then we can try to estimate its exponential tail. */

             s   += 1.0;
             sx  += (double) k + 1.0;
             sxx += ((double) k + 1.0) * ((double) k + 1.0);
             sy  += log(z / 2.0);
             sxy += ((double) k + 1.0) * log(z / 2.0);
 if(s > 3.){ sa = (s * sxy - sx * sy) / (s * sxx - sx * sx); sa = -1.0 / sa;
 sb = (sxx * sy - sx * sxy) / (s * sxx - sx * sx);}
 else {sa = 0.; sb = 0.;}
 texp[i0-1] = sa;
 if(s == 30){i0++; s = 0.; sx = 0.; sxx = 0.; sy = 0.; sxy = 0.;}
 t = (y + e) / (2.0 * y);

if(k > (int) (1000.0 * t)){
 printf(" \n\n\n...The window is 1000X larger than the autocorrelation time...\n\n\n");
 break;}
}

t = (y + e) / (2.0 * y);
sb = t * sqrt(2.0 * (2.0 * (double) k + 1.0) / (double) no);

/*  The integrated autocorrelation time is %f +- %f\n",t,sb */
/* Calculate the exponential autocorrelation time now */

/* for(i=0;i<i0;i++){s = 0.; sx = 0.;
                  for(j=i;j<i0;j++){
                     if(texp[j]<0.1)break;
                     if((j>0)&&(texp[j]>3. * texp[j-1])) break;
                       s += texp[j]; sx += 1.;}
                  if(i == 7) printf("\n");
                  if(sx < 0.1) break;
                  s = s / sx;
                  }
*/
y = sqrt( 2.0 * y * t / (double) no);
*ave = x;
*err = y;
*tau = t;
*ertau = sb;
/* fprintf(stderr,"ave = %lf err = %lf tau = %lf ertau = %lf\n",x,y,t,sb); */
free(texp);
return;
}

//This function allows me to make a graph of the polymer
void write_chain(iw, n,sigma,chain,kc)
int ***iw,n,**sigma,chain,kc;
{
fopen(outconf[kc], "w");
 int i,j,k,kd;
  int count;
  double dx,dy,ds,dsx,dsy,x1,x2,y1,y2,dx1,dy1;
  int xmin,xmax,ymin,ymax;
  if((outconf[kc] = fopen(str_conf[kc],"a"))==NULL)
  {
   printf("non riesco ad aprire il file di output numero %d !!\n",i);
   exit(1);
  }
  fprintf(outconf[kc],"%d \n",n);
  fprintf(outconf[kc]," \n");
  for(i=0;i<n;i++)
  {
     if(sigma[i][kc] >0)fprintf(outconf[kc],"1");
     else if(sigma[i][kc]==0)fprintf(outconf[kc],"0");
     else fprintf(outconf[kc],"-1");
     for(k=0;k<3;k++)fprintf(outconf[kc]," %d ", iw[i][k][chain]);
     fprintf(outconf[kc],"\n");
  }

 fprintf(outconf[kc],"\n");
 fclose(outconf[kc]);
 return;

}
/*
void write_chain_ps(iw, n,sigma,chain,kc)
int ***iw,n,**sigma,chain,kc;
{
 int i,j,k,kd;
  int count;
  double dx,dy,ds,dsx,dsy,x1,x2,y1,y2,dx1,dy1;
  int xmin,xmax,ymin,ymax;
  if((outconf[kc] = fopen(str_conf[kc],"a"))==NULL)
  {
   printf("non riesco ad aprire il file di output numero %d !!\n",i);
   exit(1);
  }
  
  fprintf(outconf[kc],"%!\n");
  fprintf(outconf[kc],"1 setlinewidth\n");
  fprintf(outconf[kc],"/cm {28.346456 mul} def\n");
  fprintf(outconf[kc],"5 cm 6 cm translate\n");
  fprintf(outconf[kc],"newpath\n");

  count = 0;
  xmin=iw[0][0][chain];
  xmax=iw[0][0][chain];
  ymin=iw[0][1][chain];
  ymax=iw[0][1][chain];

  for(j=1;j<n;j++)
  {
  if(iw[j][0][chain] < xmin)xmin=iw[j][0][chain];
  if(iw[j][0][chain] > xmax)xmax=iw[j][0][chain];
  if(iw[j][1][chain] < ymin)ymin=iw[j][1][chain];
  if(iw[j][1][chain] > ymax)ymax=iw[j][1][chain];
  }

  dx=xmax-xmin;
  dy=ymax-ymin;
  fprintf(stderr,"xmax=%d  xmin=%d ymax=%d ymin=%d dx=%lf dy=%lf\n",xmax,xmin,ymax,ymin,dx,dy);
  count = 0;
  dsx = 300./dx;
  dsy = 300./dy;

  fprintf(outconf[kc],"%.2lf %.2lf moveto\n",dsx*(iw[0][0][chain]-xmin)+10,dsx*(iw[0][1][chain]-ymin)+10);

   for(j=1;j<n;j++)
  {
  fprintf(outconf[kc],"%.2lf %.2lf lineto\n",dsx*(iw[j][0][chain]-xmin)+10,dsx*(iw[j][1][chain]-ymin)+10);
  }
     
  for(i=0;i<n;i++)
  {
     if(sigma[i][kc] > 0)fprintf(outconf[kc],"C  ");
     else fprintf(outconf[kc],"N  ");
     for(k=0;k<2;k++)fprintf(outconf[kc]," %d ", iw[i][k][chain]);
     fprintf(outconf[kc],"\n");
  }
  fprintf(outconf[kc],"stroke\n");
  fprintf(outconf[kc],"showpage\n");
 fclose(outconf[kc]);
 return;

}*/
