/* program ising */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
 
void initial(int *N,int *L,double *T,int *mcs,int *nequil,int *tsave,
             double *E,double *M,double w[],int spin[][33]);
void Metropolis(int N,int L,int spin[][33],double *E,double *M,
                double w[],double *accept);
int DeltaE(int x,int y,int L,int spin[][33]);
void data(double E,double M,double cum[]);
void output(double T,int N,int mcs,double cum[],double accept);
void save_config(int N,int L,double T,int spin[][33]);
void read_config(int N,int L,double *T,int spin[][33]);
void correl(double Ce[],double Cm[],double E,double M,
            double esave[],double msave[],int pass,int tsave);
void c_output(double Ce[],double Cm[],double cum[],int mcs,int tsave);
float rnd(int seed);
 
main()
  {
/*    Metropolis algorithm for the Ising model on a square lattice */
      int pass,mcs,nequil,N,L,tsave,spin[33][33];
      double w[17],Ce[21],Cm[21],esave[101],msave[101];
      double E,M,T,accept,cum[11];
      initial(&N,&L,&T,&mcs,&nequil,&tsave,&E,&M,w,spin);
/*    equilibrate system */
      for (pass = 1;pass <= nequil;pass++)
         Metropolis(N,L,spin,&E,&M,w,&accept);
      accept = 0;
/*    accumulate data while updating spins */
      for (pass = 1; pass <= mcs; pass++)
       {
         Metropolis(N,L,spin,&E,&M,w,&accept);
         data(E,M,cum);
         correl(Ce,Cm,E,M,esave,msave,pass,tsave);
       }
      output(T,N,mcs,cum,accept);
      c_output(Ce,Cm,cum,mcs,tsave);
   }
 
void initial(int *N,int *L,double *T,int *mcs,int *nequil,int *tsave,
             double *E,double *M,double w[],int spin[][33])
   {
      int x,y,up,right,dE,seed;
      double sum,dummy;
      printf("linear dimension of lattice = ");
      scanf("%d",L);
      printf("reduced temperature = ");
      scanf("%lf",T);
      *N = (*L)*(*L);
      printf("# MC steps per spin for equilibrium = ");
      scanf("%d",nequil);
      printf("# MC steps per spin = ");
      scanf("%d",mcs);
      printf("random number seed = ");
      scanf("%d",&seed);
/*    seed must not equal 0 to initialize rnd */
      dummy = rnd(seed); 
/*    random initial configuration */
      for (y = 1; y <= *L; y++)
        for (x = 1; x <= *L; x++)
          {
            if (rnd(0) < 0.5) 
              spin[x][y] = 1;
            else  
              spin[x][y] = -1;
            *M = *M + spin[x][y];
          }
/*    compute initial energy */
      for (y = 1; y <= *L; y++)
        {
/*       periodic boundary conditions */
         if (y==(*L)) 
           up = 1;
         else 
           up = y + 1;
         for (x = 1; x <= *L; x++)
           {
/*           conditional expression */
             if (x==(*L)) 
               right = 1; 
             else 
               right = x + 1;
             sum = spin[x][up] + spin[right][y];
             *E = *E - spin[x][y]*sum;
            }
         }
/*    compute Boltzmann probability ratios  */
      for (dE = -8; dE <= 8; dE=dE+4)
         w[dE+8] = exp(-dE/(*T));
      *tsave = 10;
   }
 
 
void Metropolis(int N,int L,int spin[][33],double *E,double *M,
                double w[],double *accept)
   {
/*    one Monte Carlo step per spin  */
      int ispin,x,y,dE;
      for (ispin=1; ispin <= N; ispin++)
         {
/*       random x and y coordinates for trial spin  */
         x = L*rnd(0) + 1;
         y = L*rnd(0) + 1;
         dE = DeltaE(x,y,L,spin);
         if (rnd(0) <= w[dE+8])
           {
            spin[x][y] = -spin[x][y];
            *accept = *accept + 1;
            *M = *M + 2*spin[x][y];
            *E = *E + dE;
           }
          }
   }
 
int DeltaE(int x,int y,int L,int spin[][33])
   {
/*    periodic boundary conditions */
      int dE,left,right,up,down;
      if (x == 1) 
        {
         left = spin[L][y];
         right = spin[2][y];
         }
      else if (x == L) 
         {
         left = spin[L-1][y];
         right = spin[1][y];
         }
      else
         {
         left = spin[x-1][y];
         right = spin[x+1][y];
         }
      if (y == 1) 
         {
         up = spin[x][2];
         down = spin[x][L];
         }
      else if (y == L) 
         {
         up = spin[x][1];
         down = spin[x][L-1];
         }
      else
         {
         up = spin[x][y+1];
         down = spin[x][y-1];
         }
      dE = 2*spin[x][y]*(left + right + up + down);
      return dE;
   }
 
void data(double E,double M,double cum[])
   {
/*    accumulate data after every Monte Carlo step per spin */
      cum[1] = cum[1] + E;
      cum[2] = cum[2] + E*E;
      cum[3] = cum[3] + M;
      cum[4] = cum[4] + M*M;
      cum[5] = cum[5] + fabs(M);
   }
 
void output(double T,int N,int mcs,double cum[],double accept)
   {
/*    average per spin */
      double norm,eave,e2ave,mave,m2ave,abs_mave;
      norm = 1.0/(mcs*N);
      accept = accept*norm;
      eave = cum[1]*norm;
      e2ave = cum[2]*norm;
      mave = cum[3]*norm;
      m2ave = cum[4]*norm;
      abs_mave = cum[5]*norm;
      printf(" \n");
      printf("Temperature = %8.4f\n",T);
      printf("acceptance probability = %8.4f\n",accept);
      printf("mean energy per spin = %8.4f\n", eave);
      printf("mean squared energy per spin = %10.4f\n", e2ave);
      printf("mean magnetization per spin = %8.4f\n", mave);
      printf("mean abs(mag) per spin = %8.4f\n", abs_mave);
      printf("mean squared mag. per spin = %10.4f\n", m2ave);
   }
 
void save_config(int N,int L,double T,int spin[][33])
   {
      int x,y;
      FILE *fp;   /* note the use of the variable type FILE */
      char filename[20];
      printf("name of file for last configuration = ");
      scanf("%s",filename);
      fp = fopen(filename,"w");
      fprintf(fp,"%f\n",T);
      for (y=1; y <= L; y++)
         for (x=1; x <= L; x++)
            fprintf(fp,"%d\n",spin[x][y]);
      fclose(fp);
   }
 
void read_config(int N,int L,double *T,int spin[][33])
   {
      int x,y;
      FILE *fp;   /* note the use of the variable type FILE */
      char filename[20];
      printf("name of file for last configuration = ");
      scanf("%s",filename);
      fp = fopen(filename,"r");
      fscanf(fp,"%f\n",T);
      for (y=1; y <= L; y++)
         for (x=1; x <= L; x++)
            fscanf(fp,"%d\n",spin[x][y]);
      fclose(fp);
   }
 
void correl(double Ce[],double Cm[],double E,double M,
            double esave[],double msave[],int pass,int tsave)
   {
/*    accumulate data for time correlation functions */
/*    save last tsave values of M and E */
      int tdiff,index0,index;
/*    index0 = array index for earliest saved time */
      index0 = ((pass - 1) % tsave) + 1;
      if (pass > tsave) 
/*      compute Ce and Cm after tsave values are saved */
        {
         index = index0;
         for (tdiff=tsave; tdiff >= 1; tdiff--)
           {
            Ce[tdiff] = Ce[tdiff] + E*esave[index];
            Cm[tdiff] = Cm[tdiff] + M*msave[index];
            index = index + 1;
            if(index > tsave) index = 1;
            }
         }
/*    save latest value in array position of earliest value */
      esave[index0] = E;
      msave[index0] = M;
    }
 
void c_output(double Ce[],double Cm[],double cum[],int mcs,int tsave)
   {
/*    compute time correlation functions */
      int tdiff;
      double ebar,e2bar,mbar,m2bar,norm;
      ebar = cum[1]/mcs;
      e2bar = cum[2]/mcs;
      Ce[0] = e2bar - ebar*ebar;
      mbar = cum[3]/mcs;
      m2bar = cum[4]/mcs;
      Cm[0] = m2bar - mbar*mbar;
      norm = 1.0/(mcs-tsave);
      printf("\n    t       Ce(t)  Cm(t)       \n\n");
      for (tdiff = 1;tdiff <= tsave;tdiff++)
        {
/*       define correlation function so that C(0) = 1 and
         C(infinity) = 0   */
         Ce[tdiff] = (Ce[tdiff]*norm - ebar*ebar)/Ce[0];
         Cm[tdiff] = (Cm[tdiff]*norm - mbar*mbar)/Cm[0];
         printf("%6d %10.5f %10.5f  \n",tdiff,Ce[tdiff],Cm[tdiff]);
        }
   }
 
float rnd(int seed)
   {
/*    linear congruential random number generator with shuffling */
/*    based on ran1 of second edition of "Numerical Recipes" */
#define min(a,b)  ((a) < (b) ? (a) : (b))
#define a 16807
#define m 2147483647  /* 2^31 - 1  */ 
#define rm 1.0/m
#define q 127773  /*  m = a*q + p  */
#define p 2836
#define n 32
#define ndiv (1 + (m-1)/n)
#define rmax (1.0 - 1.2e-7)
      static int r[n+1],r0,r1;
      int j,k;
      if (seed != 0)
/*      initialize table of random numbers  */
         {
          r1 = abs(seed);
       for (j = n+8;j>=1;j--)
         {
         k = r1/q;
         r1 = a*(r1-k*q) - p*k;
         if (r1 < 0) r1 = r1 + m;
         if (j < n) r[j] = r1;
           } 
         r0 = r[1];
         } 
/*     beginning when not initializing    */
/*     compute r1 = mod(a*r1,m) without overflows  */
       k = r1/q;
       r1 = a*(r1 - k*q) - p*k;
       if (r1 < 0) r1 = r1 + m;
       j = 1 + r0/ndiv;
       r0 = r[j];
       r[j] = r1;
       return min(rm*r0,rmax);
   }
