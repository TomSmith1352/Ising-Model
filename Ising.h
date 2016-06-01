//  Ising Model on a square lattice using Metropolis algorithm

# ifndef _ISING_H
# define _ISING_H   // make sure file is not included twice

# include <iostream>
# include <stdlib.h>
# include <math.h>

const double J = 1;      // strength of exchange interaction
const double u0 = 1;     // magnetic moment
const int up = 1;
const int down = -1;


//  data structure encapsulating the essential data of the Ising model, and all 
//  relevant operations. 
class Ising
{ 
   private:

      int * sites;         // pointer to an integer array to be allocated dynamically
      int L;               // linear dimension of lattice
      double T, M, E, B;   // temperature, magnetization, energy, external magnetic field
   
      int & spin(int index1, int index2); // retrieve spin of particular site
      void getE(void);              // calculate overall energy
      bool flip(void);             // attempt to flip a random spin

   public:

      Ising(int l=10, double t=1.0, double b=0);   // constructor with default arguments
      ~Ising()  {  if(sites) delete [] sites; }    // destructor deallocating memory

      void onestep(void);           // one Monte Carlo step per spin
      void report(ofstream & output);    // print M and E to output file stream 
      
};


// To avoid the fuss of creating another file, Ising.cpp, we included the implementation here,
// using keyword inline.

Ising::Ising(int l, double t, double b)
{ 
   L = l;
   T = t;
   B = b;
   M = 0;

   int temp = L*L;
   sites = new int[temp];
   for(int i=0; i<temp; i++) 
      {
         if (drand48()>0.5)    // this linear congruential random number generator is in <stdlib.h>
            {  sites[i] = up;
               M++;
            }
         else
            {  sites[i] = down;
               M--;
             }
       } 
    getE();
}

// indices range from -L+1 to L-1
inline int & Ising::spin(int index1, int index2)
{  
   if(index1<0) index1 += L;      // these operations guarantee periodic boundary conditions
   if(index1>=L) index1 -= L;
   if(index2<0) index2 += L;
   if(index2>=L) index2 -= L;

   int index = index1*L + index2;
   return sites[index];
}

inline void Ising::getE(void)
{  int nn_count = 0;    // nearest neighbour count

   for(int i=0; i<L; i++)
      for(int j=0; j<L; j++)
        {  if(spin(i,j)==spin(i+1, j))     // count from each spin to the right and above 
              nn_count++;
           else 
              nn_count--;
           if(spin(i, j)==spin(i, j+1))
              nn_count++;
           else
              nn_count--;
        }
    
    E = -J*nn_count - u0*B*M;
}


//  attempt to flip a random spin
inline bool Ising::flip(void)
{
    int index1 = int(drand48()*L);  // truncated
    int index2 = int(drand48()*L);  

    int delta_M = -2*spin(index1, index2);      
    int delta_nn_count = spin(index1-1, index2) + spin(index1+1, index2)
                         + spin(index1, index2-1) + spin(index1, index2+1);

    delta_nn_count = -2*spin(index1, index2)*delta_nn_count;

    double delta_E = -J*delta_nn_count -u0*B*delta_M;

    if (delta_E<=0 || drand48()<exp(-delta_E/T))
        {  spin(index1, index2) *= -1;     // spin(int, int) must return int & to modify its content
           M += delta_M;
           E += delta_E;
           return 1;
        }
    return 0;
}

inline void Ising::onestep(void)
{   int size = L*L;

    for(int i=0; i<size; i++)
       flip();
}

void Ising::report(ofstream & output)
{   output<<"\n";
    output.width(15);
    output<<M;
    output.width(15);
    output<<E;
}


# endif
