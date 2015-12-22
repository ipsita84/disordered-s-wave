#include <iostream>
#include <fstream>
#include <armadillo>
#include "gcd.h"
#include "io.h"
#include "hamil.h"

const double r = 50.0;
const unsigned int q_max = 11;
const double THRES = 1.0E-10;

using namespace std;
using namespace arma;

int main (int argc, const char * argv[])
{
    vector <double> rxvec, ryvec, rho;
    const double Pi = math::pi();
    double rrx=0, rry=0;
    const double delta = 0.1 , t=1, V=0.01 ;


    unsigned int N  = 60;//lattice size 

    for (double rx = rrx; rx <= N; rx += 1)
    {
        rxvec.push_back( double(rx) );
    }

    for (double ry = rry; ry <= N; ry += 1)
    {
        ryvec.push_back( double(ry) );
    }

    for (double rrx : rxvec)
    {
        for (double rry : ryvec)
        {   double rhosum(0);
            for (double k1y = kky; ky <=2*Pi; k1y += 2*Pi/N)
             {  for (double k1x = kkx; kx <=2*Pi; ky += 2*Pi/N)
		 double e1(0);
		 e1 = -2*t* ( cos(k1x)+ cos(k1y) ) -mu;
		 
		  { for(double k2y = kky; ky <=2*Pi; k2y += 2*Pi/N)
                     { for(double k2x = kkx; ky <=2*Pi; k2x += 2*Pi/N)
                         { double e2(0);
		           e2 = -2*t* ( cos(k2y)+ cos(k2x) ) -mu;
		           rhosum += 2*Pi*V*(-1+ e1*e2/sqrt((e1^2+delta^2)*(e2^2+delta^2)))
                                     * exp(1.0*(k2x - k1x)*rrx*ii)** exp(1.0*(k2y - 
                                                                         k1y)*rry*ii)
                                     / ( sqrt(e1^2+delta^2) + sqrt(e2^2+delta^2) );
            
                          rho.push_back( double(rhosum) );
                          }
                      }
		   }
	      }
            
        }
    }
    
    if( write_to_file("Results/result-rx-ry-rho.dat", rxvec, ryvec, rho ) )
    {
        cout << "Error in opening file to write, quitting..." << endl;
        return 1;
    }

    return 0;
}

