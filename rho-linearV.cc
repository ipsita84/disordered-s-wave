#include <iostream>
#include <fstream>
#include <armadillo>


const double THRES = 1.0E-10;

using namespace std;
using namespace arma;


//inline
//vector<double> realvec(const vector<double> & V) {
//vector<double> res; for (complex<double> val : V)
//res.push_back(val.real()); return res; }


int main (int argc, const char * argv[])
{
 
    const cx_double ii(0,1);
    const double Pi = math::pi();
    double kky=0, kkx=0;
    double rrx=1, rry=1 ;
    const double delta = 0.1 , t=1, V=-0.01, mu=0.1 ;
    double N  = 20;//lattice size 

    ofstream fout(string("result3.dat").c_str());// Opens a file for output


    for ( double rx = rrx; rx < 3*N+ THRES; rx += 1)
    {
        for ( double ry = rry; ry <3*N+THRES; ry += 1)
        {   cx_double rhosum(0);
            for (double k1y = kky; k1y <2*Pi+THRES; k1y += 2*Pi/(N))
             {  for (double k1x = kkx; k1x < 2*Pi+THRES; k1x += 2*Pi/(N))
		 { double e1(0);
		   e1 = -2*t* ( cos(k1x)+ cos(k1y) ) -mu;
		 
		   { for(double k2y = kky; k2y < 2*Pi+THRES; k2y += 2*Pi/(N))
                     { for(double k2x = kkx; k2x < 2*Pi+THRES ; k2x += 2*Pi/(N))
                         { double e2(0);
		           e2 = -2*t* ( cos(k2x)+ cos(k2y) ) -mu;
		           rhosum += 2*Pi*V*(-1+ e1*e2/sqrt((e1*e1+delta*delta)*(e2*e2    +delta*delta)))*exp(1.0*(k2x - k1x)*rx*ii)*exp(1.0*(k2y - k1y)*ry*ii)
                                     / ( sqrt(e1*e1+delta*delta) + sqrt(e2*e2+delta*delta) );
            
                          }
                      }
		    }
                  }
	      }
       
        fout << rx << '\t' << ry << '\t' << rhosum.real() << endl;
        //if (!rho.imag()) { rho has imaginary part ! };    
        }
    }

   
    fout.close();

    return 0;
}

