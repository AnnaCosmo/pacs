#include <iostream> // input output
#include <cmath> // (for sqrt)
#include <vector>
#include <tuple>
#include "readParameters.hpp"
#include "GetPot.hpp"
#include "gnuplot-iostream.hpp"// interface with gnuplot
/*!
  @file main.cpp
  @brief Temperature distribution in a 1D bar.

  @detail
    We solve  \f$ -T^{\prime\prime}(x)+act*(T(x)-T_e)=0, 0<x<L \f$ with 
    boundary conditions \f$ T(0)=To; T^\prime(L)=0\f$
    
    **************************************************
    Linear finite elements
    Iterative resolution by Gauss Siedel.
    **************************************************
    
    Example adapted by Luca Formaggia from  a code found in 
    "Simulation numerique an C++" di I. Danaila, F. Hecht e
    O. Pironneau.
*/

/*EXPLANATION:
This is the main of the solution to challenge 1.3. I implemented the 
thomas algorithm using three vector which represents first the explicit 
matrix, then the two matrixes of the LU decomposition (used for the resolution).
It's possible to chose both the type of the output and, in case, the name of 
the output file (in parameters.pot)*/ 

//! helper function
void printHelp()
{
  std::cout<<"USAGE: main [-h] [-v] -p parameterFile (default: parameters.pot)"<<std::endl;
  std::cout<<"-h this help"<<std::endl;
  std::cout<<"-v verbose output"<<std::endl;
}

//! main program
int main(int argc, char** argv)
{
  using namespace std; // avoid std::
  int status(0); // final program status
  GetPot   cl(argc, argv);
  if( cl.search(2, "-h", "--help") )
    {
      printHelp();
      return 0;
    }
  // check if we want verbosity
  bool verbose=cl.search(1,"-v");
  // Get file with parameter values
  string filename = cl.follow("parameters.pot","-p");
  cout<<"Reading parameters from "<<filename<<std::endl;
  // read parameters
  const parameters param=readParameters(filename,verbose);
  // Transfer parameters to local variables
  // I use references to save memory (not really an issue here, it is just
  // to show a possible  use of references)
  // Here I use auto (remember that you need const and & if you want constant references)
  const auto& L= param.L;  // Bar length
  const auto& a1=param.a1; // First longitudinal dimension
  const auto& a2=param.a2; //  Second longitudinal dimension
  const auto& To=param.To; // Dirichlet condition
  const auto& Te=param.Te; // External temperature (Centigrades)
  const auto& k=param.k;  // Thermal conductivity
  const auto& hc=param.hc; // Convection coefficient
  const auto&    M=param.M; // Number of grid elements
  const auto& newfilename=param.newname; //new name for the result file
  const auto& output=param.output; //tipe of output
  
  //! Precomputed coefficient for adimensional form of equation
  const auto act=2.*(a1+a2)*hc*L*L/(k*a1*a2);

  // mesh size
  const auto h=1./M;
  
  // Solution vector
  std::vector<double> theta(M+1);
  //vector used for the thomas algorithm
  std::vector<double> y(M);
  

  theta[0]=(To-Te);
  y[0]=theta[0];
    
  //three vector to implement the tridiagonal matrix
  std::vector<double> a(M,2.+h*h*act);     //main diagonal
  std::vector<double> b(M-1,-1.);   //first diagonal above the main one
  std::vector<double> c(M-1,-1.);   //first diagonal below the main one

  a[M-1]=1.;

//creating coeff. for thomas algorithm
  for(unsigned int m=0;m < M-1;++m){
     c[m]=c[m]/a[m];
     a[m+1]=a[m+1]+c[m];}


 for(unsigned int m=1;m<M;++m){
   y[m]=-c[m-1]*y[m-1];
}

 theta[M]=y[M-1]/a[M-1]; //remember that vector theta has one further component than y


 for(int m=M-1;m>0;--m){
   theta[m]=(y[m-1]-b[m-1]*theta[m+1])/a[m-1];
}




 // Analitic solution

    vector<double> thetaa(M+1);
     for(int m=0;m <= M;m++)
       thetaa[m]=Te+(To-Te)*cosh(sqrt(act)*(1-m*h))/cosh(sqrt(act));

     // writing results with format
     // x_i u_h(x_i) u(x_i) and lauch gnuplot 

     Gnuplot gp;
     std::vector<double> coor(M+1);
     std::vector<double> sol(M+1);
     std::vector<double> exact(M+1);

     
     if(output=="both"){
     cout<<"Result file: "<<newfilename<<endl; //print of the name of the file with result
     ofstream f(newfilename); 
     for(int m = 0; m<= M; m++)
       {
	 // \t writes a tab 
         f<<m*h*L<<"\t"<<Te+theta[m]<<"\t"<<thetaa[m]<<endl;
	 // An example of use of tie and tuples!
         
	 std::tie(coor[m],sol[m],exact[m])=
	   std::make_tuple(m*h*L,Te+theta[m],thetaa[m]);
       }
     // Using temporary files (another nice use of tie)
      gp<<"plot"<<gp.file1d(std::tie(coor,sol))<<
       "w lp title 'uh',"<< gp.file1d(std::tie(coor,exact))<<
      "w l title 'uex'"<<std::endl;

     f.close();
     }
     else if(output=="file"){
     cout<<"Result file: "<<newfilename<<endl; //print of the name of the file with result
     ofstream f(newfilename); 
     for(int m = 0; m<= M; m++)
       {
         // \t writes a tab 
         f<<m*h*L<<"\t"<<Te+theta[m]<<"\t"<<thetaa[m]<<endl;
        }

     f.close();
	 
     }
     else{
     for(int m = 0; m<= M; m++)
       {
	 std::tie(coor[m],sol[m],exact[m])=
	   std::make_tuple(m*h*L,Te+theta[m],thetaa[m]);
       }
     // Using temporary files (another nice use of tie)
      gp<<"plot"<<gp.file1d(std::tie(coor,sol))<<
       "w lp title 'uh',"<< gp.file1d(std::tie(coor,exact))<<
      "w l title 'uex'"<<std::endl;
}
     return status;
}
