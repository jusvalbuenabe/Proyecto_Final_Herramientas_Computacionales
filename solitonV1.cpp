#include <iostream>
#include <vector>
#include <math.h>       /* exp */
#include <cmath>     /* abs */
#include <iomanip>      // std::setprecision

typedef std::vector<std::vector<double> > Matrix; 

const double dx = 0.4;
const double dt = 0.1;
const double Tmax = 2000*dt; //Time steps
const double LX = 130*dx;
const int NT = Tmax/dt + 1;
//const int NNT = 50;
const int NX = LX/dx + 1;
//const int NNX = 50;
const double mu = 0.1;      //Mu from KdeV equation
const double eps = 0.2;     //Epsilon from KdeV eq

void get_memory(Matrix & malla);
void print(const Matrix & u);
void initial_conditions(Matrix & u);
void propagate(Matrix & u);
void init_gnuplot_contour(void);
void init_gnuplot_3D(void);

int main (void)
{
  Matrix grid;
  get_memory(grid);
  initial_conditions(grid);
  init_gnuplot_3D();
  print(grid);
  propagate(grid);
  print(grid);
  init_gnuplot_contour();
  print(grid);
  return 0;
}
void get_memory(Matrix & u)
{
  int ii = 0;
  u.resize(NX);
  for (ii = 0; ii < NX; ++ii) {
    u[ii].resize(NT);
  }
}


void init_gnuplot_contour(void)
{
  std::cout << "reset" << std::endl;
  std::cout << "unset surface" << std::endl;
  std::cout << "set xlabel \" x \" " << std::endl;
  std::cout << "set ylabel \" t \" " << std::endl;
  std::cout << "set title \" Soliton \" " << std::endl;
  // std::cout << "set terminal gif animate" << std::endl;
  //std::cout << "set out 'miprimeraanimacion.gif'" << std::endl;
  std::cout << "set contour base" << std::endl;
  std::cout << "set view map" << std::endl;
  std::cout << "set key at 0,0,0" << std::endl;
  std::cout << "set cntrparam level incremental 0.0, 0.2, 1.2" << std::endl;
  std::cout << "set pm3d" << std::endl;
}

void init_gnuplot_3D(void)
{
  std::cout << "reset" << std::endl;
  std::cout << "set contour base" << std::endl;
  std::cout << "set xlabel \" x \" " << std::endl;
  std::cout << "set ylabel \" t \" " << std::endl;
  std::cout << "set title \" Soliton \" " << std::endl;
  //std::cout << "set pm3d" << std::endl;
  std::cout << "set pm3d depthorder hidden3d 1" << std::endl;
  std::cout << "set hidden3d" << std::endl;
  std::cout << "set style fill transparent solid 0.65" << std::endl;
  std::cout << "set palette rgb 9,9,3" << std::endl;
  //std::cout << "unset colorbox" << std::endl;
  //std::cout << "set xtics 0, 10, 50" << std::endl;
  //std::cout << "set ytics 0, 0.5, 2" << std::endl;
  //std::cout << "set isosamples 5" << std::endl;
  std::cout << "unset key" << std::endl;
  //std::cout << "unset border" << std::endl;
  // std::cout << "unset tics" << std::endl;
  std::cout << "set ticslevel 0" << std::endl;
}

void print(const Matrix & u)
{
  std::cout << "splot '-' w l lw 2 " << std::endl;
  double x, y;
  int iip,jjp;
  int NNT=NX/3; //Cada Cuantas Muestras en T
  int NNX=NX; //Cada Cuantas Muestras en X
  for (int ii = 0; ii < NNX; ++ii){
    iip=ii*NX/NNX;
    x = iip*dx;
    for (int jj = 0; jj < NNT; ++jj){
      jjp=jj*NT/NNT;
      y = jjp*dt;
      std::cout << x << "  " << y << "  " <<  u[iip][jjp] << std::endl;
    }
    std::cout << std::endl;
  }
  std::cout << "e" << std::endl;
  std::cout << "pause mouse" << std::endl;
}

void initial_conditions(Matrix & u)
{
  int ii, jj, kk;
  double a1, a2, a3, fac, time; 
  
  //Initial conditions t=0
  // Initial wave form
  jj = 0;
  for ( ii = 0;  ii < NX;  ii++ ){
    u[ii][jj] = 0.5*(1.-((exp(2.*(0.2*dx*ii - 5.))-1) /(exp(2.*(0.2*dx*ii - 5.)) + 1)));
  }

  // End points
  u[0][1]   = 1.; 
  u[0][2]   = 1.;
  u[NX-1][1] = 0.;   
  u[NX-1][2] = 0.; 
  fac = mu*dt/(dx*dx*dx);              
  time = dt;
  
  // First time step
  for ( ii=1;  ii < NX-1;  ii++){
    a1 = eps*dt*(u[ii + 1][0] + u[ii][0] + u[ii-1][0]) / (dx*6.);     
    if (ii>1 && ii < NX-2){
      a2 = u[ii + 2][0] + 2.*u[ii-1][0] - 2.*u[ii + 1][0]-u[ii-2][0];
    }
    else a2 = u[ii-1][0] - u[ii + 1][0]; 
    a3 = u[ii + 1][0]-u[ii-1][0]; 
    u[ii][1] = u[ii][0] - a1*a3 - fac*a2/2.;        
  } 
  
  // Other time steps
  for ( ii = 1;  ii < NX-2;  ii++ )  {
    a1 = eps*dt*(u[ii + 1][1] + u[ii][1] + u[ii-1][1]) / (3.*dx); 
    if (ii>1 && ii < NX-3) a2 = u[ii+2][1] + 2.*u[ii-1][1] - 2.*u[ii+1][1] - u[ii-2][1];
    else a2 = u[ii-1][1] - u[ii + 1][1];  
    a3     = u[ii + 1][1] - u[ii-1][1]; 
    u[ii][2] = u[ii][0] - a1*a3 - 2.*fac*a2/3.;      
  } 
} 


void propagate(Matrix & u)
{
  int ii, jj, kk;
  double uOld, uNew, Dif;
  double a1, a2, a3, fac, time; 
  
  // Other time steps
  for ( jj=2;  jj < NT-1;  jj++ )  {                    
    // time += dt; 
    for ( ii = 1;  ii < NX-2;  ii++ )  {
      a1 = eps*dt*(u[ii + 1][jj] + u[ii][1] + u[ii-1][jj]) / (3.*dx); 
      if (ii>1 && ii < NX-2)
	a2 = u[ii+2][jj] + 2.*u[ii-1][jj] - 2.*u[ii+1][jj] - u[ii-2][jj];
      else
	a2 = u[ii-1][jj] - u[ii + 1][jj];  
      a3 = u[ii + 1][jj] - u[ii-1][jj]; 
      u[ii][jj+1] = u[ii][jj-1] - a1*a3 - 2.*fac*a2/3.;      
    }
  }
}

