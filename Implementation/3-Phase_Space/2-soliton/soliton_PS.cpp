#include <iostream>
#include <vector>
#include <math.h>       /* exp */
#include <cmath>     /* abs */
#include <iomanip>      // std::setprecision

typedef std::vector<std::vector<double> > Matrix; 

const double dx = 0.4;
const double dt = 0.1;
const double Tmax = 2000.*dt; //Time steps
const double LX = 125.*dx;
const int NT = Tmax/dt + 1;
//const int NNT = 50;
const int NX = LX/dx + 1;
//const int NNX = 50;
const double mu = 0.1;      //Mu from KdeV equation
const double eps = 0.2;     //Epsilon from KdeV eq
const double  fac = mu*dt/(dx*dx*dx);              

void get_memory(Matrix & malla);
void initial_conditions(Matrix & u, Matrix & v);
void End_Points(Matrix & u, Matrix & v);
void First_Time_Step(Matrix & u);
void propagate(Matrix & u, Matrix & v, int & jj);
void print_gif(const Matrix & u, int & jj);
void print_gif_PS(const Matrix & u, const Matrix & v, int & jj);
void init_gnuplot(void);
void init_gnuplot_PS(void);

int main (void)
{
  Matrix grid;
  Matrix Vel;
  get_memory(grid);
  get_memory(Vel);
  initial_conditions(grid,Vel);
  End_Points(grid, Vel);
  First_Time_Step(grid);
  init_gnuplot_PS();
  //init_gnuplot();
  int jj = 0;
  print_gif_PS(grid,Vel, jj);
  //print_gif_PS(Vel, jj);
  for ( jj =1; jj<NT; jj++){
    propagate(grid, Vel, jj);
    if ( jj % 80 == 0 ){
      print_gif_PS(grid,Vel, jj);
      //print_gif(Vel, jj);
    }
  }
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

void initial_conditions(Matrix & u, Matrix & v)
{
  int ii, jj;
  //Initial conditions t=0
  // Initial wave form
  jj = 0;
  for ( ii = 0;  ii < NX;  ii++ ){
    //u[ii][jj] = 0.5*(1.-(tanh(0.2*dx*ii - 5.)));
    v[ii][jj] = 0.0;
    u[ii][jj] = 0.8 * ( 1. - ( tanh( (3./12.) * dx*ii - 3.))*( tanh( (3./12.)*dx*ii - 3.) )) +  0.3 * ( 1. - ( tanh( (4.5/26.) * dx*ii-4.5))*( tanh( (4.5/26.)*dx*ii - 4.5) ));
  }
  
  for (int ii = 0; ii<NX; ii++){
    for(int jj = 1; jj<NT; jj++){
      u[ii][jj] = 0.0;
      v[ii][jj] = 0.0;
    }
  }
}

void End_Points(Matrix & u, Matrix & v)
{
  // End points
  //double a1, a2, a3; 
  for(int jj = 1; jj < NT; jj++){
    //u[0][jj]  = 1.;
    v[0][jj] = 0.0;
    u[0][jj]   = 0.8 * ( 1. - ( tanh(- 3.))*( tanh( - 3.) )) +  0.3 * ( 1. - ( tanh(-4.5))*( tanh( - 4.5) ));
    //u[NX-1][jj] = 0.;
    v[NX-1][jj] = 0.0;
    u[NX-1][jj] = 0.8 * ( 1. - ( tanh( (3./12.) * dx*(NX-1) - 3.))*( tanh( (3./12.)*dx*(NX-1) - 3.) )) +  0.3 * ( 1. - ( tanh( (4.5/26.) * dx*(NX-1)-4.5))*( tanh( (4.5/26.)*dx*(NX-1) - 4.5) ));
  }
}

void First_Time_Step(Matrix & u)
{ 
  // First time step: jj = 1 
  double a1, a2, a3; 
  int ii;
  //int jj = 1;
  for ( ii=1;  ii < NX-1;  ii++){
    a1 = eps*dt*(u[ii + 1][0] + u[ii][0] + u[ii-1][0]) / (dx*6.);     
    if (ii>1 && ii < NX-2){
      a2 = u[ii + 2][0] + 2.*u[ii-1][0] - 2.*u[ii + 1][0]-u[ii-2][0];
    }
    else a2 = u[ii-1][0] - u[ii + 1][0]; 
    a3 = u[ii + 1][0]-u[ii-1][0]; 
    u[ii][1] = u[ii][0] - a1*a3 - fac*a2/2.;        
  } 
}
 

void propagate(Matrix & u, Matrix & v, int & jj)
{
  // Other time steps 
  int ii;
  double a1, a2, a3;
  for ( ii = 1;  ii < NX-1;  ii++ )  {
    a1 = eps*dt*(u[ii + 1][jj] + u[ii][jj] + u[ii-1][jj]) / (3.*dx); 
    if (ii>1 && ii < NX-2)
      a2 = u[ii+2][jj] + 2.*u[ii-1][jj] - 2.*u[ii+1][jj] - u[ii-2][jj];
    else
      a2 = u[ii-1][jj] - u[ii + 1][jj];  
    a3 = u[ii + 1][jj] - u[ii-1][jj]; 
    u[ii][jj+1] = u[ii][jj-1] - a1*a3 - fac*a2;
    v[ii][jj]=(u[ii][jj+1]-u[ii][jj-1])/(2*dt);
  }
}

void init_gnuplot(void)
{
  std::cout << "set terminal gif animate" << std::endl;
  std::cout << "set out 'SolitonVel.gif'" << std::endl;
  std::cout << "unset key" << std::endl;
  std::cout << "set yrange [-0.02:0.02]" << std::endl;
  //std::cout << "set contour base" << std::endl;
  //std::cout << "set pm3d" << std::endl;
}

void print_gif(const Matrix & u, int & jj)
{
  std::cout << "plot '-' w l  " << std::endl;
  double x;
  for (int ii = 0; ii < NX; ++ii){
    x = ii*dx;
    std::cout << x << "  " <<  u[ii][jj] << std::endl;
    //std::cout << std::endl;
  }
  std::cout << "e" << std::endl;
  //std::cout << "pause mouse" << std::endl;
}


void init_gnuplot_PS(void)
{
  std::cout << "set terminal gif animate" << std::endl;
  std::cout << "set out 'SolitonPS.gif'" << std::endl;
  std::cout << "unset key" << std::endl;
  std::cout << "set yrange [-0.02:0.02]" << std::endl;
  std::cout << "set xrange [-1:1]" << std::endl;
  //std::cout << "set contour base" << std::endl;
  //std::cout << "set pm3d" << std::endl;
}

void print_gif_PS(const Matrix & u, const Matrix & v, int & jj)
{
  std::cout << "plot '-' w l  " << std::endl;
  // double x;
  for (int ii = 0; ii < NX; ++ii){
    //x = ii*dx;
    std::cout << u[ii][jj] << "  " <<  v[ii][jj] << std::endl;
    //std::cout << std::endl;
  }
  std::cout << "e" << std::endl;
  //std::cout << "pause mouse" << std::endl;
}
