/* From: "COMPUTATIONAL PHYSICS, 2nd Ed" 
   by RH Landau, MJ Paez, and CC Bordeianu 
   Copyright Wiley-VCH, 2007.
   Electronic Materials copyright: R Landau, Oregon State Univ, 2007;
   MJ Paez, Univ Antioquia, 2007; and CC Bordeianu, Univ Bucharest, 2007.
   Support by National Science Foundation                              
*/
// Soliton.java: Solves Kortewg-deVries Equation 

import java.io.*; 

public class Soliton {
  static double ds = 0.4;                                  // Delta x
  static double dt = 0.1;                                  // Delta t
  static int max = 2000;                                // Time steps
  static double mu = 0.1;                    // Mu from KdeV equation
  static double eps = 0.2;                   // Epsilon from KdeV eq

  public static void main(String[] argv) 
                         throws IOException, FileNotFoundException {
													 
    int i, j, k; 
    double a1, a2, a3, fac, time; 
    double u[][] = new double[131][3]; 
		
    PrintWriter w =                      // Save data in soliton.dat
          new PrintWriter(new FileOutputStream("soliton.dat"), true);
    for ( i=0;  i < 131;  i++ )                 // Initial wave form
      { u[i][0] = 0.5*(1.-((Math.exp(2*(0.2*ds*i - 5.))-1)/(Math.exp(2*(0.2*ds*i - 5.)) + 1))); }
    u[0][1]   = 1.;                                    // End points
    u[0][2]   = 1.;    
    u[130][1] = 0.;   
    u[130][2] = 0.; 
    fac = mu*dt/(ds*ds*ds);              
    time = dt; 
		                                               // First time step
    for ( i=1;  i < 130;  i++ ) {                      
      a1 = eps*dt*(u[i + 1][0] + u[i][0] + u[i-1][0]) / (ds*6.);     
      if ((i>1) && (i < 129))
        {a2 = u[i + 2][0] + 2.*u[i-1][0] - 2.*u[i + 1][0]-u[i-2][0];}
      else a2 = u[i-1][0] - u[i + 1][0]; 
      a3 = u[i + 1][0]-u[i-1][0]; 
      u[i][1] = u[i][0] - a1*a3 - fac*a2/3.;        
    }   
		                                              // Other time steps
    for ( j=1;  j < max;  j++ )  {                    
      time += dt; 
      for ( i=1;  i < 130;  i++ )  {
         a1 = eps*dt*(u[i + 1][1] + u[i][1] + u[i-1][1]) / (3.*ds); 
         if (i>1 && i < 129) a2 = u[i+2][1] + 2.*u[i-1][1] - 2.*u[i+1][1] - u[i-2][1];
         else a2 = u[i-1][1] - u[i + 1][1];  
         a3     = u[i + 1][1] - u[i-1][1]; 
         u[i][2] = u[i][0] - a1*a3 - 2.*fac*a2/3.;      
      }
      
      for (k=0; k < 131; k++) {u[k][0] = u[k][1]; u[k][1] = u[k][2];}
			                            // gnuplot format, every 200th step
      if ((j%200)==0) {              
        for ( k=0;  k < 131;  k += 2)  w.println("" + u[k][2] + "");  
        w.println( "");                     // Empty line for gnuplot
      }
    } 
    System.out.println("data stored in soliton.dat"); 
  }
} 
