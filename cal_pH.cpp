#include <iostream>
#include <cmath>
//#include <math.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <gmp.h>
#include <mpfr.h>
#include "arg.h"
#include "except.h"


using namespace std;

int main(int argc, char* argv[])
{
  int nDis = 2;
  long double a, b, c, d, f, g, h, I, i, J, j, k, r, s, t, u, 
    X1, X2, X3, am, pH, kw, pKa, kA, p=0.3333333333333333333;
  char acid_base = 'A';
  kw = 1.0116;
  pKa = 7;
  am = 0.01;

  if (nDis>0)
    cout << "\nThis program calculates pH for very dilute solutions of a weak acid or base.\n\n";

  // Provide a command line user interface
  static Config_t params[] = {

    "Kw", "Dissociation constant of water (without the 10^-14)",
    "-kw", L_DOUBLE, (caddr_t)&kw,

    "pKa", "pKa",
    "-pka", L_DOUBLE, (caddr_t)&pKa,  

    "conc", "Concentration (molarity, mol/L)",
    "-c", L_DOUBLE, (caddr_t)&am,

    "Acid/Base", "Is this acid ('A') or base ('B')",
    "-ab", CHAR, (caddr_t)&acid_base,

    "dis", "Display options; 0 - most parsimonious; 2 - most verbose",
    "-dis", INT, (caddr_t)&nDis,
	  
    0, 0, 0, NOTYPE, 0

  };
  if ( argc < 2) {
    pusage ( argv[0], params );

    if (nDis>1) {
      cout << "\n\tThe dissociation constant of water at   0 deg C is  0.1153 x 10^-14.";
      cout << "\n\tThe dissociation constant of water at  25 deg C is  1.0116 x 10^-14.";
      cout << "\n\tThe dissociation constant of water at  37 deg C is  2.2418 x 10^-14.";
      cout << "\n\tThe dissociation constant of water at  50 deg C is  5.3088 x 10^-14.";
      cout << "\n\tThe dissociation constant of water at  75 deg C is 19.4089 x 10^-14.";
      cout << "\n\tThe dissociation constant of water at 100 deg C is 54.3250 x 10^-14.\n\n";
    }

    exit (0);
  }
  if ( ppconf (argc, argv, params, true) ) {
    pusage (argv[0], params);
    exit (1);
  }

  mpfr_t gmp_tmp, gmp_a, gmp_b, gmp_c, gmp_d, gmp_f, gmp_g, gmp_h, gmp_I, gmp_i, gmp_J, gmp_j, gmp_k, gmp_r, gmp_s, gmp_t, gmp_u, 
    gmp_X1, gmp_X2, gmp_X3, gmp_am, gmp_pH, gmp_kw, gmp_pKa, gmp_kA; //, p=0.3333333333333333333;

  mpfr_inits2(128, gmp_tmp, gmp_a, gmp_b, gmp_c, gmp_d, gmp_f, gmp_g, gmp_h, gmp_I, gmp_i, gmp_J, gmp_j, gmp_k, gmp_r, gmp_s, gmp_t, gmp_u, 
	      gmp_X1, gmp_X2, gmp_X3, gmp_am, gmp_pH, gmp_kw, gmp_pKa, gmp_kA, NULL);

  mpfr_set_d (gmp_kw, kw, MPFR_RNDN);
  mpfr_mul_d (gmp_kw, gmp_kw, 1e-14, MPFR_RNDN);
  // kw *= 1e-14;

  mpfr_set_d (gmp_tmp, 10, MPFR_RNDN);
  mpfr_set_d (gmp_kA, -pKa, MPFR_RNDN);
  mpfr_pow (gmp_kA, gmp_tmp, gmp_kA, MPFR_RNDN);
  // kA = 1/(pow (10, pKa));

  switch (acid_base) {

  case 'A' :
    a = 1;
    b = kA;
    c = -(kw + kA*am);
    d = -kA*kw;
    break;

  case 'B' :
    a = 1;
    b = am+kA;
    c = - kw;
    d = -kA*kw;
    break;

  default :
    SayBye("Incorrect option\n");
  }
  
  mpfr_printf ("kA %.17Rg\n", gmp_kA);

  // f = ((3*(c/a)) - ((pow(b,2))/(pow(a,2))))/3;  
  f = ( 3*c/a - (b*b)/(a*a) ) / 3; 
  // g = ((2*(pow(b,3))/(pow(a,3)) - (((9*b)*c)/pow(a,2)) + ((27*(d/a))))) /27;
  g = ( (2*b*b*b)/(a*a*a) - (9*b*c)/(a*a) + (27*d/a) ) / 27;
  //h = (((pow(g,2))/4)+((pow(f,3)))/27);
  h = g*g*0.25 + f*f*f/27;
  
  if (h<=0 || f!=0 || g!=0) {
    I = g*g*0.25 - h;
    // i = pow(I,0.5);
    i = sqrt(I);
    //j = pow(i,p);
     j = cbrt(i);
    J = -g/(2*i);
    if (J<-1) {
      J=-1;
      printf("J is %Lf\n", J);
      LookOut("J in acos(J) is out of bounds, resetting J to -1");
    }
    k = acos(J);

    X1 = (2*j*cos(k/3)) - (b/(3*a));
    X2 = -j*( cos(k/3) + sqrt(3)*sin(k/3) ) - (b/(3*a));
    X3 = -j*( cos(k/3) - sqrt(3)*sin(k/3) ) - (b/(3*a));
    if (nDis>0)
      cout << "3 real roots; X1 = " << X1 << ", X2 = " << X2 << ", X3 = " << X3 << "\n";
  }
  else if (h>0){
    f = ((3*(c/a)) - ((pow(b,2))/(pow(a,2))))/3;
    g = ((2*(pow(b,3))/(pow(a,3)) - (((9*b)*c)/pow(a,2)) + ((27*(d/a))))) /27;
    h = (((pow(g,2))/4)+((pow(f,3)))/27);
    r = -g/2 + pow(h,0.5);
    // s = pow(r,p);
    s = cbrt(r);
    t = -g/2 - pow(h,0.5);
    u = pow(t,p);

    X1 = (s+u) - (b/(3*a));
    X2 = -(s+u)/2 -(b/(3*a)) + i*(s-u)*pow(3,0.5)/2;
    X3 = -(s+u)/2 -(b/(3*a)) - i*(s-u)*pow(3,0.5)/2;
    if (nDis>0)
      cout << "1 real root; X1 = " << X1 << ", X2 = " << X2 << ", X3 = " << X3 << "\n";
  }
  else {
    SayBye ("Numeric calculation error");
  }
	   
  pH = -log10(X1);
  cout << "pH = " << pH << "\n";

  return 1;
}
