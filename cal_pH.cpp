#include <iostream>
#include <math.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "arg.h"
#include "except.h"
#include <gsl/gsl_poly.h>

using namespace std;

int main(int argc, char* argv[])
{
  long double a, b, c, d, f, g, h, I, i, J, j, k, r, s, t, u, 
    X1, X2, X3, am, pH, kw, pKa, kA, p=0.3333333333333333333;
  char acid_base = 'A';
  kw = 1.0116;
  pKa = 7;
  am = 0.01;

  cout << "\nThis program calculates pH for very dilute solutions of a weak acid or base.\n";
  cout << "This version used GSL cubic function solver, which only accepts double precision\n";
  cout << "and thus may not be precise enough for low concentration calculation involved.\n";
  cout << "Check out the branch NoGSL_Cubic where cubic function is solved by using a manually\n";
  cout << "implemented trigonometric method with long double precision.\n\n";
  cout << "\n\tThe dissociation constant of water at   0 deg C is  0.1153 x 10^-14.";
  cout << "\n\tThe dissociation constant of water at  25 deg C is  1.0116 x 10^-14.";
  cout << "\n\tThe dissociation constant of water at  37 deg C is  2.2418 x 10^-14.";
  cout << "\n\tThe dissociation constant of water at  50 deg C is  5.3088 x 10^-14.";
  cout << "\n\tThe dissociation constant of water at  75 deg C is 19.4089 x 10^-14.";
  cout << "\n\tThe dissociation constant of water at 100 deg C is 54.3250 x 10^-14.\n\n";

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
	  
    0, 0, 0, NOTYPE, 0

  };
  if ( argc < 2) {
    pusage ( argv[0], params );
    exit (0);
  }
  if ( ppconf (argc, argv, params, true) ) {
    pusage (argv[0], params);
    exit (1);
  }
  
  kw *= 1e-14;

  kA = 1/(pow (10, pKa));

  switch (acid_base) {

  case 'A' :
    a = 1;
    b = kA;
    c = -(kw + kA*am);
    d = -kA*kw;
    break;

  case 'B' :
    a = 1;
    b = am + kA;
    c = -kw;
    d = -kA*kw;

    break;

  default :
    SayBye("Incorrect option");
  }
  
  double x0, x1, x2;
  int n_roots;
  n_roots = gsl_poly_solve_cubic (b, c, d, &x0, &x1, &x2);
  
  cout << n_roots << " roots found" << "\n";
  if (n_roots == 1) {
    cout << "X0 = " << x0 << "\t";
    pH = -log10(x0);
    cout << "pH = " << pH << "\n";
  }
  else if (n_roots == 3) {
    cout << "X0 = " << x0 << ", X1 = " << x1 << ", X2 = " << x2 << "\n";
    pH = -log10(x2);
    cout << "pH = " << pH << "\n";
  }

  return 1;
}
