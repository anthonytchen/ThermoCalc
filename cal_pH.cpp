#include <iostream>
#include <math.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "arg.h"
#include "except.h"

using namespace std;

int main(int argc, char* argv[])
{
  long double a, b, c, d, f, g, h, I, i, J, j, k, r, s, t, u, 
    X1, X2, X3, am, pH, kw, pKa, kA, p=0.3333333333333333333;
  char acid_base = 'A';
  kw = 1.0116;
  pKa = 7;
  am = 0.01;

  cout << "\nThis program calculates pH for very dilute solutions of a weak acid or base.\n\n";
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
    c = -(kw + (kA*am));
    d = -kA*kw;

    f = ((3*(c/a)) - ((pow(b,2))/(pow(a,2))))/3;
    g = ((2*(pow(b,3))/(pow(a,3)) - (((9*b)*c)/pow(a,2)) + ((27*(d/a))))) /27;
    h = (((pow(g,2))/4)+((pow(f,3)))/27);
    I = ((pow(g,2))/4 - h);
    i = pow(I,0.5);
    j = pow(i,p);
    J = -g/(2*i);
    k = acos(J);

    X1 = (2*j*cos(k/3)) - (b/(3*a));
    X2 = -j*( cos(k/3) + pow(3,0.5)*sin(k/3) ) - (b/(3*a));
    X3 = -j*( cos(k/3) - pow(3,0.5)*sin(k/3) ) - (b/(3*a));
    cout << "X1 = " << X1 << ", X2 = " << X2 << ", X3 = " << X3 << "\n";

    break;

  case 'B' :
    a = 1;
    b = kw/kA;
    c = kw*am/kA - kw;
    d = -kw*kw/kA;

    f = ((3*(c/a)) - ((pow(b,2))/(pow(a,2))))/3;
    g = ((2*(pow(b,3))/(pow(a,3)) - (((9*b)*c)/pow(a,2)) + ((27*(d/a))))) /27;
    h = (((pow(g,2))/4)+((pow(f,3)))/27);
    assert ( h > 0 );

    r = -g/2 + pow(h,0.5);
    s = pow(r,p);
    t = -g/2 - pow(h,0.5);
    u = pow(t,p);

    X1 = (s+u) - (b/(3*a));
    X2 = -(s+u)/2 -(b/(3*a)) + i*(s-u)*pow(3,0.5)/2;
    X3 = -(s+u)/2 -(b/(3*a)) - i*(s-u)*pow(3,0.5)/2;
    cout << "X1 = " << X1 << ", X2 = " << X2 << ", X3 = " << X3 << "\n";

    break;

  default :
    SayBye("Incorrect option");
  }


	   
  pH = -log10(X1);
  cout << "pH = " << pH << "\n";

  return 1;
}
