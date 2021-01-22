#include <Rcpp.h>
#include <math.h> 
#include <iostream>
using namespace Rcpp;
using namespace std;


//' Calculate dihedral angle formed by four atoms
//' 
//' @description
//' For Cartesian coordinates of atoms A-B-C-D, calculate the dihedral angle formed by viewing down the B-C axis.
//'
//' @param a length 3 vector of coordinates of atom A
//' @param b length 3 vector of coordinates of atom B
//' @param c length 3 vector of coordinates of atom C
//' @param d length 3 vector of coordinates of atom D
//' 
//' @return
//' Returns the dihedral angle (in degrees between -180 and 180).
//' 
//' @details
//' Similar to \link[bio3d]{torsion.xyz}, but with implementation in C++.
//'
//' @examples
//' torsion(c(50.051, 37.144, -4.723), c(50.044, 36.248, -3.559),
//'         c(51.296, 35.369, -3.476), c(51.930,35.119,-4.618))
//' 
//' @export
//' 
// [[Rcpp::export]]
double torsion(NumericVector a, NumericVector b, NumericVector c, NumericVector d)
{
  double xij,yij,zij,
  xkj,ykj,zkj,
  xkl,ykl,zkl,
  dxi,dyi,dzi,
  gxi,gyi,gzi,
  bi,bk,ct,
  //boi2,boj2,
  z1,z2,ap,s;
  //bioj,bjoi;

  
  
  // Calculate the vectors C,B,C
  xij =a[0]-b[0]; 
  yij =a[1]-b[1];
  zij =a[2]-b[2];
  xkj =c[0]-b[0];
  ykj =c[1]-b[1];
  zkj =c[2]-b[2];
  xkl =c[0]-d[0];
  ykl =c[1]-d[1];
  zkl =c[2]-d[2];
  
  /* Calculate the normals to the two planes n1 and n2
  this is given as the cross products:
  AB x BC
  --------- = n1
  |AB x BC|
  
  BC x CD
  --------- = n2
  |BC x CD|
  */
  dxi = yij * zkj - zij * ykj;     /* Normal to plane 1                */
  dyi = zij * xkj - xij * zkj;
  dzi = xij * ykj - yij * xkj;
  gxi = zkj * ykl - ykj * zkl;     /* Mormal to plane 2                */
  gyi = xkj * zkl - zkj * xkl;
  gzi = ykj * xkl - xkj * ykl;
  
  /* Calculate the length of the two normals                           */
  bi = dxi * dxi + dyi * dyi + dzi * dzi;
  bk = gxi * gxi + gyi * gyi + gzi * gzi;
  ct = dxi * gxi + dyi * gyi + dzi * gzi;
  
  //boi2 = 1./bi;
  //boj2 = 1./bk;
  bi   = (double)sqrt((double)bi);
  bk   = (double)sqrt((double)bk);
  
  z1   = 1./bi;
  z2   = 1./bk;
  //bioj = bi * z2;
  //bjoi = bk * z1;
  ct   = ct * z1 * z2;
  if (ct >  1.0)   ct = 1.0;
  if (ct < (-1.0)) ct = -1.0;
  ap   = acos(ct);
  
  s = xkj * (dzi * gyi - dyi * gzi)
    + ykj * (dxi * gzi - dzi * gxi)
    + zkj * (dyi * gxi - dxi * gyi);
    
    if (s < 0.0) ap = -ap;
    
    ap = (ap > 0.0) ? PI-ap : -(PI+ap);
    
    return((ap/PI)*180);
}









