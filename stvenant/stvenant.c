#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>





#define _M 2  // nbre de variables conservatives

/* gcc stvenant.c -lm */



void riemann(double *wL, double *wR, double xi, double *w);


int main(void){

  /* double hL=1; */
  /* double uL=1; */
  /* double hR=0.1737014163; */
  /* double uR=5.757005447; */

  /* double hR=1; */
  /* double uR=1; */
  /* double hL=0.1737014163; */
  /* double uL=5.757005447; */

  double hL = 2;
  double uL = 0;
  double hR = 1;
  double uR = 0;


  double wL[_M]={hL, hL * uL};
  double wR[_M]={hR, hR * uR};

  int nx = 1000;

  double xmin=-9;
  double xmax=9;

  double dx = (xmax - xmin) / nx;

  double t = 1;

  FILE *gnufile = NULL;

  gnufile = fopen("plotriem","w");

  for(int i = 0; i <= nx; i++){
    double xi = (xmin + i * dx) / t;
    double w[_M];
    riemann(wL, wR, xi, w);
    fprintf(gnufile, "%f %f %f\n",xi, w[0], w[1] / w[0]);
  }

  fclose(gnufile);

}


double fz(double hL, double hR, double uL, double uR, double hs);
double dfz(double hL, double hR, double uL, double uR, double hs);

double Z (double h1, double h2)
{
  double g = 9.81;
  return((h1 < h2 ? 0.2e1 * sqrt(g) / (sqrt(h1) + sqrt(h2)) : sqrt(g * (h1 + h2) / h1 / h2 / 0.2e1)));
}


void riemann(double *wL, double *wR, double xi, double *w){

  double g = 9.81;
  double hL = wL[0];
  double hR = wR[0];
  double uL = wL[1]/hL;
  double uR = wR[1]/hR;

  double hs = 1e-8;
  double dh = 1e8;

  //int crit = (uR - uL < 2 * sqrt(g) *(sqrt(hL) - sqrt(hR)));

  //assert(crit); // apparition du vide

  //printf("xi=%f, hL=%f uL=%f hR=%f uR=%f\n",xi, hL,uL,hR,uR);

  while(fabs(dh) > 1e-6){
    double f = fz(hL, hR, uL, uR, hs);
    double df = dfz(hL, hR, uL, uR, hs);
    dh = -f / df;
    hs += dh;
    //printf("hs=%f dh=%f f=%f df=%f\n",hs,dh,f,df);
  }

  double us = uR + (hs - hR) * Z(hs,hR);

  // left wave
  double v1m,v1p;
  if (hs < hL) { // 1-detente
    v1m = uL-sqrt(g*hL);
    v1p = us-sqrt(g*hs);
  } else { // 1-chox
    double hc = (hs+hL)/2;
    double alpha = sqrt(hs)/(sqrt(hs)+sqrt(hL));
    double uc = alpha * us + (1-alpha) * uL;
    v1m = uc - sqrt(g * hc);
    v1p = v1m;
  }

  // right wave
  double v2m,v2p;
  if (hs < hR) { // 2-detente
    v2m = us + sqrt(g*hs);
    v2p = uR + sqrt(g*hR);
  } else { // 2-choc
    double hc = (hs+hR)/2;
    double alpha = sqrt(hs)/(sqrt(hs)+sqrt(hR));
    double uc = alpha * us + (1-alpha) * uR;
    v2m = uc + sqrt(g * hc);
    v2p = v2m;
  }

  //printf("v1m=%f v1p=%f v2m=%f v2p=%f\n",v1m,v1p,v2m,v2p);

  double h,u;

  if (xi < v1m){
    u = uL;
    h = hL;
  } else if (xi < v1p) {
    u = (2 * xi + uL + 2 * sqrt(g * hL)) / 3;
    h = (u - xi) * (u - xi) / g;
  } else if (xi < v2m) {
    u = us;
    h = hs;
  } else if (xi < v2p) {
    u = (2 * xi + uR - 2 * sqrt(g * hR)) / 3;
    h = (u - xi) * (u - xi) / g;
  } else {
    u = uR;
    h = hR;
  }

  w[0] = h;
  w[1] = h * u;

}

  double fz(double hL, double hR, double uL, double uR, double hs){
  return (hs - hR) * (hs < hR ? 0.6264183906e1 / (sqrt(hs) + sqrt(hR)) : sqrt((0.4905000000e1 * hs + 0.4905000000e1 * hR) / hs / hR)) + (hs - hL) * (hs < hL ? 0.6264183906e1 / (sqrt(hs) + sqrt(hL)) : sqrt((0.4905000000e1 * hs + 0.4905000000e1 * hL) / hs / hL)) - uL + uR;
}
double dfz(double hL, double hR, double uL, double uR, double hs){
 return  (hs < hR ? 0.6264183906e1 / (sqrt(hs) + sqrt(hR)) : sqrt((0.4905000000e1 * hs + 0.4905000000e1 * hR) / hs / hR)) + (hs - hR) * (hs < hR ? -0.3132091953e1 * pow(sqrt(hs) + sqrt(hR), -0.2e1) * pow(hs, -0.1e1 / 0.2e1) : pow((0.4905000000e1 * hs + 0.4905000000e1 * hR) / hs / hR, -0.1e1 / 0.2e1) * (0.4905000000e1 / hs / hR - (0.4905000000e1 * hs + 0.4905000000e1 * hR) * pow(hs, -0.2e1) / hR) / 0.2e1) + (hs < hL ? 0.6264183906e1 / (sqrt(hs) + sqrt(hL)) : sqrt((0.4905000000e1 * hs + 0.4905000000e1 * hL) / hs / hL)) + (hs - hL) * (hs < hL ? -0.3132091953e1 * pow(sqrt(hs) + sqrt(hL), -0.2e1) * pow(hs, -0.1e1 / 0.2e1) : pow((0.4905000000e1 * hs + 0.4905000000e1 * hL) / hs / hL, -0.1e1 / 0.2e1) * (0.4905000000e1 / hs / hL - (0.4905000000e1 * hs + 0.4905000000e1 * hL) * pow(hs, -0.2e1) / hL) / 0.2e1);


}
