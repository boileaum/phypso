#include "dg_implicit.h"
#include "lobatto.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// gcc fe_implicit.c skyline.c -O -lm

int main(void) {

  //TestPLU();

  double xmin = -1;
  double xmax = 1;

  // cfl computed from the element size dx
  double cfl_dx = 0.1;
  // the actual cfl is computed from the smallest
  // distance bewteen two Gauss-Lobatto points
  //double cfl = cfl_dx * (glop(DEG, 1) - glop(DEG, 0));
  double cfl = cfl_dx;

  printf("cfl=%f\n",cfl);

  dcmplx c = (1. + I) / 2;
  //dcmplx c =1./2;

  double tmax = 0.4;

  static galerkin gal;

  gal.t = 0;
  gal_construct(&gal, xmin, xmax, cfl, tmax, c);
  int iter=0;

  while(gal.t < tmax){

    gal.t += 2* creal(gal.dt);
    // correction last time step
    if (gal.t > tmax) {
      double dt = 2* creal(gal.dt) - gal.t + tmax;
      gal.dt = dt * c;
      printf("last step\n");
      gal.t = tmax;
    }
   
    iter++;
    printf("iter=%d t=%f dt=(%f,%f)\n",iter,gal.t,
	   creal(gal.dt),cimag(gal.dt));
    gal.dt = conj(gal.dt);
    gal_step(&gal);
    /* for(int ino = 0; ino < NB_NODES; ino++){ */
    /*   printf("ino=%d w=%f\n",ino,creal(gal.wn[vindex(ino,2)])); */
    /* } */
    iter++;
    printf("iter=%d t=%f dt=(%f,%f)\n",iter,gal.t,
	   creal(gal.dt),cimag(gal.dt));
    gal.dt = conj(gal.dt);
    gal_step(&gal);

    /* for(int ino = 0; ino < NB_NODES; ino++){ */
    /*   printf("ino=%d w=%f\n",ino,creal(gal.wn[vindex(ino,2)])); */
    /* } */
    /* assert(1==2); */
  }

  gal_plot(&gal);

  for(int iv = 0; iv < M; iv++)
    printf("erreur L2 var %d=%e\n", iv,creal(gal_L2_error(&gal, iv)));

}


void source(double x, double t, dcmplx* w, dcmplx* source){

  //source[0] = 2*VMAX*x;
  source[0] = 0;
  if (M ==2) {
    //source[1] = -2*VMAX*x;
    source[1] = 0;
  }
  if (M == 3) {
    //source[1] = -2*VMAX*x;
    source[1] = 0;
    source[2] = 0;
  }

}


void bgk_project(dcmplx *w){

  dcmplx r = w[0] + w[1] + w[2];
  dcmplx q = w[2] - w[0];
  dcmplx u = q / r;

  dcmplx weq[3] = {q * (u - 1) / 2 + CSON * CSON * r / 2,
		   r * (1 - u * u - CSON * CSON),
		   q * (u + 1) / 2 + CSON * CSON * r / 2};

  w[0] = weq[0];
  w[1] = weq[1];
  w[2] = weq[2];

}


int vindex(int ipg, int iv)
{
  //return iv + ipg * M;
  return iv * NB_NODES + ipg;
}

int connec(int ie, int iloc)
{
  return (DEG + 1) * ie + iloc;
}


void gal_construct(galerkin *gal,
		   double xmin,
		   double xmax,
		   double cfl,
		   double tmax,
		   dcmplx smul)
{

  gal->xmin = xmin;
  gal->xmax = xmax;


  gal->cfl = cfl;
  gal->tmax = tmax;

  gal->smul = smul;

  gal->dx = (xmax - xmin) / NB_ELEMS;

  if (VMAX != 0) {
    gal->dt  = cfl * gal->dx / VMAX * smul;
  } else {
    gal->dt  = cfl * gal->dx * smul;
  }

  printf("dx=%f dt=(%f,%f)\n",gal->dx,creal(gal->dt),cimag(gal->dt));

  // nodes
  for(int ie = 0; ie < NB_ELEMS; ie++){
    for(int iloc = 0; iloc <= DEG; iloc++){
      double xi = glop(DEG, iloc);
      double x = gal->xmin + (ie + xi) * gal->dx;
      int ino = connec(ie, iloc);
      gal->xnode[ino] = x;
    }
  }

  // initial condition

  for(int ino = 0; ino < NB_NODES; ino++){
    double x = gal->xnode[ino];
    double t=0;
    dcmplx w[M];
    solexacte(x, t, w);
    for(int iv = 0; iv < M; iv++){
      int imem = vindex(ino,iv);
      gal->wn[imem] = w[iv];
      gal->wnm1[imem] = w[iv];
    }
  }

}

double gal_L2_error(galerkin *gal, int numvar)
{

  double gal_err = 0;
  double sav_err = 0;
  double vol = 0;
  
  FILE * ru_file = NULL;
  ru_file = fopen("fe_sav.bin", "rb" );
  dcmplx *rho_sav = NULL;
  dcmplx *u_sav = NULL;
  dcmplx *w_sav = NULL;
  double *x_sav = NULL;
  int deg_sav = 0;
  int nbelems_sav = 0;
  int nbnodes_sav = 0;

  if (ru_file){
    size_t ret;
    ret = fread(&nbelems_sav, sizeof(int), 1, ru_file);
    ret = fread(&deg_sav, sizeof(int), 1, ru_file);
    nbnodes_sav = (deg_sav + 1) * nbelems_sav;
    printf("deg=%d nbel=%d nbno=%d\n",deg_sav,nbelems_sav,nbnodes_sav),
    rho_sav = malloc(nbnodes_sav * sizeof(dcmplx));
    w_sav = malloc(nbnodes_sav * M * sizeof(dcmplx));
    u_sav = malloc(nbnodes_sav * sizeof(dcmplx));
    x_sav = malloc(nbnodes_sav * sizeof(double));
    for(int ino = 0; ino < nbnodes_sav; ino++){
      ret = fread(x_sav + ino, sizeof(double), 1, ru_file);
      for(int iv=0; iv < M; iv++) {
	ret = fread(w_sav + iv * nbnodes_sav + ino, sizeof(dcmplx), 1, ru_file);
      }
      ret = fread(rho_sav + ino, sizeof(dcmplx), 1, ru_file);
      ret = fread(u_sav + ino, sizeof(dcmplx), 1, ru_file);
      /* printf("ino=%d x=%f w=%f %f %f r=%f u=%f\n",ino, */
      /* 	     x_sav[ino], */
      /* 	     creal(w_sav[0 * nbnodes_sav + ino]), */
      /* 	     creal(w_sav[1 * nbnodes_sav + ino]), */
      /* 	     creal(w_sav[2 * nbnodes_sav + ino]), */
      /* 	     creal(rho_sav[ino]), */
      /* 	     creal(u_sav[ino])); */
    }
  } else {
    printf("no saved data for comparison\n");
  }
  
  for(int ie = 0; ie < NB_ELEMS; ie++){
    for(int kloc = 0; kloc <= DEG; kloc++){   // k: quadrature point
      int kno = connec(ie, kloc);
      double x = gal->xnode[kno];
      dcmplx w[M];
      dcmplx wex[M];
      solexacte(x, gal->t, wex);
      for(int iv = 0; iv < M; iv++){
	int kw = vindex(kno,iv);
	w[iv] =  gal->wn[kw];
      }

      dcmplx ws[M];
      if (ru_file) {
	interpolate(x_sav, w_sav,
		    deg_sav, nbnodes_sav,
		    x, ws);
      }
      
      vol += wglop(DEG, kloc) * gal->dx;
      
      double cerr = cabs(w[numvar]-wex[numvar]);
      double cerr_sav = 1e10;
      if (ru_file) cerr_sav = cabs(w[numvar]-ws[numvar]);
      //printf("%f %f\n",cabs(w[numvar]-wex[numvar]),vol);
      
      gal_err += cerr * cerr * wglop(DEG, kloc) * gal->dx;
      sav_err += cerr_sav * cerr_sav * wglop(DEG, kloc) * gal->dx;
    }
  }
  gal_err = sqrt(gal_err / vol);
  if (ru_file) {
    sav_err = sqrt(sav_err / vol);
    printf("err_sav=%e\n",sav_err);
  }

  if (ru_file){
    free(rho_sav);
    free(w_sav);
    free(u_sav);
    free(x_sav);
    fclose(ru_file);
  }


  
  return gal_err;
}

void lagrange(double *p, double *subdiv,
			 int deg, int ii, double x) {
  *p = 1;
  const int npg = deg + 1;
  for(int j = 0; j < npg; j++) {
    if (j != ii) {
      *p *= (x - subdiv[j]) / (subdiv[ii] - subdiv[j]);
    }
  }
}

// interpolate another function from a different mesh
// xi: interp points
// w_sav: stored data
// deg_sav, nbnodes_sav: quantities of saved data
// x: point where is computed the interp values
// ws: computed interpolation
void interpolate(double *xi, dcmplx *wi,
		 int deg, int nbnodes,
		 double x, dcmplx *wloc){

  double xmin = xi[0];
  double xmax = xi[nbnodes - 1];
  int nbelems = nbnodes / (deg + 1);
  assert(nbnodes % (deg + 1) == 0);

  assert(x >= xmin && x <= xmax);
  if (x == xmax) x = x - 1e-10 * (xmax-xmin);

  
  int numelem = (x-xmin)/(xmax-xmin) * nbelems;

  for(int iv = 0; iv < M; iv++) wloc[iv] = 0;
  
  for(int i = 0; i < deg + 1; i++){
    double val;
    double *sub = xi + numelem * (deg + 1);
    lagrange(&val, sub,
	     deg, i, x);
    //printf("x=%f ie=%d i=%d xi=%f val=%f\n",x,numelem,i,sub[i],val);
    int ino = numelem * (deg + 1) + i;
    for(int iv = 0; iv < M; iv++)
      wloc[iv] += wi[iv * nbnodes + ino] * val;
  }

  //printf("w=%f %f %f\n\n",creal(wloc[0]),creal(wloc[1]),creal(wloc[2]));


}


// compute and invert the local implicit matrix
void loc_assembly(galerkin *gal, dcmplx *mat, double vit){
  
  int n = DEG + 1;

  dcmplx a[n * n];
  for(int i = 0; i < n * n; i++) a[i] = 0;

  for(int i = 0; i <= DEG; i++) a[n * i + i] = 1; 

  dcmplx z = vit * gal->dt / gal->dx ;
  double wg = wglop(DEG, 0);
  
  if (vit > 0) {
     a[0] += z / wg;
  } else if (vit < 0) {
    a[n * n -1] -= z / wg;
  } else{
    assert(vit != 0);
  }
  
  for(int i = 0; i <= DEG; i++){
    for(int j = 0; j <= DEG; j++){
      a[n * i + j] += z * dlag(DEG, j, i);
    }
  }

  InvertSquare(n, a, mat); 

  // test
  /* for(int i = 0; i <= DEG; i++){ */
  /*   for(int j = 0; j <= DEG; j++){ */
  /*     dcmplx v = 0; */
  /*     for(int k =  0; k <= DEG; k++){ */
  /* 	v += a[n * i + k] * mat[n * k + j]; */
  /*     } */
  /*     v=a[n * i + j]; */
  /*     printf("i=%d j=%d aij=(%f,%f)\n",i,j, */
  /* 	     creal(v),cimag(v)); */
  /*   } */
  /* } */
  //assert(1==3);


}

void gal_step(galerkin *gal)
{


  int n = DEG + 1;

  dcmplx mat[n * n];
  // negative velocity
  int iv = 0;
  loc_assembly(gal, mat, -VMAX);
  for(int ie = NB_ELEMS - 1; ie >= 0; ie--) {
    // last glop in cell
    int ino = connec(ie, DEG);
    int iw = vindex(ino, iv);
    dcmplx wv[M];
    // boundary condition or upwind data
    if (ie < NB_ELEMS - 1) {
      wv[iv] = gal->wnm1[vindex(connec(ie + 1, 0),iv)];
    } else {
      solexacte(gal->xnode[ino],gal->t,wv);
      //printf("wR=(%f,%f)\n",creal(wv[iv]),cimag(wv[iv]));
    }
    gal->wnm1[iw] -= wv[iv] *  (-VMAX) *
      gal->dt / gal->dx / wglop(DEG, DEG);
    // local implicit scheme
    iw -= DEG; // return to first glop for matrix product
    for(int iloc = 0; iloc <= DEG; iloc++){
      gal->wn[iw + iloc] = 0;
      for(int jloc = 0; jloc <= DEG; jloc++)
	gal->wn[iw + iloc] += mat[n * iloc + jloc] * gal->wnm1[iw + jloc];
    }
    for(int iloc = 0; iloc <= DEG; iloc++) 
      gal->wnm1[iw + iloc] = gal->wn[iw + iloc];

  }

  // positive velocity
  iv = 2;
  loc_assembly(gal, mat, VMAX);
  for(int ie = 0; ie < NB_ELEMS; ie++) {
    // first glop in cell
    int ino = connec(ie, 0);
    int iw = vindex(ino, iv);
    dcmplx wv[M];
    // boundary condition or upwind data
    if (ie > 0) {
      wv[iv] = gal->wnm1[vindex(connec(ie - 1, DEG),iv)];
    } else {
      solexacte(gal->xnode[ino],gal->t,wv);
    }
    gal->wnm1[iw] += wv[iv] *  VMAX *
      gal->dt / gal->dx / wglop(DEG, 0);
    /* printf("wini=%f\n",creal(gal->wnm1[iw])); */
    /* assert(1==4); */
    // local implicit scheme
    for(int iloc = 0; iloc <= DEG; iloc++){
      gal->wn[iw + iloc] = 0;
      for(int jloc = 0; jloc <= DEG; jloc++)
	gal->wn[iw + iloc] += mat[n * iloc + jloc] * gal->wnm1[iw + jloc];
    }
    for(int iloc = 0; iloc <= DEG; iloc++) 
      gal->wnm1[iw + iloc] = gal->wn[iw + iloc];
    
  }

  for(int ino = 0; ino < NB_NODES; ino++){
    dcmplx wloc[M];
    for(int iv = 0; iv < M; iv++) wloc[iv] =gal->wn[vindex(ino,iv)]; 
    //bgk_project(wloc);
    for(int iv = 0; iv < M; iv++) gal->wn[vindex(ino,iv)]=wloc[iv]; 
  }


  // update
  for(int iw = 0; iw < WSIZE; iw++) {
    gal->wnm1[iw] = gal->wn[iw];
  }
}


void gal_plot(galerkin *gal)
{

  // commande de tracé
  // plot 'fe.dat'  using 1:9 w l, 'rho_u.dat' using 1:3 w l
  // plot 'fe.dat'  using 1:8 w l, 'rho_u.dat' using 1:2 w l
  FILE * gnufile;
  gnufile = fopen("fe.dat", "w" );

  FILE * ru_file;
  ru_file = fopen("fe.bin", "wb" );

  int temp = NB_ELEMS;
  fwrite(&temp, sizeof(int), 1, ru_file);
  temp = DEG;
  fwrite(&temp, sizeof(int), 1, ru_file);


  dcmplx wex[M];
  dcmplx wnum[M];
  for(int ino = 0; ino < NB_NODES; ino++){
    double x = gal->xnode[ino];
    fprintf(gnufile, "%f ", x);
    fwrite(&x, sizeof(double), 1, ru_file);
    double t=gal->t;
    solexacte(x, t, wex);
    for(int iv = 0; iv < M; iv++){
      int imem = vindex(ino,iv);
      wnum[iv] =  gal->wn[imem];
    }
    for(int iv = 0; iv < M; iv++){
      fprintf(gnufile, "%f ", creal(wex[iv]));
      fprintf(gnufile, "%f ", creal(wnum[iv]));
      fwrite(wnum + iv, sizeof(dcmplx), 1, ru_file);
   }
    dcmplx rho = wnum[0] + wnum[1] + wnum[2];
    dcmplx u = (wnum[2] - wnum[0]) / rho;
    fprintf(gnufile, "%f ", creal(rho));
    fprintf(gnufile, "%f ", creal(u));
    fwrite(&rho, sizeof(dcmplx), 1, ru_file);
    fwrite(&u, sizeof(dcmplx), 1, ru_file);
    fprintf(gnufile, "%f ", cimag(rho));
    fprintf(gnufile, "%f ", cimag(u));

    fprintf(gnufile, "\n");
  }
    
  fclose(gnufile);
  fclose(ru_file);


}

void solexacte(double x, double t, dcmplx* w){

  // smooth solution test case
  /* double xi = x + VMAX * t; */
  /* w[0] = 1 + cexp(- 30 * xi * xi); */
  /* xi = x; */
  /* w[1] = 1 + cexp(- 30 * xi * xi); */
  /* xi = x - VMAX * t; */
  /* w[2] = 1 + cexp(- 30 * xi * xi); */
  /* bgk_project(w); */

  /* bgk_relax(w, 1.); */

  // debug test case
  /* double xi = x + VMAX * t; */
  /* xi = x; */
  /* w[0] = (xi-2)*(xi+2); */
  /* xi = x; */
  /* w[1] = (xi-2)*(xi+2); */
  /* xi = x; */
  /* w[2] = (xi-2)*(xi+2); */
  /* w[2] = 1; */
  
  
  // Riemann problem test case
  double xi;
  xi = x + VMAX * t;
  if (xi < 0) {
    w[0] = RL * UL * (UL - 1) / 2 + CSON * CSON * RL / 2;
  } else {
    w[0] = RR * UR * (UR - 1) / 2 + CSON * CSON * RR / 2;
  }
  xi = x;
  if (xi < 0) {
    w[1] = RL * (1 - UL * UL - CSON * CSON);
  } else {
    w[1] = RR * (1 - UR * UR - CSON * CSON);
  }
  xi = x - VMAX * t;
  if (xi < 0) {
    w[2] = RL * UL * (UL + 1) / 2 + CSON * CSON * RL / 2;
  } else {
    w[2] = RR * UR * (UR + 1) / 2 + CSON * CSON * RR / 2;
  }


}

void bgk_relax(dcmplx *w, dcmplx dt){

  dcmplx r = w[0] + w[1] + w[2];
  dcmplx q = w[2] - w[0];
  dcmplx u = q / r;

  dcmplx weq[3] = {q * (u - 1) / 2 + CSON * CSON * r / 2,
		   r * (1 - u * u - CSON * CSON),
		   q * (u + 1) / 2 + CSON * CSON * r / 2};

  dcmplx dw[3] = { w[0] - weq[0],
  		   w[1] - weq[1],
  		   w[2] - weq[2]};

  //double relax = 5 / cabs(dt);
  double relax = RELAX;
  double tau = 1. / relax;
  //printf("abs dt=%f\n",10/cabs(dt));
  //double z = 0;

  dcmplx z = cexp(-relax * dt);
  z = 1. / (1 + relax * dt);
  w[0] = z * dw[0] + weq[0];
  w[1] = z * dw[1] + weq[1];
  w[2] = z * dw[2] + weq[2];

  /* dcmplx z = dt / (tau + dt / 2); */
  /* w[0] -= z * dw[0]; */
  /* w[1] -= z * dw[1]; */
  /* w[2] -= z * dw[2]; */

}

