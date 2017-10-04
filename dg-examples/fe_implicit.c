#include "fe_implicit.h"
#include "lobatto.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// gcc fe_implicit.c skyline.c -O -lm

int main(void) {

  // small sys. lin. test
  {
    dcmplx A[16] = {.2,-1,0,0,
		    -1,2,-1,0,
		    0,-1,2,-1,
		    0.1,0,-1,2};

    int sigma[4] = {0,1,2,3};

    dcmplx b[4] = {-0.8,0,0,1};

    dcmplx x[4];

    //PLU_Square(4, A, sigma);
    //PLU_Solve(4, A, sigma, b, x);
    dcmplx B[16];
    InvertSquare(4, A, B);
    PLU_Square(4, A, sigma);
    PLU_Solve(4, A, sigma, b, x);
    int n = 4;
    double res = 0;
    for(int i = 0; i < n; i++){
      for(int k = 0; k < n; k++) x[i] -= B[i*n+k]*b[k];
      //printf("b=(%f,%f)\n",creal(x[i]),cimag(x[i]));
      res += cabs(x[i]);
    }
  

    printf("test inversion, erreur=%f\n",res);
    assert(res < 1e-14);
  }


  double xmin = -1;
  double xmax = 1;

  double cfl = 0.1;

  dcmplx c = (1. + I) / 2;
  ///dcmplx c = 1. / 2;

  double tmax = 0.4;

  static galerkin gal;

  gal_construct(&gal, xmin, xmax, cfl, tmax, c);
  gal_mat_assembly(&gal);
  gal.t = 0;
  int iter=0;
  int itermax = tmax / creal(gal.dt);
  assert(itermax % 2 == 0);
  while(iter < itermax){
    
    //////// conjugate the time step and LU ////
    gal.dt = conj(gal.dt);
    for(int k = 0; k < gal.sky_flux.nmem; k++) {
      gal.sky_flux.vkgs[k] = conj(gal.sky_flux.vkgs[k]);
      gal.sky_flux.vkgi[k] = conj(gal.sky_flux.vkgi[k]);
    }
    for(int k = 0; k < WSIZE; k++)
      gal.sky_flux.vkgd[k] = conj(gal.sky_flux.vkgd[k]);
    ///////////////////////////////////////////

    gal.t += creal(gal.dt);
    iter++;
    printf("iter=%d t=%f dt=(%f,%f)\n",iter,gal.t,creal(gal.dt),cimag(gal.dt));
    gal_step(&gal);

  }

  gal_plot(&gal);

  for(int iv = 0; iv < M; iv++)
    printf("erreur L2 var %d=%f\n", iv,gal_L2_error(&gal, iv));

}

void flux(dcmplx* w, dcmplx* flux){
  flux[0] = w[0] * VMAX;
  if (M == 2) flux[1] = -w[1] * VMAX;
  if (M == 3){
    flux[1] = 0;
    flux[2] = -w[2] * VMAX;
  }
}

void dflux(dcmplx* w, dcmplx* dflux){
  dflux[0] = VMAX;
  if (M == 2) {
    dflux[1] = 0;
    dflux[2] = 0;
    dflux[3] = -VMAX;
  }
  if (M == 3) {
    dflux[0] = VMAX;
    dflux[1] = 0;
    dflux[2] = 0;

    dflux[3] = 0;
    dflux[4] = 0;
    dflux[5] = 0;

    dflux[6] = 0;
    dflux[7] = 0;
    dflux[8] = -VMAX;
  }
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

void solexacte(double x, double t, dcmplx* w){

  // smooth solution test case
  /* double xi = x - VMAX * t; */
  /* w[0] = 1+exp(- 30 * xi * xi); */
  /* xi = x; */
  /* w[1] = 1+exp(- 30 * xi * xi); */
  /* xi = x + VMAX * t; */
  /* w[2] = 1+exp(- 30 * xi * xi); */

  //Discontinuous solution test case
  double xi;
  xi = x - VMAX * t;
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
  xi = x + VMAX * t;
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
  //printf("abs dt=%f\n",10/cabs(dt));
  dcmplx z = cexp(-relax * dt);
  //double z = 0;

  w[0] = z * dw[0] + weq[0];
  w[1] = z * dw[1] + weq[1];
  w[2] = z * dw[2] + weq[2];

}


int vindex(int ipg, int iv)
{
  return iv + ipg * M;
}

int connec(int ie, int iloc)
{
  return DEG * ie + iloc;
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

  gal->dx = (xmax - xmin)/NB_ELEMS;

  if (VMAX != 0) {
    gal->dt  = cfl * gal->dx / fabs(VMAX) * smul;
  } else {
    gal->dt  = cfl * gal->dx * smul;
  }
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

  gal_mat_alloc(gal);


}

void gal_mat_alloc(galerkin *gal)
{
  // construct the matrix profile
  InitSkyline(&gal->sky_flux, NB_NODES * M);
  InitSkyline(&gal->sky_visc, NB_NODES * M);
  for(int ie = 0; ie < NB_ELEMS; ie++){
    for(int iloc = 0; iloc <= DEG; iloc++){
      int ino = connec(ie, iloc);
      for(int jloc = 0; jloc <= DEG; jloc++){
	int jno = connec(ie, jloc);
	for(int iv = 0; iv < M; iv++){
	  int iw = vindex(ino, iv);
	  for(int jv = 0; jv < M; jv++){
	    int jw = vindex(jno, jv);
	    SwitchOn(&gal->sky_flux, iw, jw);
	    SwitchOn(&gal->sky_visc, iw, jw);
	  }
	}
      }
    }
  }

  AllocateSkyline(&gal->sky_flux);
  AllocateSkyline(&gal->sky_visc);

  // construct the viscous matrix
  for(int ie = 0; ie < NB_ELEMS; ie++){
    for(int ipg = 0; ipg <= DEG; ipg++){
      for(int iloc = 0; iloc <= DEG; iloc++){
	int ino = connec(ie, iloc);
	for(int jloc = 0; jloc <= DEG; jloc++){
	  int jno = connec(ie, jloc);
	  dcmplx val = VISC / gal->dx * wglop(DEG, ipg) *
	    dlag(DEG, iloc, ipg) * dlag(DEG, jloc, ipg);
	  //printf("v=%f\n",creal(val));
	  for(int iv = 0; iv < M; iv++){
	    int iw = vindex(ino, iv);
	    int jw = vindex(jno, iv);
	    AddSkyline(&gal->sky_visc, iw, jw, val);
	  }
	}
      }
    }
  }
  
  //DisplaySkyline(&gal->sky_visc);
  
}

void gal_mat_assembly(galerkin *gal)
{

  static bool not_done = true;

  //ZeroSkyline(&gal->sky_flux);

  if (not_done){
    // copy viscosity matrix into sky_flux
    for(int k = 0; k < gal->sky_flux.nmem; k++) {
      gal->sky_flux.vkgs[k] = gal->sky_visc.vkgs[k];
      gal->sky_flux.vkgi[k] = gal->sky_visc.vkgi[k];
    }
    for(int k = 0; k < WSIZE; k++)
      gal->sky_flux.vkgd[k] = gal->sky_visc.vkgd[k];
    
    for(int ie = 0; ie < NB_ELEMS; ie++){
      for(int kloc = 0; kloc <= DEG; kloc++){   // k: quadrature point
	int kno = connec(ie, kloc);
	dcmplx w[M];
	for(int iv = 0; iv < M; iv++){
	  int kw = vindex(kno,iv);
	  w[iv] =  gal->wn[kw];
	  dcmplx val = gal->dx * wglop(DEG, kloc) / gal->dt;
	  AddSkyline(&gal->sky_flux, kw, kw, val);
	}
	dcmplx df[M * M];
	dflux(w, df);
	/* printf("jac = %f %f \n   %f %f\n\n",creal(df[0]),creal(df[1]), */
	/* 	     creal(df[2]),creal(df[3])); */
	for(int iloc = 0; iloc <= DEG; iloc++){
	  int ino = connec(ie, iloc);
	  for(int iv = 0; iv < M; iv++){
	    int iw = vindex(ino, iv);
	    for(int kv = 0; kv < M; kv++){
	      int kw = vindex(kno, kv);
	      dcmplx val = -wglop(DEG, kloc) *	    
		dlag(DEG, iloc, kloc) * df[iv * M + kv];
	      //printf("iv=%d kv=%d val=%f dlag=%f\n",iv,kv,creal(val),
	      //	   dlag(DEG, iloc, kloc));
	      AddSkyline(&gal->sky_flux, iw, kw, val);
	    }
	  }
	}
      }
    }
    //DisplaySkyline(&gal->sky_flux);
    //assert(1==2);
    // boundary conditions
    for(int iv = 0; iv < M; iv++){
      int ino = connec(0,0);
      int iw = vindex(ino,iv);
      AddSkyline(&gal->sky_flux, iw, iw, 1e20);
      ino = connec(NB_ELEMS - 1, DEG);
      iw = vindex(ino,iv);
      AddSkyline(&gal->sky_flux, iw, iw, 1e20);
    }
    FactoLU(&gal->sky_flux);  
    not_done = false;
    printf("first factolu finished...\n");
  }
}

void gal_res_assembly(galerkin *gal)
{

  // apply viscous term
  MatVectSkyline(&gal->sky_visc, gal->wn, gal->res);

  // change sign
  for(int iw = 0; iw < WSIZE; iw++) gal->res[iw] = -gal->res[iw];


  for(int ie = 0; ie < NB_ELEMS; ie++){
    for(int kloc = 0; kloc <= DEG; kloc++){   // k: quadrature point
      int kno = connec(ie, kloc);
      dcmplx w[M], wm1[M];
      for(int iv = 0; iv < M; iv++){
	int kw = vindex(kno,iv);
	w[iv] =  gal->wn[kw];
	wm1[iv] =  gal->wnm1[kw];
      }
      dcmplx s[M];
      source(gal->xnode[kno], gal->t, w, s);      
      for(int iv = 0; iv < M; iv++){
	int kw = vindex(kno,iv);
	gal->res[kw] -= gal->dx * wglop(DEG, kloc) *
	  (w[iv] - gal->wnm1[kw]) / gal->dt;
	gal->res[kw] += gal->dx * wglop(DEG, kloc) * s[iv];
	//printf("%f\n",creal(w[iv] - gal->wnm1[kw]));
      }
      dcmplx f[M];
      flux(w, f);
      for(int iloc = 0; iloc <= DEG; iloc++){
	int ino = connec(ie, iloc);
	for(int iv = 0; iv < M; iv++){
	  int iw = vindex(ino, iv);
	  dcmplx val = -wglop(DEG, kloc) *	    
	    dlag(DEG, iloc, kloc) * f[iv];
	  gal->res[iw] -= val;
	}
      }
    }
  }
  /* for(int iw = 0; iw < WSIZE; iw++){ */
  /*   printf("res=%f\n",creal(gal->res[iw])); */
  /* } */
}

double gal_L2_error(galerkin *gal, int numvar)
{

  double gal_err = 0;
  double vol = 0;
  
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
      
      vol += wglop(DEG, kloc) * gal->dx;
      
      double cerr = cabs(w[numvar]-wex[numvar]);
      //printf("%f %f\n",cabs(w[numvar]-wex[numvar]),vol);
      
      gal_err += cerr * cerr * wglop(DEG, kloc) * gal->dx;
    }
  }
  gal_err = sqrt(gal_err / vol);
  return gal_err;
}


// non linear iteration solve (result stored into dw)
void gal_mat_solve(galerkin *gal)
{
  for(int iter = 0; iter < NEWTON; iter++){
    gal_res_assembly(gal);
    gal_mat_assembly(gal);
    FastSolveSkyline(&gal->sky_flux,gal->res,gal->dw);
    double crit = 0;
    for(int iw = 0; iw < WSIZE; iw++) {
      gal->wn[iw] += gal->dw[iw];
      double crit2 = cabs(gal->dw[iw]);
      crit = crit > crit2 ? crit : crit2;
    }
    //printf("iter newton %d res=%f\n",iter, crit);
  }
      
}

void gal_step(galerkin *gal)
{
  gal_mat_solve(gal);

  for(int ino = 0; ino < NB_NODES; ino++){
    dcmplx *w = gal->wn + vindex(ino, 0);
    //bgk_relax(w, gal->dt); // warning: assume gathered w_i
  }


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

  dcmplx wex[M];
  dcmplx wnum[M];
  for(int ino = 0; ino < NB_NODES; ino++){
    double x = gal->xnode[ino];
    fprintf(gnufile, "%f ", x);
    double t=gal->t;
    solexacte(x, t, wex);
    for(int iv = 0; iv < M; iv++){
      int imem = vindex(ino,iv);
      wnum[iv] =  gal->wn[imem];
    }
    for(int iv = 0; iv < M; iv++){
      fprintf(gnufile, "%f ", creal(wex[iv]));
      fprintf(gnufile, "%f ", creal(wnum[iv]));
    }
    dcmplx rho = wnum[0] + wnum[1] + wnum[2];
    dcmplx u = (wnum[2] - wnum[0]) / rho;
    fprintf(gnufile, "%f ", creal(rho));
    fprintf(gnufile, "%f ", creal(u));
    fprintf(gnufile, "%f ", cimag(rho));
    fprintf(gnufile, "%f ", cimag(u));

    fprintf(gnufile, "\n");
  }
    
  fclose(gnufile);


}




    
