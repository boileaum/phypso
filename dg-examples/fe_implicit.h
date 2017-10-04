#ifndef _FE_IMPLICIT_H
#define _FE_IMPLICIT_H
#include "skyline.h"


///////////////// physical data ////////////////////

// nb of kinetic equations
#define M 3

// sound speed
#define CSON (0.6)

// relaxation parameter
#define RELAX (10000)

// Riemann data
#define RL 2
#define UL -0.

#define RR 1
#define UR 0.

// estimated maximal wave speed
#define VMAX (-1.)

// viscosity
#define VISC (1e-10)
//#define VISC 1

// physical flux
void flux(dcmplx* w, dcmplx* flux);

// jacobian of the physical flux
void dflux(dcmplx* w, dcmplx* dflux);

// exact or reference solution 
void solexacte(double x, double t, dcmplx* w);

///////////// interp. data ////////////////////

// number of finite elements
#define NB_ELEMS (100)
// polynomial order
#define DEG (3)


// nb of finite element nodes
#define NB_NODES (DEG * NB_ELEMS + 1)

// nb of unknowns on the mesh
#define WSIZE (M * NB_NODES)

// max newton iteration
#define NEWTON 1

/////////// main container ////////////

// a struct for managing all this...
typedef struct galerkin{

  dcmplx wn[WSIZE];
  dcmplx wnm1[WSIZE];
  dcmplx dw[WSIZE];
  dcmplx res[WSIZE];

  // nodes positions
  double xnode[NB_NODES];
  
  // mesh bounds
  double xmin;
  double xmax;
  double dx;

  // final time, time step and cfl
  double tmax, cfl, t;
  dcmplx dt;
  
  // time step multiplicator (1/2 or (1+I)/2)
  dcmplx smul;

  // skyline matrices
  // convective part
  Skyline sky_flux;

  // viscous part
  Skyline sky_visc;

} galerkin;

// constructor: parameters and init. condition.
void gal_construct(galerkin *gal,
		   double xmin,
		   double xmax,
		   double cfl,
		   double tmax,
		   dcmplx smult);

// data access functions
// memory location of var iv at gauss point ipg
int vindex(int ipg, int iv);

// connectivity array
// return global node index of local node iloc in elem ie
int connec(int ie, int iloc);


// destructor: free memory
void gal_destruct(galerkin *gal);

// matrix allocation and init.
void gal_mat_alloc(galerkin *gal);

// non linear iteration assembly (matrix part)
void gal_mat_assembly(galerkin *gal);

// non linear iteration assembly (residual part)
void gal_res_assembly(galerkin *gal);

// non linear iteration solve (result stored into dw)
void gal_mat_solve(galerkin *gal);

// advance one time step
void gal_step(galerkin *gal);

// outputs
void gal_plot(galerkin *gal);

// error measure
double gal_L2_error(galerkin *gal, int numvar);

// BGK projection
void bgk_relax(dcmplx *w, dcmplx dt);

#endif
