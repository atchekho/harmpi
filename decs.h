/***********************************************************************************
    Copyright 2006 Charles F. Gammie, Jonathan C. McKinney, Scott C. Noble, 
                   Gabor Toth, and Luca Del Zanna

                        HARM  version 1.0   (released May 1, 2006)

    This file is part of HARM.  HARM is a program that solves hyperbolic 
    partial differential equations in conservative form using high-resolution
    shock-capturing techniques.  This version of HARM has been configured to 
    solve the relativistic magnetohydrodynamic equations of motion on a 
    stationary black hole spacetime in Kerr-Schild coordinates to evolve
    an accretion disk model. 

    You are morally obligated to cite the following two papers in his/her 
    scientific literature that results from use of any part of HARM:

    [1] Gammie, C. F., McKinney, J. C., \& Toth, G.\ 2003, 
        Astrophysical Journal, 589, 444.

    [2] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006, 
        Astrophysical Journal, 641, 626.

   
    Further, we strongly encourage you to obtain the latest version of 
    HARM directly from our distribution website:
    http://rainman.astro.uiuc.edu/codelib/


    HARM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    HARM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HARM; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/*************************************************************************
      COMPILE-TIME PARAMETERS : 
*************************************************************************/

//which problem
#define MONOPOLE_PROBLEM_1D 1
#define MONOPOLE_PROBLEM_2D 2
#define TORUS_PROBLEM 2
#define BONDI_PROBLEM_1D 3
#define BONDI_PROBLEM_2D 4

#define WHICHPROBLEM MONOPOLE_PROBLEM_1D


/** here are the few things that we change frequently **/

#if WHICHPROBLEM == MONOPOLE_PROBLEM_1D
#define N1       (2*386)        /* number of physical zones in X1-direction */
#define N2       (1)        /* number of physical zones in X2-direction */
#elif WHICHPROBLEM == MONOPOLE_PROBLEM_2D
#define N1       (2*386)        /* number of physical zones in X1-direction */
#define N2       (256)        /* number of physical zones in X2-direction */
#elif WHICHPROBLEM == TORUS_PROBLEM
#define N1       (256)        /* number of physical zones in X1-direction */
#define N2       (256)        /* number of physical zones in X2-direction */
#elif WHICHPROBLEM == BONDI_PROBLEM_1D
#define N1       (256)        /* number of physical zones in X1-direction */
#define N2       (1)        /* number of physical zones in X2-direction */
#elif WHICHPROBLEM == BONDI_PROBLEM_2D
#define N1       (256)        /* number of physical zones in X1-direction */
#define N2       (256)        /* number of physical zones in X2-direction */
#endif


#define NPR        (8)        /* number of primitive variables */
#define NDIM       (4)        /* number of total dimensions.  Never changes */
#define NPG        (4)        /* number of positions on grid for grid functions */
#define COMPDIM    (2)        /* number of non-trivial spatial dimensions used in computation */

#define NIMG       (4)        /* Number of types of images to make, kind of */


/* whether or not to use Font's  adiabatic/isothermal prim. var. inversion method: */
#define DO_FONT_FIX (1)

/* how many cells near the poles to stabilize, choose 0 for no stabilization */
#define POLEFIX (2)

/* whether or not to rescale primitive variables before interpolating them for flux/BC's */
#define RESCALE     (0)

/** FIXUP PARAMETERS, magnitudes of rho and u, respectively, in the floor : **/
#define RHOMIN	(1.e-4)
#define UUMIN	(1.e-6)
#define RHOMINLIMIT (1.e-20)
#define UUMINLIMIT  (1.e-20)
#define POWRHO (4)

#define FLOORFACTOR (1.)
#define BSQORHOMAX (50.*FLOORFACTOR)
#define BSQOUMAX (250.*FLOORFACTOR)
#define UORHOMAX (50.*FLOORFACTOR)


/* A numerical convenience to represent a small non-zero quantity compared to unity:*/
#define SMALL	(1.e-20)

/* Max. value of gamma, the lorentz factor */
#define GAMMAMAX (50.)

/* maximum fractional increase in timestep per timestep */
#define SAFE	(1.3)


#define COORDSINGFIX 1
// whether to move polar axis to a bit larger theta
// theta value where singularity is displaced to
#define SINGSMALL (1.E-20)


/* I/O format strings used herein : */
#define FMT_DBL_OUT "%28.18e"
#define FMT_INT_OUT "%10d"



/*************************************************************************
    MNEMONICS SECTION 
*************************************************************************/
/* boundary condition mnemonics */
#define OUTFLOW	(0)
#define SYMM	(1)
#define ASYMM	(2)
#define FIXED	(3)

/* mnemonics for primitive vars; conserved vars */
#define RHO	(0)	
#define UU	(1)
#define U1	(2)
#define U2	(3)
#define U3	(4)
#define B1	(5)
#define B2	(6)
#define B3	(7)

/* mnemonics for dimensional indices */
#define TT	(0)     
#define RR	(1)
#define TH	(2)
#define PH	(3)

/* mnemonics for centering of grid functions */
#define FACE1	(0)	
#define FACE2	(1)
#define CORN	(2)
#define CENT	(3)

/* mnemonics for slope limiter */
#define MC	(0)
#define VANL	(1)
#define MINM	(2)

/* mnemonics for diagnostic calls */
#define INIT_OUT	(0)
#define DUMP_OUT	(1)
#define IMAGE_OUT	(2)
#define LOG_OUT		(3)
#define FINAL_OUT	(4)

/* Directional Mnemonics */
// -------------> r
// |         3    
// |        1-0   
// |         2    
// v            
// theta      
#define X1UP    (0)
#define X1DN    (1)
#define X2UP    (2)
#define X2DN    (3)


/* failure modes */
#define FAIL_UTOPRIM        (1)
#define FAIL_VCHAR_DISCR    (2)
#define FAIL_COEFF_NEG	    (3)
#define FAIL_COEFF_SUP	    (4)
#define FAIL_GAMMA          (5)
#define FAIL_METRIC         (6)


/* For rescale() operations: */
#define FORWARD 1
#define REVERSE 2

/*************************************************************************
    GLOBAL ARRAY SECTION 
*************************************************************************/
extern double  a_p[N1+4][N2+4][NPR] ;	/* space for primitive vars */
extern double  a_dq[N1+4][N2+4][NPR] ;  /* slopes */
extern double  a_F1[N1+4][N2+4][NPR] ;	/* fluxes */
extern double  a_F2[N1+4][N2+4][NPR] ;	/* fluxes */
extern double  a_ph[N1+4][N2+4][NPR] ;	/* half-step primitives */
extern int     a_pflag[N1+4][N2+4];	/* identifies failure points */

/* for debug */
extern double psave[N1][N2][NPR] ;
extern double fimage[NIMG][N1*N2];
extern int    failimage[5][N1*N2];

/* grid functions */
extern double a_conn[N1+4][N2+4][NDIM][NDIM][NDIM] ;
extern double a_gcon[N1+4][N2+4][NPG][NDIM][NDIM] ;
extern double a_gcov[N1+4][N2+4][NPG][NDIM][NDIM] ;
extern double a_gdet[N1+4][N2+4][NPG] ;

extern double (*   p)[N2+4][NPR] ;
extern double (*  dq)[N2+4][NPR] ;
extern double (*  F1)[N2+4][NPR] ;
extern double (*  F2)[N2+4][NPR] ;
extern double (*  ph)[N2+4][NPR] ;
extern int    (*  pflag)[N2+4];
extern double (* conn)[N2+4][NDIM][NDIM][NDIM] ;
extern double (* gcon)[N2+4][NPG][NDIM][NDIM] ;
extern double (* gcov)[N2+4][NPG][NDIM][NDIM] ;
extern double (* gdet)[N2+4][NPG] ;


#if(DO_FONT_FIX)
extern double Katm[N1];
#endif


/*************************************************************************
    GLOBAL VARIABLES SECTION 
*************************************************************************/
/* physics parameters */
extern double a ;
extern double gam ;

/* numerical parameters */
extern double Rin,Rout,hslope,R0,fractheta ;
extern double cour ;
extern double dV,dx[NPR],startx[NPR] ;
extern double dt ;
extern double t,tf ;
extern double x1curr,x2curr ;
extern int nstep ;

/* output parameters */
extern double DTd ;
extern double DTl ;
extern double DTi ;
extern int    DTr ;
extern int    dump_cnt ;
extern int    image_cnt ;
extern int    rdump_cnt ;
extern int    nstroke ;

/* global flags */
extern int failed ;
extern int lim ;
extern double defcon ;

/* diagnostics */
extern double mdot ;
extern double edot ;
extern double ldot ;

/* set global variables that indicate current local metric, etc. */
extern int icurr,jcurr,pcurr ;
struct of_geom {
	double gcon[NDIM][NDIM] ;
	double gcov[NDIM][NDIM] ;
	double g ;
} ;

struct of_state {
	double ucon[NDIM] ;
	double ucov[NDIM] ;
	double bcon[NDIM] ;
	double bcov[NDIM] ;
} ;

/*************************************************************************
    MACROS
*************************************************************************/
/* loop over all active zones */
#define ZLOOP for(i=0;i<N1;i++)for(j=0;j<N2;j++)

/* loop over all active zones */
#define IMAGELOOP for(j=0;j<N2;j++)for(i=0;i<N1;i++)

/* specialty loop */
extern int istart,istop,jstart,jstop ;
#define ZSLOOP(istart,istop,jstart,jstop) \
        for(i=istart;i<=istop;i++)\
	for(j=jstart;j<=jstop;j++)

/* loop over Primitive variables */
#define PLOOP  for(k=0;k<NPR;k++)
/* loop over all Dimensions; second rank loop */
#define DLOOP  for(j=0;j<NDIM;j++) for(k=0;k<NDIM;k++)
/* loop over all Dimensions; first rank loop */
#define DLOOPA for(j=0;j<NDIM;j++)
/* loop over all Space dimensions; second rank loop */
#define SLOOP  for(j=1;j<NDIM;j++) for(k=1;k<NDIM;k++)
/* loop over all Space dimensions; first rank loop */
#define SLOOPA for(j=1;j<NDIM;j++)


extern double fval1,fval2;
#define MY_MIN(fval1,fval2) ( ((fval1) < (fval2)) ? (fval1) : (fval2))
#define MY_MAX(fval1,fval2) ( ((fval1) > (fval2)) ? (fval1) : (fval2))

#define delta(i,j) ( (i == j) ? 1. : 0.)
#define dot(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]) 


/*************************************************************************
    FUNCTION DECLARATIONS 
*************************************************************************/
double bl_gdet_func(double r, double th) ;
double bsq_calc(double *pr, struct of_geom *geom) ;
int    gamma_calc(double *pr, struct of_geom *geom, double *gamma) ;
double gdet_func(double lgcov[][NDIM]) ;
double mink(int j, int k) ; 
double ranc(int seed) ;
double slope_lim(double y1, double y2, double y3) ;

int restart_init(void) ;

void area_map(int i, int j, double prim[][N2+4][NPR]) ;
void bcon_calc(double *pr, double *ucon, double *ucov, double *bcon) ;
void blgset(int i, int j, struct of_geom *geom);
void bl_coord(double *X, double *r, double *th) ;
void bl_gcon_func(double r, double th, double gcov[][NDIM]) ;
void bl_gcov_func(double r, double th, double gcov[][NDIM]) ;
void bound_prim(double pr[][N2+4][NPR]) ;
void conn_func(double *X, struct of_geom *geom, double lconn[][NDIM][NDIM]) ;
void coord(int i, int j, int loc, double *X) ;
void diag(int call_code) ;
void diag_flux(double F1[][N2+4][NPR], double F2[][N2+4][NPR]) ;
void dump(FILE *fp) ;
void gdump(FILE *fp) ;

void fail(int fail_type) ;
void fixup(double (* pv)[N2+4][NPR]) ;
void fixup1zone( int i, int j, double prim[NPR] ) ;
void fixup_utoprim( double (*pv)[N2 + 4][NPR] )  ;
void set_Katm(void);
int  get_G_ATM( double *g_tmp );
void fix_flux(double F1[][N2+4][NPR], double F2[][N2+4][NPR]) ;
void flux_ct(double F1[][N2+4][NPR],double F2[][N2+4][NPR]) ;
void gaussj(double **tmp, int n, double **b, int m) ;
void gcon_func(double lgcov[][NDIM], double lgcon[][NDIM]) ;
void gcov_func(double *X, double lgcov[][NDIM]) ;
void get_geometry(int i, int j, int loc, struct of_geom *geom) ;
void get_state(double *pr, struct of_geom *geom, struct of_state *q) ;
void image_all( int image_count ) ;
void init(void) ;
void lower(double *a, struct of_geom *geom, double *b) ;
void ludcmp(double **a, int n, int *indx, double *d) ;
void mhd_calc(double *pr, int dir, struct of_state *q, double *mhd)  ;
void misc_source(double *ph, int ii, int jj, struct of_geom *geom,
		                struct of_state *q, double *dU, double Dt) ;
void primtoflux(double *pa, struct of_state *q, int dir, struct of_geom *geom,
			double *fl) ;
void primtoU(double *p, struct of_state *q, struct of_geom *geom, double *U);
void raise(double *v1, struct of_geom *geom, double *v2) ;
void rescale(double *pr, int which, int dir, int ii, int jj, int face, 
			struct of_geom *geom) ;
void restart_write(void) ;
void restart_read(FILE *fp) ;
void set_arrays(void) ;
void set_grid(void) ;
void set_points(void) ;
void step_ch(void) ;
void source(double *pa, struct of_geom *geom, int ii, int jj, double *Ua,double Dt) ;
void timestep(void) ;
void u_to_v(double *pr, int i, int j) ;
void ucon_calc(double *pr, struct of_geom *geom, double *ucon) ;
void usrfun(double *pr,int n,double *beta,double **alpha) ;
void Utoprim(double *Ua, struct of_geom *geom, double *pa) ;
int Utoprim_2d(double U[NPR], double gcov[NDIM][NDIM], double gcon[NDIM][NDIM], 
               double gdet, double prim[NPR]);
int Utoprim_1dvsq2fix1(double U[NPR], double gcov[NDIM][NDIM], double gcon[NDIM][NDIM], 
               double gdet, double prim[NPR], double K );
int Utoprim_1dfix1(double U[NPR], double gcov[NDIM][NDIM], double gcon[NDIM][NDIM], 
               double gdet, double prim[NPR], double K );

void vchar(double *pr, struct of_state *q, struct of_geom *geom,
		int dir, double *cmax, double *cmin) ;

int invert_matrix( double A[][NDIM], double Ainv[][NDIM] );
int LU_decompose( double A[][NDIM], int permute[] );
void LU_substitution( double A[][NDIM], double B[], int permute[] );


