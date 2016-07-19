//Modified by Alexander Tchekhovskoy: MPI+3D
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

#include "decs.h"

/***********************************************************************************************/
/***********************************************************************************************
 primtoflux():
 ---------
 --  calculate fluxes in direction dir,
 
 ***********************************************************************************************/

void primtoflux(double *pr, struct of_state *q, int dir,
		struct of_geom *geom, double *flux)
{
	int j,k, m ;
	double mhd[NDIM], EE4[NDIM], EE5[NDIM], EUAD[NDIM], EEVAR[NDIM];

	/* particle number flux */
	flux[RHO] = pr[RHO]*q->ucon[dir] ;

	mhd_calc(pr, dir, q, mhd) ;

	/* MHD stress-energy tensor w/ first index up, 
	 * second index down. */
	flux[UU] = mhd[0] + flux[RHO] ;
	flux[U1] = mhd[1] ;
	flux[U2] = mhd[2] ;
	flux[U3] = mhd[3] ;

	/* dual of Maxwell tensor */
	flux[B1]  = q->bcon[1]*q->ucon[dir] - q->bcon[dir]*q->ucon[1] ;
	flux[B2]  = q->bcon[2]*q->ucon[dir] - q->bcon[dir]*q->ucon[2] ;
	flux[B3]  = q->bcon[3]*q->ucon[dir] - q->bcon[dir]*q->ucon[3] ;
#if(DOKTOT)
        flux[KTOT] = flux[RHO]*pr[KTOT] ;
#endif
  
  PLOOP flux[m] *= geom->g ;
  
  
}

/* calculate "conserved" quantities; provided strictly for
 * historical reasons */
void primtoU(double *pr, struct of_state *q, struct of_geom *geom, double *U)
{
  
  primtoflux(pr,q,0,geom, U) ;
  return ;
}

/* calculate magnetic field four-vector */
void bcon_calc(double *pr, double *ucon, double *ucov, double *bcon)
{
  int j ;
  
  bcon[TT] = pr[B1]*ucov[1] + pr[B2]*ucov[2] + pr[B3]*ucov[3] ;
  for(j=1;j<4;j++)
    bcon[j] = (pr[B1-1+j] + bcon[TT]*ucon[j])/ucon[TT] ;
  
  return ;
}

/* MHD stress tensor, with first index up, second index down */
void mhd_calc(double *pr, int dir, struct of_state *q, double *mhd)
{
  int j ;
  double r,u,P,w,bsq,eta,ptot ;
  
  r = pr[RHO] ;
  u = pr[UU] ;
  P = (gam - 1.)*u ;
  w = P + r + u ;
  bsq = dot(q->bcon,q->bcov) ;
  eta = w + bsq ;
  ptot = P + 0.5*bsq;
  
  /* single row of mhd stress tensor,
   * first index up, second index down */
  DLOOPA mhd[j] = eta*q->ucon[dir]*q->ucov[j]
  + ptot*delta(dir,j) - q->bcon[dir]*q->bcov[j] ;
  
}

/* add in source terms to equations of motion */
void source(double *ph, struct of_geom *geom, int ii, int jj, int kk, double *dU,
            double Dt)
{
  double mhd[NDIM][NDIM], EE4[NDIM][NDIM],EE5[NDIM][NDIM], EUAD[NDIM][NDIM], EEVAR[NDIM][NDIM];
  int j,k,m ;
  struct of_state q ;
  
  get_state(ph, geom, &q) ;
  mhd_calc(ph, 0, &q, mhd[0]) ;
  mhd_calc(ph, 1, &q, mhd[1]) ;
  mhd_calc(ph, 2, &q, mhd[2]) ;
  mhd_calc(ph, 3, &q, mhd[3]) ;
  
  
  
  
  /* contract mhd stress tensor with connection */
  PLOOP dU[m] = 0. ;
  DLOOP {
    dU[UU] += mhd[j][k]*conn[ii][jj][kk][k][0][j] ;
    dU[U1] += mhd[j][k]*conn[ii][jj][kk][k][1][j] ;
    dU[U2] += mhd[j][k]*conn[ii][jj][kk][k][2][j] ;
    dU[U3] += mhd[j][k]*conn[ii][jj][kk][k][3][j] ;

  }


  
  
  
  
  //misc_source(ph, ii, jj, geom, &q, dU, Dt) ;
  
  PLOOP dU[m] *= geom->g ;
  
  /* done! */
}

/* returns b^2 (i.e., twice magnetic pressure) */
double bsq_calc(double *pr, struct of_geom *geom)
{
  struct of_state q ;
  
  get_state(pr,geom,&q) ;
  return( dot(q.bcon,q.bcov) ) ;
}

/* find ucon, ucov, bcon, bcov from primitive variables */
void get_state(double *pr, struct of_geom *geom, struct of_state *q)
{
  
  /* get ucon */
  ucon_calc(pr, geom, q->ucon) ;
  lower(q->ucon, geom, q->ucov) ;
  bcon_calc(pr, q->ucon, q->ucov, q->bcon) ;
  lower(q->bcon, geom, q->bcov) ;
  
  return ;
}

/* find relative 4-velocity from 4-velocity (both in code coords) */
void ucon_to_utcon(double *ucon,struct of_geom *geom, double *utcon)
{
  double alpha, beta[NDIM], gamma;
  int j;
  
  /* now solve for v-- we can use the same u^t because
   * it didn't change under KS -> KS' */
  alpha = 1./sqrt(-geom->gcon[TT][TT]) ;
  SLOOPA beta[j] = geom->gcon[TT][j]*alpha*alpha ;
  gamma = alpha*ucon[TT] ;
  
 
  utcon[0] = 0;
  SLOOPA utcon[j] = ucon[j] + gamma*beta[j]/alpha ;
}

/* find contravariant four-velocity */
void ucon_calc(double *pr, struct of_geom *geom, double *ucon)
{
  double alpha,gamma ;
  double beta[NDIM] ;
  int j ;
  
  alpha = 1./sqrt(-geom->gcon[TT][TT]) ;
  SLOOPA beta[j] = geom->gcon[TT][j]*alpha*alpha ;
  
  if( gamma_calc(pr,geom,&gamma) ) {
    fflush(stderr);
    fprintf(stderr,"\nucon_calc(): gamma failure \n");
    fflush(stderr);
    fail(FAIL_GAMMA);
  }
  
  ucon[TT] = gamma/alpha ;
  SLOOPA ucon[j] = pr[U1+j-1] - gamma*beta[j]/alpha ;
  
  return ;
}

void ut_calc_3vel(double *vcon, struct of_geom *geom, double *ut)
{
  double AA, BB, CC, DD, one_over_alpha_sq;
  //compute the Lorentz factor based on contravariant 3-velocity
  AA =     geom->gcov[TT][TT] ;
  BB = 2.*(geom->gcov[TT][1]*vcon[1] +
           geom->gcov[TT][2]*vcon[2] +
           geom->gcov[TT][3]*vcon[3]) ;
  CC = geom->gcov[1][1]*vcon[1]*vcon[1] +
  geom->gcov[2][2]*vcon[2]*vcon[2] +
  geom->gcov[3][3]*vcon[3]*vcon[3] +
  2.*(geom->gcov[1][2]*vcon[1]*vcon[2] +
      geom->gcov[1][3]*vcon[1]*vcon[3] +
      geom->gcov[2][3]*vcon[2]*vcon[3]) ;
  
  DD = -1./(AA+BB+CC);
  
  one_over_alpha_sq = -geom->gcon[TT][TT];

  if(DD<one_over_alpha_sq) {
    DD = one_over_alpha_sq;
  }
  
  *ut = sqrt(DD);

}

/* find gamma-factor wrt normal observer */
int gamma_calc(double *pr, struct of_geom *geom, double *gamma)
{
  double qsq ;
  
  qsq =     geom->gcov[1][1]*pr[U1]*pr[U1]
  + geom->gcov[2][2]*pr[U2]*pr[U2]
  + geom->gcov[3][3]*pr[U3]*pr[U3]
  + 2.*(geom->gcov[1][2]*pr[U1]*pr[U2]
        + geom->gcov[1][3]*pr[U1]*pr[U3]
        + geom->gcov[2][3]*pr[U2]*pr[U3]) ;
  
  if( qsq < 0. ){
    if( fabs(qsq) > 1.E-10 ){ // then assume not just machine precision
      fprintf(stderr,"gamma_calc():  failed: i,j,k,qsq = %d %d %d %28.18e \n", icurr,jcurr,kcurr,qsq);
      fprintf(stderr,"v[1-3] = %28.18e %28.18e %28.18e  \n",pr[U1],pr[U2],pr[U3]);
      *gamma = 1.;
      return (1);
    }
    else qsq=1.E-10; // set floor
  }
  
  *gamma = sqrt(1. + qsq) ;
  
  return(0) ;
}


/*
 * VCHAR():
 *
 * calculate components of magnetosonic velocity
 * corresponding to primitive variables p
 *
 * cfg 7-10-01
 *
 */

void vchar(double *pr, struct of_state *q, struct of_geom *geom, int js,
           double *vmax, double *vmin)
{
  double discr,vp,vm,bsq,EE,EF,va2,cs2,cms2,rho,u ;
  double Acov[NDIM],Bcov[NDIM],Acon[NDIM],Bcon[NDIM] ;
  double Asq,Bsq,Au,Bu,AB,Au2,Bu2,AuBu,A,B,C ;
  int j ;
  
  
  DLOOPA Acov[j] = 0. ;
  Acov[js] = 1. ;
  raise(Acov,geom,Acon) ;
  
  DLOOPA Bcov[j] = 0. ;
  Bcov[TT] = 1. ;
  raise(Bcov,geom,Bcon) ;
  
  /* find fast magnetosonic speed */
  bsq = dot(q->bcon,q->bcov) ;
  rho = pr[RHO] ;
  u = pr[UU] ;
  EF = rho + gam*u ;
  EE = bsq + EF ;
  va2 = bsq/EE ;
  cs2 = gam*(gam - 1.)*u/EF ;
  
  //	if(cs2 < 0.) cs2 = SMALL ;
  //	if(cs2 > 1.) cs2 = 1. ;
  //	if(va2 < 0.) va2 = SMALL ;
  //	if(va2 > 1.) va2 = 1. ;
  
  cms2 = cs2 + va2 - cs2*va2 ;	/* and there it is... */
  
  //cms2 *= 1.1 ;
  
  /* check on it! */
  if(cms2 < 0.) {
    fail(FAIL_COEFF_NEG) ;
    cms2 = SMALL ;
  }
  if(cms2 > 1.) {
    fail(FAIL_COEFF_SUP) ;
    cms2 = 1. ;
  }
  
  /* now require that speed of wave measured by observer
   q->ucon is cms2 */
  Asq = dot(Acon,Acov) ;
  Bsq = dot(Bcon,Bcov) ;
  Au =  dot(Acov,q->ucon) ;
  Bu =  dot(Bcov,q->ucon) ;
  AB =  dot(Acon,Bcov) ;
  Au2 = Au*Au ;
  Bu2 = Bu*Bu ;
  AuBu = Au*Bu ;
  
  A =      Bu2  - (Bsq + Bu2)*cms2 ;
  B = 2.*( AuBu - (AB + AuBu)*cms2 ) ;
  C =      Au2  - (Asq + Au2)*cms2 ;
  
  discr = B*B - 4.*A*C ;
  if((discr<0.0)&&(discr>-1.e-10)) discr=0.0;
  else if(discr < -1.e-10) {
    fprintf(stderr,"\n\t %g %g %g %g %g\n",A,B,C,discr,cms2) ;
    fprintf(stderr,"\n\t q->ucon: %g %g %g %g\n",q->ucon[0],q->ucon[1],
            q->ucon[2],q->ucon[3]) ;
    fprintf(stderr,"\n\t q->bcon: %g %g %g %g\n",q->bcon[0],q->bcon[1],
            q->bcon[2],q->bcon[3]) ;
    fprintf(stderr,"\n\t Acon: %g %g %g %g\n",Acon[0],Acon[1],
            Acon[2],Acon[3]) ;
    fprintf(stderr,"\n\t Bcon: %g %g %g %g\n",Bcon[0],Bcon[1],
            Bcon[2],Bcon[3]) ;
    fail(FAIL_VCHAR_DISCR) ;
    discr = 0. ;
  }
  
  discr = sqrt(discr) ;
  vp = -(-B + discr)/(2.*A) ;
  vm = -(-B - discr)/(2.*A) ;
  
  if(vp > vm) {
    *vmax = vp ;
    *vmin = vm ;
  }
  else {
    *vmax = vm ;
    *vmin = vp ;
  }
  
  return ;
}


/* Add any additional source terms (e.g. cooling functions) */

void misc_source(double *ph, double *phxp1, double *phxm1, double *phyp1, double *phym1,int ii, int jj, int kk, struct of_geom *geom,
                 struct of_state *q, double *dU, double Dt)
{
  
  double dudx1,dudx2;
  
  
  /* This is merely an example and does not represent any physical source term that I can think of */
  /* Place your calculation for the extra source terms here */
  //  dU[RHO]  += ph[RHO] ;
  //  dU[UU ]  += ph[UU ] ;
  //  dU[U1 ]  += ph[U1 ] ;
  //  dU[U2 ]  += ph[U2 ] ;
  //  dU[U3 ]  += ph[U3 ] ;
  
  
  
  
  
  
}
