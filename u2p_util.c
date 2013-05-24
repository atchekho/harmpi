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



#include "u2p_defs.h"


void primtoU_g( FTYPE prim[], FTYPE gcov[][NDIM], FTYPE gcon[][NDIM], FTYPE gdet, 
		FTYPE U[] );
void ucon_calc_g(FTYPE prim[],FTYPE gcov[][NDIM],FTYPE gcon[][4],FTYPE ucon[]);
void raise_g(FTYPE vcov[], FTYPE gcon[][NDIM], FTYPE vcon[]);
void lower_g(FTYPE vcon[], FTYPE gcov[][NDIM], FTYPE vcov[]);
void ncov_calc(FTYPE gcon[][NDIM],FTYPE ncov[]) ;
void bcon_calc_g(FTYPE prim[],FTYPE ucon[],FTYPE ucov[],FTYPE ncov[],FTYPE bcon[]); 
FTYPE pressure_rho0_u(FTYPE rho0, FTYPE u);
FTYPE pressure_rho0_w(FTYPE rho0, FTYPE w);


/**********************************************************************/
/******************************************************************
   primtoU_g(): 

       -- calculates the conserved variables from the primitive variables 
            and the metric;
       -- assumes that the conserved and primitive variables are defined ala HARM:

              /  rho u^t           \
         U =  |  T^t_\mu + rho u^t |  sqrt(-det(g_{\mu\nu}))
              \   B^i              /

             /    rho        \
	 P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /

**************************************************************************/

void primtoU_g(
	       FTYPE prim[NPR],       /* primitive variables */
	       FTYPE gcov[NDIM][NDIM],    /* covariant (index dn) form of metric */
	       FTYPE gcon[NDIM][NDIM],    /* contravariant (index up) form of metric */
	       FTYPE gdet,                /* sqrt of -1 times det(g_{\mu \nu}) */
	       FTYPE U[NPR]           /* matrix of derivatives */
	       ) {
  int i,j ;
  FTYPE rho0 ;
  static FTYPE ucon[NDIM],ucov[NDIM],bcon[NDIM],bcov[NDIM],ncov[NDIM] ;
  FTYPE gamma,n_dot_b,bsq,u,p,w, alpha ;

	
  /* Calculate auxiliary quantities: */
  alpha = 1.0/sqrt(-gcon[0][0]);

  ucon_calc_g(prim,gcov,gcon,ucon) ;
  lower_g(ucon,gcov,ucov) ;
  ncov_calc(gcon,ncov) ;

  gamma = -ncov[0]*ucon[0] ;

  bcon_calc_g(prim,ucon,ucov,ncov,bcon) ;
  lower_g(bcon,gcov,bcov) ;

  n_dot_b = 0. ;
  for(i=0;i<4;i++) n_dot_b += ncov[i]*bcon[i] ;
  bsq = 0. ;
  for(i=0;i<4;i++) bsq += bcov[i]*bcon[i] ;

  rho0 = prim[RHO] ;
  u = prim[UU] ;
  p = pressure_rho0_u(rho0,u) ;
  w = rho0 + u + p ;

  // Now set the conserved variables themselves, using HARM's definition:
  U[RHO] = ucon[0]*rho0 ;

  for( i = 0; i < 4; i++) {
    U[QCOV0+i] = gamma*(w + bsq)*ucov[i] 
      - (p + bsq/2.)*ncov[i] 
      + n_dot_b*bcov[i] ;

    U[QCOV0+i] /= alpha;
  }

  U[QCOV0] = U[QCOV0] + U[RHO];
  U[BCON1] = prim[BCON1] ;
  U[BCON2] = prim[BCON2] ;
  U[BCON3] = prim[BCON3] ;

  for(i = 0; i < NPR; i++ ) {
    U[i] *= gdet;
  }

  return ;
}

/**********************************************************************/
/******************************************************************

  ucon_calc_g(): 
    
       -- calculates the contravariant (up) components of the four-velocity
          given the primitive variables, of which the velocity is 
          \tilde{u}^i = \gamma v^j  where v^j is the velocity of the 
          flow w.r.t a normal observer to the coordinates;

       -- also requires the metric and inverse metric;

       -- assumes:

             /    rho        \
	 P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /

******************************************************************/
void ucon_calc_g(FTYPE prim[NPR],FTYPE gcov[NDIM][NDIM],FTYPE gcon[NDIM][NDIM],
		 FTYPE ucon[NDIM])
{
  FTYPE u_tilde_con[4] ;
  FTYPE u_tilde_sq ;
  FTYPE gamma,lapse ;
  int i,j ;
	
  u_tilde_con[0] = 0. ;
  u_tilde_con[1] = prim[UTCON1] ;
  u_tilde_con[2] = prim[UTCON2] ;
  u_tilde_con[3] = prim[UTCON3] ;

  u_tilde_sq = 0. ;
  for(i=0;i<NDIM;i++)
    for(j=0;j<NDIM;j++)
      u_tilde_sq += gcov[i][j]*u_tilde_con[i]*u_tilde_con[j] ;
  u_tilde_sq = fabs(u_tilde_sq) ;

  gamma = sqrt(1. + u_tilde_sq) ;

  lapse = sqrt(-1./gcon[0][0]) ;

  for(i=0;i<NDIM;i++) ucon[i] = u_tilde_con[i] - lapse*gamma*gcon[0][i] ;

  return ;
}

/**********************************************************************/
/******************************************************************
   
    raise_g():
 
         -- calculates the contravariant form of a covariant tensor, 
            using the inverse of the metric;

******************************************************************/
void raise_g(FTYPE vcov[NDIM], FTYPE gcon[NDIM][NDIM], FTYPE vcon[NDIM])
{
  int i,j;

  for(i=0;i<NDIM;i++) {
    vcon[i] = 0. ;
    for(j=0;j<NDIM;j++) 
      vcon[i] += gcon[i][j]*vcov[j] ;
  }

  return ;
}

/**********************************************************************/
/******************************************************************

     lower_g():
  
          -- calculates the ocvariant form of a contravariant tensor 
             using the metric;

******************************************************************/
void lower_g(FTYPE vcon[NDIM], FTYPE gcov[NDIM][NDIM], FTYPE vcov[NDIM])
{
  int i,j;

  for(i=0;i<NDIM;i++) {
    vcov[i] = 0. ;
    for(j=0;j<NDIM;j++) 
      vcov[i] += gcov[i][j]*vcon[j] ;
  }

  return ;
}

/**********************************************************************/
/******************************************************************

     ncov_calc(): 

         -- calculates the covariant form of the normal vector to our 
            spacelike hypersurfaces ala the ADM formalism.

         -- requires the inverse metric;

******************************************************************/
void ncov_calc(FTYPE gcon[NDIM][NDIM],FTYPE ncov[NDIM]) 
{
  FTYPE lapse ;
  int i;

  lapse = sqrt(-1./gcon[0][0]) ;

  ncov[0] = -lapse ;
  for( i = 1; i < NDIM; i++) { 
    ncov[i] = 0. ;
  }

  return ;
}

/**********************************************************************/
/******************************************************************

    bcon_calc_g(): 
  
        -- using the primitive variables, contra-/co-variant 4-vel., 
           and covariant normal vector, calculate the contravariant 
           form of the magnetic 4-vector b^\mu (the small "b" in HARM);
       -- assumes:

             /    rho        \
	 P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /


******************************************************************/
void bcon_calc_g(FTYPE prim[NPR],FTYPE ucon[NDIM],FTYPE ucov[NDIM],
		 FTYPE ncov[NDIM],FTYPE bcon[NDIM]) 
{
  static FTYPE Bcon[NDIM] ;
  FTYPE u_dot_B ;
  FTYPE gamma ;
  int i ;

  // Bcon = \mathcal{B}^\mu  of the paper:
  Bcon[0] = 0. ;
  for(i=1;i<NDIM;i++) Bcon[i] = -ncov[0] * prim[BCON1+i-1] ;

  u_dot_B = 0. ;
  for(i=0;i<NDIM;i++) u_dot_B += ucov[i]*Bcon[i] ;

  gamma = -ucon[0]*ncov[0] ;
  for(i=0;i<NDIM;i++) bcon[i] = (Bcon[i] + ucon[i]*u_dot_B)/gamma ;
}


/**********************************************************************/
/******************************************************************

    gamma_calc_g(): 
  
        -- using the primitive variables, contra-/co-variant 4-vel., 
           and covariant normal vector, calculate the contravariant 
           form of the magnetic 4-vector b^\mu (the small "b" in HARM);

       -- assumes:

             /    rho        \
	 P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /

******************************************************************/
int gamma_calc_g(FTYPE *pr, FTYPE gcov[NDIM][NDIM], FTYPE *gamma)
{
  FTYPE utsq ;
  int j,k;

  utsq =    gcov[1][1]*pr[UTCON1]*pr[UTCON1]
    + gcov[2][2]*pr[UTCON2]*pr[UTCON2]
    + gcov[3][3]*pr[UTCON3]*pr[UTCON3]
    + 2.*(  gcov[1][2]*pr[UTCON1]*pr[UTCON2]
	    + gcov[1][3]*pr[UTCON1]*pr[UTCON3]
	    + gcov[2][3]*pr[UTCON2]*pr[UTCON3]) ;

  if(utsq<0.0){
    if(fabs(utsq)>1.E-10){ // then assume not just machine precision
      return (1);
    }
    else utsq=fabs(utsq); // set floor
  }

  *gamma = sqrt(1. + utsq) ;

  return(0) ;
}


/**************************************************/
/**************************************************
  The following functions assume a Gamma-law EOS:
***************************************************/

/* 
pressure as a function of rho0 and u 
this is used by primtoU and Utoprim_?D
*/
FTYPE pressure_rho0_u(FTYPE rho0, FTYPE u)
{
  return((GAMMA - 1.)*u) ;
}


  
/* 
pressure as a function of rho0 and w = rho0 + u + p 
this is used by primtoU and Utoprim_1D
*/
FTYPE pressure_rho0_w(FTYPE rho0, FTYPE w)
{
  return((GAMMA-1.)*(w - rho0)/GAMMA) ;
}


