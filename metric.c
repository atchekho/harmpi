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

/***************************************************************************/
/***************************************************************************
    coord():
    -------
       -- given the indices i,j and location in the cell, return with 
          the values of X1,X2 there;  
       -- the locations are defined by : 
           -----------------------
           |                     |
           |                     |
           |FACE1   CENT         |
           |                     |
           |CORN    FACE2        |
           ----------------------
***************************************************************************/
void coord(int i, int j, int loc, double *X)
{
        if(loc == FACE1) {
                X[1] = startx[1] + i*dx[1] ;
                X[2] = startx[2] + (j + 0.5)*dx[2] ;
        }
        else if(loc == FACE2) {
                X[1] = startx[1] + (i + 0.5)*dx[1] ;
                X[2] = startx[2] + j*dx[2] ;
        }
        else if(loc == CENT) {
                X[1] = startx[1] + (i + 0.5)*dx[1] ;
                X[2] = startx[2] + (j + 0.5)*dx[2] ;
        }
        else {
                X[1] = startx[1] + i*dx[1] ;
                X[2] = startx[2] + j*dx[2] ;
        }

        return ;
}

/* assumes gcov has been set first; returns determinant */
double gdet_func(double gcov[][NDIM]) 
{
  int i,j,k;
  int permute[NDIM]; 
  double gcovtmp[NDIM][NDIM];
  double detg;

  for( i = 0 ; i < NDIM*NDIM ; i++ ) {  gcovtmp[0][i] = gcov[0][i]; }
  if( LU_decompose( gcovtmp,  permute ) != 0  ) { 
    fprintf(stderr, "gdet_func(): singular matrix encountered! \n");
    fail(FAIL_METRIC);
  }
  detg = 1.;
  DLOOPA detg *= gcovtmp[j][j];
  return( sqrt(fabs(detg)) );

}

/* invert gcov to get gcon */
void gcon_func(double gcov[][NDIM], double gcon[][NDIM])
{
  invert_matrix( gcov, gcon );
}

/***************************************************************************/
/***************************************************************************
  conn_func():
  -----------

   -- this gives the connection coefficient
	\Gamma^{i}_{j,k} = conn[..][i][j][k]
   --  where i = {1,2,3,4} corresponds to {t,r,theta,phi}

***************************************************************************/

/* Sets the spatial discretization in numerical derivatives : */
#define DELTA 1.e-5

/* NOTE: parameter hides global variable */
void conn_func(double *X, struct of_geom *geom, double conn[][NDIM][NDIM])
{
	int i,j,k,l ;
	double tmp[NDIM][NDIM][NDIM] ;
	double Xh[NDIM],Xl[NDIM] ;
	double gh[NDIM][NDIM] ;
	double gl[NDIM][NDIM] ;

	for(k=0;k<NDIM;k++) {
		for(l=0;l<NDIM;l++) Xh[l] = X[l] ;
		for(l=0;l<NDIM;l++) Xl[l] = X[l] ;
		Xh[k] += DELTA ;
		Xl[k] -= DELTA ;
		gcov_func(Xh,gh) ;
		gcov_func(Xl,gl) ;

		for(i=0;i<NDIM;i++)
		for(j=0;j<NDIM;j++) 
			conn[i][j][k] = (gh[i][j] - gl[i][j])/(Xh[k] - Xl[k]) ;
	}

	/* now rearrange to find \Gamma_{ijk} */
	for(i=0;i<NDIM;i++)
	for(j=0;j<NDIM;j++)
	for(k=0;k<NDIM;k++) 
		tmp[i][j][k] = 0.5*(conn[j][i][k] + conn[k][i][j] - conn[k][j][i]) ;

	/* finally, raise index */
	for(i=0;i<NDIM;i++)
	for(j=0;j<NDIM;j++)
	for(k=0;k<NDIM;k++)  {
		conn[i][j][k] = 0. ;
		for(l=0;l<NDIM;l++) conn[i][j][k] += geom->gcon[i][l]*tmp[l][j][k] ;
	}

	/* done! */
}

/* Lowers a contravariant rank-1 tensor to a covariant one */
void lower(double *ucon, struct of_geom *geom, double *ucov)
{

	ucov[0] = geom->gcov[0][0]*ucon[0] 
		+ geom->gcov[0][1]*ucon[1] 
		+ geom->gcov[0][2]*ucon[2] 
		+ geom->gcov[0][3]*ucon[3] ;
	ucov[1] = geom->gcov[1][0]*ucon[0] 
		+ geom->gcov[1][1]*ucon[1] 
		+ geom->gcov[1][2]*ucon[2] 
		+ geom->gcov[1][3]*ucon[3] ;
	ucov[2] = geom->gcov[2][0]*ucon[0] 
		+ geom->gcov[2][1]*ucon[1] 
		+ geom->gcov[2][2]*ucon[2] 
		+ geom->gcov[2][3]*ucon[3] ;
	ucov[3] = geom->gcov[3][0]*ucon[0] 
		+ geom->gcov[3][1]*ucon[1] 
		+ geom->gcov[3][2]*ucon[2] 
		+ geom->gcov[3][3]*ucon[3] ;

        return ;
}

/* Raises a covariant rank-1 tensor to a contravariant one */
void raise(double *ucov, struct of_geom *geom, double *ucon)
{

	ucon[0] = geom->gcon[0][0]*ucov[0] 
		+ geom->gcon[0][1]*ucov[1] 
		+ geom->gcon[0][2]*ucov[2] 
		+ geom->gcon[0][3]*ucov[3] ;
	ucon[1] = geom->gcon[1][0]*ucov[0] 
		+ geom->gcon[1][1]*ucov[1] 
		+ geom->gcon[1][2]*ucov[2] 
		+ geom->gcon[1][3]*ucov[3] ;
	ucon[2] = geom->gcon[2][0]*ucov[0] 
		+ geom->gcon[2][1]*ucov[1] 
		+ geom->gcon[2][2]*ucov[2] 
		+ geom->gcon[2][3]*ucov[3] ;
	ucon[3] = geom->gcon[3][0]*ucov[0] 
		+ geom->gcon[3][1]*ucov[1] 
		+ geom->gcon[3][2]*ucov[2] 
		+ geom->gcon[3][3]*ucov[3] ;

        return ;
}

/* load local geometry into structure geom */
void get_geometry(int ii, int jj, int kk, struct of_geom *geom)
{
	int j,k ;

	//-new DLOOP geom->gcov[j][k] = gcov[ii][jj][kk][j][k] ;
	//-new DLOOP geom->gcon[j][k] = gcon[ii][jj][kk][j][k] ;
	for(j=0;j<=NDIM*NDIM-1;j++){
	  geom->gcon[0][j] = gcon[ii][jj][kk][0][j];
	  geom->gcov[0][j] = gcov[ii][jj][kk][0][j];
	}
	geom->g = gdet[ii][jj][kk] ;
	icurr = ii ;
	jcurr = jj ;
	pcurr = kk ;
}

#undef DELTA

/* Minkowski metric; signature +2 */
double mink(int i, int j)
{
	if(i == j) {
		if(i == 0) return(-1.) ;
		else return(1.) ;
	}
	else return(0.) ;
}

/* Boyer-Lindquist ("bl") metric functions */
void blgset(int i, int j, struct of_geom *geom)
{
	double r,th,X[NDIM] ;

	coord(i,j,CENT,X) ;
	bl_coord(X,&r,&th) ;

	if(th < 0) th *= -1. ;
	if(th > M_PI) th = 2.*M_PI - th ;

	geom->g = bl_gdet_func(r,th) ;
	bl_gcov_func(r,th,geom->gcov) ;
	bl_gcon_func(r,th,geom->gcon) ;
}

double bl_gdet_func(double r, double th)
{
	double a2,r2 ;

	a2 = a*a ;
	r2 = r*r ;
	return( 
		r*r*fabs(sin(th))*(1. + 0.5*(a2/r2)*(1. + cos(2.*th)))
	) ;
}

void bl_gcov_func(double r, double th, double gcov[][NDIM])
{
	int j,k ;
	double sth,cth,s2,a2,r2,DD,mu ;

	DLOOP gcov[j][k] = 0. ;

	sth = fabs(sin(th)) ;
	s2 = sth*sth ;
	cth = cos(th) ;
	a2 = a*a ;
	r2 = r*r ;
	DD = 1. - 2./r + a2/r2 ;
	mu = 1. + a2*cth*cth/r2 ;
	
	gcov[TT][TT] = -(1. - 2./(r*mu)) ;
	gcov[TT][3] = -2.*a*s2/(r*mu) ;
	gcov[3][TT] = gcov[TT][3] ;
	gcov[1][1] = mu/DD ;
	gcov[2][2] = r2*mu ;
	gcov[3][3] = r2*sth*sth*(1. + a2/r2 + 2.*a2*s2/(r2*r*mu)) ;

}

void bl_gcon_func(double r, double th, double gcon[][NDIM])
{
	int j,k ;
	double sth,cth,a2,r2,r3,DD,mu ;

	DLOOP gcon[j][k] = 0. ;

	sth = sin(th) ;
	cth = cos(th) ;

#if(COORDSINGFIX)
	if (fabs(sth) < SINGSMALL) {
	  if(sth>=0) sth=SINGSMALL;
	  if(sth<0) sth=-SINGSMALL;
	}
#endif

	a2 = a*a ;
	r2 = r*r ;
	r3 = r2*r ;
	DD = 1. - 2./r + a2/r2 ;
	mu = 1. + a2*cth*cth/r2 ;

	gcon[TT][TT] = -1. - 2.*(1. + a2/r2)/(r*DD*mu) ;
	gcon[TT][3] = -2.*a/(r3*DD*mu) ;
	gcon[3][TT] = gcon[TT][3] ;
	gcon[1][1] = DD/mu ;
	gcon[2][2] = 1./(r2*mu) ;
	gcon[3][3] = (1. - 2./(r*mu))/(r2*sth*sth*DD) ;


}

