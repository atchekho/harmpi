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

/** 
 *
 * this file contains all the coordinate dependent
 * parts of the code, except the initial and boundary
 * conditions 
 *
 **/

/* should return boyer-lindquist coordinte of point */
void bl_coord(double *X, double *r, double *th)
{

	*r = exp(X[1]) + R0 ;
	*th = M_PI*X[2] + ((1. - hslope)/2.)*sin(2.*M_PI*X[2]) ;

	// avoid singularity at polar axis
#if(COORDSINGFIX)
	if(fabs(*th)<SINGSMALL){
	  if((*th)>=0) *th =  SINGSMALL;
	  if((*th)<0)  *th = -SINGSMALL;
	}
	if(fabs(M_PI - (*th)) < SINGSMALL){
	  if((*th)>=M_PI) *th = M_PI+SINGSMALL;
	  if((*th)<M_PI)  *th = M_PI-SINGSMALL;
	}
#endif
	
	return ;
}

/* insert metric here */
void gcov_func(double *X, double gcov[][NDIM])
{
        int j,k ;
        double sth,cth,s2,rho2 ;
	double r,th ;
	double tfac,rfac,hfac,pfac ;

        DLOOP gcov[j][k] = 0. ;

	bl_coord(X,&r,&th) ;

        cth = cos(th) ;

	//-orig        sth = fabs(sin(th)) ;
        //-orig        if(sth < SMALL) sth = SMALL ;

        sth = sin(th) ;

        s2 = sth*sth ;
        rho2 = r*r + a*a*cth*cth ;

	tfac = 1. ;
	rfac = r - R0 ;
	hfac = M_PI + (1. - hslope)*M_PI*cos(2.*M_PI*X[2]) ;
	pfac = 1. ;

        gcov[TT][TT] = (-1. + 2.*r/rho2)      * tfac*tfac ;
        gcov[TT][1] = (2.*r/rho2)             * tfac*rfac ;
        gcov[TT][3] = (-2.*a*r*s2/rho2)       * tfac*pfac ;

        gcov[1][TT] = gcov[TT][1] ;
        gcov[1][1] = (1. + 2.*r/rho2)         * rfac*rfac ;
        gcov[1][3] = (-a*s2*(1. + 2.*r/rho2)) * rfac*pfac ;

        gcov[2][2] = rho2                     * hfac*hfac ;

        gcov[3][TT] = gcov[TT][3] ;
        gcov[3][1]  = gcov[1][3] ;
        gcov[3][3] = s2*(rho2 + a*a*s2*(1. + 2.*r/rho2)) * pfac*pfac ;
}

/* some grid location, dxs */
void set_points()
{
        startx[1] = log(Rin - R0) ;
        startx[2] = 0. ;
        dx[1] = log((Rout - R0)/(Rin - R0))/N1 ;
        dx[2] = 1./N2 ;
}

void fix_flux(double F1[][N2+4][NPR], double F2[][N2+4][NPR])
{
	int i,j,k ;

	for(i=-1;i<N1;i++) {

		F1[i][-1][B2] = -F1[i][0][B2] ;
		F1[i][N2][B2] = -F1[i][N2-1][B2] ;

		PLOOP F2[i][N2][k] = 0. ;
		PLOOP F2[i][0][k] = 0. ;
	}

	for(j=-1;j<N2;j++) {
	  if(F1[0 ][j][RHO]>0) F1[0 ][j][RHO]=0;
	  if(F1[N1][j][RHO]<0) F1[N1][j][RHO]=0;
	}

	return;
}


void rescale(double *pr, int which, int dir, int ii, int jj, int face, struct of_geom *geom)
{
	double scale[NPR], r, th, X[NDIM];
	int k ;

	coord(ii, jj, face, X);
	bl_coord(X, &r, &th);

	if (dir == 1) {
		// optimized for pole
		scale[RHO] = pow(r, 1.5);
		scale[UU] = scale[RHO] * r;
		scale[U1] = scale[RHO];
		scale[U2] = 1.0;
		scale[U3] = r * r;
		scale[B1] = r * r;
		scale[B2] = r * r;
		scale[B3] = r * r;
	} else if (dir == 2) {
		scale[RHO] = 1.0;
		scale[UU]  = 1.0;
		scale[U1] = 1.0;
		scale[U2] = 1.0;
		scale[U3] = 1.0;
		scale[B1] = 1.0;
		scale[B2] = 1.0;
		scale[B3] = 1.0;
	} 

	if (which == FORWARD) {	// rescale before interpolation
		PLOOP pr[k] *= scale[k];
	} else if (which == REVERSE) {	// unrescale after interpolation
		PLOOP pr[k] /= scale[k];
	} else {
		fprintf(stderr,"no such rescale type!\n") ;
		exit(100) ;
	}
}
