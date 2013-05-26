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

/* bound array containing entire set of primitive variables */

void bound_prim( double prim[][N2+4][NPR] )
{
        int i,j,k,jref ;
	void inflow_check(double *pr, int ii, int jj, int type );
        struct of_geom geom ;

        /* inner r boundary condition: u, gdet extrapolation */
        for(j=0;j<N2;j++) {
#if( RESCALE )
		get_geometry(0,j,CENT,&geom) ;
		rescale(prim[0][j],FORWARD, 1, 0,j,CENT,&geom) ;
#endif
                PLOOP prim[-1][j][k] = prim[0][j][k] ;
                PLOOP prim[-2][j][k] = prim[0][j][k] ;
                pflag[-1][j] = pflag[0][j] ;
                pflag[-2][j] = pflag[0][j] ;

#if( RESCALE )
		get_geometry(0,j,CENT,&geom) ;
		rescale(prim[0][j],REVERSE, 1, 0,j,CENT,&geom) ;
		get_geometry(-1,j,CENT,&geom) ;
		rescale(prim[-1][j],REVERSE, 1, -1,j,CENT,&geom) ;
		get_geometry(-2,j,CENT,&geom) ;
		rescale(prim[-2][j],REVERSE, 1, -2,j,CENT,&geom) ;
#endif

        }

        /* outer r BC: outflow */
        for(j=0;j<N2;j++) {
#if( RESCALE )
		get_geometry(N1-1,j,CENT,&geom) ;
		rescale(prim[N1-1][j],FORWARD, 1, N1-1,j,CENT,&geom) ;
#endif

		PLOOP prim[N1  ][j][k] = prim[N1-1][j][k] ;
		PLOOP prim[N1+1][j][k] = prim[N1-1][j][k] ;
		pflag[N1  ][j] = pflag[N1-1][j] ;
		pflag[N1+1][j] = pflag[N1-1][j] ;

#if( RESCALE )
		get_geometry(N1-1,j,CENT,&geom) ;
		rescale(prim[N1-1][j],REVERSE, 1, N1-1,j,CENT,&geom) ;
		get_geometry(N1,j,CENT,&geom) ;
		rescale(prim[N1][j],REVERSE, 1, N1,j,CENT,&geom) ;
		get_geometry(N1+1,j,CENT,&geom) ;
		rescale(prim[N1+1][j],REVERSE, 1, N1+1,j,CENT,&geom) ;
#endif
        }

        /* make sure there is no inflow at the inner boundary */
	for(i=-2;i<=-1;i++)  for(j=-2;j<N2+2;j++) { 
	  inflow_check(prim[-1][j],i,j,0) ;
	  inflow_check(prim[-2][j],i,j,0) ;
	}

        /* make sure there is no inflow at the outer boundary */
        for(i=N1;i<=N1+1;i++)  for(j=-2;j<N2+2;j++) { 
	  inflow_check(prim[N1  ][j],i,j,1) ;
	  inflow_check(prim[N1+1][j],i,j,1) ;
	}

#if(POLEFIX)
	//copy all densities and B^phi in; interpolate linearly transverse velocity
	jref = POLEFIX;
        for(i=-2;i<N1+2;i++) {
	  for(j=0;j<jref;j++) {
	    PLOOP {
	      if(k==B1 || k==B2) 
		//don't touch magnetic fields
		continue;
	      else if(k==U2) {
		//linear interpolation of transverse velocity (both poles)
		prim[i][j][k] = (j+0.5)/(jref+0.5) * prim[i][jref][k];
		prim[i][N2-1-j][k] = (j+0.5)/(jref+0.5) * prim[i][N2-1-jref][k];
	      }
	      else {
		//everything else copy (both poles)
		prim[i][j][k] = prim[i][jref][k];
		prim[i][N2-1-j][k] = prim[i][N2-1-jref][k];
	      }
	    }
	  }
	}	
#endif
        /* polar BCs */
        for(i=-2;i<=N1+1;i++) { 
	  PLOOP {
                prim[i][-1  ][k] = prim[i][0   ][k] ;
                prim[i][-2  ][k] = prim[i][1   ][k] ;
                prim[i][N2  ][k] = prim[i][N2-1][k] ;
                prim[i][N2+1][k] = prim[i][N2-2][k] ;
	  }
	  pflag[i][-1  ] = pflag[i][0   ] ;
	  pflag[i][-2  ] = pflag[i][1   ] ;
	  pflag[i][N2  ] = pflag[i][N2-1] ;
	  pflag[i][N2+1] = pflag[i][N2-2] ;
        }

        /* make sure b and u are antisymmetric at the poles */
        for(i=-2;i<N1+2;i++) {
                for(j=-2;j<0;j++) {
                        prim[i][j][U2] *= -1. ;
                        prim[i][j][B2] *= -1. ;
                }
                for(j=N2;j<N2+2;j++) {
                        prim[i][j][U2] *= -1. ;
                        prim[i][j][B2] *= -1. ;
                }
        }
}

void inflow_check(double *pr, int ii, int jj, int type )
{
        struct of_geom geom ;
        double ucon[NDIM] ;
        int j,k ;
        double alpha,beta1,gamma,vsq ;

        get_geometry(ii,jj,CENT,&geom) ;
        ucon_calc(pr, &geom, ucon) ;

        if( ((ucon[1] > 0.) && (type==0)) || ((ucon[1] < 0.) && (type==1)) ) { 
                /* find gamma and remove it from primitives */
	  if( gamma_calc(pr,&geom,&gamma) ) { 
	    fflush(stderr);
	    fprintf(stderr,"\ninflow_check(): gamma failure \n");
	    fflush(stderr);
	    fail(FAIL_GAMMA);
	  }
	    pr[U1] /= gamma ;
	    pr[U2] /= gamma ;
	    pr[U3] /= gamma ;
	    alpha = 1./sqrt(-geom.gcon[0][0]) ;
	    beta1 = geom.gcon[0][1]*alpha*alpha ;

	    /* reset radial velocity so radial 4-velocity
	     * is zero */
	    pr[U1] = beta1/alpha ;

	    /* now find new gamma and put it back in */
	    vsq = 0. ;
	    SLOOP vsq += geom.gcov[j][k]*pr[U1+j-1]*pr[U1+k-1] ;
	    if( fabs(vsq) < 1.e-13 )  vsq = 1.e-13;
	    if( vsq >= 1. ) { 
	      vsq = 1. - 1./(GAMMAMAX*GAMMAMAX) ;
	    }
	    gamma = 1./sqrt(1. - vsq) ;
	    pr[U1] *= gamma ;
	    pr[U2] *= gamma ;
	    pr[U3] *= gamma ;

	    /* done */
	  }
	  else
	    return ;

	}

