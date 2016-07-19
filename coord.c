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

/** 
 *
 * this file contains all the coordinate dependent
 * parts of the code, except the initial and boundary
 * conditions 
 *
 **/

static void vofx_gammiecoords(double *X, double *V);  //original gammie coordinates
static void vofx_sjetcoords( double *X, double *V );  //original coordinates
static void vofx_cylindrified( double *Xin, void (*vofx)(double*, double*), double *Vout ); //coordinate "cylindrifier"

/* should return boyer-lindquist coordinte of point */
void bl_coord_vec(double *X, double *V)
{

    if (BL == 1){
      //use original coordinates
      vofx_gammiecoords( X, V );
    }
    else{
        //use uniform coordinates
        V[0] = X[0];
        V[1] = X[1];
        V[2] = X[2];
        V[3] = X[3];
    }

	// avoid singularity at polar axis
#if(COORDSINGFIX && BL)
	if(fabs(V[2])<SINGSMALL){
	  if((V[2])>=0) V[2] =  SINGSMALL;
	  if((V[2])<0)  V[2] = -SINGSMALL;
	}
	if(fabs(M_PI - V[2]) < SINGSMALL){
	  if((V[2])>=M_PI) V[2] = M_PI+SINGSMALL;
	  if((V[2])<M_PI)  V[2] = M_PI-SINGSMALL;
	}
#endif
	
	return ;
}

//in this version, -1<=X[2]<=1 maps onto 0<=theta<=pi
//(in the original version, 0<=X[2]<=1 maps onto 0<=theta<=pi
void vofx_gammiecoords(double *X, double *V)
{
  double theexp = X[1];
  
  if( X[1] > x1br ) {
    theexp += cpow2 * pow(X[1]-x1br,npow2);  //hyperexponential for X[1] > x1br
    
  }
  
  V[0] = X[0];
  V[1] = exp(theexp) + R0 ;
  V[2] = M_PI_2*(1.0+X[2]) + ((1. - hslope)/2.)*sin(M_PI*(1.0+X[2])) ;
  V[3] = X[3] ;

}
void bl_coord(double *X, double *r, double *th, double *phi)
{
  double V[NDIM];
  bl_coord_vec(X,V);
  *r = V[1];
  *th = V[2];
  *phi = V[3];
  return ;
}


/* insert metric here */
void gcov_func(double *X, double gcovp[][NDIM])
{
  int i,j,k,l ;
  double sth,cth,s2,rho2 ;
  double r,th,phi ;
  double tfac,rfac,hfac,pfac ;
  double gcov[NDIM][NDIM];
  double dxdxp[NDIM][NDIM];
  
  DLOOP gcov[j][k] = 0. ;
  DLOOP gcovp[j][k] = 0.;
  
  bl_coord(X,&r,&th,&phi) ;
  
  cth = cos(th) ;
  
  //-orig        sth = fabs(sin(th)) ;
  //-orig        if(sth < SMALL) sth = SMALL ;
  
  sth = sin(th) ;
  
  s2 = sth*sth ;
  rho2 = r*r + a*a*cth*cth ;
  
  //    tfac = 1. ;
  //    if (r<rbr){
  //	rfac = r - R0 ;
  //    }
  //    else{
  //        rfac = (r-R0)*(1.+npow2*cpow2*pow((-x1br+X[1]),npow2-1.));
  //    }
  //    hfac = M_PI + (1. - hslope)*M_PI*cos(2.*M_PI*X[2]) ;
  //    pfac = 1. ;
  
  //compute Jacobian x1,x2,x3 -> r,th,phi
  dxdxp_func(X, dxdxp);
  
  if(GR==1){
    gcov[TT][TT] = (-1. + 2.*r/rho2);
    gcov[TT][1] = (2.*r/rho2);
    gcov[TT][3] = (-2.*a*r*s2/rho2);
    
    gcov[1][TT] = gcov[TT][1] ;
    gcov[1][1] = (1. + 2.*r/rho2);
    gcov[1][3] = (-a*s2*(1. + 2.*r/rho2));
    
    gcov[2][2] = rho2;
    
    gcov[3][TT] = gcov[TT][3] ;
    gcov[3][1]  = gcov[1][3] ;
    gcov[3][3] = s2*(rho2 + a*a*s2*(1. + 2.*r/rho2));
    
    //convert to code coordinates
    for(i=0;i<NDIM;i++){
      for(j=0;j<NDIM;j++){
        gcovp[i][j] = 0;
        for(k=0;k<NDIM;k++) {
          for(l=0;l<NDIM;l++){
            gcovp[i][j] += gcov[k][l]*dxdxp[k][i]*dxdxp[l][j];
          }
        }
      }
    }
  }
  else {
    gcovp[TT][TT] = -1. ;
    gcovp[1][1] = 1. ;
    gcovp[2][2] = 1. ;
    gcovp[3][3] = 1. ;
  }
}

/* some grid location, dxs */
void set_points()
{
  int j, iter;
  int dim;
  double lenx[NDIM];
  double x1max, x1max0, dxmax;
  
  //first, define the full range of coordinates, then modify them for MPI case
  if (BL == 0){
    startx[1] = 0. ;
    startx[2] = 0. ;
    startx[3] = 0. ;
    lenx[1] = 1.;
    lenx[2] = 1.;
    lenx[3] = 1.;
  }
  else if(BL == 1){
    
    const double RELACC = 1e-14;
    const int ITERMAX = 50;
    
    x1br = log( rbr - R0 );
    
    if( Rout < rbr ) {
      x1max = log(Rout-R0);
    }
    else {
      x1max0 = 1;
      x1max = 2;
      
      //find the root via iterations
      for( iter = 0; iter < ITERMAX; iter++ ) {
        if( fabs((x1max - x1max0)/x1max) < RELACC ) {
          break;
        }
        x1max0 = x1max;
        dxmax= (pow( (log(Rout-R0) - x1max0)/cpow2, 1./npow2 ) + x1br) - x1max0;
        
        // need a slight damping factor
        double dampingfactor=0.5;
        x1max = x1max0 + dampingfactor*dxmax;
        if (x1max> log(Rout-R0)){x1max = log(Rout-R0);}
      }
      
      if( iter == ITERMAX ) {
        if(mpi_rank==MASTER) {
          printf( "Error: iteration procedure for finding x1max has not converged: x1max = %g, dx1max/x1max = %g, iter = %d\n",
                 x1max, (x1max-x1max0)/x1max, iter );
          printf( "Error: iteration procedure for finding x1max has not converged: rbr= %g, x1br = %g, log(Rout-R0) = %g\n",
                 rbr, x1br, log(Rout-R0) );
        }
        exit(1);
      }
      else {
        if(mpi_rank==MASTER) printf( "x1max = %g (dx1max/x1max = %g, itno = %d)\n", x1max, (x1max-x1max0)/x1max, iter );
      }
    }
    startx[1] = log(Rin - R0) ;  //minimum values
    startx[2] = -1.+(1.-fractheta) ;   //minimum values
    startx[3] = 0. ;   //minimum values
    lenx[1] = x1max - startx[1];
    lenx[2] = 2.*fractheta;
    lenx[3] = 2.*M_PI*fracphi;
  }
  
  dx[1] = lenx[1]/mpi_ntot[1] ; //radial
  dx[2] = lenx[2]/mpi_ntot[2] ; //theta
  dx[3] = lenx[3]/mpi_ntot[3] ; //phi

}

void fix_flux(double F1[][N2M][N3M][NPR], double F2[][N2M][N3M][NPR], double F3[][N2M][N3M][NPR])
{
    int i,j,k,m ;
    //NOTE: MPI checks below always true if MPI is disabled
    if(BL){
//#pragma omp parallel for schedule(static,MY_MAX((N1+1)*(N3+1)/nthreads,1)) collapse(2) default(none) shared(F1,F2,nthreads) private(i,k,m)
        for(i=-D1;i<N1+D1;i++) {
          for(k=-D3;k<N3+D3;k++) {

            if(mpi_coords[2] == 0) {
              F1[i][-1][k][B2] = -F1[i][0][k][B2] ;
              F3[i][-1][k][B2] = -F3[i][0][k][B2] ;
            }
            if(mpi_coords[2] == mpi_dims[2]-1) {
              F1[i][N2][k][B2] = -F1[i][N2-1][k][B2] ;
              F3[i][N2][k][B2] = -F3[i][N2-1][k][B2] ;
            }

            if(INFLOW ==  0){
                if(mpi_coords[2] == mpi_dims[2]-1) PLOOP F2[i][N2][k][m] = 0. ;
                if(mpi_coords[2] == 0) PLOOP F2[i][0][k][m] = 0. ;
            }
          }
        }
    }

    
    if(INFLOW == 0 ){
//#pragma omp parallel for schedule(static,MY_MAX((N2+1)*(N3+1)/nthreads,1)) collapse(2) default(none) shared(F1,nthreads) private(j,k)
        for(j=-D2;j<N2+D2;j++) {
          for(k=-D3;k<N3+D3;k++) {
            if(mpi_coords[1] == 0) if(F1[0 ][j][k][RHO]>0) F1[0 ][j][k][RHO]=0;
            if(mpi_coords[1] == mpi_dims[1]-1) if(F1[N1][j][k][RHO]<0) F1[N1][j][k][RHO]=0;
          }
        }
    }

    return;
}


void rescale(double *pr, int which, int dir, int ii, int jj, int kk, int face, struct of_geom *geom)
{
  double scale[NPR], r, th, phi, X[NDIM];
  int m ;
  
  coord(ii, jj, kk, face, X);
  bl_coord(X, &r, &th, &phi);
  
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
  else if (dir == 3) {
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
    PLOOP pr[m] *= scale[m];
  } else if (which == REVERSE) {	// unrescale after interpolation
    PLOOP pr[m] /= scale[m];
  } else {
    if(MASTER==mpi_rank) fprintf(stderr,"no such rescale type!\n") ;
    exit(100) ;
  }
}

