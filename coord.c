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
#if(!DOCYLINDRIFYCOORDS)
      //use original coordinates
      vofx_gammiecoords( X, V );
#else
      //apply cylindrification to original coordinates
      //[this internally calls vofx_sjetcoords()]
      vofx_cylindrified( X, vofx_gammiecoords, V );
#endif
      
    }
  else if (BL == 2){
#if(!DOCYLINDRIFYCOORDS)
    //use original coordinates
    vofx_sjetcoords( X, V );
#else
    //apply cylindrification to original coordinates
    //[this internally calls vofx_sjetcoords()]
    vofx_cylindrified( X, vofx_sjetcoords, V );
#endif
  }
  
    else{
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
  else if(BL == 1 || BL == 2){
    
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

//////////////////////////////////////////////////////////////////////////////////////////
//
//  DISK-JET COORDINATES
//
//////////////////////////////////////////////////////////////////////////////////////////

void vofx_sjetcoords( double *X, double *V )
{
  //for SJETCOORDS
  double theexp;
  double Ftrgen( double x, double xa, double xb, double ya, double yb );
  double limlin( double x, double x0, double dx, double y0 );
  double minlin( double x, double x0, double dx, double y0 );
  double mins( double f1, double f2, double df );
  double maxs( double f1, double f2, double df );
  double thetaofx2(double x2, double ror0nu);
  double  fac, faker, ror0nu;
  double fakerdisk, fakerjet;
  double rbeforedisk, rinsidedisk, rinsidediskmax, rafterdisk;
  double ror0nudisk, ror0nujet, thetadisk, thetajet;
  
  V[0] = X[0];
  
  theexp = X[1];
  
  if( X[1] > x1br ) {
    theexp += cpow2 * pow(X[1]-x1br,npow2);
  }
  V[1] = R0+exp(theexp);
  
#if(0) //JON's method
  myhslope=2.0-Qjet*pow(V[1]/r1jet,-njet*(0.5+1.0/M_PI*atan(V[1]/r0grid-rsjet/r0grid)));
  
  if(X[2]<0.5){
    V[2] = M_PI * X[2] + ((1. - myhslope) / 2.) * mysin(2. * M_PI * X[2]);
  }
  else{
    V[2] = M_PI - (M_PI * (1.0-X[2])) + ((1. - myhslope) / 2.) * (-mysin(2. * M_PI * (1.0-X[2])));
  }
#elif(1) //SASHA's
  fac = Ftrgen( fabs(X[2]), global_fracdisk, 1-global_fracjet, 0, 1 );
  
  rbeforedisk = mins( V[1], global_r0disk, 0.5*global_r0disk );
  rinsidedisk = 1.;
  rinsidediskmax = 1.;
  rafterdisk = maxs( 1, 1 + (V[1]-global_rdiskend)*global_r0jet/(global_rjetend*global_r0disk*rinsidediskmax), 0.5*global_rdiskend*global_r0jet/(global_rjetend*global_r0disk*rinsidediskmax) );
  
  fakerdisk = rbeforedisk * rinsidedisk * rafterdisk;
  
  fakerjet = mins( V[1], global_r0jet, 0.5*global_r0jet ) * maxs( 1, V[1]/global_rjetend, 0.5 );
  
  ror0nudisk = pow( (fakerdisk - global_rsjet*Rin)/global_r0grid, global_jetnu/2 );
  ror0nujet = pow( (fakerjet - global_rsjet*Rin)/global_r0grid, global_jetnu/2 );
  thetadisk = thetaofx2( X[2], ror0nudisk );
  thetajet = thetaofx2( X[2], ror0nujet );
  V[2] = fac*thetajet + (1 - fac)*thetadisk;
#else
  V[2] = M_PI_2 * (1.0+ X[2]);
#endif
  
  // default is uniform \phi grid
  V[3]=X[3];
}

double thetaofx2(double x2, double ror0nu)
{
  double theta;
  if( x2 < -0.5 ) {
    theta = 0       + atan( tan((x2+1)*M_PI_2)/ror0nu );
  }
  else if( x2 >  0.5 ) {
    theta = M_PI    + atan( tan((x2-1)*M_PI_2)/ror0nu );
  }
  else {
    theta = M_PI_2 + atan( tan(x2*M_PI_2)*ror0nu );
  }
  return(theta);
}


//////////////////////////////////////////////////////////////////////////////////////////
//
//  CYLINDRIFICATION
//
//////////////////////////////////////////////////////////////////////////////////////////

//smooth step function:
// Ftr = 0 if x < 0, Ftr = 1 if x > 1 and smoothly interps. in btw.
double Ftr( double x )
{
  double res;
  
  if( x <= 0 ) {
    res = 0;
  }
  else if( x >= 1 ) {
    res = 1;
  }
  else {
    res = (64 + cos(5*M_PI*x) + 70*sin((M_PI*(-1 + 2*x))/2.) + 5*sin((3*M_PI*(-1 + 2*x))/2.))/128.;
  }
  
  return( res );
}

double Ftrgenlin( double x, double xa, double xb, double ya, double yb )
{
  double Ftr( double x );
  double res;
  
  res = (x*ya)/xa + (-((x*ya)/xa) + ((x - xb)*(1 - yb))/(1 - xb) + yb)*Ftr((x - xa)/(-xa + xb));
  
  return( res );
}

//goes from ya to yb as x goes from xa to xb
double Ftrgen( double x, double xa, double xb, double ya, double yb )
{
  double Ftr( double x );
  double res;
  
  res = ya + (yb-ya)*Ftr( (x-xa)/(xb-xa) );
  
  return( res );
}

double Fangle( double x )
{
  double res;
  
  if( x <= -1 ) {
    res = 0;
  }
  else if( x >= 1 ) {
    res = x;
  }
  else {
    res = (1 + x + (-140*sin((M_PI*(1 + x))/2.) + (10*sin((3*M_PI*(1 + x))/2.))/3. + (2*sin((5*M_PI*(1 + x))/2.))/5.)/(64.*M_PI))/2.;
  }
  
  return( res );
  
}

double limlin( double x, double x0, double dx, double y0 )
{
  double Fangle( double x );
  return( y0 - dx * Fangle(-(x-x0)/dx) );
}

double minlin( double x, double x0, double dx, double y0 )
{
  double Fangle( double x );
  return( y0 + dx * Fangle((x-x0)/dx) );
}

double mins( double f1, double f2, double df )
{
  double limlin( double x, double x0, double dx, double y0 );
  return( limlin(f1, f2, df, f2) );
}

double maxs( double f1, double f2, double df )
{
  double mins( double f1, double f2, double df );
  return( -mins(-f1, -f2, df) );
}

//=mins if dir < 0
//=maxs if dir >= 0
double minmaxs( double f1, double f2, double df, double dir )
{
  double mins( double f1, double f2, double df );
  double maxs( double f1, double f2, double df );
  if( dir>=0 ) {
    return( maxs(f1, f2, df) );
  }
  
  return( mins(f1, f2, df) );
}

static double sinth0( double *X0, double *X, void (*vofx)(double*, double*) );
static double sinth1in( double *X0, double *X, void (*vofx)(double*, double*) );
static double th2in( double *X0, double *X, void (*vofx)(double*, double*) );
static void to1stquadrant( double *Xin, double *Xout, int *ismirrored );
static double func1( double *X0, double *X,  void (*vofx)(double*, double*) );
static double func2( double *X0, double *X,  void (*vofx)(double*, double*) );

//Converts copies Xin to Xout and converts
//but sets Xout[2] to lie in the 1st quadrant, i.e. Xout[2] \in [-1,0])
//if the point had to be mirrored
void to1stquadrant( double *Xin, double *Xout, int *ismirrored )
{
  double ntimes;
  int j;
  
  DLOOPA Xout[j] = Xin[j];
  
  //bring the angle variables to -2..2 (for X) and -2pi..2pi (for V)
  ntimes = floor( (Xin[2]+2.0)/4.0 );
  //this forces -2 < Xout[2] < 2
  Xout[2] -= 4 * ntimes;
  
  *ismirrored = 0;
  
  if( Xout[2] > 0. ) {
    Xout[2] = -Xout[2];
    *ismirrored = 1-*ismirrored;
  }
  
  //now force -1 < Xout[2] < 0
  if( Xout[2] < -1. ) {
    Xout[2] = -2. - Xout[2];
    *ismirrored = 1-*ismirrored;
  }
}

double sinth0( double *X0, double *X, void (*vofx)(double*, double*) )
{
  double V0[NDIM];
  double Vc0[NDIM];
  double Xc0[NDIM];
  int j;
  
  //X1 = {0, X[1], X0[1], 0}
  DLOOPA Xc0[j] = X[j];
  Xc0[2] = X0[2];
  
  vofx( Xc0, Vc0 );
  vofx( X0, V0 );
  
  
  return( V0[1] * sin(V0[2]) / Vc0[1] );
}

double sinth1in( double *X0, double *X, void (*vofx)(double*, double*) )
{
  double V[NDIM];
  double V0[NDIM];
  double V0c[NDIM];
  double X0c[NDIM];
  int j;
  
  //X1 = {0, X[1], X0[1], 0}
  DLOOPA X0c[j] = X0[j];
  X0c[2] = X[2];
  
  vofx( X, V );
  vofx( X0c, V0c );
  vofx( X0, V0 );
  
  return( V0[1] * sin(V0c[2]) / V[1] );
}


double th2in( double *X0, double *X, void (*vofx)(double*, double*) )
{
  double V[NDIM];
  double V0[NDIM];
  double Vc0[NDIM];
  double Xc0[NDIM];
  double Xcmid[NDIM];
  double Vcmid[NDIM];
  int j;
  double res;
  double th0;
  
  DLOOPA Xc0[j] = X[j];
  Xc0[2] = X0[2];
  vofx( Xc0, Vc0 );
  
  DLOOPA Xcmid[j] = X[j];
  Xcmid[2] = 0;
  vofx( Xcmid, Vcmid );
  
  vofx( X0, V0 );
  vofx( X, V );
  
  th0 = asin( sinth0(X0, X, vofx) );
  
  res = (V[2] - Vc0[2])/(Vcmid[2] - Vc0[2]) * (Vcmid[2]-th0) + th0;
  
  return( res );
}

//Adjusts V[2]=theta so that a few innermost cells around the pole
//become cylindrical
//ASSUMES: poles are at
//            X[2] = -1 and +1, which correspond to
//            V[2] = 0 and pi
void vofx_cylindrified( double *Xin, void (*vofx)(double*, double*), double *Vout )
{
  double npiovertwos;
  double X[NDIM], V[NDIM];
  double Vin[NDIM];
  double X0[NDIM], V0[NDIM];
  double Xtr[NDIM], Vtr[NDIM];
  double f1, f2, dftr;
  double sinth, th;
  int j, ismirrored;
  
  vofx( Xin, Vin );
  
  // BRING INPUT TO 1ST QUADRANT:  X[2] \in [-1 and 0]
  to1stquadrant( Xin, X, &ismirrored );
  vofx( X, V );
  
  //initialize X0: cylindrify region
  //X[1] < X0[1] && X[2] < X0[2] (value of X0[3] not used)
  X0[0] = Xin[0];
  X0[1] = global_x10;
  X0[2] = global_x20;
  X0[3] = 0;
  vofx( X0, V0 );
  
  //{0, roughly midpoint between grid origin and x10, -1, 0}
  DLOOPA Xtr[j] = X[j];
  Xtr[1] = log( 0.5*( exp(X0[1])+exp(startx[1]) ) );   //always bound to be between startx[1] and X0[1]
  vofx( Xtr, Vtr );
  
  f1 = func1( X0, X, vofx );
  f2 = func2( X0, X, vofx );
  dftr = func2( X0, Xtr, vofx ) - func1( X0, Xtr, vofx );
  
  // Compute new theta
  sinth = maxs( V[1]*f1, V[1]*f2, Vtr[1]*fabs(dftr)+SMALL ) / V[1];
  
  th = asin( sinth );
  
  //initialize Vout with the original values
  DLOOPA Vout[j] = Vin[j];
  
  //apply change in theta in the original quadrant
  if( 0 == ismirrored ) {
    Vout[2] = Vin[2] + (th - V[2]);
  }
  else {
    //if mirrrored, flip the sign
    Vout[2] = Vin[2] - (th - V[2]);
  }
}

double func1( double *X0, double *X,  void (*vofx)(double*, double*) )
{
  double V[NDIM];
  
  vofx( X, V );
  
  return( sin(V[2]) );
}

double func2( double *X0, double *X,  void (*vofx)(double*, double*) )
{
  double V[NDIM];
  double Xca[NDIM];
  double func2;
  int j;
  double sth1in, sth2in, sth1inaxis, sth2inaxis;
  
  //{0, X[1], -1, 0}
  DLOOPA Xca[j] = X[j];
  Xca[2] = -1;
  
  vofx( X, V );
  
  sth1in = sinth1in( X0, X, vofx );
  sth2in = sin( th2in(X0, X, vofx) );
  
  sth1inaxis = sinth1in( X0, Xca, vofx );
  sth2inaxis = sin( th2in(X0, Xca, vofx) );
  
  func2 = minmaxs( sth1in, sth2in, fabs(sth2inaxis-sth1inaxis)+SMALL, X[1] - X0[1] );
  
  return( func2 );
}

