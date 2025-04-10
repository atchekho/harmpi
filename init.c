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

/*
 *
 * generates initial conditions for a fishbone & moncrief disk 
 * with exterior at minimum values for density & internal energy.
 *
 * cfg 8-10-01
 *
 */

#include "decs.h"
#include <float.h>


void coord_transform(double *pr,int i, int j,int k) ;
double compute_Amax( double (*A)[N2+D2][N3+D3] );
double compute_B_from_A( double (*A)[N2+D2][N3+D3], double (*p)[N2M][N3M][NPR] );
double normalize_B_by_maxima_ratio(double beta_target, double (*p)[N2M][N3M][NPR], double *norm_value);
double normalize_B_by_beta(double beta_target, double (*p)[N2M][N3M][NPR], double rmax, double *norm_value);

/////////////////////
//magnetic field geometry and normalization
#define NORMALFIELD (0)

#define WHICHFIELD NORMALFIELD

#define NORMALIZE_FIELD_BY_MAX_RATIO (1)
#define NORMALIZE_FIELD_BY_BETAMIN (2)
#define WHICH_FIELD_NORMALIZATION NORMALIZE_FIELD_BY_BETAMIN
//end magnetic field
//////////////////////

//////////////////////
//torus density normalization
#define THINTORUS_NORMALIZE_DENSITY (1)
#define DOAUTOCOMPUTEENK0 (1)

#define NORMALIZE_BY_TORUS_MASS (1)
#define NORMALIZE_BY_DENSITY_MAX (2)

#define DENSITY_NORMALIZATION NORMALIZE_BY_DENSITY_MAX
//torus density normalization
//////////////////////

double rmax = 0.;
double rhomax = 1.;

void init()
{
  void init_bondi(void);
  void init_torus(void);
  void init_sndwave(void);
  void init_entwave(void);
  void init_monopole(double Rout_val);

  switch( WHICHPROBLEM ) {
  case MONOPOLE_PROBLEM_1D:
  case MONOPOLE_PROBLEM_2D:
    init_monopole(1e3);
    break;
  case BZ_MONOPOLE_2D:
    init_monopole(100.);
    break;
  case TORUS_PROBLEM:
    init_torus();
    break;
  case SNDWAVE_TEST :
    init_sndwave() ;
    break ;
  case ENTWAVE_TEST :
    init_entwave() ;
    break ;
  case BONDI_PROBLEM_1D:
  case BONDI_PROBLEM_2D:
    init_bondi();
    break;
  }

}

void init_torus()
{
  int i,j,k ;
  double r,th,phi,sth,cth ;
  double ur,uh,up,u,rho ;
  double X[NDIM] ;
  struct of_geom geom ;

  /* for disk interior */
  double l,rin,lnh,expm2chi,up1 ;
  double DD,AA,SS,thin,sthin,cthin,DDin,AAin,SSin ;
  double kappa,hm1 ;

  /* for magnetic field */
  double A[N1+D1][N2+D2][N3+D3] ;
  double rho_av,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ;
  double lfish_calc(double rmax) ;
  
  int iglob, jglob, kglob;
  double rancval;
  
  double amax, aphipow;
  
  double Amin, Amax, cutoff_frac = 0.01;
  
  /* some physics parameters */
  gam = 5./3. ;

  /* disk parameters (use fishbone.m to select new solutions) */
  a = 0.9 ;
  rin = 6. ;
  rmax = 13. ;
  l = lfish_calc(rmax) ;

  kappa =1.e-3;
  beta = 100. ;

  /* some numerical parameters */
  lim = MC ;
  failed = 0 ;	/* start slow */
  cour = .8 ;
  dt = 1.e-5 ;
  R0 = 0.0 ;
  Rin = 0.87*(1. + sqrt(1. - a*a)) ;  //.98
  Rout = 1e5;
  rbr = 400.;
  npow2=4.0; //power exponent
  cpow2=1.0; //exponent prefactor (the larger it is, the more hyperexponentiation is)


  t = 0. ;
  hslope = 0.3 ;

  if(N2!=1) {
    //2D problem, use full pi-wedge in theta
    fractheta = 1.;
  }
  else{
    //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
    fractheta = 1.e-2;
  }
  
  fracphi = 1.;

  //cylindrification parameters
  global_x10 = 3.5;  //radial distance in MCOORD until which the innermost angular cell is cylinrdical
  global_x20 = -1. + 1./mpi_ntot[2];     //This restricts grid cylindrification to the one
  //single grid cell closest to the pole (other cells virtually unaffeced, so there evolution is accurate).
  //This trick minimizes the resulting pole deresolution and relaxes the time step.
  //The innermost grid cell is evolved inaccurately whether you resolve it or not, and it will be fixed
  //by POLEFIX (see bounds.c).
  
  set_arrays() ;
  set_grid() ;

  get_phys_coord(5,0,0,&r,&th,&phi) ;
  if(MASTER==mpi_rank) {
    fprintf(stderr,"r[5]: %g\n",r) ;
    fprintf(stderr,"r[5]/rhor: %g",r/(1. + sqrt(1. - a*a))) ;
    if( r > 1. + sqrt(1. - a*a) ) {
      fprintf(stderr, ": INSUFFICIENT RESOLUTION, ADD MORE CELLS INSIDE THE HORIZON\n" );
    }
    else {
      fprintf(stderr, "\n");
    }
  }

  /* output choices */
  tf = 10000.0 ;

  DTd = 10.; /* dumping frequency, in units of M */
  DTl = 10. ;	/* logfile frequency, in units of M */
  DTi = 10. ; 	/* image file frequ., in units of M */
  DTr = 10. ; /* restart file frequ., in units of M */
  DTr01 = 100. ; /* restart file frequ., in timesteps */

  /* start diagnostic counters */
  dump_cnt = 0 ;
  image_cnt = 0 ;
  rdump_cnt = 0 ;
  rdump01_cnt = 0 ;
  defcon = 1. ;

  rhomax = 0. ;
  umax = 0. ;
  //ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
  for(iglob=0;iglob<mpi_ntot[1];iglob++) {
    for(jglob=0;jglob<mpi_ntot[2];jglob++) {
      for(kglob=0;kglob<mpi_ntot[3];kglob++) {
        
        rancval = ranc(0);
        i = iglob-mpi_startn[1];
        j = jglob-mpi_startn[2];
        k = kglob-mpi_startn[3];
        if(i<0 ||
           j<0 ||
           k<0 ||
           i>=N1 ||
           j>=N2 ||
           k>=N3){
          continue;
        }
        get_phys_coord(i,j,k,&r,&th,&phi) ;

        sth = sin(th) ;
        cth = cos(th) ;

        /* calculate lnh */
        DD = r*r - 2.*r + a*a ;
        AA = (r*r + a*a)*(r*r + a*a) - DD*a*a*sth*sth ;
        SS = r*r + a*a*cth*cth ;

        thin = M_PI/2. ;
        sthin = sin(thin) ;
        cthin = cos(thin) ;
        DDin = rin*rin - 2.*rin + a*a ;
        AAin = (rin*rin + a*a)*(rin*rin + a*a) 
                - DDin*a*a*sthin*sthin ;
        SSin = rin*rin + a*a*cthin*cthin ;

        if(r >= rin) {
          lnh = 0.5*log((1. + sqrt(1. + 4.*(l*l*SS*SS)*DD/
                  (AA*sth*AA*sth)))/(SS*DD/AA)) 
                  - 0.5*sqrt(1. + 4.*(l*l*SS*SS)*DD/(AA*AA*sth*sth))
                  - 2.*a*r*l/AA 
                  - (0.5*log((1. + sqrt(1. + 4.*(l*l*SSin*SSin)*DDin/
                  (AAin*AAin*sthin*sthin)))/(SSin*DDin/AAin)) 
                  - 0.5*sqrt(1. + 4.*(l*l*SSin*SSin)*DDin/
                          (AAin*AAin*sthin*sthin)) 
                  - 2.*a*rin*l/AAin ) ;
        }
        else
          lnh = 1. ;


        /* regions outside torus */
        if(lnh < 0. || r < rin) {
          //reset density and internal energy to zero outside torus
          rho = 0.; //1.e-7*RHOMIN ;
          u = 0.; //1.e-7*UUMIN ;

          /* these values are demonstrably physical
             for all values of a and r */
          /*
          ur = -1./(r*r) ;
          uh = 0. ;
          up = 0. ;
          */

          ur = 0. ;
          uh = 0. ;
          up = 0. ;

          /*
          get_geometry(i,j,CENT,&geom) ;
          ur = geom.gcon[0][1]/geom.gcon[0][0] ;
          uh = geom.gcon[0][2]/geom.gcon[0][0] ;
          up = geom.gcon[0][3]/geom.gcon[0][0] ;
          */

          p[i][j][k][RHO] = rho ;
          p[i][j][k][UU] = u ;
          p[i][j][k][U1] = ur ;
          p[i][j][k][U2] = uh ;
          p[i][j][k][U3] = up ;
        }
        /* region inside magnetized torus; u^i is calculated in
         * Boyer-Lindquist coordinates, as per Fishbone & Moncrief,
         * so it needs to be transformed at the end */
        else { 
          hm1 = exp(lnh) - 1. ;
          rho = pow(hm1*(gam - 1.)/(kappa*gam),
                                  1./(gam - 1.)) ; 
          u = kappa*pow(rho,gam)/(gam - 1.) ;
          ur = 0. ;
          uh = 0. ;

          /* calculate u^phi */
          expm2chi = SS*SS*DD/(AA*AA*sth*sth) ;
          up1 = sqrt((-1. + sqrt(1. + 4.*l*l*expm2chi))/2.) ;
          up = 2.*a*r*sqrt(1. + up1*up1)/sqrt(AA*SS*DD) +
                  sqrt(SS/AA)*up1/sth ;


          p[i][j][k][RHO] = rho ;
          if(rho > rhomax) rhomax = rho ;
          p[i][j][k][UU] = u*(1. + 4.e-2*(rancval-0.5)) ;
          if(u > umax && r > rin) umax = u ;
          p[i][j][k][U1] = ur ;
          p[i][j][k][U2] = uh ;

          p[i][j][k][U3] = up ;

          /* convert from 4-vel in BL coords to relative 4-vel in code coords */
          coord_transform(p[i][j][k],i,j,k) ;
        }

        p[i][j][k][B1] = 0. ;
        p[i][j][k][B2] = 0. ;
        p[i][j][k][B3] = 0. ;
      }
    }
  }

#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&rhomax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&umax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
  
  /* Normalize the densities so that max(rho) = 1 */
  if(MASTER==mpi_rank) fprintf(stderr,"rhomax: %g\n",rhomax) ;
  ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
          p[i][j][k][RHO] /= rhomax ;
          p[i][j][k][UU]  /= rhomax ;
  }
  umax /= rhomax ;
  kappa *= pow(rhomax,gam-1);
  global_kappa = kappa;
  rhomax = 1. ;

  bound_prim(p) ;

  if (WHICHFIELD == NORMALFIELD) {
    aphipow = 0.;
  }
  else {
    fprintf(stderr, "Unknown field type: %d\n", (int)WHICHFIELD);
    exit(321);
  }

  /* first find corner-centered vector potential */
  ZSLOOP(0,N1-1+D1,0,N2-1+D2,0,N3-1+D3) A[i][j][k] = 0. ;
  ZSLOOP(0,N1-1+D1,0,N2-1+D2,0,N3-1+D3) {
          /* radial field version */
          /*
          coord(i,j,k,CORN,X) ;
          bl_coord(X,&r,&th,&phi) ;
         
          A[i][j][k] = (1-cos(th)) ;
          */

    
          /* vertical field version */
          /*
          coord(i,j,k,CORN,X) ;
          bl_coord(X,&r,&th,&phi) ;

          A[i][j][k] = r*r*sin(th)*sin(th) ;
          */
    
    

          /* field-in-disk version */
          /* flux_ct */
    
      //cannot use get_phys_coords() here because it can only provide coords at CENT
      coord(i,j,k,CORN,X) ;
      bl_coord(X,&r,&th,&phi) ;


          rho_av = 0.25*(
                  p[i][j][k][RHO] +
                  p[i-1][j][k][RHO] +
                  p[i][j-1][k][RHO] +
                  p[i-1][j-1][k][RHO]) ;

          q = pow(r,aphipow)*rho_av/rhomax ;
          if (WHICHFIELD == NORMALFIELD) {
            q -= 0.2;
          }
          if(q > 0.) A[i][j][k] = q ;

  }
  
  fixup(p) ;

  /* now differentiate to find cell-centered B,
     and begin normalization */
  
  bsq_max = compute_B_from_A(A,p);
  
  if(WHICHFIELD == NORMALFIELD) {
    if(MASTER==mpi_rank)
      fprintf(stderr,"initial bsq_max: %g\n",bsq_max) ;

    /* finally, normalize to set field strength */
    beta_act =(gam - 1.)*umax/(0.5*bsq_max) ;

    if(MASTER==mpi_rank)
      fprintf(stderr,"initial beta: %g (should be %g)\n",beta_act,beta) ;
    
    if(WHICH_FIELD_NORMALIZATION == NORMALIZE_FIELD_BY_BETAMIN) {
      beta_act = normalize_B_by_beta(beta, p, 10*rmax, &norm);
    }
    else if(WHICH_FIELD_NORMALIZATION == NORMALIZE_FIELD_BY_MAX_RATIO) {
      beta_act = normalize_B_by_maxima_ratio(beta, p, &norm);
    }
    else {
      if(i_am_the_master) {
        fprintf(stderr, "Unknown magnetic field normalization %d\n",
                WHICH_FIELD_NORMALIZATION);
        #ifdef MPI
        MPI_Finalize();
	#endif
        exit(2345);
      }
    }

    if(MASTER==mpi_rank)
      fprintf(stderr,"final beta: %g (should be %g); normalization factor: %g\n",beta_act,beta,norm) ;
  }

    
  /* enforce boundary conditions */
  fixup(p) ;
  bound_prim(p) ;




#if( DO_FONT_FIX )
  set_Katm();
#endif 


}

//note that only axisymmetric A is supported
double compute_Amax( double (*A)[N2+D2][N3+D3] )
{
  double Amax = 0.;
  int i, j, k;
  struct of_geom geom;
  
  ZLOOP {
    if(A[i][j][k] > Amax) Amax = A[i][j][k];
  }
  
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&Amax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
  
  return(Amax);
}


//note that only axisymmetric A is supported
double compute_B_from_A( double (*A)[N2+D2][N3+D3], double (*p)[N2M][N3M][NPR] )
{
  double bsq_max = 0., bsq_ij ;
  int i, j, k;
  struct of_geom geom;
  
  ZLOOP {
    get_geometry(i,j,k,CENT,&geom) ;
    
    /* flux-ct */
    p[i][j][k][B1] = -(A[i][j][k] - A[i][j+1][k]
                       + A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[2]*geom.g) ;
    p[i][j][k][B2] = (A[i][j][k] + A[i][j+1][k]
                      - A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[1]*geom.g) ;
    
    p[i][j][k][B3] = 0. ;
    
    bsq_ij = bsq_calc(p[i][j][k],&geom) ;
    if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
  }
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&bsq_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif

  return(bsq_max);
}

double normalize_B_by_maxima_ratio(double beta_target, double (*p)[N2M][N3M][NPR], double *norm_value)
{
  double beta_act, bsq_ij, u_ij, umax = 0., bsq_max = 0.;
  double norm;
  int i, j, k;
  struct of_geom geom;
  
  ZLOOP {
    get_geometry(i,j,k,CENT,&geom) ;
    bsq_ij = bsq_calc(p[i][j][k],&geom) ;
    if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
    u_ij = p[i][j][k][UU];
    if(u_ij > umax) umax = u_ij;
  }
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&umax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&bsq_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif

  /* finally, normalize to set field strength */
  beta_act =(gam - 1.)*umax/(0.5*bsq_max) ;
  
  norm = sqrt(beta_act/beta_target) ;
  bsq_max = 0. ;
  ZLOOP {
    p[i][j][k][B1] *= norm ;
    p[i][j][k][B2] *= norm ;
    
    get_geometry(i,j,k,CENT,&geom) ;
    bsq_ij = bsq_calc(p[i][j][k],&geom) ;
    if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
  }
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&bsq_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
  
  beta_act = (gam - 1.)*umax/(0.5*bsq_max) ;

  if(norm_value) {
    *norm_value = norm;
  }
  return(beta_act);
}

//normalize the magnetic field using the values inside r < rmax
double normalize_B_by_beta(double beta_target, double (*p)[N2M][N3M][NPR], double rmax, double *norm_value)
{
  double beta_min = 1e100, beta_ij, beta_act, bsq_ij, u_ij, umax = 0., bsq_max = 0.;
  double norm;
  int i, j, k;
  struct of_geom geom;
  double X[NDIM], r, th, ph;
  
  ZLOOP {
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th, &ph);
    if (r>rmax) {
      continue;
    }
    get_geometry(i,j,k,CENT,&geom) ;
    bsq_ij = bsq_calc(p[i][j][k],&geom) ;
    u_ij = p[i][j][k][UU];
    beta_ij = (gam - 1.)*u_ij/(0.5*(bsq_ij+SMALL)) ;
    if(beta_ij < beta_min) beta_min = beta_ij ;
  }
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&beta_min,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif
  
  /* finally, normalize to set field strength */
  beta_act = beta_min;
  
  norm = sqrt(beta_act/beta_target) ;
  beta_min = 1e100;
  ZLOOP {
    p[i][j][k][B1] *= norm ;
    p[i][j][k][B2] *= norm ;
    p[i][j][k][B3] *= norm ;
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th, &ph);
    get_geometry(i,j,k,CENT,&geom) ;
    bsq_ij = bsq_calc(p[i][j][k],&geom) ;
    u_ij = p[i][j][k][UU];
    beta_ij = (gam - 1.)*u_ij/(0.5*(bsq_ij+SMALL)) ;
    if(r<rmax && beta_ij < beta_min) beta_min = beta_ij ;
  }
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&beta_min,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif
  
  beta_act = beta_min;

  if(norm_value) {
    *norm_value = norm;
  }

  return(beta_act);
}


void init_bondi()
{
	int i,j,k ;
	double r,th,phi,sth,cth ;
	double ur,uh,up,u,rho ;
	double X[NDIM] ;
	struct of_geom geom ;
	double rhor;

	/* for disk interior */
	double l,rin,lnh,expm2chi,up1 ;
	double DD,AA,SS,thin,sthin,cthin,DDin,AAin,SSin ;
	double kappa,hm1 ;

	/* for magnetic field */
	double A[N1+1][N2+1][N3+1] ;
	double rho_av,rhomax,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ;
	double rmax, lfish_calc(double rmax) ;

	/* some physics parameters */
	gam = 4./3. ;

	/* black hole parameters */
        a = 0.9375 ;

	kappa = 1.e-3 ;

	/* radius of the inner edge of the initial density distribution */
	rin = 10.;

        /* some numerical parameters */
        lim = MC ;
        failed = 0 ;	/* start slow */
        cour = 0.9 ;
        dt = 1.e-5 ;
	rhor = (1. + sqrt(1. - a*a)) ;
	R0 = -2*rhor ;
        Rin = 0.5*rhor ;
        Rout = 1e3 ;
        rbr = Rout*10.;
        npow2=4.0; //power exponent
        cpow2=1.0; //exponent prefactor (the larger it is, the more hyperexponentiation is)


        t = 0. ;
        hslope = 1.0 ; //uniform angular grid

	if(N2!=1) {
	  //2D problem, use full pi-wedge in theta
	  fractheta = 1.;
	}
	else{
	  //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
	  fractheta = 1.e-2;
	}
        fracphi = 1.;

        set_arrays() ;
        set_grid() ;

	coord(-2,0,0,CENT,X) ;
	bl_coord(X,&r,&th,&phi) ;
	fprintf(stderr,"rmin: %g\n",r) ;
	fprintf(stderr,"rmin/rm: %g\n",r/(1. + sqrt(1. - a*a))) ;

        /* output choices */
	tf = Rout ;

	DTd = 2. ;	/* dumping frequency, in units of M */
	DTl = 2. ;	/* logfile frequency, in units of M */
	DTi = 2. ; 	/* image file frequ., in units of M */
        DTr = 50 ; /* restart file frequ., in units of M */
        DTr01 = 1000 ; /* restart file frequ., in timesteps */

	/* start diagnostic counters */
	dump_cnt = 0 ;
	image_cnt = 0 ;
	rdump_cnt = 0 ;
        rdump01_cnt = 0 ;
	defcon = 1. ;

	rhomax = 0. ;
	umax = 0. ;
	ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
		coord(i,j,k,CENT,X) ;
		bl_coord(X,&r,&th,&phi) ;

		sth = sin(th) ;
		cth = cos(th) ;

		/* regions outside uniform density distribution */
		if(r < rin) {
			rho = 1.e-7*RHOMIN ;
                        u = 1.e-7*UUMIN ;

			/* these values are demonstrably physical
			   for all values of a and r */
			/*
                        ur = -1./(r*r) ;
                        uh = 0. ;
			up = 0. ;
			*/

			ur = 0. ;
			uh = 0. ;
			up = 0. ;

			/*
			get_geometry(i,j,CENT,&geom) ;
                        ur = geom.gcon[0][1]/geom.gcon[0][0] ;
                        uh = geom.gcon[0][2]/geom.gcon[0][0] ;
                        up = geom.gcon[0][3]/geom.gcon[0][0] ;
			*/

			p[i][j][k][RHO] = rho ;
			p[i][j][k][UU] = u ;
			p[i][j][k][U1] = ur ;
			p[i][j][k][U2] = uh ;
			p[i][j][k][U3] = up ;
		}
		/* region inside initial uniform density */
		else { 
		  rho = 1.;
		  u = kappa*pow(rho,gam)/(gam - 1.) ;
		  ur = 0. ;
		  uh = 0. ;


		  p[i][j][k][RHO] = rho ;
		  if(rho > rhomax) rhomax = rho ;
		  p[i][j][k][UU] = u;
		  if(u > umax && r > rin) umax = u ;
		  p[i][j][k][U1] = ur ;
		  p[i][j][k][U2] = uh ;
		  p[i][j][k][U3] = up ;
		  
		  /* convert from 4-vel to 3-vel */
		  coord_transform(p[i][j][k],i,j,k) ;
		}

		p[i][j][k][B1] = 0. ;
		p[i][j][k][B2] = 0. ;
		p[i][j][k][B3] = 0. ;

	}

	fixup(p) ;
	bound_prim(p) ;
    

#if(0) //disable for now
	/* first find corner-centered vector potential */
	ZSLOOP(0,N1,0,N2,0,N3) A[i][j][k] = 0. ;
        ZSLOOP(0,N1,0,N2,0,N3) {
                /* vertical field version */
                /*
                coord(i,j,l,CORN,X) ;
                bl_coord(X,&r,&th,&phi) ;

                A[i][j][k] = 0.5*r*sin(th) ;
                */

                /* field-in-disk version */
		/* flux_ct */
                rho_av = 0.25*(
                        p[i][j][RHO] +
                        p[i-1][j][RHO] +
                        p[i][j-1][RHO] +
                        p[i-1][j-1][RHO]) ;

                q = rho_av/rhomax - 0.2 ;
                if(q > 0.) A[i][j][k] = q ;

        }

	/* now differentiate to find cell-centered B,
	   and begin normalization */
	bsq_max = 0. ;
	ZLOOP {
		get_geometry(i,j,k,CENT,&geom) ;

		/* flux-ct */
		p[i][j][B1] = -(A[i][j][k] - A[i][j+1][k]
				+ A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[2]*geom.g) ;
		p[i][j][B2] = (A[i][j][k] + A[i][j+1][k]
				- A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[1]*geom.g) ;

		p[i][j][B3] = 0. ;

		bsq_ij = bsq_calc(p[i][j][k],&geom) ;
		if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
	}
	fprintf(stderr,"initial bsq_max: %g\n",bsq_max) ;

	/* finally, normalize to set field strength */
	beta_act = (gam - 1.)*umax/(0.5*bsq_max) ;
	fprintf(stderr,"initial beta: %g (should be %g)\n",beta_act,beta) ;
	norm = sqrt(beta_act/beta) ;
	bsq_max = 0. ;
	ZLOOP {
		p[i][j][k][B1] *= norm ;
		p[i][j][k][B2] *= norm ;

		get_geometry(i,j,k,CENT,&geom) ;
		bsq_ij = bsq_calc(p[i][j][k],&geom) ;
		if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
	}
	beta_act = (gam - 1.)*umax/(0.5*bsq_max) ;
	fprintf(stderr,"final beta: %g (should be %g)\n",beta_act,beta) ;

	/* enforce boundary conditions */
	fixup(p) ;
	bound_prim(p) ;
    
#endif

    

    
#if( DO_FONT_FIX ) 
	set_Katm();
#endif 


}

void init_monopole(double Rout_val)
{
	int i,j,k ;
	double r,th,phi,sth,cth ;
	double ur,uh,up,u,rho ;
	double X[NDIM] ;
	struct of_geom geom ;
	double rhor;

	/* for disk interior */
	double l,rin,lnh,expm2chi,up1 ;
	double DD,AA,SS,thin,sthin,cthin,DDin,AAin,SSin ;
	double kappa,hm1 ;

	/* for magnetic field */
	double A[N1+1][N2+1] ;
	double rho_av,rhomax,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ;
	double rmax, lfish_calc(double rmax) ;

	/* some physics parameters */
	gam = 4./3. ;

	/* disk parameters (use fishbone.m to select new solutions) */
        a = 0.9375 ;
        rin = 6. ;
        rmax = 12. ;
        l = lfish_calc(rmax) ;

	kappa = 1.e-3 ;
	beta = 1.e2 ;

        /* some numerical parameters */
        lim = MC ;
        failed = 0 ;	/* start slow */
        cour = 0.9 ;
        dt = 1.e-5 ;
	rhor = (1. + sqrt(1. - a*a)) ;
	R0 = -4*rhor;
        Rin = 0.7*rhor ;
        Rout = Rout_val ;
        rbr = Rout*10.;
    npow2=4.0; //power exponent
    cpow2=1.0; //exponent prefactor (the larger it is, the more hyperexponentiation is)

        t = 0. ;
        hslope = 1. ;

	if(N2!=1) {
	  //2D problem, use full pi-wedge in theta
	  fractheta = 1.;
	}
	else{
	  //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
	  fractheta = 1.e-2;
	}

        fracphi = 1.;

        set_arrays() ;
        set_grid() ;

	coord(-2,0,0,CENT,X) ;
	bl_coord(X,&r,&th,&phi) ;
	fprintf(stderr,"rmin: %g\n",r) ;
	fprintf(stderr,"rmin/rm: %g\n",r/(1. + sqrt(1. - a*a))) ;

        /* output choices */
	tf = 2*Rout ;

	DTd = 1. ;	/* dumping frequency, in units of M */
	DTl = 50. ;	/* logfile frequency, in units of M */
	DTi = 50. ; 	/* image file frequ., in units of M */
        DTr = 1. ; /* restart file frequ., in units of M */
	DTr01 = 1000 ; 	/* restart file frequ., in timesteps */

	/* start diagnostic counters */
	dump_cnt = 0 ;
	image_cnt = 0 ;
	rdump_cnt = 0 ;
        rdump01_cnt = 0 ;
	defcon = 1. ;

	rhomax = 0. ;
	umax = 0. ;
	ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
	  coord(i,j,k,CENT,X) ;
	  bl_coord(X,&r,&th,&phi) ;

	  sth = sin(th) ;
	  cth = cos(th) ;

	  /* rho = 1.e-7*RHOMIN ; */
	  /* u = 1.e-7*UUMIN ; */

	  /* rho = pow(r,-4.)/BSQORHOMAX; */
	  /* u = pow(r,-4.*gam)/BSQOUMAX; */

	  rho = RHOMINLIMIT+(r/10./rhor)/pow(r,4)/BSQORHOMAX;
	  u = UUMINLIMIT+(r/10./rhor)/pow(r,4)/BSQORHOMAX;

	  /* these values are demonstrably physical
	     for all values of a and r */
	  /*
	    ur = -1./(r*r) ;
	    uh = 0. ;
	    up = 0. ;
	  */

	  ur = 0. ;
	  uh = 0. ;
	  up = 0. ;

	  /*
	    get_geometry(i,j,CENT,&geom) ;
	    ur = geom.gcon[0][1]/geom.gcon[0][0] ;
	    uh = geom.gcon[0][2]/geom.gcon[0][0] ;
	    up = geom.gcon[0][3]/geom.gcon[0][0] ;
	  */

	  p[i][j][k][RHO] = rho ;
	  p[i][j][k][UU] = u ;
	  p[i][j][k][U1] = ur ;
	  p[i][j][k][U2] = uh ;
	  p[i][j][k][U3] = up ;
	  p[i][j][k][B1] = 0. ;
	  p[i][j][k][B2] = 0. ;
	  p[i][j][k][B3] = 0. ;
	}

	rhomax = 1. ;
	fixup(p) ;
	bound_prim(p) ;

        //leave A[][] a 2D array for now, which means that magnetic field will be axisymmetric
	/* first find corner-centered vector potential */
	ZSLOOP(0,N1,0,N2,0,0) A[i][j] = 0. ;
        ZSLOOP(0,N1,0,N2,0,0) {
#if(0)
                /* vertical field version */
                coord(i,j,k,CORN,X) ;
                bl_coord(X,&r,&th,&phi) ;
                A[i][j] = 0.5*pow(r*sin(th),2);
#elif(1)
                /* radial (monopolar) field version */
                coord(i,j,k,CORN,X) ;
                bl_coord(X,&r,&th,&phi) ;
                A[i][j] = (1-cos(th)) ;
#endif

        }

	/* now differentiate to find cell-centered B,
	   and begin normalization */
	bsq_max = 0. ;
	ZLOOP {
		get_geometry(i,j,k,CENT,&geom) ;

		/* flux-ct */
		p[i][j][k][B1] = -(A[i][j] - A[i][j+1]
				+ A[i+1][j] - A[i+1][j+1])/(2.*dx[2]*geom.g) ;
		p[i][j][k][B2] = (A[i][j] + A[i][j+1]
				- A[i+1][j] - A[i+1][j+1])/(2.*dx[1]*geom.g) ;

		p[i][j][k][B3] = 0. ;

		bsq_ij = bsq_calc(p[i][j][k],&geom) ;
		if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
	}
	fprintf(stderr,"initial bsq_max: %g\n",bsq_max) ;

	/* enforce boundary conditions */
	fixup(p) ;
	bound_prim(p) ;
    
    



#if( DO_FONT_FIX )
	set_Katm();
#endif 


}

void init_entwave()
{
  int i,j,k ;
  double x,y,z,sth,cth ;
  double ur,uh,up,u,rho ;
  double X[NDIM] ;
  struct of_geom geom ;
  double rhor;
  
  double myrho, myu, mycs, myv;
  double delta_rho;
  double cosa, sina;
  double delta_ampl = 1.e-1; //amplitude of the wave
  double k_vec_x = 2 * M_PI;  //wavevector
  double k_vec_y = 2 * M_PI;
  double k_vec_len = sqrt( k_vec_x * k_vec_x + k_vec_y * k_vec_y );
  double tfac = 1.e3; //factor by which to reduce velocity
  

  /* some physics parameters */
  gam = 5./3. ;
  
  /* some numerical parameters */
  lim = VANL ;
  failed = 0 ;	/* start slow */
  cour = 0.9 ;
  dt = 1.e-5 ;
  
  t = 0. ;
  hslope = 1. ;
  
  if(N2!=1) {
    //2D problem, use full pi-wedge in theta
    fractheta = 1.;
  }
  else{
    //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
    fractheta = 1.e-2;
  }
  
  fracphi = 1.;

  set_arrays() ;
  set_grid() ;
  
  
  myrho = 1.;
  myu = 4.*myrho / (gam * (gam-1));  //so that mycs is unity
  
  mycs = sqrt(gam * (gam-1) * myu / myrho);  //background sound speed
  myv = 1.;  //velocity with a magnitude of 1
  
  /* output choices */

    tf = tfac;///mycs;
  
  DTd = tf/10. ;	/* dumping frequency, in units of M */
  DTl = tf/10. ;	/* logfile frequency, in units of M */
  DTi = tf/10. ; 	/* image file frequ., in units of M */
  DTr = tf/10. ; /* restart file frequ., in units of M */
  DTr01 = 1000 ; 	/* restart file frequ., in timesteps */
  
  /* start diagnostic counters */
  dump_cnt = 0 ;
  image_cnt = 0 ;
  rdump_cnt = 0 ;
  rdump01_cnt = 0 ;
  defcon = 1. ;
  

  ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
    coord(i,j,k,CENT,X) ;
    bl_coord(X,&x,&y,&z) ;
    
    //applying the perturbations
    delta_rho = delta_ampl * cos( k_vec_x * x + k_vec_y * y );
    
  //  p[i][j][k][RHO] = myrho + delta_rho;
      if(x<.2){
          p[i][j][k][RHO] = 1.;
          p[i][j][k][U2] = 0.;

      }
      else if (x>.8){
          p[i][j][k][RHO] = 1.;
          p[i][j][k][U2] = 0.;
      }
      else{
          p[i][j][k][RHO] = 1e4;
          p[i][j][k][U2] = 0.;

      }
    p[i][j][k][UU] = myu/(tfac*tfac);
    p[i][j][k][U1] = myv/tfac;
    //p[i][j][k][U2] = myv/tfac;
    //p[i][j][k][U2] = myv/tfac;//cos(k_vec_x*x);
    p[i][j][k][U3] = 0 ;
    p[i][j][k][B1] = 0. ;
    p[i][j][k][B2] = 0. ;
    p[i][j][k][B3] = 0. ;
  }
  
  /* enforce boundary conditions */
  
  
  fixup(p) ;
  bound_prim(p) ;
  
  
  
  
  
#if( DO_FONT_FIX )
  set_Katm();
#endif
  
  
}

void init_sndwave()
{
  int i,j,k ;
  double x,y,z,sth,cth ;
  double ur,uh,up,u,rho ;
  double X[NDIM] ;
  struct of_geom geom ;
  
  double myrho, myu, mycs, myv;
  double delta_rho;
  double cosa, sina;
  double delta_ampl = 1e-5; //amplitude of the wave
  double k_vec_x = 2 * M_PI;  //wavevector
  double k_vec_y = 0;
  double k_vec_len = sqrt( k_vec_x * k_vec_x + k_vec_y * k_vec_y );
  double tfac = 1e4; //factor by which to reduce velocity
  
  
  /* some physics parameters */
  gam = 5./3. ;
  
  /* some numerical parameters */
  lim = VANL ;
  failed = 0 ;	/* start slow */
  cour = 0.9 ;
  dt = 1.e-5 ;
  
  t = 0. ;
  hslope = 1. ;
  
  if(N2!=1) {
    //2D problem, use full pi-wedge in theta
    fractheta = 1.;
  }
  else{
    //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
    fractheta = 1.e-2;
  }
  
  fracphi = 1.;
  
  set_arrays() ;
  set_grid() ;
  
  
  myrho = 1.;
  myu = myrho / (gam * (gam-1));  //so that mycs is unity
  
  mycs = sqrt(gam * (gam-1) * myu / myrho);  //background sound speed
  
  /* output choices */
  
  tf = tfac/mycs;
  
  DTd = tf/10. ;	/* dumping frequency, in units of M */
  DTl = tf/10. ;	/* logfile frequency, in units of M */
  DTi = tf/10. ; 	/* image file frequ., in units of M */
  DTr = tf/10. ; /* restart file frequ., in units of M */
  DTr01 = 1000 ; /* restart file frequ., in timesteps */
  
  /* start diagnostic counters */
  dump_cnt = 0 ;
  image_cnt = 0 ;
  rdump_cnt = 0 ;
  rdump01_cnt = 0 ;
  defcon = 1. ;
  
  ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
    coord(i,j,k,CENT,X) ;
    bl_coord(X,&x,&y,&z) ;
    
    //applying the perturbations
    delta_rho = delta_ampl * cos( k_vec_x * x + k_vec_y * y );
    
    p[i][j][k][RHO] = myrho + delta_rho;
    p[i][j][k][UU] = (myu + gam * myu * delta_rho / myrho)/(tfac*tfac);
    p[i][j][k][U1] = (delta_rho/myrho * mycs * k_vec_x / k_vec_len)/tfac;
    p[i][j][k][U2] = (delta_rho/myrho * mycs * k_vec_y / k_vec_len)/tfac;
    p[i][j][k][U3] = 0 ;
    p[i][j][k][B1] = 0. ;
    p[i][j][k][B2] = 0. ;
    p[i][j][k][B3] = 0. ;
  }
  
  /* enforce boundary conditions */
  
  
  fixup(p) ;
  bound_prim(p) ;
  
  
  
  
  
#if( DO_FONT_FIX )
  set_Katm();
#endif
  
  
}


/* this version starts w/ BL 4-velocity and
 * converts to relative 4-velocities in modified
 * Kerr-Schild coordinates */

void coord_transform(double *pr,int ii, int jj, int kk)
{
  double X[NDIM],r,th,phi,ucon[NDIM],uconp[NDIM],trans[NDIM][NDIM],tmp[NDIM] ;
  double AA,BB,CC,discr ;
  double utconp[NDIM], dxdxp[NDIM][NDIM], dxpdx[NDIM][NDIM] ;
  struct of_geom geom ;
  struct of_state q ;
  int i,j,k,m ;

  coord(ii,jj,kk,CENT,X) ;
  bl_coord(X,&r,&th,&phi) ;
  blgset(ii,jj,kk,&geom) ;

  ucon[1] = pr[U1] ;
  ucon[2] = pr[U2] ;
  ucon[3] = pr[U3] ;

  AA =     geom.gcov[TT][TT] ;
  BB = 2.*(geom.gcov[TT][1]*ucon[1] +
           geom.gcov[TT][2]*ucon[2] +
           geom.gcov[TT][3]*ucon[3]) ;
  CC = 1. +
          geom.gcov[1][1]*ucon[1]*ucon[1] +
          geom.gcov[2][2]*ucon[2]*ucon[2] +
          geom.gcov[3][3]*ucon[3]*ucon[3] +
      2.*(geom.gcov[1][2]*ucon[1]*ucon[2] +
          geom.gcov[1][3]*ucon[1]*ucon[3] +
          geom.gcov[2][3]*ucon[2]*ucon[3]) ;

  discr = BB*BB - 4.*AA*CC ;
  ucon[TT] = (-BB - sqrt(discr))/(2.*AA) ;
  /* now we've got ucon in BL coords */

  /* transform to Kerr-Schild */
  /* make transform matrix */
  DLOOP trans[j][k] = 0. ;
  DLOOPA trans[j][j] = 1. ;
  trans[0][1] = 2.*r/(r*r - 2.*r + a*a) ;
  trans[3][1] = a/(r*r - 2.*r + a*a) ;

  /* transform ucon */
  DLOOPA tmp[j] = 0. ;
  DLOOP tmp[j] += trans[j][k]*ucon[k] ;
  DLOOPA ucon[j] = tmp[j] ;
  /* now we've got ucon in KS coords */

  /* transform to KS' coords */
  /* dr^\mu/dx^\nu jacobian, where x^\nu are internal coords */
  dxdxp_func(X, dxdxp);
  /* dx^\mu/dr^\nu jacobian */
  invert_matrix(dxdxp, dxpdx);
  
  for(i=0;i<NDIM;i++) {
    uconp[i] = 0;
    for(j=0;j<NDIM;j++){
      uconp[i] += dxpdx[i][j]*ucon[j];
    }
  }
  //old way of doing things for Gammie coords
  //ucon[1] *= (1./(r - R0)) ;
  //ucon[2] *= (1./(M_PI + (1. - hslope)*M_PI*cos(2.*M_PI*X[2]))) ;
  //ucon[3] *= 1.; //!!!ATCH: no need to transform since will use phi = X[3]

  get_geometry(ii, jj, kk, CENT, &geom);
  
  /* now solve for relative 4-velocity that is used internally in the code:
   * we can use the same u^t because it didn't change under KS -> KS' */
  ucon_to_utcon(uconp,&geom,utconp);
  
  pr[U1] = utconp[1] ;
  pr[U2] = utconp[2] ;
  pr[U3] = utconp[3] ;

  /* done! */
}


double lfish_calc(double r)
{
	return(
   ((pow(a,2) - 2.*a*sqrt(r) + pow(r,2))*
      ((-2.*a*r*(pow(a,2) - 2.*a*sqrt(r) + pow(r,2)))/
         sqrt(2.*a*sqrt(r) + (-3. + r)*r) +
        ((a + (-2. + r)*sqrt(r))*(pow(r,3) + pow(a,2)*(2. + r)))/
         sqrt(1 + (2.*a)/pow(r,1.5) - 3./r)))/
    (pow(r,3)*sqrt(2.*a*sqrt(r) + (-3. + r)*r)*(pow(a,2) + (-2. + r)*r))
	) ;
}


