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

/**
 *
 * this contains the generic piece of code for advancing
 * the primitive variables
 *
 **/


#include "decs.h"

/** algorithmic choices **/

/* use local lax-friedrichs or HLL flux:  these are relative weights on each numerical flux */
#define HLLF  (0.0)
#define LAXF  (1.0)

/** end algorithmic choices **/


double advance( double pi[][N2M][N3M][NPR], double pb[][N2M][N3M][NPR],
               double Dt, double pf[][N2M][N3M][NPR], double *ndt1, double *ndt2, double *ndt3) ;

double fluxcalc( double pr[][N2M][N3M][NPR], double F[][N2M][N3M][NPR],
                int dir ) ;
void   flux_cd(double F1[][N2M][N3M][NPR], double F2[][N2M][N3M][NPR]) ;



/***********************************************************************************************/
/***********************************************************************************************
 step_ch():
 ---------
 -- handles the sequence of making the time step, the fixup of unphysical values,
 and the setting of boundary conditions;
 
 -- also sets the dynamically changing time step size;
 
 ***********************************************************************************************/
void step_ch(double *ndt1, double *ndt2, double *ndt3)
{
  double ndt;
    int i,j,k,m;
#ifdef MPI
  double mpi_buf[3];
#endif
  
    
  if(MASTER==mpi_rank) fprintf(stderr,"h") ;
  advance(p, p, 0.5*dt, ph, ndt1, ndt2, ndt3) ;   /* time step primitive variables to the half step */
  
    
	fixup(ph) ;         /* Set updated densities to floor, set limit for gamma */
	bound_prim(ph) ;    /* Set boundary conditions for primitive variables, flag bad ghost zones */
	fixup_utoprim(ph);  /* Fix the failure points using interpolation and updated ghost zone values */
    
        //ZLOOP eVol(p,p,ph,0.5*dt,i,j,k,0,1); /* Account for floor heat */
	bound_prim(ph) ;    /* Reset boundary conditions with fixed up points */
    
    
    
  /* Repeat and rinse for the full time (aka corrector) step:  */
  
  if(MASTER==mpi_rank) fprintf(stderr,"f") ;
#pragma omp parallel for schedule(static,N1*N2*N3/nthreads) collapse(3) default(none) shared(psave,p,nthreads) private(i,j,k,m)
  ZLOOP PLOOP psave[i][j][k][m] = p[i][j][k][m] ;
  advance(p, ph, dt,    p, ndt1, ndt2, ndt3) ;

  
  fixup(p) ;
  bound_prim(p) ;
  fixup_utoprim(p);

  bound_prim(p) ;
  
  
  
  /* Determine next time increment based on current characteristic speeds: */
  if(dt < 1.e-9 && MASTER==mpi_rank) {
    fprintf(stderr,"timestep too small: dt = %g\n", dt) ;
    exit(11) ;
  }
  
  
  /* increment time */
  
  t += dt ;
  
  /* set next timestep */
#ifdef MPI
  //find the minimum ndt across all MPI processes
  mpi_buf[0] = *ndt1;
  mpi_buf[1] = *ndt2;
  mpi_buf[2] = *ndt3;
  MPI_Allreduce(MPI_IN_PLACE,mpi_buf,3,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  *ndt1 = mpi_buf[0];
  *ndt2 = mpi_buf[1];
  *ndt3 = mpi_buf[2];
#endif
  ndt = defcon * 1./(1./(*ndt1) + 1./(*ndt2) + 1./(*ndt3)) ;
  if(ndt > SAFE*dt) ndt = SAFE*dt ;
  dt = ndt ;
  if(t + dt > tf) dt = tf - t ;  /* but don't step beyond end of run */
  
#ifdef MPI
  mpi_buf[0] = t;
  mpi_buf[1] = dt;
  MPI_Bcast(mpi_buf,2,MPI_DOUBLE,MASTER,MPI_COMM_WORLD);
  t = mpi_buf[0];
  dt = mpi_buf[1];
#endif
  
  /* done! */
}



/***********************************************************************************************/
/***********************************************************************************************
 advance():
 ---------
 -- responsible for what happens during a time step update, including the flux calculation,
 the constrained transport calculation (aka flux_ct()), the finite difference
 form of the time integral, and the calculation of the primitive variables from the
 update conserved variables;
 -- also handles the "fix_flux()" call that sets the boundary condition on the fluxes;
 
 ***********************************************************************************************/
double advance(
               double pi[][N2M][N3M][NPR],
               double pb[][N2M][N3M][NPR],
               double Dt,
               double pf[][N2M][N3M][NPR],
               double *ndt1,
               double *ndt2,
               double *ndt3
               )
{
  int i,j,k,m, dummy;
  double ndt,U[NPR],dU[NPR], usave;
  struct of_geom geom ;
  struct of_state q ;
  double kappa;
  
  double Uo[NPR],po[NPR] ;
  
  int was_floor_activated = 0;
  
  
#pragma omp parallel for schedule(static,N1*N2*N3/nthreads) collapse(3) default(none) shared(pf,pi,nthreads) private(i,j,k,m)
  ZLOOP PLOOP pf[i][j][k][m] = pi[i][j][k][m] ;        /* needed for Utoprim */
  
  if(MASTER==mpi_rank) fprintf(stderr,"0") ;
  *ndt1 = fluxcalc(pb, F1, 1) ;
  *ndt2 = fluxcalc(pb, F2, 2) ;
  *ndt3 = fluxcalc(pb, F3, 3) ;
  
  fix_flux(F1,F2,F3) ;
  
  flux_ct(F1,F2,F3) ;
  
  /* evaluate diagnostics based on fluxes */
  diag_flux(F1,F2) ;
  
  if(MASTER==mpi_rank) fprintf(stderr,"1") ;
  /** now update pi to pf **/
  
#pragma omp parallel for schedule(static,N1*N2*N3/nthreads) collapse(3) default(none) \
    shared(pb,pi,pf,pflag,F1,F2,F3,dx,nthreads,failimage,Katm,Dt) \
    private(i,j,k,m,q,dU,geom,U,usave)
  ZLOOP {

    get_geometry(i,j,k,CENT,&geom) ;

    source(pb[i][j][k],&geom,i,j,k,dU,Dt) ;
    get_state(pb[i][j][k],&geom,&q) ;
    //misc_source(pb[i][j][k],pb[i+1][j][k],pb[i-1][j][k], pb[i][j+1][k],pb[i][j-1][k],i, j, k, &geom,&q,dU, Dt);

    get_state(pi[i][j][k],&geom,&q) ;
    primtoU(pi[i][j][k],&q,&geom,U) ;

      
    PLOOP {
      U[m] += Dt*(
#if( N1 != 1 )
		  - (F1[i+1][j][k][m] - F1[i][j][k][m])/dx[1]
#endif
#if( N2 != 1 )
		  - (F2[i][j+1][k][m] - F2[i][j][k][m])/dx[2]
#endif
#if( N3 != 1 )
		  - (F3[i][j][k+1][m] - F3[i][j][k][m])/dx[3]
#endif
		  + dU[m]

		  ) ;

    }

    
    
    pflag[i][j][k] = Utoprim_2d(U, geom.gcov, geom.gcon, geom.g, pf[i][j][k]);
    
    
    if( pflag[i][j][k] ) {
      failimage[i][j][k][0]++ ;
      was_floor_activated = 1;
    }
    else {
      was_floor_activated = 0;
    }

    
#if( DO_FONT_FIX || DOKTOT )
    if (DOKTOT) {
      kappa = pf[i][j][k][KTOT];
    }
    else {
      kappa = Katm[i];
    }
    if( pflag[i][j][k] ) {
      
      pflag[i][j][k] = Utoprim_1dvsq2fix1(U, geom.gcov, geom.gcon, geom.g, pf[i][j][k], kappa );
      if( pflag[i][j][k] ) {
        was_floor_activated = 2;
	failimage[i][j][k][1]++ ;
	pflag[i][j][k] = Utoprim_1dfix1(U, geom.gcov, geom.gcon, geom.g, pf[i][j][k], kappa);
        if( pflag[i][j][k] ) {
          was_floor_activated = 3;
          failimage[i][j][k][2]++ ;
        }
      }
    }
#endif
    
#if(DOKTOT)
    compute_ktot(pi,pb,pf,Dt,i,j,k,was_floor_activated,0);
#endif
    
    
  }
  
    
  
  
  if(MASTER==mpi_rank)  fprintf(stderr,"2") ;
  
  return(0) ;
}

#if(DOKTOT)
void compute_ktot(double pi[][N2M][N3M][NPR],double prh[][N2M][N3M][NPR], double pr[][N2M][N3M][NPR], int i, int j, int k, double Dt, int was_floor_activated, int is_after_fixup)
{
    
    double kappatot;

    kappatot =(gam-1.)*pr[i][j][k][UU]*pow(pr[i][j][k][RHO],-gam);
    
    //reset to actual internal energy/Entropy
    pr[i][j][k][KTOT] += (kappatot-pr[i][j][k][KTOT]);
    
    return;
}
/* this function sets the entropy
 * based on the initial conditions */
void init_entropy(){
    int i, j, k;
    ZSLOOP(-N1G,N1+N1G-1,-N2G,N2G+N2-1,-N3G,N3G+N3-1){
        p[i][j][k][KTOT] = (gam-1.)*p[i][j][k][UU]*pow(p[i][j][k][RHO],-gam);
    }
}
#endif

/***********************************************************************************************/
/***********************************************************************************************
 fluxcalc():
 ---------
 -- sets the numerical fluxes, avaluated at the cell boundaries using the slope limiter
 slope_lim();
 
 -- only has HLL and Lax-Friedrichs  approximate Riemann solvers implemented;
 
 ***********************************************************************************************/
double fluxcalc(
                double pr[][N2M][N3M][NPR],
                double F[][N2M][N3M][NPR],
                int dir
                )
{
  int i,j,k,m,idel,jdel,kdel,face ;
  double p_l[NPR],p_r[NPR],F_l[NPR],F_r[NPR],U_l[NPR],U_r[NPR] ;
  double cmax_l,cmax_r,cmin_l,cmin_r,cmax,cmin,ndt,dtij ;
  double ctop, s_l, s_r ;
  double cmax_lw,cmax_rw,cmin_lw,cmin_rw,cmaxw,cminw ;
  double ctopw, dummy;
  double  tmp;
  
  struct of_geom geom ;
  struct of_state state_l,state_r ;
  struct of_state state_lw,state_rw ;
  
  void rescale(double *pr, int which, int dir, int ii, int jj, int kk, int face, struct of_geom *geom) ;
  double bsq ;
  
  if( (N1==1 && dir == 1) || (N2==1 && dir == 2) || (N3==1 && dir == 3)) {
    //if current direction is trivial, reset fluxes in the current direction to zero so no evolution due to this direction
#pragma omp parallel for schedule(static,(N1+2*D1)*(N2+2*D2)*(N3+2*D3)/nthreads) collapse(3) default(none) shared(dq,F,nthreads) private(i,j,k,m)
    ZSLOOP(-D1,N1-1+D1,-D2,N2-1+D2,-D3,N3-1+D3) PLOOP {
      dq[i][j][k][m] = 0.;
      F[i][j][k][m] = 0.;
    }
    ndt = 1.e9;
    return(ndt) ;
  }

  if     (dir == 1) {idel = 1; jdel = 0; kdel = 0; face = FACE1;}
  else if(dir == 2) {idel = 0; jdel = 1; kdel = 0; face = FACE2;}
  else if(dir == 3) {idel = 0; jdel = 0; kdel = 1; face = FACE3;}
  else { exit(10); }
  
#if(RESCALE)
  /** evaluate slopes of primitive variables **/
  /* first rescale */
#pragma omp parallel for schedule(static,N1M*N2M*N3M/nthreads) default(none) collapse(3) shared(pr,dir,nthreads) private(i,j,k,m,geom)
  ZSLOOP(-N1G,N1+N1G-1,-N2G,N2+N2G-1,-N3G,N3+N3G-1) {
    get_geometry(i,j,k,CENT,&geom) ;
    rescale(pr[i][j][k],FORWARD, dir, i,j,k,CENT,&geom) ;
  }
#endif
  /* then evaluate slopes on active grid plus 1-cell-wide layer of ghost cells */
#pragma omp parallel for schedule(static,(N1+2*D1)*(N2+2*D2)*(N3+2*D3)/nthreads) collapse(3) default(none) shared(dq,pr,idel,jdel,kdel,nthreads) private(i,j,k,m)
 ZSLOOP(-D1,N1-1+D1,-D2,N2-1+D2,-D3,N3-1+D3){
    PLOOP {
   
      dq[i][j][k][m] = slope_lim(
                              pr[i-idel][j-jdel][k-kdel][m],
                              pr[i][j][k][m],
                              pr[i+idel][j+jdel][k+kdel][m]
                              ) ;
    }
  
  }
  
  
  ndt = 1.e9 ;
#pragma omp parallel for schedule(static,(D1*(N1+!idel)+1)*(D2*(N2+!jdel)+1)*(D3*(N3+!kdel)+1)/nthreads) \
    collapse(3) \
    shared(F,dq,idel,jdel,kdel,nthreads,dir,pr,face,dx,cour,mpi_startn,mpi_ntot) \
    private(i,j,k,m,geom,p_l,p_r,state_l,state_r,F_l,F_r,U_l,U_r,cmin_l,cmin_r,cmin,cmax,ctop,dtij,tmp,cmax_l,cmax_r) \
    default(none) \
    reduction(min:ndt)
  ZSLOOP(-D1*!idel,D1*N1,-D2*!jdel,D2*N2,-D3*!kdel,D3*N3){ //-jdel,N1,-idel,N2) {
    
    /* this avoids problems on the pole */
    
    
    if(dir == 2 && (j+mpi_startn[2] == 0 || j+mpi_startn[2] == mpi_ntot[2]) && BL) {
      //FIXFLUXMARK: fix_flux-like zero out of fluxes in x2-dir
      PLOOP F[i][j][k][m] = 0. ;
    }
    else {
      
      PLOOP {
        p_l[m] = pr[i-idel][j-jdel][k-kdel][m]
        + 0.5*dq[i-idel][j-jdel][k-kdel][m] ;
        p_r[m] = pr[i][j][k][m]
        - 0.5*dq[i][j][k][m] ;
      }
      
      get_geometry(i,j,k,face,&geom) ;
      
#if(RESCALE)
      rescale(p_l,REVERSE,dir,i,j,k,face,&geom) ;
      rescale(p_r,REVERSE,dir,i,j,k,face,&geom) ;
#endif
      get_state(p_l,&geom,&state_l) ;
      get_state(p_r,&geom,&state_r) ;
      
      primtoflux(p_l,&state_l,dir,&geom,F_l) ;
      primtoflux(p_r,&state_r,dir,&geom,F_r) ;
      
      primtoflux(p_l,&state_l,TT, &geom,U_l) ;
      primtoflux(p_r,&state_r,TT, &geom,U_r) ;
      
      vchar(p_l,&state_l,&geom,dir,&cmax_l,&cmin_l) ;
      vchar(p_r,&state_r,&geom,dir,&cmax_r,&cmin_r) ;
      
      cmax = fabs(MY_MAX(MY_MAX(0., cmax_l),  cmax_r)) ;
      cmin = fabs(MY_MAX(MY_MAX(0.,-cmin_l), -cmin_r)) ;
      ctop = MY_MAX(cmax,cmin) ;
      
      PLOOP {
          F[i][j][k][m] =
          HLLF*(
                (cmax*F_l[m] + cmin*F_r[m]
                 - cmax*cmin*(U_r[m] - U_l[m]))/
                (cmax + cmin + SMALL)
                ) +
          LAXF*(
                0.5*(F_l[m] + F_r[m]
                     - ctop*(U_r[m] - U_l[m]))
                
                ) ;
        }
        
//                        }
//                    }
      
      /* evaluate restriction on timestep */
      cmax = MY_MAX(cmax,cmin) ;
      dtij = cour*dx[dir]/cmax ;
      if(dtij < ndt) ndt = dtij ;
        
    }
  }
  
#if(RESCALE)
#pragma omp parallel for schedule(static,N1M*N2M*N3M/nthreads) default(none) collapse(3) shared(pr,dir,nthreads) private(i,j,k,m,geom)
  ZSLOOP(-N1G,N1+N1G-1,-N2G,N2+N2G-1,-N3G,N3+N3G) {
    get_geometry(i,j,k,CENT,&geom) ;
    rescale(pr[i][j][k],REVERSE,dir,i,j,k,CENT,&geom) ;
  }
#endif
  
  
  return(ndt) ;
  
}

/***********************************************************************************************/
/***********************************************************************************************
 flux_ct():
 ---------
 -- performs the flux-averaging used to preserve the del.B = 0 constraint (see Toth 2000);
 
 ***********************************************************************************************/
void flux_ct(double F1[][N2M][N3M][NPR], double F2[][N2M][N3M][NPR], double F3[][N2M][N3M][NPR])
{
  int i,j,k ;

  /* calculate EMFs */
#define DOE1 (N2>1 && N3>1)
#define DOE2 (N1>1 && N3>1)
#define DOE3 (N1>1 && N2>1)
  
  /* Toth approach: just average */
#pragma omp parallel for schedule(static,(N1+D1)*(N2+D2)*(N3+D3)/nthreads) collapse(3) default(none) shared(emf,F1,F2,F3,nthreads) private(i,j,k)
  ZSLOOP(0,N1-1+D1,0,N2-1+D2,0,N3-1+D3) {
//E_i = e_{ijk} F_j B_k
#if(DOE1)
    //E1
    emf[1][i][j][k] = 0.25*((F2[i][j][k][B3] + F2[i][j][k-1][B3])
                          - (F3[i][j][k][B2] + F3[i][j-1][k][B2])
                            ) ;
#endif

#if(DOE2)
    //E2
    emf[2][i][j][k] = 0.25*((F3[i][j][k][B1] + F3[i-1][j][k][B1])
                          - (F1[i][j][k][B3] + F1[i][j][k-1][B3])
                            ) ;
#endif

#if(DOE3)
    //E3
    emf[3][i][j][k] = 0.25*((F1[i][j][k][B2] + F1[i][j-1][k][B2])
                          - (F2[i][j][k][B1] + F2[i-1][j][k][B1])
                           ) ;
#endif
  }
  /* rewrite EMFs as fluxes, after Toth */
#pragma omp parallel for schedule(static,(N1+D1)*N2*N3/nthreads) collapse(3) default(none) shared(emf,F1,F2,F3,nthreads) private(i,j,k)
  ZSLOOP(0,N1-1+D1,0,N2-1,0,N3-1) {
#if(N1>1)
    F1[i][j][k][B1] = 0. ;
#endif
#if(DOE3)
    F1[i][j][k][B2] =  0.5*(emf[3][i][j][k] + emf[3][i][j+1][k]) ;
#endif
#if(DOE2)
    F1[i][j][k][B3] = -0.5*(emf[2][i][j][k] + emf[2][i][j][k+1]) ;
#endif
  }

#pragma omp parallel for schedule(static,N1*(N2+D2)*N3/nthreads) collapse(3) default(none) shared(emf,F1,F2,F3,nthreads) private(i,j,k)
  ZSLOOP(0,N1-1,0,N2-1+D2,0,N3-1) {
#if(DOE3)
    F2[i][j][k][B1] = -0.5*(emf[3][i][j][k] + emf[3][i+1][j][k]) ;
#endif
#if(DOE1)
    F2[i][j][k][B3] =  0.5*(emf[1][i][j][k] + emf[1][i][j][k+1]) ;
#endif
#if(N2>1)
    F2[i][j][k][B2] = 0. ;
#endif
  }

#pragma omp parallel for schedule(static,N1*N2*(N3+D3)/nthreads) collapse(3) default(none) shared(emf,F1,F2,F3,nthreads) private(i,j,k)
  ZSLOOP(0,N1-1,0,N2-1,0,N3-1+D3) {
#if(DOE2)
    F3[i][j][k][B1] =  0.5*(emf[2][i][j][k] + emf[2][i+1][j][k]) ;
#endif
#if(DOE1)
    F3[i][j][k][B2] = -0.5*(emf[1][i][j][k] + emf[1][i][j+1][k]) ;
#endif
#if(N3>1)
    F3[i][j][k][B3] = 0. ;
#endif
  }
}
