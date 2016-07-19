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

#define FLOOP for(m=0;m<B1;m++)

/* apply floors to density, internal energy */

void fixup(double (* pv)[N2M][N3M][NPR]) 
{
  int i,j,k ;

#pragma omp parallel for schedule(static,N1*N2*N3/nthreads) collapse(3) default(none) shared(pv,nthreads) private(i,j,k)
  ZLOOP {
    fixup1zone( i, j, k, pv[i][j][k] );
  }

}

//SASMARK: time-dependent floor specific to Rodrigo's torus
//returns fixed-density floors (independent of magnetic field)
void get_rho_u_floor( double r, double th, double phi, double *rho_floor, double *u_floor )
{
  double rhoflr, uuflr;
  double uuscal,rhoscal;
  
  if (BL == 0){
    rhoflr =1e-6;
    uuflr  = pow(rhoflr,gam) ;
  }
  else{
    rhoscal = pow(r,-POWRHO) ;
    //rhoscal = pow(r,-2) ;
    uuscal = pow(rhoscal,gam); //rhoscal/r ;
    rhoflr = RHOMIN*rhoscal; //this is Rodrigo's rhot
    uuflr  = UUMIN*uuscal;
    
    if( rhoflr < RHOMINLIMIT ) rhoflr = RHOMINLIMIT;
    if( uuflr  < UUMINLIMIT  ) uuflr  = UUMINLIMIT;
  }

  *rho_floor = rhoflr;
  *u_floor = uuflr;

}

void fixup1zone( int i, int j, int k, double pv[NPR] )
{
  void get_rho_u_floor( double r, double th, double phi, double *rho_floor, double *u_floor );
  double r,th,phi,X[NDIM], rhoflr,uuflr ;
  double bsq, ughat, agame, thex,gfac;
  double f, gamma, alpha, aold, thold;
  struct of_geom geom ;
  double pv_prefloor[NPR], dpv[NPR];
  double U_prefloor[NPR], dU[NPR], U[NPR];
  struct of_state q ;
  int flag;
  int m;
  double uel4;
  int dofloor;
  double kappa;
#if(DRIFT_FLOOR)
  double betapar, betasq, betasqmax, ucondr[NDIM], ucovdr[NDIM];
  double udrsq, udotB, Bcov[NDIM], Bcon[NDIM], Bsq, B, wold;
  double QdotB, wnew, x, vpar, one_over_ucondr_t, vcon[NDIM], ut;
  double udotBnew, QdotBnew, reldiff, ucon[NDIM], utcon[NDIM];
  double vparold, ut2, trans;
#endif

  
  //make a backup
  PLOOP pv_prefloor[m] = pv[m];

  coord(i,j,k,CENT,X) ;
  get_phys_coord(i,j,k,&r,&th,&phi) ;

  get_rho_u_floor(r,th,phi,&rhoflr,&uuflr);
  
  //compute the square of fluid frame magnetic field (twice magnetic pressure)
  get_geometry(i,j,k,CENT,&geom) ;
  get_state(pv,&geom,&q) ;
  bsq = dot(q.bcon,q.bcov) ;
  
  //tie floors to the local values of magnetic field and internal energy density
#if(1)
  if( rhoflr < bsq / BSQORHOMAX ) rhoflr = bsq / BSQORHOMAX;
  if( uuflr < bsq / BSQOUMAX ) uuflr = bsq / BSQOUMAX;
  if( rhoflr < pv[UU] / UORHOMAX ) rhoflr = pv[UU] / UORHOMAX;
#endif


  dofloor = 0;
  /* floor on density and internal energy density (momentum *not* conserved) */
    if(pv[RHO] < rhoflr ){
        pv[RHO] = rhoflr;
        dofloor = 1;
    }
  
    if(pv[UU]  < uuflr  )   {
      
        pv[UU]  = uuflr;
        dofloor = 1;
    }

#if( DRIFT_FLOOR )
  if(dofloor && (trans=10.*bsq/MY_MIN(pv[RHO],pv[UU])-1.) > 0. ) {
    //Apply floors in the drift frame. Since the drift frame
    //is only well-defined in the regions where magnetic field is non-negligible
    //only use it when trans > 0 <=> B^2/4pi > 0.1 min(rho c^2, ug)
    //transitional variable
    if (trans>1.) {
      trans = 1.;
    }
    //set velocity to drift velocity
    betapar = -q.bcon[0]/((bsq+SMALL)*q.ucon[0]);
    betasq = betapar*betapar*bsq;
    betasqmax = 1.-1./(GAMMAMAX*GAMMAMAX);
    if(betasq>betasqmax) {
      betasq = betasqmax;
    }
    gamma = 1./sqrt(1-betasq);
    for(m = 0; m < NDIM; m++) {
      ucondr[m] = gamma*(q.ucon[m]+betapar*q.bcon[m]);
    }
    //the below prints out -1 and thereby verifies velocity normalization
    //lower(ucondr, &geom, ucovdr);
    //udrsq = dot(ucondr,ucovdr);
    //printf("(%d,%d,%d): udrsq = %g\n",i,j,k,udrsq);

    Bcon[0] = 0.;
    for(m = 1; m < NDIM; m++) {
      Bcon[m] = pv[B1-1+m];
    }
    
    lower(Bcon, &geom, Bcov);
    udotB = dot(q.ucon,Bcov);
    Bsq = dot(Bcon,Bcov);
    B = sqrt(Bsq);

    //check if u^t is correctly computed from udrift^t and vpar
    //vparold = udotB/(B*q.ucon[0]);
    //ut2 = sqrt(ucondr[0]*ucondr[0]/(1.-ucondr[0]*ucondr[0]*vparold*vparold));
    //reldiff = 2.*(q.ucon[0]-ut2)/(fabs(q.ucon[0])+fabs(ut2));
    //if( reldiff > 1e-15 ) {
    //  printf("(%d,%d,%d): utcode = %g, utmy = %g, reldiff = %g\n",
    //         i+mpi_startn[1],j+mpi_startn[2],k+mpi_startn[3],
    //         q.ucon[0], ut2, reldiff);
    //}
    
    
    //enthalpy before the floors
    wold = pv_prefloor[RHO]+pv_prefloor[UU]*gam;

    //B^\mu Q_\mu = (B^\mu u_\mu) (\rho+u+p) u^t (eq. (26) divided by alpha; Noble et al. 2006)
    QdotB = udotB*wold*q.ucon[0];
    
    //now, we apply the floors to w and recompute vpar, parallel velocity along
    //the magnetic field, which is reduced due to the floor addition in the drift frame

    //enthalpy after the floors
    wnew = pv[RHO]+pv[UU]*gam;
    //wnew = wold;
    
    x = 2.*QdotB/(B*wnew*ucondr[0]+SMALL);
    
    //new parallel velocity
    vpar = x/( ucondr[0]*(1.+sqrt(1.+x*x)) );

//    reldiff = 2.*(vpar-vparold)/(fabs(vpar)+fabs(vparold));
//    if( reldiff > 1e-13 ) {
//      printf("(%d,%d,%d): vpar = %g, vparold = %g, reldiff = %g\n",
//             i+mpi_startn[1],j+mpi_startn[2],k+mpi_startn[3],
//             vpar, vparold, reldiff);
//    }

    
    one_over_ucondr_t = 1./ucondr[0];
    
    //new contravariant 3-velocity, v^i
    vcon[0] = 1.;
    for (m=1; m<NDIM; m++) {
      //parallel (to B) plus perpendicular (to B) velocities
      vcon[m] = vpar*Bcon[m]/(B+SMALL) + ucondr[m]*one_over_ucondr_t;
    }
    
    //compute u^t corresponding to the new v^i
    ut_calc_3vel(vcon, &geom, &ut);
    
    for(m=0;m<NDIM;m++) {
      ucon[m] = ut*vcon[m];
    }
    ucon_to_utcon(ucon, &geom, utcon);
                  
    //now convert 3-vel to relative 4-velocity and put it into pv[U1..U3]
    //\tilde u^i = u^t(v^i-g^{ti}/g^{tt})
    for(m=1;m<NDIM;m++) {
      pv[m+UU] = utcon[m]*trans + pv_prefloor[m+UU]*(1.-trans);
    }
    
    //
    // CHECKS
    //
    
//    reldiff = 2.*(q.ucon[0]-ut)/(fabs(q.ucon[0])+fabs(ut));
//    if( reldiff > 1e-15 ) {
//      printf("(%d,%d,%d): utcode = %g, utmy = %g, reldiff = %g\n",
//             i+mpi_startn[1],j+mpi_startn[2],k+mpi_startn[3],
//             q.ucon[0], ut, reldiff);
//    }
    
    
//    for(m=U1;m<=U3;m++) {
//      reldiff = 2.*(pv[m]-pv_prefloor[m])/(fabs(pv[m])+fabs(pv_prefloor[m]));
//      if( reldiff > 1e-14 ) {
//        printf("(%d,%d,%d): pv[%d] = %g, pv[%d] = %g, diff = %g\n",
//               i+mpi_startn[1],j+mpi_startn[2],k+mpi_startn[3],
//               m, pv[m], m, pv_prefloor[m], reldiff);
//      }
//    }

#if(0) //check whether QdotB is unchanged by the floor
    //check if QdotB stayed the same
    get_state(pv,&geom,&q) ;
    udotBnew = dot(q.ucon,Bcov);
    //enthalpy before the floors
    wnew = pv[RHO]+pv[UU]*gam;
    
    //B^\mu Q_\mu = (B^\mu u_\mu) (\rho+u+p) u^t (eq. (26) divided by alpha; Noble et al. 2006)
    QdotBnew = udotBnew*wnew*q.ucon[0];
    
    reldiff = 2.*fabs(QdotBnew - QdotB)/(fabs(QdotBnew)+fabs(QdotB)+SMALL);
    if( reldiff > 1e-10 ) {
      printf("(%d,%d,%d): QdotBold = %g, QdotBnew = %g, diff = %g\n",
             i+mpi_startn[1],j+mpi_startn[2],k+mpi_startn[3],
             QdotB, QdotBnew, reldiff);
    }
#endif //end check on QdotB

#endif
}

#if(DOKTOT)
    if(pv_prefloor[RHO] >0){
        pv[KTOT] = pv[KTOT]* pow(pv[RHO]/pv_prefloor[RHO],-gam);
    }

#endif
    
     /* limit gamma wrt normal observer */

  if( gamma_calc(pv,&geom,&gamma) ) {
    /* Treat gamma failure here as "fixable" for fixup_utoprim() */
    pflag[i][j][k] = -333;
    failimage[i][j][k][3]++ ;
#if(DOFLR)
      //treat no gamma solution as the floor activation
      pv[FLR] = 1.;
#endif
  }
  else { 
    if(gamma > GAMMAMAX) {
      f = sqrt(
	       (GAMMAMAX*GAMMAMAX - 1.)/
	       (gamma*gamma - 1.)
	       ) ;
      pv[U1] *= f ;	
      pv[U2] *= f ;	
      pv[U3] *= f ;	
#if(DOFLR)
      //treat limiting gamma as the floor activation
      pv[FLR] = 1.;
#endif
    }
  }



    


  return;
}


/**************************************************************************************
 INTERPOLATION STENCILS:  
 ------------------------
   -- let the stencils be characterized by the following numbering convention:

           1 2 3 
           8 x 4      where x is the point at which we are interpolating 
           7 6 5
*******************************************************************************************/

//!!!ATCH: leave the fixup stencils unchanged, i.e., never couple cells at different k for fixups (for now)

/* 12345678 */
#define AVG8(pr,i,j,k,m)  \
        (0.125*(pr[i-1][j+1][k][m]+pr[i][j+1][k][m]+pr[i+1][j+1][k][m]+pr[i+1][j][k][m]+pr[i+1][j-1][k][m]+pr[i][j-1][k][m]+pr[i-1][j-1][k][m]+pr[i-1][j][k][m]))

/* 2468  */
#define AVG4_1(pr,i,j,k,m) (0.25*(pr[i][j+1][k][m]+pr[i][j-1][k][m]+pr[i-1][j][k][m]+pr[i+1][j][k][m]))

/* 1357  */
#define AVG4_2(pr,i,j,k,m) (0.25*(pr[i+1][j+1][k][m]+pr[i+1][j-1][k][m]+pr[i-1][j+1][k][m]+pr[i-1][j-1][k][m]))

/* + shaped,  Linear interpolation in X1 or X2 directions using only neighbors in these direction */
/* 48  */
#define AVG2_X1(pr,i,j,k,m) (0.5*(pr[i-1][j  ][k][m]+pr[i+1][j  ][k][m]))
/* 26  */
#define AVG2_X2(pr,i,j,k,m) (0.5*(pr[i  ][j-1][k][m]+pr[i  ][j+1][k][m]))

/* x shaped,  Linear interpolation diagonally along both X1 and X2 directions "corner" neighbors */
/* 37  */
#define AVG2_1_X1X2(pr,i,j,k,m) (0.5*(pr[i-1][j-1][k][m]+pr[i+1][j+1][k][m]))
/* 15  */
#define AVG2_2_X1X2(pr,i,j,k,m) (0.5*(pr[i-1][j+1][k][m]+pr[i+1][j-1][k][m]))

/*******************************************************************************************
  fixup_utoprim(): 

    -- figures out (w/ pflag[]) which stencil to use to interpolate bad point from neighbors;

    -- here we use the following numbering scheme for the neighboring cells to i,j:  

                      1  2  3 
                      8  x  4        where "x" is the (i,j) cell or the cell to be interpolated
                      7  6  5

 *******************************************************************************************/


void fixup_utoprim( double (*pv)[N2M][N3M][NPR] )  
{
  int i, j, k, m;
  static int pf[9];
    

  /* Flip the logic of the pflag[] so that it now indicates which cells are good  */
#pragma omp parallel for schedule(static,N1M*N2M*N3M/nthreads) collapse(3) \
    default(none) \
    shared(pflag,pv,nthreads) \
    private(i,j,k,m,pf)
  ZSLOOP(-N1G,(N1+N1G-1),-N2G,(N2+N2G-1),-N3G,(N3+N3G-1)) { pflag[i][j][k] = !pflag[i][j][k] ; }

  /* Fix the interior points first */
#pragma omp parallel for schedule(static,N1*N2*N3/nthreads) collapse(3) \
default(none) \
shared(pflag,pv,nthreads,failimage,game,game4,game5,gam) \
private(i,j,k,m,pf,uel4,uel5,ugold,ueldis,uels,ughat,ughats)
  ZSLOOP(0,(N1-1),0,(N2-1),0,(N3-1)) {
    if( pflag[i][j][k] == 0 ) {
      pf[1] = pflag[i-1][j+1][k];   pf[2] = pflag[i][j+1][k];  pf[3] = pflag[i+1][j+1][k];
      pf[8] = pflag[i-1][j  ][k];                              pf[4] = pflag[i+1][j  ][k];
      pf[7] = pflag[i-1][j-1][k];   pf[6] = pflag[i][j-1][k];  pf[5] = pflag[i+1][j-1][k];

      /* Now the pf's  are true if they represent good points */

//      if(      pf[1]&&pf[2]&&pf[3]&&pf[4]&&pf[5]&&pf[6]&&pf[7]&&pf[8] ){ FLOOP pv[i][j][k] = AVG8(            pv,i,j,k)                   ; }
//      else if(        pf[2]&&       pf[4]&&       pf[6]&&       pf[8] ){ FLOOP pv[i][j][k] = AVG4_1(          pv,i,j,k)                   ; }
//      else if( pf[1]&&       pf[3]&&       pf[5]&&       pf[7]        ){ FLOOP pv[i][j][k] = AVG4_2(          pv,i,j,k)                   ; }
//      else if(               pf[3]&&pf[4]&&              pf[7]&&pf[8] ){ FLOOP pv[i][j][k] = 0.5*(AVG2_1_X1X2(pv,i,j,k)+AVG2_X1(pv,i,j,k)); }
//      else if(        pf[2]&&pf[3]&&              pf[6]&&pf[7]        ){ FLOOP pv[i][j][k] = 0.5*(AVG2_1_X1X2(pv,i,j,k)+AVG2_X2(pv,i,j,k)); }
//      else if( pf[1]&&              pf[4]&&pf[5]&&              pf[8] ){ FLOOP pv[i][j][k] = 0.5*(AVG2_2_X1X2(pv,i,j,k)+AVG2_X1(pv,i,j,k)); }
//      else if( pf[1]&&pf[2]&&              pf[5]&&pf[6]               ){ FLOOP pv[i][j][k] = 0.5*(AVG2_2_X1X2(pv,i,j,k)+AVG2_X2(pv,i,j,k)); }
//      else if(               pf[3]&&                     pf[7]        ){ FLOOP pv[i][j][k] = AVG2_1_X1X2(     pv,i,j,k)                   ; }
//      else if( pf[1]&&                     pf[5]                      ){ FLOOP pv[i][j][k] = AVG2_2_X1X2(     pv,i,j,k)                   ; }
//      else if(        pf[2]&&                     pf[6]               ){ FLOOP pv[i][j][k] = AVG2_X2(         pv,i,j,k)                   ; }
//      else if(                      pf[4]&&                     pf[8] ){ FLOOP pv[i][j][k] = AVG2_X1(         pv,i,j,k)                   ; }

// Old way:
      if(             pf[2]&&       pf[4]&&       pf[6]&&       pf[8] ){
          FLOOP pv[i][j][k][m] = AVG4_1(          pv,i,j,k,m)                   ;
      }
      else if( pf[1]&&       pf[3]&&       pf[5]&&       pf[7]        ){
          FLOOP pv[i][j][k][m] = AVG4_2(          pv,i,j,k,m)                   ;
      }
      else{ 
//	fflush(stderr);
//	fprintf(stderr,"fixup_utoprim()1: no good stencils1: i j pflag = %d %d : \n %4d %4d %4d \n %4d %4d %4d \n %4d %4d %4d \n\n", 
//		i,j,pflag[i-1][j+1],pflag[i][j+1],pflag[i+1][j+1],pflag[i-1][j],pflag[i][j],pflag[i+1][j],
//		pflag[i-1][j-1],pflag[i][j-1],pflag[i+1][j-1]);
//	fflush(stderr);
	failimage[i][j][k][4]++ ;


        for( m = RHO; m <= UU; m++ ) {
            
            pv[i][j][k][m] = 0.5*( AVG4_1(pv,i,j,k,m) + AVG4_2(pv,i,j,k,m) );
        
        }

	pv[i][j][k][U1] = pv[i][j][k][U2] = pv[i][j][k][U3] = 0.;
      }
        

      pflag[i][j][k] = 0;                  /* The cell has been fixed so we can use it for interpolation elsewhere */
      fixup1zone( i, j, k, pv[i][j][k] ) ;  /* Floor and limit gamma the interpolated value */
    }
      
      
      // if (pv[i][j][k][KEL5] != pv[i][j][k][KTOT]){
      //     fprintf(stderr, "In utoprim fixup! %d %d %d \n current k5, ktot %g %g \n prefloor k5, ktot %g %g \n", i,j,k, pv[i][j][k][KEL5], pv[i][j][k][KTOT],k5init, ktotinit);
      //     exit(1);
      // }
  }

  return;
}

#if( DO_FONT_FIX ) 
/***********************************************************************
   set_Katm():

       -- sets the EOS constant used for Font's fix. 

       -- see utoprim_1dfix1.c and utoprim_1dvsq2fix1.c  for more
           information. 

       -- uses the initial floor values of rho/u determined by fixup1zone()

       -- we assume here that Constant X1,r is independent of theta,X2

***********************************************************************/
void set_Katm( void )
{
  int i, j, k, m, G_type ;
  double prim[NPR], G_tmp, K, *Katm_fullgrid;
  double r,th,phi,X[NDIM];

  G_type = get_G_ATM( &G_tmp );

  fflush(stdout);
  if(MASTER==mpi_rank) fprintf(stdout,"G_tmp = %26.20e \n", G_tmp );
  fflush(stdout);

  j = 0;
  k = 0;
  for( i = 0 ; i < N1; i++ ) {
    
    PLOOP prim[m] = 0.;
    prim[RHO] = prim[UU] = -1.;

    coord(i,j,k,CENT,X) ;
    get_phys_coord(i,j,k,&r,&th,&phi) ;

    fixup1zone( i, j, k, prim );
    
    Katm[i] = (gam - 1.) * prim[UU] / pow( prim[RHO], G_tmp ) ;
    
    //fflush(stdout);
    //if(mpi_rank==2) fprintf(stdout,"Katm[%d] = %26.20e, r=%26.20e, x1=%26.20e \n", i+mpi_startn[1], Katm[i], r, X[1] );
    //fflush(stdout);
    
  }
#ifdef MPI
  //allocate 1D array that spans the full radial grid used by all MPI processes
  Katm_fullgrid = calloc(mpi_ntot[1],sizeof(double));
  
  //copy the jglobal = kglobal = 0 values that were filled correctly above into Katm_fullgrid[]
  if(mpi_coords[2]==0 && mpi_coords[3]==0){
    for (i=0; i<N1; i++) {
      Katm_fullgrid[i+mpi_startn[1]] = Katm[i];
    }
  }
  //exchange the info between the MPI processes
  MPI_Allreduce(MPI_IN_PLACE,Katm_fullgrid,mpi_ntot[1],MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  //copy back the relevant piece from Katm_fullgrid[] to Katm[]
  for (i=0; i<N1; i++) {
    Katm[i] = Katm_fullgrid[i+mpi_startn[1]];
  }
  free(Katm_fullgrid);
#endif
  return;
}
#endif


#undef FLOOP 
