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

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************

utoprim_1dfix1.c: 
---------------

  -- uses eq. (27) of Noble  et al. or the "momentum equation" and ignores
        the energy equation (29) in order to use the additional EOS;
  
    Uses the 1D_W method: 
       -- solves for one independent variable (W) via a 1D
          Newton-Raphson method 
       -- can be used (in principle) with a general equation of state. 

  -- Currently returns with an error state (>0) if a negative rest-mass
      density or internal energy density is calculated.  You may want 
      to change this aspect of the code so that it still calculates the 
      velocity and so that you can floor the densities.  If you want to 
      change this aspect of the code please comment out the "return(retval)"
      statement after "retval = 5;" statement in Utoprim_new_body();

******************************************************************************/


#include "u2p_util.h"

#define NEWT_DIM 1

#define LTRACE 0

/* these variables need to be shared between the functions
   Utoprim_1D, residual, and utsq */
static FTYPE Bsq,QdotBsq,Qtsq,Qdotn,D, K_atm ;
static FTYPE W_for_gnr2, rho_for_gnr2, W_for_gnr2_old, rho_for_gnr2_old;
#pragma omp threadprivate(Bsq,QdotBsq,Qtsq,Qdotn,D, K_atm, W_for_gnr2, rho_for_gnr2, W_for_gnr2_old, rho_for_gnr2_old)

// Declarations: 
static FTYPE vsq_calc(FTYPE W);
static FTYPE u_of_p(FTYPE p);
static FTYPE pressure_of_rho(FTYPE rho0);
static int Utoprim_new_body(FTYPE U[], FTYPE gcov[NDIM][NDIM],  FTYPE gcon[NDIM][NDIM], FTYPE gdet,  FTYPE prim[]);
static void func_1d_orig1(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n);
static void func_1d_orig2(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n);
static int general_newton_raphson( FTYPE x[], int n, void (*funcd) (FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], FTYPE *, FTYPE *, int) );
static void func_gnr2_rho(FTYPE x[], FTYPE dx[], FTYPE resid[], FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n);
static int gnr2( FTYPE x[], int n, void (*funcd) (FTYPE [], FTYPE [], FTYPE [], FTYPE [][NEWT_DIM], FTYPE *, FTYPE *, int) );

/**********************************************************************/
/******************************************************************

  Utoprim_1dfix1():

  -- Driver for new prim. var. solver.  The driver just translates
     between the two sets of definitions for U and P.  The user may 
     wish to alter the translation as they see fit.  

     It assumes that on input/output:


              /  rho u^t           \
         U =  |  T^t_t   + rho u^t |  sqrt(-det(g_{\mu\nu}))
              |  T^t_\mu           |  
              \   B^i              /


             /    rho        \
	 P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /


     ala HARM. 

   Arguments:
       U[NPR]    = conserved variables (current values on input/output);
       gcov[NDIM][NDIM] = covariant form of the metric ;
       gcon[NDIM][NDIM] = contravariant form of the metric ;
       gdet             = sqrt( - determinant of the metric) ;
       prim[NPR] = primitive variables (guess on input, calculated values on 
                                        output if there are no problems);
  
   -- NOTE: for those using this routine for special relativistic MHD and are
            unfamiliar with metrics, merely set 
              gcov = gcon = diag(-1,1,1,1)  and gdet = 1.  ;

   -- different from Utoprim_2d() in that it also assumes an EOS of 
            P = K rho^\Gamma  or   P = K rho    
        depending on the value of G_ATM flag in u2p_defs.h

******************************************************************/

int Utoprim_1dfix1(FTYPE U[NPR], FTYPE gcov[NDIM][NDIM], FTYPE gcon[NDIM][NDIM], 
	       FTYPE gdet, FTYPE prim[NPR], FTYPE K )
{

  FTYPE U_tmp[NPR], prim_tmp[NPR];
  int i, j, ret; 
  FTYPE alpha;


  /* First update the primitive B-fields */
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] / gdet ;

  if( U[0] <= 0. ) {
    return(-100);
  }

  K_atm = K ; 

  /* Set the geometry variables: */
  alpha = 1.0/sqrt(-gcon[0][0]);
  
  /* Transform the CONSERVED variables into the new system */
  U_tmp[RHO] = alpha * U[RHO] / gdet;
  U_tmp[UU]  = alpha * (U[UU] - U[RHO])  / gdet ;


    
  for( i = UTCON1; i <= UTCON3; i++ ) {
    U_tmp[i] = alpha * U[i] / gdet ;
  }
  for( i = BCON1; i <= BCON3; i++ ) {
    U_tmp[i] = alpha * U[i] / gdet ;
  }

  /* Transform the PRIMITIVE variables into the new system */
  for( i = 0; i < BCON1; i++ ) {
    prim_tmp[i] = prim[i];
  }
  for( i = BCON1; i <= BCON3; i++ ) {
    prim_tmp[i] = alpha*prim[i];
  }

  ret = Utoprim_new_body(U_tmp, gcov, gcon, gdet, prim_tmp);

  /* Transform new primitive variables back if there was no problem : */ 
  if( ret == 0 ) { 
    for( i = 0; i < BCON1; i++ ) {
      prim[i] = prim_tmp[i];
    }
  }

#if( DOKTOT )
  prim[KTOT] = U[KTOT]/U[RHO];
#endif

  

  
  
  return( ret ) ;

}


/**********************************************************************/
/**********************************************************************************

  Utoprim_new_body():

     -- Attempt an inversion from U to prim using the initial guess prim.

     -- This is the main routine that calculates auxiliary quantities for the 
        Newton-Raphson routine. 

  -- assumes that 
             /  rho gamma        \
         U = |  alpha T^t_\mu    |
             \  alpha B^i        /



               /    rho        \
	prim = |    uu         |
               | \tilde{u}^i   |
               \  alpha B^i   /


return:  (i*100 + j)  where 
         i = 0 ->  Newton-Raphson solver either was not called (yet or not used) 
                   or returned successfully;
             1 ->  Newton-Raphson solver did not converge to a solution with the 
                   given tolerances;
             2 ->  Newton-Raphson procedure encountered a numerical divergence 
                   (occurrence of "nan" or "+/-inf" ;
	     
         j = 0 -> success 
             1 -> failure: some sort of failure in Newton-Raphson; 
             2 -> failure: utsq<0 w/ initial p[] guess;
	     3 -> failure: W<0 or W>W_TOO_BIG
             4 -> failure: v^2 > 1 
             5 -> failure: rho,uu <= 0 ;

**********************************************************************************/

static int Utoprim_new_body(FTYPE U[NPR], FTYPE gcov[NDIM][NDIM], 
			    FTYPE gcon[NDIM][NDIM], FTYPE gdet,  FTYPE prim[NPR])
{

  FTYPE x_1d[1];
  FTYPE QdotB,Bcon[NDIM],Bcov[NDIM],Qcov[NDIM],Qcon[NDIM],ncov[NDIM],ncon[NDIM],Qsq,Qtcon[NDIM];
  FTYPE rho0,u,p,w,gammasq,gamma,gtmp,W_last,W,utsq,vsq,tmpdiff ;
    FTYPE alpha, ucovt, utsqp1, aco, bco, cco, pevar, agame, the, utt;
  int i,j, retval, i_increase ;
    double dummy;

  // Assume ok initially:
  retval = 0;

  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;

  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:
  Bcon[0] = 0. ;
  for(i=1;i<4;i++) Bcon[i] = U[BCON1+i-1] ;

  lower_g(Bcon,gcov,Bcov) ;

  for(i=0;i<4;i++) Qcov[i] = U[QCOV0+i] ;
  raise_g(Qcov,gcon,Qcon) ;


  Bsq = 0. ;
  for(i=1;i<4;i++) Bsq += Bcon[i]*Bcov[i] ;

  QdotB = 0. ;
  for(i=0;i<4;i++) QdotB += Qcov[i]*Bcon[i] ;
  QdotBsq = QdotB*QdotB ;
  
  ncov_calc(gcon,ncov) ;
  raise_g(ncov,gcon,ncon);

  Qdotn = Qcon[0]*ncov[0] ;

  Qsq = 0. ;
  for(i=0;i<4;i++) Qsq += Qcov[i]*Qcon[i] ;

  Qtsq = Qsq + Qdotn*Qdotn ;

  D = U[RHO] ;

  /* calculate W from last timestep and use for guess */
  utsq = 0. ;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++) utsq += gcov[i][j]*prim[UTCON1+i-1]*prim[UTCON1+j-1] ;


  if( (utsq < 0.) && (fabs(utsq) < 1.0e-13) ) { 
    utsq = fabs(utsq);
  }
  if(utsq < 0. || utsq > UTSQ_TOO_BIG) {
    retval = 2;
    return(retval) ;
  }

  gammasq = 1. + utsq ;
  gamma  = sqrt(gammasq);
	
  // Always calculate rho from D and gamma so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(vsq) . 
  rho0 = D / gamma ;
  p = pressure_of_rho( rho0 );
  u = u_of_p(p);
  w = rho0 + u + p ;


#if(LTRACE)
  if( rho0 <= 0. ) { 
    fprintf(stderr,"beg neg rho fix1, rho,D,gamma = %26.20e %26.20e %26.20e \n", rho0, D, gamma); fflush(stderr);
    rho0 = fabs(rho0);
  }
#endif 
    
  W_last = w*gammasq ;

  // Make sure that W is large enough so that v^2 < 1 : 
  i_increase = 0;
  while( (( W_last*W_last*W_last * ( W_last + 2.*Bsq ) 
	    - QdotBsq*(2.*W_last + Bsq) ) <= W_last*W_last*(Qtsq-Bsq*Bsq))
	 && (i_increase < 10) ) {
    W_last *= 10.;
    i_increase++;
  }
  
  W_for_gnr2 = W_for_gnr2_old = W_last;
  rho_for_gnr2 = rho_for_gnr2_old = rho0;

  // Calculate W: 
  x_1d[0] = W_last;
  
#if( USE_ISENTROPIC )   
  retval = general_newton_raphson( x_1d, 1, func_1d_orig1);
#else
  retval = general_newton_raphson( x_1d, 1, func_1d_orig2);
#endif

  W = x_1d[0];

  /* Problem with solver, so return denoting error before doing anything further */
  if( (retval != 0) || (W == FAIL_VAL) ) {
    retval = retval*100+1;
#if(LTRACE)
    fprintf(stderr,"fix1: retval, W, rho = %d %26.20e %26.20e \n", retval,W, rho_for_gnr2);fflush(stderr);
#endif
    return(retval);
  }
  else{
    if(W <= 0. || W > W_TOO_BIG) {
      retval = 3;
#if(LTRACE)
      fprintf(stderr,"fix1: retval, W, rho = %d %26.20e %26.20e \n", retval,W, rho_for_gnr2);fflush(stderr);
#endif
      return(retval) ;
    }
  }


  // Calculate v^2 : 
  vsq = vsq_calc(W) ;
  if( vsq >= 1. ) {
    retval = 4;
#if(LTRACE)
    fprintf(stderr,"fix1: retval, W, rho,vsq = %d %26.20e %26.20e %26.20e \n", retval,W, rho_for_gnr2,vsq);fflush(stderr);
#endif
    return(retval) ;
  }

  // Recover the primitive variables from the scalars and conserved variables:
  gtmp = sqrt(1. - vsq);
  gamma = 1./gtmp ;
  rho0 = D * gtmp;

  w = W * (1. - vsq) ;

  p = pressure_of_rho( rho0 );
  u = u_of_p(p);

  // User may want to handle this case differently, e.g. do NOT return upon 
  // a negative rho/u, calculate v^i so that rho/u can be floored by other routine:
  if( (rho0 <= 0.) || (u <= 0.) ) { 
    retval = 5;
#if(LTRACE)
    fprintf(stderr,"fix1: retval, W, rho,vsq,u = %d %26.20e %26.20e %26.20e %26.20e \n", retval,W, rho0,vsq,u);fflush(stderr);
#endif
    return(retval) ;
  }

   prim[RHO] = rho0 ;
   prim[UU] = u ;
    
  for(i=1;i<4;i++)  Qtcon[i] = Qcon[i] + ncon[i] * Qdotn;
  for(i=1;i<4;i++) prim[UTCON1+i-1] = gamma/(W+Bsq) * ( Qtcon[i] + QdotB*Bcon[i]/W ) ;
	
    
  /* set field components */
  for(i = BCON1; i <= BCON3; i++) prim[i] = U[i] ;


  /* done! */
  return(retval) ;

}


/**********************************************************************/
/****************************************************************************
   vsq_calc(): 
    
      -- evaluate v^2 (spatial, normalized velocity) from 
            W = \gamma^2 w 

****************************************************************************/
static FTYPE vsq_calc(FTYPE W)
{
	FTYPE Wsq,Xsq;
	
	Wsq = W*W ;
	Xsq = (Bsq + W) * (Bsq + W);

	return(  ( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq) );
}

/**********************************************************************/
/****************************************************************************
   dvsq_dW(): 
    
      -- evaluate the partial derivative of v^2 w.r.t. W

****************************************************************************/
static FTYPE dvsq_dW(FTYPE W)
{
	FTYPE W3,X3,X;
	
	X = Bsq + W;
	W3 = W*W*W ;
	X3 = X*X*X;

	return( -2.*( Qtsq/X3  +  QdotBsq * (3*W*X + Bsq*Bsq) / ( W3 * X3 )  )  );
}


/**********************************************************************/
/************************************************************

  general_newton_raphson(): 

    -- performs Newton-Rapshon method on an arbitrary system.

    -- inspired in part by Num. Rec.'s routine newt();

    Arguements: 

       -- x[]   = set of independent variables to solve for;
       -- n     = number of independent variables and residuals;
       -- funcd = name of function that calculates residuals, etc.;

*****************************************************************/
static int general_newton_raphson( FTYPE x[], int n, 
				   void (*funcd) (FTYPE [], FTYPE [], FTYPE [], 
						  FTYPE [][NEWT_DIM], FTYPE *, 
						  FTYPE *, int) )
{
  FTYPE f, df, dx[NEWT_DIM], x_old[NEWT_DIM], resid[NEWT_DIM], 
    jac[NEWT_DIM][NEWT_DIM];
  FTYPE errx, x_orig[NEWT_DIM];
  int    n_iter, id, jd, i_extra, doing_extra;
  FTYPE dW,dvsq,vsq_old,vsq,W,W_old;

  int   keep_iterating, i_increase;


  // Initialize various parameters and variables:
  errx = 1. ; 
  df = f = 1.;
  i_extra = doing_extra = 0;
  for( id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;

  vsq_old = vsq = W = W_old = 0.;

  n_iter = 0;


  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    (*funcd) (x, dx, resid, jac, &f, &df, n);  /* returns with new dx, f, df */

    /* Save old values before calculating the new: */
    errx = 0.;
    for( id = 0; id < n ; id++) {
      x_old[id] = x[id] ;
    }

    /* don't use line search : */
    for( id = 0; id < n ; id++) {
      x[id] += dx[id]  ;
    }

    //METHOD specific:
    i_increase = 0;
    while( (( x[0]*x[0]*x[0] * ( x[0] + 2.*Bsq ) - 
	      QdotBsq*(2.*x[0] + Bsq) ) <= x[0]*x[0]*(Qtsq-Bsq*Bsq))
	   && (i_increase < 10) ) {
      x[0] -= (1.*i_increase) * dx[0] / 10. ;
      i_increase++;
    }

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/

    /* For the new criterion, always look at error in "W" : */
    // METHOD specific:
    errx  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);


    /****************************************/
    /* Make sure that the new x[] is physical : */
    /****************************************/
    x[0] = fabs(x[0]);


    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/
    
    if( (fabs(errx) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++ ;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0)) || 
	(i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)


  /*  Check for bad untrapped divergences : */
  if( (finite(f)==0) || (finite(df)==0) || (finite(x[0])==0)  ) {
#if(LTRACE)
    fprintf(stderr,"\ngnr not finite, f,df,x_o,x,W_o,W,rho_o,rho = %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e \n",
	    f,df,x[0],x_old[0],W_for_gnr2_old,W_for_gnr2,rho_for_gnr2_old,rho_for_gnr2); fflush(stderr); 
#endif
    return(2);
  }


  if( fabs(errx) > MIN_NEWT_TOL){
#if(LTRACE)
    fprintf(stderr," totalcount = %d   0   %d  %26.20e \n",n_iter,i_extra,errx); fflush(stderr); 
#endif
    return(1);
  }
  if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    //fprintf(stderr," totalcount = %d   1   %d  %26.20e \n",n_iter,i_extra,errx); fflush(stderr);
    return(0);
  }
  if( fabs(errx) <= NEWT_TOL ){
    //fprintf(stderr," totalcount = %d   2   %d  %26.20e \n",n_iter,i_extra,errx); fflush(stderr); 
    return(0);
  }

  return(0);

}


/***********************************************************/
/********************************************************************** 

  gnr2()

    -- used to calculate rho from W

*****************************************************************/
static int gnr2( FTYPE x[], int n, 
			    void (*funcd) (FTYPE [], FTYPE [], FTYPE [], 
					 FTYPE [][NEWT_DIM],FTYPE *,FTYPE *,int) )
{
  FTYPE f, df, dx[NEWT_DIM], x_old[NEWT_DIM], resid[NEWT_DIM], 
    jac[NEWT_DIM][NEWT_DIM];
  FTYPE errx, x_orig[NEWT_DIM];
  int    n_iter, id,jd, i_extra, doing_extra;
  FTYPE dW,dvsq,vsq_old,vsq,W,W_old;

  int   keep_iterating;


  // Initialize various parameters and variables:
  errx = 1. ; 
  df = f = 1.;
  i_extra = doing_extra = 0;
  for( id = 0; id < n ; id++)  x_old[id] = x_orig[id] = x[id] ;

  n_iter = 0;


  /* Start the Newton-Raphson iterations : */
  keep_iterating = 1;
  while( keep_iterating ) { 

    (*funcd) (x, dx, resid, jac, &f, &df, n);  /* returns with new dx, f, df */

    /* Save old values before calculating the new: */
    errx = 0.;
    for( id = 0; id < n ; id++) {
      x_old[id] = x[id] ;
    }

    /* Make the newton step: */
    for( id = 0; id < n ; id++) {
      x[id] += dx[id] ;
    }

    /* Calculate the convergence criterion */
    for( id = 0; id < n ; id++) {
      errx  += (x[id]==0.) ?  fabs(dx[id]) : fabs(dx[id]/x[id]);
    }
    errx /= 1.*n;

    /* Make sure that the new x[] is physical : */
    // METHOD specific:
    x[0] = fabs(x[0]);


    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*   before stopping                                                         */
    if( (fabs(errx) <= NEWT_TOL2) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 0;
    }

    if( doing_extra == 1 ) i_extra++ ;

    // See if we've done the extra iterations, or have done too many iterations:
    if( ((fabs(errx) <= NEWT_TOL2)&&(doing_extra == 0)) || 
	(i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }  


  /*  Check for bad untrapped divergences : */
  if( (finite(f)==0) || (finite(df)==0) || (finite(x[0])==0)  ) {
#if(LTRACE)
    fprintf(stderr,"\ngnr2 not finite, f,df,x_o,x,W_o,W,rho_o,rho = %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e %26.20e \n",
	    f,df,x[0],x_old[0],W_for_gnr2_old,W_for_gnr2,rho_for_gnr2_old,rho_for_gnr2); fflush(stderr); 
#endif
    return(2);
  }

  // Return in different ways depending on whether a solution was found:
  if( fabs(errx) > MIN_NEWT_TOL){
#if(LTRACE)
    fprintf(stderr," totalcount2 = %d   0   %d  %26.20e \n",n_iter,i_extra,errx); fflush(stderr); 
#endif
    return(1);
  }
  if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    //fprintf(stderr," totalcount2 = %d   1   %d  %26.20e \n",n_iter,i_extra,errx); fflush(stderr); 
    return(0);
  }
  if( fabs(errx) <= NEWT_TOL ){
    //fprintf(stderr," totalcount2 = %d   2   %d  %26.20e \n",n_iter,i_extra,errx); fflush(stderr); 
    return(0);
  }

  return(0);

}



/**********************************************************************/
/*********************************************************************************
   func_1d_orig[1,2](): 

        -- calculates the residuals, and Newton step for general_newton_raphson();
        -- for this method, x=W here;

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
 *********************************************************************************/
//isentropic version:   eq.  (27)
static void func_1d_orig1(FTYPE x[], FTYPE dx[], FTYPE resid[], 
			 FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n)
{
  int retval, ntries;
  FTYPE  Dc, t1, t10,  t2 ,  t21,  t23,  t26,  t29,  t3 ,  t30;
  FTYPE  t32,  t33,  t34,  t38,  t5 ,  t51,  t67, t8, W, x_rho[1], rho, rho_g ;

 
  W  = x[0];
  W_for_gnr2_old = W_for_gnr2;
  W_for_gnr2 = W;

  // get rho from NR:
  rho_g = x_rho[0] = rho_for_gnr2;
  
  ntries = 0;
  while (  (retval = gnr2( x_rho, 1, func_gnr2_rho)) &&  ( ntries++ < 10 )  ) { 
    rho_g *= 10.;
    x_rho[0] = rho_g;
  }

#if(LTRACE)
  if( x_rho[0] <= 0. ) { 
    fprintf(stderr,"gnr2 neg rho = %d ,rho_n,rho,rho_o,W,W_o = %26.20e %26.20e %26.20e %26.20e %26.20e \n", retval, x_rho[0], rho_for_gnr2, rho_for_gnr2_old, x[0], W_for_gnr2_old);
    fflush(stderr);
  }

  if( retval ) { 
    fprintf(stderr,"gnr2 retval = %d ,rho_n,rho,rho_o,W,W_o = %26.20e %26.20e %26.20e %26.20e %26.20e \n", retval, x_rho[0], rho_for_gnr2, rho_for_gnr2_old, x[0], W_for_gnr2_old);
    fflush(stderr);
  }
#endif 

  rho_for_gnr2_old = rho_for_gnr2; 
  rho = rho_for_gnr2 = x_rho[0];

  Dc = D;
  t1 = Dc*Dc;
  t2 = QdotBsq*t1;
  t3 = t2*Bsq;
  t5 = Bsq*Bsq;
  t8 = t1*Bsq;
  t10 = t1*W;
  t21 = W*W;
  t23 = rho*rho;
  t26 = 1/t1;
  resid[0] = (t3+(2.0*t2+((Qtsq-t5)*t1
			  +(-2.0*t8-t10)*W)*W)*W+(t5+(2.0*Bsq+W)*W)*t21*t23)*t26/t21;
  t29 = t1*t1;
  t30 = QdotBsq*t29;
  t32 = GAMMA*K_atm;
  t33 = pow(rho,1.0*GAMMA);
  t34 = t32*t33;
  t38 = t23 * t33;
  t51 = GAMMA*t1*K_atm*t33;
  t67 = t21*W;
  jac[0][0] = -2.0*(t30*Bsq*t34+(t30*t34
			   +((-t38*Bsq*t32+Bsq*GAMMA*t1*K_atm*t33)*t1
			     +(-t38*GAMMA*K_atm+t51)*t1*W)*t21)*W
	      +((-t3+(-t2+(-t8-t10)*t21)*W)*W+(-t5-Bsq*W)*t67*t23)*t23)*t26/(t51-W*t23)/t67;

  dx[0] = -resid[0]/jac[0][0];

  *f = 0.5*resid[0]*resid[0];
  *df = -2. * (*f);

  return;

}

//isothermal version:  eq. (27)
static void func_1d_orig2(FTYPE x[], FTYPE dx[], FTYPE resid[], 
			  FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n)
{
  double Dc ,   t1 ,   t10,   t2 ,   t21,   t23,   t26,   t3 ,   t5 ,   t8, W, rho ;

  W  = x[0];
  Dc = D;
  t1 = Dc*Dc;
  t2 = QdotBsq*t1;
  t3 = t2*Bsq;
  t5 = Bsq*Bsq;
  t8 = t1*Bsq;
  t10 = t1*W;
  t21 = W*W;
  rho = t1 * ( 1. + GAMMA*K_atm/(GAMMA-1.) ) / W;
  t23 = rho*rho;
  t26 = 1/t1;
  resid[0] = (t3+(2.0*t2+((Qtsq-t5)*t1+(-2.0*t8-t10)*W)*W)*W+(t5+(2.0*Bsq+W)*W)*t21*t23)*t26/t21;
  jac[0][0] = -2.0*(t3+(t2+(t8+t10)*t21)*W+(t5+Bsq*W)*t21*t23)*t26/t21/W;

  dx[0] = -resid[0]/jac[0][0];

  *f = 0.5*resid[0]*resid[0];
  *df = -2. * (*f);

  return;
}



/**********************************************************************/
/*********************************************************************************
   func_gnr2_rho():

        -- residual/jacobian routine to calculate rho from W via the polytrope:

        W  =  ( 1 + GAMMA * K_atm * rho^(GAMMA-1)/(GAMMA-1) ) D^2 / rho

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
 *********************************************************************************/
// for the isentropic version:   eq.  (27)
static void func_gnr2_rho(FTYPE x[], FTYPE dx[], FTYPE resid[], 
			 FTYPE jac[][NEWT_DIM], FTYPE *f, FTYPE *df, int n)
{

  FTYPE A, B, C, rho, W, B0;
  
  A = D*D;
  B0 = A * GAMMA * K_atm ;
  B  =  B0 / (GAMMA - 1.);
  rho = x[0];
  W = W_for_gnr2;
  C = pow( rho, GAMMA - 1. );


  resid[0] = rho*W - A - B*C ;
		
  jac[0][0] = W  -  B0 * C / rho;

  dx[0] = -resid[0]/jac[0][0];
  *f = 0.5*resid[0]*resid[0];
  *df = -2. * (*f);

  return;

}

int get_G_ATM( FTYPE *g_tmp )  
{
  *g_tmp = G_ATM ;

  return( USE_ISENTROPIC );

}


/**********************************************************************
 ********************************************************************** 
   
 The following routines specify the equation of state.  All routines 
  above here should be indpendent of EOS.  If the user wishes 
  to use another equation of state, the below functions must be replaced 
  by equivalent routines based upon the new EOS. 

 **********************************************************************
**********************************************************************/

/* 
pressure as a function of rho0 and w = rho0 + u + p 
this is used by primtoU and Utoprim_1D
*/
static FTYPE pressure_of_rho(FTYPE rho0)
{

  return( K_atm * pow( rho0, G_ATM )  );

}

/* 
internal energy density as a function of the pressure
*/
static FTYPE u_of_p(FTYPE p)
{
  
  return( p / (GAMMA - 1.) ) ;

}


/****************************************************************************** 
             END   OF   UTOPRIM_1D.C
 ******************************************************************************/

