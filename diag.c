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

/* all diagnostics subroutine */

void append_rank(char *name)
{
  char suff[MAXLEN];
  //file suffixes that will be different by MPI rank
  sprintf(suff, "_%04d", mpi_rank);
  strcat(name,suff);
}

void diag(call_code)
int call_code ;
{
  int di, dj, dk;
  char efnam[MAXLEN], dfnam[MAXLEN],ifnam[MAXLEN] ;
  int i,j,k,m ;
  FILE *dump_file;
  double U[NPR],pp,e,rmed,divb,divbmax,divbmax_global,e_fin,m_fin,gamma ;
  struct of_geom geom ;
  struct of_state q ;
  int imax,jmax,kmax ;
  int istart,istop,jstart,jstop,kstart,kstop;
  
  static double e_init,m_init ;
  static FILE *ener_file = NULL;

#if(DOENER)
  if(call_code==INIT_OUT || NULL == ener_file) {
    /* set things up */
    strcpy(efnam, "ener.out");
    append_rank(efnam);
    ener_file = fopen(efnam,"a") ;
    if(ener_file==NULL) {
      fprintf(stderr,"error opening energy output file, %s\n", efnam) ;
      exit(1) ;
    }
  }
#endif
  
  /* calculate conserved quantities */
  if(call_code==INIT_OUT ||
     call_code==LOG_OUT ||
     call_code==DIVB_OUT ||
     (call_code==FINAL_OUT &&
     !failed)) {
    pp = 0. ;
    e = 0. ;
    rmed = 0. ;
    divbmax = 0. ;
    imax = 0 ;
    jmax = 0 ;
    kmax = 0;
    di = (N1>1);
    dj = (N2>1);
    dk = (N3>1);
    //all cells not immediately adjacent to physical boundaries
    istart = is_physical_bc(1, 0);
    istop  = N1 - 1 - is_physical_bc(1, 1);
    jstart = 0*is_physical_bc(2, 0);
    jstop  = N2 - 0*is_physical_bc(2, 1);
    kstart = 0*is_physical_bc(3, 0);
    kstop  = N3 - 0*is_physical_bc(3, 1);
    ZSLOOP(istart,istop,jstart,jstop,kstart,kstop) {
      get_geometry(i,j,k,CENT,&geom) ;
      get_state(p[i][j][k],&geom,&q) ;
      primtoU(p[i][j][k],&q,&geom,U) ;
      
      rmed += U[RHO]*dV ;
      pp += U[U3]*dV ;
      e += U[UU]*dV ;

      /* flux-ct defn */
      divb = fabs(
#if(N1>1)
                  0.25*(
                        + p[i  ][j   ][k   ][B1]*gdet[i  ][j   ][k   ][CENT]
                        + p[i  ][j   ][k-dk][B1]*gdet[i  ][j   ][k-dk][CENT]
                        + p[i  ][j-dj][k   ][B1]*gdet[i  ][j-dj][k   ][CENT]
                        + p[i  ][j-dj][k-dk][B1]*gdet[i  ][j-dj][k-dk][CENT]
                        - p[i-1][j   ][k   ][B1]*gdet[i-1][j   ][k   ][CENT]
                        - p[i-1][j   ][k-dk][B1]*gdet[i-1][j   ][k-dk][CENT]
                        - p[i-1][j-dj][k   ][B1]*gdet[i-1][j-dj][k   ][CENT]
                        - p[i-1][j-dj][k-dk][B1]*gdet[i-1][j-dj][k-dk][CENT]
                        )/dx[1]
#endif
#if(N2>1)
                 +0.25*(
                        + p[i   ][j  ][k   ][B2]*gdet[i   ][j  ][k   ][CENT]
                        + p[i   ][j  ][k-dk][B2]*gdet[i   ][j  ][k-dk][CENT]
                        + p[i-di][j  ][k   ][B2]*gdet[i-di][j  ][k   ][CENT]
                        + p[i-di][j  ][k-dk][B2]*gdet[i-di][j  ][k-dk][CENT]
                        - p[i   ][j-1][k   ][B2]*gdet[i   ][j-1][k   ][CENT]
                        - p[i   ][j-1][k-dk][B2]*gdet[i   ][j-1][k-dk][CENT]
                        - p[i-di][j-1][k   ][B2]*gdet[i-di][j-1][k   ][CENT]
                        - p[i-di][j-1][k-dk][B2]*gdet[i-di][j-1][k-dk][CENT]
                       )/dx[2]
#endif
#if(N3>1)
                 +0.25*(
                        + p[i   ][j   ][k  ][B3]*gdet[i   ][j   ][k  ][CENT]
                        + p[i-di][j   ][k  ][B3]*gdet[i-di][j   ][k  ][CENT]
                        + p[i   ][j-dj][k  ][B3]*gdet[i   ][j-dj][k  ][CENT]
                        + p[i-di][j-dj][k  ][B3]*gdet[i-di][j-dj][k  ][CENT]
                        - p[i   ][j   ][k-1][B3]*gdet[i   ][j   ][k-1][CENT]
                        - p[i-di][j   ][k-1][B3]*gdet[i-di][j   ][k-1][CENT]
                        - p[i   ][j-dj][k-1][B3]*gdet[i   ][j-dj][k-1][CENT]
                        - p[i-di][j-dj][k-1][B3]*gdet[i-di][j-dj][k-1][CENT]
                   )/dx[3]
#endif
                  );
      if(divb > divbmax) {
        imax = i+mpi_startn[1] ;
        jmax = j+mpi_startn[2] ;
        kmax = k+mpi_startn[3] ;
        divbmax = divb ;
      }
    }
  }
  
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(&divbmax,&divbmax_global,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&e,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&rmed,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&pp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
  divbmax_global = divbmax;
#endif
  
  if(call_code == INIT_OUT) {
    e_init = e ;
    m_init = rmed ;
  }
  
  if(MASTER == mpi_rank && call_code == FINAL_OUT) {
    e_fin = e ;
    m_fin = rmed ;
    fprintf(stderr,"\n\nEnergy: ini,fin,del: %g %g %g\n",
            e_init,e_fin,(e_fin-e_init)/e_init) ;
    fprintf(stderr,"mass: ini,fin,del: %g %g %g\n",
            m_init,m_fin,(m_fin-m_init)/m_init) ;
  }
  
  if(call_code == INIT_OUT ||
     call_code == LOG_OUT ||
     call_code == FINAL_OUT ||
     call_code == DIVB_OUT) {
    if (divbmax==divbmax_global) {
      fprintf(stderr,"LOG      t=%g \t divbmax: %d %d %d %g\n",
            t,imax,jmax,kmax,divbmax) ;
    }
  }

#if(DOENER)
  if(call_code == INIT_OUT ||
     call_code == LOG_OUT ||
     call_code == FINAL_OUT) {
    fprintf(ener_file,"%10.5g %10.5g %10.5g %10.5g %15.8g %15.8g ",
            t,rmed,pp,e,p[N1/2][N2/2][N3/2][UU]*pow(p[N1/2][N2/2][N3/2][RHO],-gam),
            p[N1/2][N2/2][N3/2][UU]) ;
    fprintf(ener_file,"%15.8g %15.8g %15.8g ",mdot,edot,ldot) ;
    fprintf(ener_file,"\n") ;
    fflush(ener_file) ;
  }
#endif
  
  /* gdump only at code start */
  if(call_code == INIT_OUT) {
    /* make grid dump file */
    
    gdump(0) ;
    gdump2(0) ;
  }
  
  /* dump at regular intervals */
  if(call_code == INIT_OUT ||
     call_code == DUMP_OUT ||
     call_code == FINAL_OUT) {
    /* make regular dump file */
    dump(dump_cnt,0) ;
    dump_cnt++ ;
  }

  /* dump at regular intervals */
  if(call_code == DIVB_OUT) {
    /* make regular dump file */
    dump(-dump_cnt,0) ;
    //dump_cnt++ ;
  }
  
  /* rdump at regular intervals */
  if(call_code == INIT_OUT ||
     call_code == RDUMP_OUT ||
     call_code == FINAL_OUT) {
    /* make rdump file */
    //increment rdump count *before* writing the restart dump file
    //so that get the correct counter value if restart from the file
    rdump_cnt++ ;
    restart_write(rdump_cnt-1);
    
  }

  /* image dump at regular intervals */
  if(call_code == IMAGE_OUT ||
     call_code == INIT_OUT ||
     call_code == FINAL_OUT) {
    
    fdump( image_cnt );
    
    image_cnt++ ;
  }
}


/** some diagnostic routines **/


void fail(int fail_type)
{

	failed = 1 ;

	fprintf(stderr,"\n\nfail: %d %d %d %d\n",icurr,jcurr,kcurr,fail_type) ;

	area_map(icurr,jcurr,kcurr,p) ;
	
	fprintf(stderr,"fail\n") ;

	diag(FINAL_OUT) ;

#ifdef MPI
        MPI_Finalize();
#endif
	/* for diagnostic purposes */
	exit(0) ;
}



/* map out region around failure point */
void area_map(int i, int j, int k, double prim[][N2M][N3M][NPR])
{
	int m ;

	fprintf(stderr,"area map\n") ;

	PLOOP {
		fprintf(stderr,"variable %d \n",m) ;
		fprintf(stderr,"i = \t %12d %12d %12d\n",i-1,i,i+1) ;
		fprintf(stderr,"j = %d \t %12.5g %12.5g %12.5g\n",
				j+1,
				prim[i-1][j+1][k][m],
				prim[i][j+1][k][m],
				prim[i+1][j+1][k][m]) ;
		fprintf(stderr,"j = %d \t %12.5g %12.5g %12.5g\n",
				j,
				prim[i-1][j][k][m],
				prim[i][j][k][m],
				prim[i+1][j][k][m]) ;
		fprintf(stderr,"j = %d \t %12.5g %12.5g %12.5g\n",
				j-1,
				prim[i-1][j-1][k][m],
				prim[i][j-1][k][m],
				prim[i+1][j-1][k][m]) ;
	}

	/* print out other diagnostics here */

}

/* evaluate fluxed based diagnostics; put results in
 * global variables */

void diag_flux(double F1[][N2M][N3M][NPR], double F2[][N2M][N3M][NPR])
{
	int j,k ;

        mdot = edot = ldot = 0. ;
//#pragma omp parallel for schedule(static,MY_MAX(N2*N3/nthreads,1)) collapse(2) \
//    default(none) \
//    shared(F1,F2,dx,nthreads) \
//    private(j,k) \
//    reduction(+:mdot) reduction(-:edot) reduction(+:ldot)
        for(j=0;j<N2;j++) {
          for(k=0;k<N3;k++) {
                mdot += F1[0][j][k][RHO]*dx[2]*dx[3] ;
                edot -= (F1[0][j][k][UU] - F1[0][j][k][RHO])*dx[2]*dx[3] ;
                ldot += F1[0][j][k][U3] *dx[2]*dx[3] ;
          }
        }
}

