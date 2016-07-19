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


#include <complex.h>

////////////////////////////////////////////////////////
//
//  MPI section
//
////////////////////////////////////////////////////////

#ifdef MPI
///////////////////////////////////////////////
// how to write dumps and gdumps in parallel
void *mpi_file_buffer = NULL;
dumptype *dump_buffer = NULL;
gdumptype *gdump_buffer = NULL;
gdump2type *gdump2_buffer = NULL;
rdumptype *rdump_buffer = NULL;
fdumptype *fdump_buffer = NULL;
//
////////////////////////////////////////////////
int mpi_numtasks;
MPI_Comm mpi_cartcomm;
int mpi_reorder;
int mpi_nbrs[NDIM][2];
//primitives
MPI_Request mpi_reqs_send[NDIM][2];
MPI_Request mpi_reqs_recv[NDIM][2];
MPI_Status mpi_stat_send[NDIM][2];
MPI_Status mpi_stat_recv[NDIM][2];
double mpi_buf_send[NDIM][2][(NMAX+2*NG)*(NMAX+2*NG)*NG*NPR];
double mpi_buf_recv[NDIM][2][(NMAX+2*NG)*(NMAX+2*NG)*NG*NPR];
//pflag
MPI_Request mpi_reqs_send_pflag[NDIM][2];
MPI_Request mpi_reqs_recv_pflag[NDIM][2];
MPI_Status mpi_stat_send_pflag[NDIM][2];
MPI_Status mpi_stat_recv_pflag[NDIM][2];
int mpi_buf_send_pflag[NDIM][2][(NMAX+2*NG)*(NMAX+2*NG)*NG*NPR];
int mpi_buf_recv_pflag[NDIM][2][(NMAX+2*NG)*(NMAX+2*NG)*NG*NPR];
MPI_Datatype gdump_file_type, gdump_cell_type;
MPI_Datatype gdump2_file_type, gdump2_cell_type;
MPI_Datatype dump_file_type, dump_cell_type;
MPI_Datatype rdump_file_type, rdump_cell_type;
MPI_Datatype fdump_file_type, fdump_cell_type;
#endif
int mpi_rank;
int mpi_dims[NDIM];
int mpi_coords[NDIM];
int mpi_ntot[NDIM];
int mpi_ntile[NDIM];
int mpi_startn[NDIM];
int mpi_periods[NDIM];
int i_am_the_master;

/*************************************************************************
 GLOBAL ARRAYS SECTION
 *************************************************************************/


double  a_p[N1M][N2M][N3M][NPR] ;	/* space for primitive vars */
double  a_dq[N1M][N2M][N3M][NPR] ;  /* slopes */
double  a_dsource[N1M][N2M][N3M] ;  /* sources */
double  a_duscon[N1M][N2M][N3M] ;  /* sources */
double  a_sour[N1M][N2M][N3M] ;  /* sources */
double  a_F1[N1M][N2M][N3M][NPR] ;	/* fluxes */
double  a_F2[N1M][N2M][N3M][NPR] ;	/* fluxes */
double  a_F3[N1M][N2M][N3M][NPR] ;	/* fluxes */

double  a_ph[N1M][N2M][N3M][NPR] ;	/* half-step primitives */
double  a_psave[N1M][N2M][N3M][NPR] ;	/* initial primitives */
double  a_pbound[N1M][N2M][N3M][NPR] ;	/* boundary conditions */


int     a_pflag[N1M][N2M][N3M];	/* identifies failure points */


/* for debug */
//double fimage[NIMG][N1*N2*N3];
long long    failimage[N1][N2][N3][NIMG];


/* grid functions */
double a_conn[N1M][N2M][N3M][NDIM][NDIM][NDIM] ;
double a_gcon[N1M][N2M][N3M][NPG][NDIM][NDIM] ;
double a_gcov[N1M][N2M][N3M][NPG][NDIM][NDIM] ;
double a_gdet[N1M][N2M][N3M][NPG] ;
double a_phys_coords[NDIM][N1M][N2M][N3M];
double emf[NDIM][N1+D1][N2+D2][N3+D3] ; //OPTMARK: could reduce NDIM to 1 in 2D and eliminate completely in 1D

double (*   p)[N2M][N3M][NPR] ;
double (*   psave)[N2M][N3M][NPR];
double (*   pbound)[N2M][N3M][NPR];
double (*  dq)[N2M][N3M][NPR] ;
double (*  dsource)[N2M][N3M] ;
double (*  duscon)[N2M][N3M] ;
double (*  sour)[N2M][N3M] ;
double (*  F1)[N2M][N3M][NPR] ;
double (*  F2)[N2M][N3M][NPR] ;
double (*  F3)[N2M][N3M][NPR] ;
double (*  ph)[N2M][N3M][NPR] ;
int    (*  pflag)[N2M][N3M];
double (* conn)[N2M][N3M][NDIM][NDIM][NDIM] ;
double (* gcon)[N2M][N3M][NPG][NDIM][NDIM] ;
double (* gcov)[N2M][N3M][NPG][NDIM][NDIM] ;
double (* gdet)[N2M][N3M][NPG] ;
double (* phys_coords)[N1M][N2M][N3M];








#if(DO_FONT_FIX)
double Katm[N1];
#endif

/*************************************************************************
    GLOBAL VARIABLES SECTION 
*************************************************************************/
/* physics parameters */
double a ;
double gam ;
double game;
double game4 = 4./3.;
double game5 = 5./3.;
double mrat = 1836.152672;
double kelmin;
double qosc = 1.e-7;
double fel0;
double felfloor;



/* numerical parameters */
double Rin,Rout,hslope,R0,fractheta,fracphi ;
double x1br, cpow2,npow2, rbr;
double cour ;
double dV,dx[NDIM],startx[NDIM] ;
double dt ;
double t,tf ;
double rcurr,hcurr ;
int istart,istop,jstart,jstop ;
int icurr,jcurr,kcurr,pcurr,ihere,jhere,phere ;
double dminarg1,dminarg2 ;
int nstep ;
double fval1,fval2;
int threadid = 0; //openmp id
int nthreads = 1; //openmp number of threads
//#pragma omp threadprivate(threadid)
//#pragma omp threadshared(nthreads)

/* output parameters */
double DTd ;
double DTl ;
double DTi ;
double DTr ;
int    DTr01 ;
int    dump_cnt ;
int    image_cnt ;
int    rdump_cnt ;
int    rdump01_cnt ;
int    nstroke ;

/* global flags */
int    failed ;
int    lim ;
double defcon ;

/* diagnostics */
double mdot = 0. ;
double edot = 0. ;
double ldot = 0. ;

/* current local position */
//int icurr,jcurr,kcurr,pcurr ; //!!!ATCH: duplicated above, commented this out

double tdump,trdump,timage,tlog ;

////////////////////////////////
//SJETCOORDS
////////////////////////////////
double global_fracdisk;
double global_fracjet;
double global_jetnu;
double global_rsjet;
double global_r0grid;
double global_r0jet;
double global_rjetend;
double global_r0disk;
double global_rdiskend;
double global_x10;
double global_x20;

//TORUS
double global_kappa;
double aphipow;
