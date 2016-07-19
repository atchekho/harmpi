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

#define WG (WRITEGHOST!=0)

size_t write_to_dump( int is_dry_run, FILE *fp, dumptype *buf, dumptype val )
{
  if (is_dry_run) {
    return(1L);
  }
  if (fp) {
    return fwrite( &val, sizeof(dumptype), 1, fp );
  }
  if (buf) {
    *buf = val;
    return(1L);
  }
  return(0L);
}

size_t write_to_gdump( int is_dry_run, FILE *fp, gdumptype *buf, gdumptype val )
{
    if (is_dry_run) {
        return(1L);
    }
    if (fp) {
        return fwrite( &val, sizeof(gdumptype), 1, fp );
    }
    if (buf) {
        *buf = val;
        return(1L);
    }
    return(0L);
}

size_t write_to_gdump2( int is_dry_run, FILE *fp, gdump2type *buf, gdump2type val )
{
    if (is_dry_run) {
        return(1L);
    }
    if (fp) {
        return fwrite( &val, sizeof(gdump2type), 1, fp );
    }
    if (buf) {
        *buf = val;
        return(1L);
    }
    return(0L);
}

size_t write_to_rdump( int is_dry_run, FILE *fp, rdumptype *buf, rdumptype val )
{
    if (is_dry_run) {
        return(1L);
    }
    if (fp) {
        return fwrite( &val, sizeof(rdumptype), 1, fp );
    }
    if (buf) {
        *buf = val;
        return(1L);
    }
    return(0L);
}


size_t dump(int dumpno, int is_dry_run)
{
  int i,j,k,m ;
  double divb ;
  double X[NDIM] ;
  double r,th,phi,vmin,vmax ;
  struct of_geom geom ;
  struct of_state q ;
  int di, dj, dk;
  size_t n;
  char file_name[MAXLEN];
  FILE *fp;
  int is_write;
  long long offset;
  
  if (!is_dry_run){
    /***************************************************************
     Write header information
     ***************************************************************/
    sprintf(file_name,"dumps/dump%03d",dumpno) ;
#if MPI && !DO_PARALLEL_WRITE
      append_rank(file_name);
#endif
    if(i_am_the_master) fprintf(stderr,"DUMP     file=%s\n",file_name) ;
  
    if (!DO_PARALLEL_WRITE || i_am_the_master) {
        fp = fopen(file_name,"wb") ;
        
        if(fp==NULL) {
          fprintf(stderr,"error opening dump file\n") ;
          exit(2) ;
        }

        fprintf(fp, FMT_DBL_OUT, t        );
        //output per-dumpfile resolution
#if MPI && DO_PARALLEL_WRITE
        fprintf(fp, FMT_INT_OUT, mpi_ntot[1]);
        fprintf(fp, FMT_INT_OUT, mpi_ntot[2]);
        fprintf(fp, FMT_INT_OUT, mpi_ntot[3]);
#else
        fprintf(fp, FMT_INT_OUT, N1);
        fprintf(fp, FMT_INT_OUT, N2);
        fprintf(fp, FMT_INT_OUT, N3);
#endif
        //output total resolution
        fprintf(fp, FMT_INT_OUT, mpi_ntot[1]);
        fprintf(fp, FMT_INT_OUT, mpi_ntot[2]);
        fprintf(fp, FMT_INT_OUT, mpi_ntot[3]);
        //output numbers of ghost cells written to the dumps
#if MPI && DO_PARALLEL_WRITE
        fprintf(fp, FMT_INT_OUT, 0);
        fprintf(fp, FMT_INT_OUT, 0);
        fprintf(fp, FMT_INT_OUT, 0);
#else
        fprintf(fp, FMT_INT_OUT, WG*N1G);
        fprintf(fp, FMT_INT_OUT, WG*N2G);
        fprintf(fp, FMT_INT_OUT, WG*N3G);
#endif
        fprintf(fp, FMT_DBL_OUT, startx[1]);
        fprintf(fp, FMT_DBL_OUT, startx[2]);
        fprintf(fp, FMT_DBL_OUT, startx[3]);
        fprintf(fp, FMT_DBL_OUT, dx[1]    );
        fprintf(fp, FMT_DBL_OUT, dx[2]    );
        fprintf(fp, FMT_DBL_OUT, dx[3]    );
        fprintf(fp, FMT_DBL_OUT, tf       );
        fprintf(fp, FMT_INT_OUT, nstep    );
        fprintf(fp, FMT_DBL_OUT, a        );
        fprintf(fp, FMT_DBL_OUT, gam      );
        fprintf(fp, FMT_DBL_OUT, cour     );
        fprintf(fp, FMT_DBL_OUT, DTd      );
        fprintf(fp, FMT_DBL_OUT, DTl      );
        fprintf(fp, FMT_DBL_OUT, DTi      );
        fprintf(fp, FMT_DBL_OUT, DTr      );
        fprintf(fp, FMT_INT_OUT, DTr01    );
        fprintf(fp, FMT_INT_OUT, dump_cnt );
        fprintf(fp, FMT_INT_OUT, image_cnt);
        fprintf(fp, FMT_INT_OUT, rdump_cnt);
        fprintf(fp, FMT_INT_OUT, rdump01_cnt);
        fprintf(fp, FMT_DBL_OUT, dt       );
        fprintf(fp, FMT_INT_OUT, lim      );
        fprintf(fp, FMT_INT_OUT, failed   );
        fprintf(fp, FMT_DBL_OUT, Rin      );
        fprintf(fp, FMT_DBL_OUT, Rout     );
        fprintf(fp, FMT_DBL_OUT, hslope   );
        fprintf(fp, FMT_DBL_OUT, R0       );
        fprintf(fp, FMT_INT_OUT, (int)NPR );
        fprintf(fp, FMT_INT_OUT, (int)DOKTOT );
        fprintf(fp, FMT_DBL_OUT, fractheta);
        fprintf(fp, FMT_DBL_OUT, fracphi  );
        fprintf(fp, FMT_DBL_OUT, rbr      );
        fprintf(fp, FMT_DBL_OUT, npow2    );
        fprintf(fp, FMT_DBL_OUT, cpow2    );
        fprintf(fp, FMT_INT_OUT, (int)BL    );
      
        fprintf(fp,"\n") ;
#if MPI && DO_PARALLEL_WRITE
        fclose(fp) ;
        fp = NULL;
#endif
    }
    else {
      fp = NULL;
    }
    /***************************************************************
     Done writing header information
     ***************************************************************/
  }

  //needed for divb calc below
  di = (N1>1);
  dj = (N2>1);
  dk = (N3>1);
  
  n = 0L;
  ZSLOOP(0-WG*N1G,N1-1+WG*N1G,0-WG*N2G,N2-1+WG*N2G,0-WG*N3G,N3-1+WG*N3G) {
    coord(i,j,k,CENT,X) ;
    get_phys_coord(i,j,k,&r,&th,&phi) ;
    
    //output global i, j, k
    n += write_to_dump( is_dry_run, fp, dump_buffer+n, i + mpi_startn[1] );
    n += write_to_dump( is_dry_run, fp, dump_buffer+n, j + mpi_startn[2] );
    n += write_to_dump( is_dry_run, fp, dump_buffer+n, k + mpi_startn[3] );
    n += write_to_dump( is_dry_run, fp, dump_buffer+n, X[1] );
    n += write_to_dump( is_dry_run, fp, dump_buffer+n, X[2] );
    n += write_to_dump( is_dry_run, fp, dump_buffer+n, X[3] );
    n += write_to_dump( is_dry_run, fp, dump_buffer+n, r );
    n += write_to_dump( is_dry_run, fp, dump_buffer+n, th );
    n += write_to_dump( is_dry_run, fp, dump_buffer+n, phi);
    
    PLOOP n += write_to_dump(is_dry_run, fp, dump_buffer+n, p[i][j][k][m]);
      
    /* divb flux-ct defn; corner-centered.  Use
     only interior corners */
    if((N1 == 1 || mpi_periods[1] || (i + mpi_startn[1] > 0 && i + mpi_startn[1] < mpi_ntot[1]-1 && i > -N1G && i < N1+N1G-1)) &&
       (N2 == 1 || mpi_periods[2] || (j + mpi_startn[2] >= 0 && j + mpi_startn[2] <= mpi_ntot[2] && j > -N2G && j < N2+N2G-1)) &&
       (N3 == 1 || mpi_periods[3] || (k + mpi_startn[3] >= 0 && k + mpi_startn[3] <= mpi_ntot[3] && k > -N3G && k < N3+N3G-1)) ){
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
    }
    else divb = 0. ;
    
    n += write_to_dump(is_dry_run, fp, dump_buffer+n, divb     );
    
    
    if(!failed) {
      get_geometry(i,j,k, CENT,&geom) ;
      get_state(p[i][j][k],&geom,&q) ;
      
      for(m=0;m<NDIM;m++) n += write_to_dump(is_dry_run, fp, dump_buffer+n, q.ucon[m]) ;
      for(m=0;m<NDIM;m++) n += write_to_dump(is_dry_run, fp, dump_buffer+n, q.ucov[m]) ;
      for(m=0;m<NDIM;m++) n += write_to_dump(is_dry_run, fp, dump_buffer+n, q.bcon[m]) ;
      for(m=0;m<NDIM;m++) n += write_to_dump(is_dry_run, fp, dump_buffer+n, q.bcov[m]) ;
      
      vchar(p[i][j][k],&q,&geom,1,&vmax,&vmin) ;
      n += write_to_dump(is_dry_run, fp, dump_buffer+n, vmin );
      n += write_to_dump(is_dry_run, fp, dump_buffer+n, vmax );
      
      vchar(p[i][j][k],&q,&geom,2,&vmax,&vmin) ;
      n += write_to_dump(is_dry_run, fp, dump_buffer+n, vmin );
      n += write_to_dump(is_dry_run, fp, dump_buffer+n, vmax );
      
      vchar(p[i][j][k],&q,&geom,3,&vmax,&vmin) ;
      n += write_to_dump(is_dry_run, fp, dump_buffer+n, vmin );
      n += write_to_dump(is_dry_run, fp, dump_buffer+n, vmax );
      
      n += write_to_dump(is_dry_run, fp, dump_buffer+n, geom.g );
    }
    if (is_dry_run) {
      return(n);
    }
  }
#if(!DO_PARALLEL_WRITE)
  fclose(fp);
#endif
#if MPI && DO_PARALLEL_WRITE
  is_write = 1; //write, not read
  offset = -1; //append
  parallel_readwrite(file_name,dump_buffer,DUMP_FILE, is_write, offset);
  return(n);
#endif
  return(0L);
}

size_t gdump(int is_dry_run)
{
  int i,j,k,l,m ;
  double divb ;
  double X[NDIM] ;
  double r,th,phi,vmin,vmax ;
  struct of_geom geom ;
  struct of_state q ;
  double dxdxp[NDIM][NDIM];
  size_t n;
  char file_name[MAXLEN];
  FILE *fp;
  int is_write;
  long long offset;

  if (!is_dry_run) {
    /***************************************************************
     Write header information
     ***************************************************************/
    sprintf(file_name,"dumps/gdump") ;
#if MPI && !DO_PARALLEL_WRITE
      append_rank(file_name);
#endif
    if(i_am_the_master) fprintf(stderr,"GDUMP    file=%s\n",file_name) ;

    if (!DO_PARALLEL_WRITE || i_am_the_master) {

        fp = fopen(file_name,"wb") ;
        
        if(fp==NULL) {
          fprintf(stderr,"error opening grid dump file\n") ;
          exit(2) ;
        }

        fprintf(fp, FMT_DBL_OUT, t        );
        //output per-dumpfile resolution
#if MPI && DO_PARALLEL_WRITE
        fprintf(fp, FMT_INT_OUT, mpi_ntot[1]);
        fprintf(fp, FMT_INT_OUT, mpi_ntot[2]);
        fprintf(fp, FMT_INT_OUT, mpi_ntot[3]);
#else
        fprintf(fp, FMT_INT_OUT, N1);
        fprintf(fp, FMT_INT_OUT, N2);
        fprintf(fp, FMT_INT_OUT, N3);
#endif
        //output total resolution
        fprintf(fp, FMT_INT_OUT, mpi_ntot[1]);
        fprintf(fp, FMT_INT_OUT, mpi_ntot[2]);
        fprintf(fp, FMT_INT_OUT, mpi_ntot[3]);
        //output numbers of ghost cells written to the dumps
#if MPI && DO_PARALLEL_WRITE
        fprintf(fp, FMT_INT_OUT, 0);
        fprintf(fp, FMT_INT_OUT, 0);
        fprintf(fp, FMT_INT_OUT, 0);
#else
        fprintf(fp, FMT_INT_OUT, WG*N1G);
        fprintf(fp, FMT_INT_OUT, WG*N2G);
        fprintf(fp, FMT_INT_OUT, WG*N3G);
#endif
        fprintf(fp, FMT_DBL_OUT, startx[1]);
        fprintf(fp, FMT_DBL_OUT, startx[2]);
        fprintf(fp, FMT_DBL_OUT, startx[3]);
        fprintf(fp, FMT_DBL_OUT, dx[1]    );
        fprintf(fp, FMT_DBL_OUT, dx[2]    );
        fprintf(fp, FMT_DBL_OUT, dx[3]    );
        fprintf(fp, FMT_DBL_OUT, tf       );
        fprintf(fp, FMT_INT_OUT, nstep    );
        fprintf(fp, FMT_DBL_OUT, a        );
        fprintf(fp, FMT_DBL_OUT, gam      );
        fprintf(fp, FMT_DBL_OUT, cour     );
        fprintf(fp, FMT_DBL_OUT, DTd      );
        fprintf(fp, FMT_DBL_OUT, DTl      );
        fprintf(fp, FMT_DBL_OUT, DTi      );
        fprintf(fp, FMT_DBL_OUT, DTr      );
        fprintf(fp, FMT_INT_OUT, DTr01    );
        fprintf(fp, FMT_INT_OUT, dump_cnt );
        fprintf(fp, FMT_INT_OUT, image_cnt);
        fprintf(fp, FMT_INT_OUT, rdump_cnt);
        fprintf(fp, FMT_INT_OUT, rdump01_cnt);
        fprintf(fp, FMT_DBL_OUT, dt       );
        fprintf(fp, FMT_INT_OUT, lim      );
        fprintf(fp, FMT_INT_OUT, failed   );
        fprintf(fp, FMT_DBL_OUT, Rin      );
        fprintf(fp, FMT_DBL_OUT, Rout     );
        fprintf(fp, FMT_DBL_OUT, hslope   );
        fprintf(fp, FMT_DBL_OUT, R0       );
        fprintf(fp, FMT_INT_OUT, (int)NPR );
        fprintf(fp, FMT_INT_OUT, (int)DOKTOT );
        fprintf(fp, FMT_DBL_OUT, fractheta);
        fprintf(fp, FMT_DBL_OUT, fracphi  );
        fprintf(fp, FMT_DBL_OUT, rbr      );
        fprintf(fp, FMT_DBL_OUT, npow2    );
        fprintf(fp, FMT_DBL_OUT, cpow2    );
        fprintf(fp, FMT_INT_OUT, (int)BL    );
        
        fprintf(fp,"\n") ;
#if MPI && DO_PARALLEL_WRITE
        fclose(fp) ;
        fp = NULL;
#endif
     }
     else {
       fp = NULL;
     }
     /***************************************************************
     Done writing header information
     ***************************************************************/
  }
		
  n = 0L;
  ZSLOOP(0-WG*N1G,N1-1+WG*N1G,0-WG*N2G,N2-1+WG*N2G,0-WG*N3G,N3-1+WG*N3G) {
    coord(i,j,k,CENT,X) ;
    get_phys_coord(i,j,k,&r,&th,&phi) ;
    
    //output global i, j, k
    n += write_to_gdump( is_dry_run, fp, gdump_buffer+n, i + mpi_startn[1] );
    n += write_to_gdump( is_dry_run, fp, gdump_buffer+n, j + mpi_startn[2] );
    n += write_to_gdump( is_dry_run, fp, gdump_buffer+n, k + mpi_startn[3] );
    n += write_to_gdump( is_dry_run, fp, gdump_buffer+n, X[1]       );
    n += write_to_gdump( is_dry_run, fp, gdump_buffer+n, X[2]       );
    n += write_to_gdump( is_dry_run, fp, gdump_buffer+n, X[3]       );
    n += write_to_gdump( is_dry_run, fp, gdump_buffer+n, r          );
    n += write_to_gdump( is_dry_run, fp, gdump_buffer+n, th         );
    n += write_to_gdump( is_dry_run, fp, gdump_buffer+n, phi        );
    
    if(!failed) {
      get_geometry(i,j,k,CENT,&geom) ;
      
      //g_{ml}
      for(m=0;m<NDIM;m++)
        for(l=0;l<NDIM;l++)
          n += write_to_gdump( is_dry_run, fp, gdump_buffer+n, geom.gcov[m][l]) ;
      
      //g^{ml}
      for(m=0;m<NDIM;m++)
        for(l=0;l<NDIM;l++)
          n += write_to_gdump( is_dry_run, fp, gdump_buffer+n, geom.gcon[m][l]) ;
      
      //(-deg(g))**0.5
      n += write_to_gdump( is_dry_run, fp, gdump_buffer+n, geom.g );
      
      //dr^i/dx^j
      dxdxp_func(X, dxdxp);
      for(m=0;m<NDIM;m++) {
        for(l=0;l<NDIM;l++) {
          n += write_to_gdump( is_dry_run, fp, gdump_buffer+n, dxdxp[m][l]) ;
        }
      }
    }
    if (is_dry_run) {
      return(n);
    }
  }
#if(!DO_PARALLEL_WRITE)
  fclose(fp);
#endif

#if MPI && DO_PARALLEL_WRITE
  is_write = 1; //write, not read
  offset = -1; //append
  parallel_readwrite(file_name,gdump_buffer,GDUMP_FILE,is_write,offset);
  return(n);
#endif
  return(0L);
}

size_t gdump2(int is_dry_run)
{
    int i,j,k,l,m ;
    double divb ;
    double X[NDIM], Xloc[NDIM] ;
    double r,th,phi,vmin,vmax ;
    double rloc,thloc,philoc;
    int whichloc;
    struct of_geom geom ;
    struct of_state q ;
    double dxdxp[NDIM][NDIM];
    size_t n;
    char file_name[MAXLEN];
    FILE *fp;
    int is_write;
    long long offset;
    
    if (!is_dry_run) {
        /***************************************************************
         Write header information
         ***************************************************************/
        sprintf(file_name,"dumps/gdump2") ;
#if MPI && !DO_PARALLEL_WRITE
        append_rank(file_name);
#endif
        if(i_am_the_master) fprintf(stderr,"GDUMP2    file=%s\n",file_name) ;
        
        if (!DO_PARALLEL_WRITE || i_am_the_master) {
            
            fp = fopen(file_name,"wb") ;
            
            if(fp==NULL) {
                fprintf(stderr,"error opening grid dump file\n") ;
                exit(2) ;
            }
            
            fprintf(fp, FMT_DBL_OUT, t        );
            //output per-dumpfile resolution
#if MPI && DO_PARALLEL_WRITE
            fprintf(fp, FMT_INT_OUT, mpi_ntot[1]);
            fprintf(fp, FMT_INT_OUT, mpi_ntot[2]);
            fprintf(fp, FMT_INT_OUT, mpi_ntot[3]);
#else
            fprintf(fp, FMT_INT_OUT, N1);
            fprintf(fp, FMT_INT_OUT, N2);
            fprintf(fp, FMT_INT_OUT, N3);
#endif
            //output total resolution
            fprintf(fp, FMT_INT_OUT, mpi_ntot[1]);
            fprintf(fp, FMT_INT_OUT, mpi_ntot[2]);
            fprintf(fp, FMT_INT_OUT, mpi_ntot[3]);
            //output numbers of ghost cells written to the dumps
#if MPI && DO_PARALLEL_WRITE
            fprintf(fp, FMT_INT_OUT, 0);
            fprintf(fp, FMT_INT_OUT, 0);
            fprintf(fp, FMT_INT_OUT, 0);
#else
            fprintf(fp, FMT_INT_OUT, WG*N1G);
            fprintf(fp, FMT_INT_OUT, WG*N2G);
            fprintf(fp, FMT_INT_OUT, WG*N3G);
#endif
            fprintf(fp, FMT_DBL_OUT, startx[1]);
            fprintf(fp, FMT_DBL_OUT, startx[2]);
            fprintf(fp, FMT_DBL_OUT, startx[3]);
            fprintf(fp, FMT_DBL_OUT, dx[1]    );
            fprintf(fp, FMT_DBL_OUT, dx[2]    );
            fprintf(fp, FMT_DBL_OUT, dx[3]    );
            fprintf(fp, FMT_DBL_OUT, tf       );
            fprintf(fp, FMT_INT_OUT, nstep    );
            fprintf(fp, FMT_DBL_OUT, a        );
            fprintf(fp, FMT_DBL_OUT, gam      );
            fprintf(fp, FMT_DBL_OUT, cour     );
            fprintf(fp, FMT_DBL_OUT, DTd      );
            fprintf(fp, FMT_DBL_OUT, DTl      );
            fprintf(fp, FMT_DBL_OUT, DTi      );
            fprintf(fp, FMT_DBL_OUT, DTr      );
            fprintf(fp, FMT_INT_OUT, DTr01    );
            fprintf(fp, FMT_INT_OUT, dump_cnt );
            fprintf(fp, FMT_INT_OUT, image_cnt);
            fprintf(fp, FMT_INT_OUT, rdump_cnt);
            fprintf(fp, FMT_INT_OUT, rdump01_cnt);
            fprintf(fp, FMT_DBL_OUT, dt       );
            fprintf(fp, FMT_INT_OUT, lim      );
            fprintf(fp, FMT_INT_OUT, failed   );
            fprintf(fp, FMT_DBL_OUT, Rin      );
            fprintf(fp, FMT_DBL_OUT, Rout     );
            fprintf(fp, FMT_DBL_OUT, hslope   );
            fprintf(fp, FMT_DBL_OUT, R0       );
            fprintf(fp, FMT_INT_OUT, (int)NPR );
            fprintf(fp, FMT_INT_OUT, (int)DOKTOT );
            fprintf(fp, FMT_DBL_OUT, fractheta);
            fprintf(fp, FMT_DBL_OUT, fracphi  );
            fprintf(fp, FMT_DBL_OUT, rbr      );
            fprintf(fp, FMT_DBL_OUT, npow2    );
            fprintf(fp, FMT_DBL_OUT, cpow2    );
            fprintf(fp, FMT_INT_OUT, (int)BL    );
            
            fprintf(fp,"\n") ;
#if MPI && DO_PARALLEL_WRITE
            fclose(fp) ;
            fp = NULL;
#endif
        }
        else {
            fp = NULL;
        }
        /***************************************************************
         Done writing header information
         ***************************************************************/
    }
    
    n = 0L;
    ZSLOOP(0-WG*N1G,N1-1+WG*N1G,0-WG*N2G,N2-1+WG*N2G,0-WG*N3G,N3-1+WG*N3G) {
        coord(i,j,k,CENT,X) ;
        get_phys_coord(i,j,k,&r,&th,&phi) ;
        
        //output global i, j, k
        n += write_to_gdump2( is_dry_run, fp, gdump2_buffer+n, i + mpi_startn[1] );
        n += write_to_gdump2( is_dry_run, fp, gdump2_buffer+n, j + mpi_startn[2] );
        n += write_to_gdump2( is_dry_run, fp, gdump2_buffer+n, k + mpi_startn[3] );
        n += write_to_gdump2( is_dry_run, fp, gdump2_buffer+n, X[1]       );
        n += write_to_gdump2( is_dry_run, fp, gdump2_buffer+n, X[2]       );
        n += write_to_gdump2( is_dry_run, fp, gdump2_buffer+n, X[3]       );

        //coordinates of all locations within the cell
        for(whichloc=FACE1; whichloc<=EDGE3; whichloc++) {
            coord(i,j,k,whichloc,Xloc) ;
            bl_coord(Xloc,&rloc,&thloc,&philoc) ;
            n += write_to_gdump2( is_dry_run, fp, gdump2_buffer+n, rloc          );
            n += write_to_gdump2( is_dry_run, fp, gdump2_buffer+n, thloc         );
            n += write_to_gdump2( is_dry_run, fp, gdump2_buffer+n, philoc        );
        }
        
        if(!failed) {
            get_geometry(i,j,k,CENT,&geom) ;
            
            //(-deg(g))**0.5
            n += write_to_gdump2( is_dry_run, fp, gdump2_buffer+n, geom.g );
            
//          don't write this out to save space
//            //dr^i/dx^j
//            dxdxp_func(X, dxdxp);
//            for(m=0;m<NDIM;m++) {
//                for(l=0;l<NDIM;l++) {
//                    n += write_to_gdump2( is_dry_run, fp, gdump2_buffer+n, dxdxp[m][l]) ;
//                }
//            }
        }
        if (is_dry_run) {
            return(n);
        }
    }
#if(!DO_PARALLEL_WRITE)
    fclose(fp);
#endif
    
#if MPI && DO_PARALLEL_WRITE
    is_write = 1; //write, not read
    offset = -1; //append
    parallel_readwrite(file_name,gdump2_buffer,GDUMP2_FILE,is_write,offset);
    return(n);
#endif
    return(0L);
}


void fdump(int dumpno)
{
  FILE *fp ;
  int idum,i,j,k,m ;
  char file_name[100];
  size_t ntowrite, nhavewritten;
  size_t n;
  long long offset;
  int is_write;
  
  sprintf(file_name,"dumps/fdump%03d",dumpno) ;
  
#if MPI && !DO_PARALLEL_WRITE
  append_rank(file_name);
#endif
  if(i_am_the_master) fprintf(stderr,"FDUMP file=%s\n", file_name) ;
  
  if (!DO_PARALLEL_WRITE || i_am_the_master) {
    /*************************************************************
     Write the header of the restart file:
     *************************************************************/
    fp = fopen(file_name,"wb") ;
    
    if(fp==NULL) {
      fprintf(stderr,"error opening fail dump file\n") ;
      exit(2) ;
    }
    //output per-dumpfile resolution
    fprintf(fp, FMT_DBL_OUT, t        );
    //output per-dumpfile resolution
#if MPI && DO_PARALLEL_WRITE
    fprintf(fp, FMT_INT_OUT, mpi_ntot[1]);
    fprintf(fp, FMT_INT_OUT, mpi_ntot[2]);
    fprintf(fp, FMT_INT_OUT, mpi_ntot[3]);
#else
    fprintf(fp, FMT_INT_OUT, N1);
    fprintf(fp, FMT_INT_OUT, N2);
    fprintf(fp, FMT_INT_OUT, N3);
#endif
    //output total resolution
    fprintf(fp, FMT_INT_OUT, mpi_ntot[1]);
    fprintf(fp, FMT_INT_OUT, mpi_ntot[2]);
    fprintf(fp, FMT_INT_OUT, mpi_ntot[3]);
    //output numbers of ghost cells written to the dumps
    fprintf(fp, FMT_INT_OUT, 0);
    fprintf(fp, FMT_INT_OUT, 0);
    fprintf(fp, FMT_INT_OUT, 0);
    fprintf(fp, FMT_DBL_OUT, startx[1]);
    fprintf(fp, FMT_DBL_OUT, startx[2]);
    fprintf(fp, FMT_DBL_OUT, startx[3]);
    fprintf(fp, FMT_DBL_OUT, dx[1]    );
    fprintf(fp, FMT_DBL_OUT, dx[2]    );
    fprintf(fp, FMT_DBL_OUT, dx[3]    );
    fprintf(fp, FMT_DBL_OUT, tf       );
    fprintf(fp, FMT_INT_OUT, nstep    );
    fprintf(fp, FMT_DBL_OUT, a        );
    fprintf(fp, FMT_DBL_OUT, gam      );
    fprintf(fp, FMT_DBL_OUT, cour     );
    fprintf(fp, FMT_DBL_OUT, DTd      );
    fprintf(fp, FMT_DBL_OUT, DTl      );
    fprintf(fp, FMT_DBL_OUT, DTi      );
    fprintf(fp, FMT_DBL_OUT, DTr      );
    fprintf(fp, FMT_INT_OUT, DTr01    );
    fprintf(fp, FMT_INT_OUT, dump_cnt );
    fprintf(fp, FMT_INT_OUT, image_cnt);
    fprintf(fp, FMT_INT_OUT, rdump_cnt);
    fprintf(fp, FMT_INT_OUT, rdump01_cnt);
    fprintf(fp, FMT_DBL_OUT, dt       );
    fprintf(fp, FMT_INT_OUT, lim      );
    fprintf(fp, FMT_INT_OUT, failed   );
    fprintf(fp, FMT_DBL_OUT, Rin      );
    fprintf(fp, FMT_DBL_OUT, Rout     );
    fprintf(fp, FMT_DBL_OUT, hslope   );
    fprintf(fp, FMT_DBL_OUT, R0       );
    fprintf(fp, FMT_INT_OUT, (int)NPR );
    fprintf(fp, FMT_INT_OUT, (int)DOKTOT );
    fprintf(fp, FMT_DBL_OUT, fractheta);
    fprintf(fp, FMT_DBL_OUT, fracphi  );
    fprintf(fp, FMT_DBL_OUT, rbr      );
    fprintf(fp, FMT_DBL_OUT, npow2    );
    fprintf(fp, FMT_DBL_OUT, cpow2    );
    fprintf(fp, FMT_INT_OUT, (int)BL    );
      
    fprintf(fp,"\n") ;
#if MPI && DO_PARALLEL_WRITE
    fclose(fp) ;
    fp = NULL;
#endif
  }
  else {
    fp = NULL;
  }

#if MPI && DO_PARALLEL_WRITE
  is_write = 1; //write, not read
  offset = -1;  //append
//  n = 0L;
//  ZLOOP for(m=0;m<NIMG;m++) {
//    fdump_buffer[n] = failimage[i][j][k][m];
//    n++;
//  }
  parallel_readwrite(file_name,failimage,FDUMP_FILE,is_write,offset);
#else
  /*************************************************************
   Write the body of the restart file:
   *************************************************************/
  ntowrite = 1L*N1*N2*N3*NIMG;
  nhavewritten = fwrite( failimage, sizeof(fdumptype), ntowrite, fp );
  
  if( nhavewritten != ntowrite ) {
    fprintf(stderr, "MPI rank %d: Could not write to fdump file: nhavewritten = %ld, ntowrite = %ld, wrong format or corrupted file?\n", mpi_rank, (long)nhavewritten, (long)ntowrite);
    fflush(stderr);
#if MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif
    exit(1);
  }
  
  if( ferror(fp) ) {
    fprintf(stderr, "MPI rank %d: Error writing to fdump file\n", mpi_rank);
    fflush(stderr);
#if MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif
    exit(1);
  }
  
  //added this to ensure the entirety of the file is pushed out to the disk
  //to avoid a corrupted file upon code kill
  fflush(fp) ;
  
  fclose(fp) ;
#endif
  
  /* Reset array after every image dump: */
  ZLOOP for( m = 0 ; m < NIMG; m++ ) {
    failimage[i][j][k][m] = 0L;
  }

  return;
}


