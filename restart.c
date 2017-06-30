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


/* restart functions; restart_init and restart_dump */

#include "decs.h"



/***********************************************************************/
/***********************************************************************
  restart_write():
     -- writes current state of primitive variables to the 
        checkpointing/restart file. 
     -- uses ASCII text format ;
     -- when changing this routine, be sure to make analogous changes 
        in restart_read();
************************************************************************/
void restart_write(int dumpno)
{
  FILE *fp ;
  int idum,i,j,k,m ;
  char file_name[100];
  size_t ntowrite, nhavewritten;
  size_t n;
  long long offset;
  int is_write;

  if (dumpno < 0) {
      sprintf(file_name,"dumps/rdump%d",rdump01_cnt%2) ;
      //increment the counter so that upon restart get the incremented value
      rdump01_cnt++ ;
  }
  else {
      sprintf(file_name,"dumps/rdump%03d",dumpno) ;
  }
#if MPI && !DO_PARALLEL_WRITE
      append_rank(file_name);
#endif
  if(i_am_the_master) fprintf(stderr,"RESTART file=%s\n", file_name) ;
   
  if (!DO_PARALLEL_WRITE || i_am_the_master) {
      /*************************************************************
       Write the header of the restart file:
       *************************************************************/
      fp = fopen(file_name,"wb") ;
      
      if(fp==NULL) {
          fprintf(stderr,"error opening restart dump file\n") ;
          exit(2) ;
      }
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
      fprintf(fp, FMT_INT_OUT, mpi_ntot[1]);
      fprintf(fp, FMT_INT_OUT, mpi_ntot[2]);
      fprintf(fp, FMT_INT_OUT, mpi_ntot[3]);
      //output numbers of ghost cells written to the dumps
#if MPI && DO_PARALLEL_WRITE
      fprintf(fp, FMT_INT_OUT, 0);
      fprintf(fp, FMT_INT_OUT, 0);
      fprintf(fp, FMT_INT_OUT, 0);
#else
      fprintf(fp, FMT_INT_OUT, N1G);
      fprintf(fp, FMT_INT_OUT, N2G);
      fprintf(fp, FMT_INT_OUT, N3G);
#endif
      fprintf(fp, FMT_INT_OUT, mpi_startn[1]);
      fprintf(fp, FMT_INT_OUT, mpi_startn[2]);
      fprintf(fp, FMT_INT_OUT, mpi_startn[3]);
      fprintf(fp, FMT_DBL_OUT, t        );
      fprintf(fp, FMT_DBL_OUT, tf       );
      fprintf(fp, FMT_INT_OUT, nstep    );
      fprintf(fp, FMT_DBL_OUT, a        );
      fprintf(fp, FMT_DBL_OUT, gam      );
      fprintf(fp, FMT_DBL_OUT, game     );
      fprintf(fp, FMT_DBL_OUT, game4    );
      fprintf(fp, FMT_DBL_OUT, game5    );
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
      fprintf(fp, FMT_DBL_OUT, fractheta);
      fprintf(fp, FMT_DBL_OUT, fracphi  );
      fprintf(fp, FMT_DBL_OUT, rbr      );
      fprintf(fp, FMT_DBL_OUT, npow2    );
      fprintf(fp, FMT_DBL_OUT, cpow2    );
      fprintf(fp, FMT_DBL_OUT, global_x10);
      fprintf(fp, FMT_DBL_OUT, global_x20);
      fprintf(fp, FMT_DBL_OUT, mrat      );
      fprintf(fp, FMT_DBL_OUT, fel0      );
      fprintf(fp, FMT_DBL_OUT, felfloor  );
      fprintf(fp, FMT_DBL_OUT, tdump     );
      fprintf(fp, FMT_DBL_OUT, trdump    );
      fprintf(fp, FMT_DBL_OUT, timage    );
      fprintf(fp, FMT_DBL_OUT, tlog      );
      fprintf(fp, FMT_DBL_OUT, global_fracdisk  );
      fprintf(fp, FMT_DBL_OUT, global_fracjet   );
      fprintf(fp, FMT_DBL_OUT, global_r0disk    );
      fprintf(fp, FMT_DBL_OUT, global_rdiskend  );
      fprintf(fp, FMT_DBL_OUT, global_r0jet     );
      fprintf(fp, FMT_DBL_OUT, global_rjetend   );
      fprintf(fp, FMT_DBL_OUT, global_jetnu     );
      fprintf(fp, FMT_DBL_OUT, global_rsjet     );
      fprintf(fp, FMT_DBL_OUT, global_r0grid    );

      fprintf(fp," \n");
#if MPI && DO_PARALLEL_WRITE
          fclose(fp) ;
          fp = NULL;
#endif
  }
  else {
    fp = NULL;
  }

  n = 0L;
#if MPI && DO_PARALLEL_WRITE
  ZLOOP {
    PLOOP {
      n+= write_to_rdump(0, NULL, rdump_buffer+n, p[i][j][k][m]);
    }
  }
  is_write = 1; //write, not read
  offset = -1;  //append
  parallel_readwrite(file_name,rdump_buffer,RDUMP_FILE,is_write,offset);
#else
  /*************************************************************
   Write the body of the restart file:
   *************************************************************/
  ntowrite = 1L*N1M*N2M*N3M*NPR;
  nhavewritten = fwrite( a_p, sizeof(double), ntowrite, fp );

  if( nhavewritten != ntowrite ) {
    fprintf(stderr, "MPI rank %d: Could not write to restart file: nhavewritten = %ld, ntowrite = %ld, wrong format or corrupted file?\n", mpi_rank, (long)nhavewritten, (long)ntowrite);
    fflush(stderr);
#if MPI
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif
    exit(1);
  }

  
  if( ferror(fp) ) {
      fprintf(stderr, "MPI rank %d: Error writing to restart dump file\n", mpi_rank);
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
  return;
}

/***********************************************************************/
/***********************************************************************
  restart_init():
     -- main driver for setting initial conditions from a checkpoint 
        or restart file. 
     -- determines if there are any restart files to use and then 
         lets the  user choose if there are more than one file. 
     -- then calls initializes run with restart data;
************************************************************************/
int restart_init()
{
  FILE *fp, *fp1, *fp0 ;
  char ans[100], suff[100],fname0[100],fname1[100] ;
  //int strncmp(char *s1, char *s2, int n) ;
  int i,j,k,m ;
  int dumpno;
  int restart_success;

  /* set up global arrays */
  set_arrays() ;

  /********************************************************************
   Attempt to restart from a checkpoint file.  If successful, then
   we need to read in data, assign grid functions and define the grid: 
  ********************************************************************/
  //first attempt to restart from rdump0
  dumpno = 0;
  restart_success = restart_read(dumpno) ;
  if(restart_success) {
    if(i_am_the_master) {
      fprintf(stderr, "Successfully restarted from rdump0\n");
    }
  }
  else {
    dumpno = 1;
    restart_success = restart_read(dumpno) ;
    if(restart_success) {
      if(i_am_the_master) {
        fprintf(stderr, "Successfully restarted from rdump1\n");
      }
    }
    else {
      if(i_am_the_master) {
        fprintf(stderr,"No restart file\n") ;
      }
      return(0);
    }
  }

  /* set half-step primitives everywhere */
  ZSLOOP(-N1G,N1+N1G-1,-N2G,N2+N2G-1,-N3G,N3+N3G-1) PLOOP ph[i][j][k][m] = p[i][j][k][m] ;

  /* set metric functions */
  set_grid() ;

  /* bound */
  bound_prim(p) ;
  bound_prim(ph) ;


#if( DO_FONT_FIX ) 
  set_Katm();
#endif 

  /***********************************************************************
    Make any changes to parameters in restart file  here: 
      e.g., cour = 0.4 , change in limiter...
  ************************************************************************/
  //lim = MC ;
  //cour = 0.9 ;
  //lim = VANL ;
  //tf = 4000. ;
  tf = 20000.;

  if(i_am_the_master) fprintf(stderr,"done with restart init.\n") ;

  /* done! */
  return(1) ;

}

/***********************************************************************/
/***********************************************************************
  restart_read():
     -- reads in data from the restart file, which is specified in 
         restart_init() but is usually named "dumps/rdump[0,1]" 
     -- returns 1 on success, 0 on failure
************************************************************************/
int restart_read(int dumpno)
{
  int fscanf_and_bcast(FILE *fp, char *format, ...);

  FILE *fp ;
  int idum,i,j,k,m ;
  char file_name[100];
  size_t ntoread, nhaveread;
  int test;
  char s[MAXLEN];
  size_t n;
  long offset;
  int is_write;
  int file_exists;

  sprintf(file_name,"dumps/rdump%d",dumpno) ;
#if MPI && !DO_PARALLEL_WRITE
    append_rank(file_name);
#endif

  /*************************************************************
   Write the header of the restart file:
   *************************************************************/
  if(i_am_the_master || !DO_PARALLEL_WRITE) {
    fp = fopen(file_name,"rb") ;
    file_exists = (fp!=NULL);
  }
  else {
    fp = NULL;
  }

#if MPI
#if DO_PARALLEL_WRITE
  MPI_Bcast(&file_exists,1,MPI_INT,MASTER,MPI_COMM_WORLD);
#else
  MPI_Allreduce(MPI_IN_PLACE,&file_exists,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
#endif
#endif //MPI
  
  if(!file_exists) {
    if(i_am_the_master)
      fprintf(stderr,"Error opening %s file, does not exist?\n", file_name) ;
    return(0) ;
  }
  
  /*************************************************************
          READ the header of the restart file: 
  *************************************************************/
  fscanf_and_bcast(fp, "%d", &idum  );
  if(!DO_PARALLEL_WRITE && idum != N1) {
    if(i_am_the_master) {
      fprintf(stderr,"error reading restart file; N1 differs: read: %d, N1 = %d\n",
              idum, (int)N1) ;
      fflush(stderr);
    }
    exit(3) ;
  }
  else if(DO_PARALLEL_WRITE && idum != mpi_ntot[1]) {
    if(i_am_the_master) {
      fprintf(stderr,"error reading restart file; N1 differs: read: %d, ntot1 = %d\n",
              idum, mpi_ntot[1]) ;
      fflush(stderr);
    }
    exit(3) ;
  }
  fscanf_and_bcast(fp, "%d", &idum  );
  if(!DO_PARALLEL_WRITE && idum != N2) {
    if(i_am_the_master) {
      fprintf(stderr,"error reading restart file; N2 differs\n") ;
      fflush(stderr);
    }
    exit(4) ;
  }
  else if(DO_PARALLEL_WRITE && idum != mpi_ntot[2]) {
    if(i_am_the_master) {
      fprintf(stderr,"error reading restart file; N2 differs: read: %d, ntot2 = %d\n",
              idum, mpi_ntot[2]) ;
      fflush(stderr);
    }
    exit(4) ;
  }

  fscanf_and_bcast(fp, "%d", &idum  );
  if(!DO_PARALLEL_WRITE && idum != N3) {
    if(i_am_the_master) {
      fprintf(stderr,"error reading restart file; N3 differs\n") ;
      fflush(stderr);
    }
    exit(4) ;
  }
  else if(DO_PARALLEL_WRITE && idum != mpi_ntot[3]) {
    if(i_am_the_master) {
      fprintf(stderr,"error reading restart file; N3 differs: read: %d, ntot3 = %d\n",
              idum, mpi_ntot[3]) ;
      fflush(stderr);
    }
    exit(4) ;
  }

  fscanf_and_bcast(fp, "%d", &idum  );
  if(idum != mpi_ntot[1]) {
    if(i_am_the_master) {
      fprintf(stderr,"error reading restart file; ntot1 differs: read %d, ntot1 = %d\n",
              idum, mpi_ntot[1]) ;
      fflush(stderr);
    }
    exit(3) ;
  }
  fscanf_and_bcast(fp, "%d", &idum  );
  if(idum != mpi_ntot[2]) {
    if(i_am_the_master) {
      fprintf(stderr,"error reading restart file; ntot2 differs\n") ;
      fflush(stderr);
    }
    exit(4) ;
  }
  
  fscanf_and_bcast(fp, "%d", &idum  );
  if(idum != mpi_ntot[3]) {
    if(i_am_the_master) {
      fprintf(stderr,"error reading restart file; ntot3 differs\n") ;
      fflush(stderr);
    }
    exit(4) ;
  }
  fscanf_and_bcast(fp, "%d", &idum  );
  if(!DO_PARALLEL_WRITE && idum != N1G) {
    if(i_am_the_master) {
      fprintf(stderr,"error reading restart file; N1G differs\n") ;
      fflush(stderr);
    }
    exit(3) ;
  }
  else if(DO_PARALLEL_WRITE && idum != 0) {
    if(i_am_the_master) {
      fprintf(stderr,"error reading restart file; N1G differs: read: %d, N1G = %d\n",
              idum, (int)0) ;
      fflush(stderr);
    }
    exit(4) ;
  }
  fscanf_and_bcast(fp, "%d", &idum  );
  if(!DO_PARALLEL_WRITE && idum != N2G) {
    if(i_am_the_master) {
      fprintf(stderr,"error reading restart file; N2G differs\n") ;
      fflush(stderr);
    }
    exit(4) ;
  }
  else if(DO_PARALLEL_WRITE && idum != 0) {
    if(i_am_the_master) {
      fprintf(stderr,"error reading restart file; N2G differs: read: %d, N2G = %d\n",
              idum, (int)0) ;
      fflush(stderr);
    }
    exit(4) ;
  }
  
  fscanf_and_bcast(fp, "%d", &idum  );
  if(!DO_PARALLEL_WRITE && idum != N3G) {
    if(i_am_the_master) {
      fprintf(stderr,"error reading restart file; N3G differs\n") ;
      fflush(stderr);
    }
    exit(4) ;
  }
  else if(DO_PARALLEL_WRITE && idum != 0) {
    if(i_am_the_master) {
      fprintf(stderr,"error reading restart file; N3G differs: read: %d, N3G = %d\n",
              idum, (int)0) ;
      fflush(stderr);
    }
    exit(4) ;
  }

  fscanf_and_bcast(fp, "%d", &idum  );
  if((i_am_the_master || !DO_PARALLEL_WRITE) && idum != mpi_startn[1]) {
    fprintf(stderr,"error reading restart file; mpi_startn[1] differs\n") ;
    fflush(stderr);
    exit(3) ;
  }
  fscanf_and_bcast(fp, "%d", &idum  );
  if((i_am_the_master || !DO_PARALLEL_WRITE) && idum != mpi_startn[2]) {
    fprintf(stderr,"error reading restart file; mpi_startn[2] differs\n") ;
    fflush(stderr);
    exit(4) ;
  }
  
  fscanf_and_bcast(fp, "%d", &idum  );
  if((i_am_the_master || !DO_PARALLEL_WRITE) && idum != mpi_startn[3]) {
    fprintf(stderr,"error reading restart file; mpi_startn[3] differs\n") ;
    fflush(stderr);
    exit(4) ;
  }
  
  fscanf_and_bcast(fp, "%lf", &t        );
  fscanf_and_bcast(fp, "%lf", &tf       );
  fscanf_and_bcast(fp, "%d",  &nstep    );
  fscanf_and_bcast(fp, "%lf", &a        );
  fscanf_and_bcast(fp, "%lf", &gam      );
  fscanf_and_bcast(fp, "%lf", &game     );
  fscanf_and_bcast(fp, "%lf", &game4    );
  fscanf_and_bcast(fp, "%lf", &game5    );
  fscanf_and_bcast(fp, "%lf", &cour     );
  fscanf_and_bcast(fp, "%lf", &DTd      );
  fscanf_and_bcast(fp, "%lf", &DTl      );
  fscanf_and_bcast(fp, "%lf", &DTi      );
  fscanf_and_bcast(fp, "%lf", &DTr      );
  fscanf_and_bcast(fp, "%d",  &DTr01    );
  fscanf_and_bcast(fp, "%d",  &dump_cnt );
  fscanf_and_bcast(fp, "%d",  &image_cnt);
  fscanf_and_bcast(fp, "%d",  &rdump_cnt);
  fscanf_and_bcast(fp, "%d",  &rdump01_cnt);
  fscanf_and_bcast(fp, "%lf", &dt       );
  fscanf_and_bcast(fp, "%d",  &lim      );
  fscanf_and_bcast(fp, "%d",  &failed   );
  fscanf_and_bcast(fp, "%lf", &Rin      );
  fscanf_and_bcast(fp, "%lf", &Rout     );
  fscanf_and_bcast(fp, "%lf", &hslope   );
  fscanf_and_bcast(fp, "%lf", &R0       );
  fscanf_and_bcast(fp, "%lf", &fractheta);
  fscanf_and_bcast(fp, "%lf", &fracphi  );
  fscanf_and_bcast(fp, "%lf", &rbr      );
  fscanf_and_bcast(fp, "%lf", &npow2    );
  fscanf_and_bcast(fp, "%lf", &cpow2    );
  fscanf_and_bcast(fp, "%lf", &global_x10    );
  fscanf_and_bcast(fp, "%lf", &global_x20    );
  fscanf_and_bcast(fp, "%lf", &tdump    );
  fscanf_and_bcast(fp, "%lf", &trdump   );
  fscanf_and_bcast(fp, "%lf", &timage   );
  fscanf_and_bcast(fp, "%lf", &tlog     );
  fscanf_and_bcast(fp, "%lf", &global_fracdisk   );
  fscanf_and_bcast(fp, "%lf", &global_fracjet    );
  fscanf_and_bcast(fp, "%lf", &global_r0disk     );
  fscanf_and_bcast(fp, "%lf", &global_rdiskend   );
  fscanf_and_bcast(fp, "%lf", &global_r0jet      );
  fscanf_and_bcast(fp, "%lf", &global_rjetend    );
  fscanf_and_bcast(fp, "%lf", &global_jetnu      );
  fscanf_and_bcast(fp, "%lf", &global_rsjet      );
  fscanf_and_bcast(fp, "%lf", &global_r0grid      );
  
  //skip to the start of the data by
  //reading the rest of the line including the newline ('\n') character
  if(fp) fgets(s, MAXLEN, fp);

#if MPI && DO_PARALLEL_WRITE
  if (i_am_the_master) {
    //record the location of the first byte beyond the header
    offset = ftell(fp);
    fclose(fp) ;
    fp = NULL;
  }
#endif
  /*************************************************************
   READ the body of the restart file:
   *************************************************************/
  
#if MPI && DO_PARALLEL_WRITE 
  {
    MPI_Bcast(&offset,1,MPI_LONG,MASTER,MPI_COMM_WORLD);
    is_write = 0; //read, not write
    parallel_readwrite(file_name,rdump_buffer,RDUMP_FILE,is_write,offset);
    n = 0L;
    ZLOOP {
      PLOOP {
        p[i][j][k][m]=rdump_buffer[n];
        n++;
      }
                
    }
  }
#else 
  {
    ntoread = 1L*N1M*N2M*N3M*NPR;
    nhaveread = fread( a_p, sizeof(double), ntoread, fp );
    
    if( nhaveread != ntoread || feof(fp) ) {
      fprintf(stderr, "MPI rank %d: Restart file ended prematurely: nhaveread = %ld, ntoread = %ld, wrong format or corrupted file?\n", mpi_rank, (long)nhaveread, (long)ntoread);
      fflush(stderr);
#if MPI
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
#endif
      exit(1);
    }

    
    
    if( ferror(fp) ) {
      fprintf(stderr, "MPI rank %d: Error reading from restart dump file\n", mpi_rank);
      fflush(stderr);
#if MPI
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
#endif
      exit(1);
    }
    
    test = fgetc( fp );
    
    if( EOF != test || !feof(fp) ) {
      fprintf(stderr, "MPI rank %d: Restart has more data than expected, wrong format or corrupted file?\n", mpi_rank);
      fflush(stderr);
#if MPI
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,1);
#endif
      exit(1);
    }
  }
#endif

  return(1) ;
}

//Reads the in the values from the string and
//broadcasts the values from master to the rest
//NOTES: fp should be non-Null in the master thread
//       fp is ignored in other threads
int fscanf_and_bcast(FILE *fp, char *format, ...)
{
  va_list args, args_copy;
  va_start (args, format);
  va_copy (args_copy, args);
  double *pdouble;
  int *pint;
  int ret;
  if (i_am_the_master) {
    if(fp) {
      ret = vfscanf(fp, format, args_copy);
    }
    else {
      fprintf(stderr, "fscanf_and_bcast(): fp = NULL, yet I am the master\n");
      fflush(stderr);
      ret = 0;
    }
  }
#if MPI
  if (!strcmp(format,"%lf")) {
      pdouble = va_arg(args,double*);
      MPI_Bcast(pdouble,1,MPI_DOUBLE,MASTER,MPI_COMM_WORLD);
      MPI_Bcast(&ret,1,MPI_INT,MASTER,MPI_COMM_WORLD);
  }
  else if (!strcmp(format,"%d")) {
      pint = va_arg(args,int*);
      MPI_Bcast(pint,1,MPI_INT,MASTER,MPI_COMM_WORLD);
      MPI_Bcast(&ret,1,MPI_INT,MASTER,MPI_COMM_WORLD);
  }
  else {
    fprintf(stderr,"Unknown format specifier, %s\n", format);
    fflush(stderr);
    ret = 0;
  }
#else
  fprintf(stderr,"Cannot communicate since MPI is not defined\n");
  fflush(stderr);
  ret = 0;
#endif
  va_end (args);
  va_end(args_copy);
  return(ret);
}

#undef FMT_DBL_OUT
#undef FMT_INT_OUT
