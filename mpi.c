//
//  mpi.c
//  HARM2D
//
//  Created by Alexander Tchekhovskoy on 2/22/15.
//  Copyright (c) 2015. All rights reserved.
//

#include "mpi.h"
#include "decs.h"
#include <unistd.h>
#include <errno.h>

static size_t dump_buffer_size, gdump_buffer_size, gdump2_buffer_size, rdump_buffer_size, fdump_buffer_size;

//initializes MPI communicator and sets up MPI-related bookkeeping
//if compiled without an MPI library, then sets the code up for serial run
int mpi_init(int argc,char *argv[])
{
  int dim;
#ifdef MPI
  int ntiles=1;
  int error = 0;
#endif
  
  ///////////////////////////
  //
  //both MPI and serial
  //
  ///////////////////////////
  
  //are we using periodic boundary conditions?
  mpi_periods[1] = 0;
  mpi_periods[2] = 0;
  mpi_periods[3] = 0;
#if PERIODIC && !BL
  mpi_periods[1] = 1;
  mpi_periods[2] = 1;
  mpi_periods[3] = 1;
#endif
#if BL
  mpi_periods[3] = 1;
#endif
  

#ifdef MPI
  ntiles=1;
  error = 0;

  //some MPI settings
  mpi_reorder = 0; //not supported anyway
  
  //check if the supplied arguments contain the numbers of tiles in each direction
  if (argc < MPI_NDIM+1) {
    printf("Usage: %s ntiles1 ntiles2 ntiles3\n", argv[0]);
    exit(1);
  }

  //initialize MPI and figure out the number of tasks
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_numtasks);

  //convert the arguments to the executable into dimensions
  //SASMARK: no checking for errors
  for( dim = 1; dim < NDIM; dim++){
    mpi_dims[dim] = atoi(argv[dim]);
    ntiles *= mpi_dims[dim];
  }

  mpi_ntile[1] = N1;
  mpi_ntile[2] = N2;
  mpi_ntile[3] = N3;
  
  for(dim=1;dim<NDIM;dim++) {
    mpi_ntot[dim]=mpi_ntile[dim]*mpi_dims[dim];
    if(mpi_dims[dim]>1 && mpi_ntile[dim]==1){
      MPI_Finalize();
      fprintf(stderr,"Cannot have N%d = 1 for nonunity Ntile%d = %d\n", dim, dim, mpi_ntile[dim]);
      exit(1);
    }
  }
  
  if (ntiles != mpi_numtasks) {
    MPI_Finalize();
    //the following gets printed by each MPI task, would be good to avoid this
    printf("The number of tiles, %d, does not match the number of tasks, %d\n",
           ntiles, mpi_numtasks);
    fflush(stdout);
    exit(1);
  }

  //create MPI communicator
  //NOTE: adding +1 to mpi_ndims and mpi_periods so that dir goes from 1 through 3
  MPI_Cart_create(MPI_COMM_WORLD, MPI_NDIM, mpi_dims+1, mpi_periods+1, mpi_reorder, &mpi_cartcomm);
  //figure out the rank of the current MPI process
  MPI_Comm_rank(mpi_cartcomm, &mpi_rank);
  i_am_the_master = (MASTER == mpi_rank);
  //convert the rank to the Cartesian coordinates of the current MPI process
  //NOTE: adding +1 to mpi_coords so that dir goes from 1 through 3
  MPI_Cart_coords(mpi_cartcomm, mpi_rank, MPI_NDIM, mpi_coords+1);

  //print out diagnostics for the main MPI process
  if( 0 == mpi_rank ) {
    printf("CPU geometry: %d %d %d (total number of tasks: %d)\n", mpi_dims[1], mpi_dims[2], mpi_dims[3], mpi_numtasks);
    printf("Tile size: %d %d %d\n", N1, N2, N3);
    printf("Total resolution: %d %d %d\n", mpi_ntot[1], mpi_ntot[2], mpi_ntot[3]);
  }
  
  printf( "MPI process %4d (%4d,%4d,%4d) has PID = %d\n", mpi_rank, mpi_coords[1], mpi_coords[2], mpi_coords[3], getpid() );
  fflush(stdout);
  //sleep(20);

  //compute and store neighbor ranks: 3 dimensions times 2 directions per dimension
  for( dim = 1; dim < NDIM; dim++) {
    MPI_Cart_shift(mpi_cartcomm, dim-1, 1, &mpi_nbrs[dim][0], &mpi_nbrs[dim][1]);
    //starting value of index
    mpi_startn[dim] = mpi_coords[dim]*mpi_ntile[dim];
  }
  
#else
  //only one process without MPI
  for(dim=1;dim<NDIM;dim++){
    mpi_coords[dim]=0;
    mpi_startn[dim]=0;
    mpi_dims[dim]=1;
  }
  mpi_ntot[1]=N1;
  mpi_ntot[2]=N2;
  mpi_ntot[3]=N3;
  mpi_rank = MASTER;
  i_am_the_master = 1;
#endif  //end #ifdef MPI

  return(0);
}

//when called with stage = 0, initializes rdump_buffer
//when called with stage = 1, initializes dump_buffer and gdump_buffer
void initialize_parallel_write(int stage)
{
#if MPI && DO_PARALLEL_WRITE
    size_t nvars_dump, nvars_gdump, nvars_gdump2, nvars_rdump, nvars_fdump;
    size_t max_buffer_size_bytes, dump_buffer_size_bytes, gdump_buffer_size_bytes, gdump2_buffer_size_bytes, rdump_buffer_size_bytes, fdump_buffer_size_bytes;
    int array_of_distribs[NDIM], array_of_dargs[NDIM];
    int dim;
    int is_dry_run = 1;


    //figure out the amount of memory needed to hold each dump type
    if (stage) {
      //check if various dumps fit into the above-allocated buffer
      nvars_dump = dump(0, is_dry_run);
      nvars_gdump = gdump(is_dry_run);
      nvars_gdump2 = gdump2(is_dry_run);
      nvars_fdump = NIMG;
      
      dump_buffer_size = nvars_dump*N1*N2*N3;
      gdump_buffer_size = nvars_gdump*N1*N2*N3;
      gdump2_buffer_size = nvars_gdump2*N1*N2*N3;
      fdump_buffer_size = nvars_fdump*N1*N2*N3;
      
      dump_buffer_size_bytes = dump_buffer_size*sizeof(dumptype);
      gdump_buffer_size_bytes = gdump_buffer_size*sizeof(gdumptype);
      gdump2_buffer_size_bytes = gdump2_buffer_size*sizeof(gdump2type);
      fdump_buffer_size_bytes = fdump_buffer_size*sizeof(fdumptype);
      if (i_am_the_master) {
        printf("dump   size = %.2lg GB\n", dump_buffer_size_bytes*mpi_dims[1]*mpi_dims[2]*mpi_dims[3]/(1024.*1024.*1024.));
        printf("gdump  size = %.2lg GB\n", gdump_buffer_size_bytes*mpi_dims[1]*mpi_dims[2]*mpi_dims[3]/(1024.*1024.*1024.));
        printf("gdump2 size = %.2lg GB\n", gdump2_buffer_size_bytes*mpi_dims[1]*mpi_dims[2]*mpi_dims[3]/(1024.*1024.*1024.));
        printf("fdump  size = %.2lg GB\n", fdump_buffer_size_bytes*mpi_dims[1]*mpi_dims[2]*mpi_dims[3]/(1024.*1024.*1024.));
      }
      max_buffer_size_bytes = MY_MAX(dump_buffer_size_bytes,gdump_buffer_size_bytes);
      max_buffer_size_bytes = MY_MAX(max_buffer_size_bytes,gdump2_buffer_size_bytes);
      max_buffer_size_bytes = MY_MAX(max_buffer_size_bytes,fdump_buffer_size_bytes);
    }
    else {
      nvars_rdump = NPR;
      rdump_buffer_size = nvars_rdump*N1*N2*N3;
      rdump_buffer_size_bytes = rdump_buffer_size*sizeof(rdumptype);
      max_buffer_size_bytes = rdump_buffer_size_bytes;
    }
  
    //if already allocated, free memory
    if(mpi_file_buffer) {
      free(mpi_file_buffer);
      mpi_file_buffer = NULL;
    }
    //and then allocate anew to make sure the largest of the dumps fits byte-wise
    mpi_file_buffer = (void*) malloc(max_buffer_size_bytes);
    if (!mpi_file_buffer) {
      fprintf(stderr,"Rank %d could not allocate %ld bytes for holding mpi_file_buffer", mpi_rank, max_buffer_size_bytes);
      MPI_Abort(MPI_COMM_WORLD, errno);
    }
    //all arrays now will share the same memory
    //this possible because different types of dumps are written out in sequence
    dump_buffer = (dumptype*)mpi_file_buffer;
    gdump_buffer = (gdumptype*)mpi_file_buffer;
    gdump2_buffer = (gdump2type*)mpi_file_buffer;
    rdump_buffer = (rdumptype*)mpi_file_buffer;
    fdump_buffer = (fdumptype*)failimage; //can write directly since already contiguous array (because no ghost cells)

    //create MPI file types for each of dump types
  
    //initialize MPI arrays
    for (dim=0; dim<NDIM; dim++) {
        array_of_distribs[dim] = MPI_DISTRIBUTE_BLOCK;
        array_of_dargs[dim] = MPI_DISTRIBUTE_DFLT_DARG;
    }

    if (stage) {
      //create new cell and file types for GDUMP
      MPI_Type_contiguous(nvars_gdump, MPI_GDUMP_TYPE, &gdump_cell_type);
      MPI_Type_commit(&gdump_cell_type);
      MPI_Type_create_darray(mpi_numtasks, mpi_rank, 3,
                             mpi_ntot+1, array_of_distribs+1,
                             array_of_dargs+1, mpi_dims+1, MPI_ORDER_C,
                             gdump_cell_type, &gdump_file_type);
      MPI_Type_commit(&gdump_file_type);

      //create new cell and file types for GDUMP2
      MPI_Type_contiguous(nvars_gdump2, MPI_GDUMP2_TYPE, &gdump2_cell_type);
      MPI_Type_commit(&gdump2_cell_type);
      MPI_Type_create_darray(mpi_numtasks, mpi_rank, 3,
                             mpi_ntot+1, array_of_distribs+1,
                             array_of_dargs+1, mpi_dims+1, MPI_ORDER_C,
                             gdump2_cell_type, &gdump2_file_type);
      MPI_Type_commit(&gdump2_file_type);

      //create new cell and types for DUMP
      MPI_Type_contiguous(nvars_dump, MPI_DUMP_TYPE, &dump_cell_type);
      MPI_Type_commit(&dump_cell_type);
      MPI_Type_create_darray(mpi_numtasks, mpi_rank, 3,
                             mpi_ntot+1, array_of_distribs+1,
                             array_of_dargs+1, mpi_dims+1, MPI_ORDER_C,
                             dump_cell_type, &dump_file_type);
      MPI_Type_commit(&dump_file_type);

      //create new cell and types for FDUMP
      MPI_Type_contiguous(nvars_fdump, MPI_FDUMP_TYPE, &fdump_cell_type);
      MPI_Type_commit(&fdump_cell_type);
      MPI_Type_create_darray(mpi_numtasks, mpi_rank, 3,
                             mpi_ntot+1, array_of_distribs+1,
                             array_of_dargs+1, mpi_dims+1, MPI_ORDER_C,
                             fdump_cell_type, &fdump_file_type);
      MPI_Type_commit(&fdump_file_type);

      //if(i_am_the_master) fprintf(stderr, "dump_buffer_size = %ld bytes, nvars_dump = %ld\n", (long int)dump_buffer_size, (long int)nvars_dump);
      //if(i_am_the_master) fprintf(stderr, "gdump_buffer_size = %ld bytes, nvars_gdump = %ld\n", (long int)gdump_buffer_size, (long int)nvars_gdump);

    }
    else {
      //create new cell and types for RDUMP
      MPI_Type_contiguous(nvars_rdump, MPI_RDUMP_TYPE, &rdump_cell_type);
      MPI_Type_commit(&rdump_cell_type);
      MPI_Type_create_darray(mpi_numtasks, mpi_rank, 3,
                             mpi_ntot+1, array_of_distribs+1,
                             array_of_dargs+1, mpi_dims+1, MPI_ORDER_C,
                             rdump_cell_type, &rdump_file_type);
      MPI_Type_commit(&rdump_file_type);
      //if(i_am_the_master) fprintf(stderr, "rdump_buffer_size = %ld bytes, nvars_rdump = %ld\n", (long int)rdump_buffer_size, (long int)nvars_rdump);
    }
#endif
}

void de_initialize_parallel_write()
{
#if MPI && DO_PARALLEL_WRITE
    if (mpi_file_buffer) {
        free(mpi_file_buffer);
        mpi_file_buffer = NULL;
    }
#endif
}

void parallel_readwrite(char *file_name, void *dump_buffer,
                            int type_of_file, int is_write, long long offset)
{
#if MPI && DO_PARALLEL_WRITE
  MPI_File fh;
  MPI_Status status;
  MPI_Datatype mpi_elementary_type, mpi_file_type;
  int file_open_error, file_write_error ;
  int error_string_length;
  char error_string[BUFSIZ];
  MPI_Offset file_size;
  int count;
  void *mpi_buffer;
  size_t mpi_buffer_size;
  int mode;
  MPI_Offset mpi_offset;
    
  MPI_Barrier(MPI_COMM_WORLD);
  
  if (is_write) {
    mode = MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND;
  }
  else {
    mode = MPI_MODE_RDONLY;
  }
  file_open_error = MPI_File_open(MPI_COMM_WORLD, file_name,
                                  mode,
                                  MPI_INFO_NULL, &fh);
  if (file_open_error != MPI_SUCCESS) {
    MPI_Error_string(file_open_error, error_string,
                     &error_string_length);
    fprintf(stderr, "parallel_readwrite(): error opening file: %3d: %s\n", mpi_rank, error_string);
    MPI_Abort(MPI_COMM_WORLD, file_open_error);
    
    /* It is still OK to abort, because we have failed to
     open the file. */
    
  }
  else {
    
//    if (i_am_the_master)
//      chmod(file_name, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if (offset < 0L) {
      if(is_write) {
        MPI_File_get_position(fh, &mpi_offset);
        offset = mpi_offset;
      }
      else {
        offset = 0L;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //differentiate data type and buffers involved based on file type
    if( DUMP_FILE == type_of_file ) {
        mpi_elementary_type = MPI_DUMP_TYPE;
        mpi_file_type = dump_file_type;
        mpi_buffer = (void*)dump_buffer;
        mpi_buffer_size = dump_buffer_size;
    }
    else if( GDUMP_FILE == type_of_file){
        mpi_elementary_type = MPI_GDUMP_TYPE;
        mpi_file_type = gdump_file_type;
        mpi_buffer = (void*)gdump_buffer;
        mpi_buffer_size = gdump_buffer_size;
    }
    else if( GDUMP2_FILE == type_of_file){
      mpi_elementary_type = MPI_GDUMP2_TYPE;
      mpi_file_type = gdump2_file_type;
      mpi_buffer = (void*)gdump2_buffer;
      mpi_buffer_size = gdump2_buffer_size;
    }
    else if( RDUMP_FILE == type_of_file){
        mpi_elementary_type = MPI_RDUMP_TYPE;
        mpi_file_type = rdump_file_type;
        mpi_buffer = (void*)rdump_buffer;
        mpi_buffer_size = rdump_buffer_size;
    }
    else if( FDUMP_FILE == type_of_file){
      mpi_elementary_type = MPI_FDUMP_TYPE;
      mpi_file_type = fdump_file_type;
      mpi_buffer = (void*)fdump_buffer;
      mpi_buffer_size = fdump_buffer_size;
    }
    else {
        if(i_am_the_master)
            fprintf(stderr, "Unknown file type %d\n", type_of_file);
        MPI_File_close(&fh);
        MPI_Finalize();
        exit(2);
    }
    MPI_File_set_view(fh, offset, mpi_elementary_type, mpi_file_type, "native", MPI_INFO_NULL);
    if (is_write) {
      file_write_error =
      MPI_File_write_all(fh, mpi_buffer, mpi_buffer_size, mpi_elementary_type,
                         &status);
    }
    else {
      file_write_error =
      MPI_File_read_all(fh, mpi_buffer, mpi_buffer_size, mpi_elementary_type,
                         &status);
    }
    if (file_write_error != MPI_SUCCESS) {
      MPI_Error_string(file_write_error, error_string,
                       &error_string_length);
      fprintf(stderr, "parallel_readwrite(): error %s file: %3d: %s\n",
              (is_write)?("writing"):("reading"), mpi_rank, error_string);
      MPI_File_close(&fh);
      //if (i_am_the_master) MPI_File_delete(file_name, MPI_INFO_NULL);
      MPI_Finalize();
      exit(1);
    }
//    MPI_Get_count(&status, MPI_FLOAT, &count);
//    MPI_File_get_size(fh, &file_size);
//    if(1) {
//      printf("%3d: wrote %d floats, expected to write %lld floats\n", mpi_rank, count, (long long int)dump_buffer_size);
//      printf("%3d: file size is %lld bytes, header-related offset is %lld\n", mpi_rank, file_size, offset);
//    }

    MPI_File_close(&fh);
  }
#endif
}









































