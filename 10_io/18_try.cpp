#include "hdf5.h"
#include "stdlib.h"
#include <cassert>
#include <cstdio>
#include <chrono>

using namespace std;
#define H5FILE_NAME     "SDS_pat.h5"
#define DATASETNAME 	"IntArray" 
#define NX     10000                      /* dataset dimensions */
#define NY     10000
#define RANK   2
#define RANK1  1

int main (int argc, char **argv)
{
	/*
	 * HDF5 APIs definitions
	 */ 	
	hid_t       file_id, dset_id;         /* file and dataset identifiers */
	hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
	hsize_t     dimsf[2];                 /* dataset dimensions */
	hsize_t     dimsm[1];                 /* dataset dimensions */
	int         *data;                    /* pointer to data buffer to write */
	hsize_t	count[2];	          /* hyperslab selection parameters */
	hsize_t	stride[2];
	hsize_t	offset[2];
	hsize_t     block[2];
	hid_t	plist_id;                 /* property list identifier */
	int         i;
	herr_t	status;

	/*
	 * MPI variables
	 */
	int mpi_size, mpi_rank;
	MPI_Comm comm  = MPI_COMM_WORLD;
	MPI_Info info  = MPI_INFO_NULL;

	/*
	 * Initialize MPI
	 */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(comm, &mpi_size);
	MPI_Comm_rank(comm, &mpi_rank);  
	/*
	 * Exit if number of processes is not 4. 
	 */
	if (mpi_size != 4) {
		printf("This example to set up to use only 2 processes \n");
		printf("Quitting...\n");
		return 0;
	}

	/* 
	 * Set up file access property list with parallel I/O access
	 */
	plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, comm, info);

	/*
	 * Create a new file collectively and release property list identifier.
	 */
	file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
	H5Pclose(plist_id);


	/*
	 * Create the dataspace for the dataset.
	 */
	dimsf[0] = NX;
	dimsf[1] = NY;
	dimsm[0] = NX*NY/4; 
	block[0] = 2500;
	block[1] = 2500;  
	filespace = H5Screate_simple(RANK, dimsf, NULL); 
	memspace  = H5Screate_simple(RANK1, dimsm, NULL); 

	/*
	 * Create the dataset with default properties and close filespace.
	 */
	dset_id = H5Dcreate(file_id, DATASETNAME, H5T_NATIVE_INT, filespace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(filespace);

	/* 
	 * Each process defines dataset in memory and writes it to the hyperslab
	 * in the file.
	 */
	count[0] = 2;
	count[1] = 2;
	stride[0] = 5000;
	stride[1] = 5000;
	if(mpi_rank == 0) {
		offset[0] = 0;
		offset[1] = 0;
	} 
	if(mpi_rank == 1) {
		offset[0] = 2500;
		offset[1] = 0;
	}
	if(mpi_rank == 2) {
		offset[0] = 0;
		offset[1] = 2500;
	} 
	if(mpi_rank == 3) {
		offset[0] = 2500;
		offset[1] = 2500;
	} 

	/*
	 * Select hyperslab in the file.
	 */
	filespace = H5Dget_space(dset_id);
	status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, stride, count, block);

	/*
	 * Initialize data buffer 
	 */
	data = (int *) malloc(sizeof(int)*dimsm[0]);
	for (i=0; i < (int)dimsm[0]; i++) {
		data[i] = mpi_rank + 1;
	}

	/*
	 * Create property list for collective dataset write.
	 */
	plist_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	auto tic = chrono::steady_clock::now();

	status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
			plist_id, data);

	auto toc = chrono::steady_clock::now();
	free(data);

	/*
	 * Close/release resources.
	 */
	H5Dclose(dset_id);
	H5Sclose(filespace);
	H5Sclose(memspace);
	H5Pclose(plist_id);
	H5Fclose(file_id);
	double time = chrono::duration<double>(toc - tic).count();
	printf("N=%d: %lf s (%lf GB/s)\n",NX*NY,time,4*NX*NY/time/1e9); 
	MPI_Finalize();

	return 0;
}     
