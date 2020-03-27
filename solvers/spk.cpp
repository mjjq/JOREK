// Script for solving linear system Ax=B using STRUMPACK
// function spk can be called from the external code
// phase indicates the step (solver initialization, 
// reordering and factorization, solving, and finalizing)
// I. Holod 31/01/2020
#ifdef USE_STRUMPACK
#include <iostream>
#include <string>
#include "hdf5.h"
#include <math.h>
#include "mkl_spblas.h"

#include "StrumpackSparseSolverMPIDist.hpp"
#include "sparse/CSRMatrix.hpp"
#include <chrono>

#include <sys/sysinfo.h>
#include <unistd.h>
using namespace strumpack;

extern "C" void convert2csr(int *indx, int *n, int *nnz, int **irn, int **jcn, double **val);
int* distribute(int n, int P);
extern "C" void spk(void) {}

//===============================================================================//
extern "C" void spk_init(StrumpackSparseSolverMPIDist<double,int>** spss_,MPI_Fint* comm_) {

  StrumpackSparseSolverMPIDist<double,int>* spss= *spss_;
  MPI_Comm comm=MPI_Comm_f2c(*comm_);
  int thread_level,rank,P;
  double eps=1e-36, epsr=1.e-12, rnorm;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &P);
  MPI_Query_thread(&thread_level);
  if (thread_level != MPI_THREAD_FUNNELED && rank == 0)
    std::cout << "MPI implementation does not support MPI_THREAD_FUNNELED"
              << std::endl;

  *spss_= new StrumpackSparseSolverMPIDist<double,int>(comm);
  spss = *spss_;

  spss->options().set_matching(MatchingJob::MAX_DIAGONAL_PRODUCT_SCALING);
  spss->options().set_reordering_method(ReorderingStrategy::METIS);    
  spss->options().enable_METIS_NodeNDP();

//  spss->options().set_reordering_method(ReorderingStrategy::PARMETIS);
//  spss->options().set_reordering_method(ReorderingStrategy::SCOTCH);    
//
//  spss->options().set_Krylov_solver(KrylovSolver::PREC_GMRES);
  spss->options().set_Krylov_solver(KrylovSolver::DIRECT);
//  spss->options().set_Krylov_solver(KrylovSolver::REFINE);
  spss->options().set_rel_tol(epsr);
  spss->options().set_abs_tol(eps);  
  spss->options().set_maxit(200);  
  spss->options().set_gmres_restart(50);    
  spss->options().set_verbose(true);  
//  spss->options().enable_HSS();
//  spss->options().use_HSS();
//  spss->options().set_HSS_min_sep_size(432);
//  spss->options().HSS_options().set_rel_tol(1e-6);  
//  spss->options().HSS_options().set_abs_tol(1e-10);  
//  spss->options().enable_BLR();
//  spss->options().use_BLR();  
//  spss->options().BLR_options().set_rel_tol(1e-4);  
//  spss->options().BLR_options().set_abs_tol(1e-8);    

  return;
}
//===============================================================================//
extern "C" void spk_set_mat(int* n_,int* nnz_,int** irn_,int** jcn_,double** val_,
		StrumpackSparseSolverMPIDist<double,int>** spss_,MPI_Fint* comm_,bool* upd_) {
// set and factorize (distributed) matrix

  int n=*n_;
  int nnz=*nnz_;
  int* irn=*irn_;  
  int* jcn=*jcn_;  
  double* val=*val_;
  bool upd=*upd_;

  StrumpackSparseSolverMPIDist<double,int>* spss= *spss_;   

  MPI_Comm comm=MPI_Comm_f2c(*comm_);
  int rank,P;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &P);	
  std::chrono::steady_clock::time_point t0, t1;  

  int* dist = distribute(n, P);
  int n_local = dist[rank+1]-dist[rank];

  if (!rank) std::cout<<"n = "<< n <<" nnz = "<< nnz << " irn[nnz] = "<< irn[nnz-1] << std::endl;  
  
  int indx = 1; // Fortran indexing
  
  convert2csr(&indx, &n, &nnz, &irn, &jcn, &val);

  // splitting the matrix
  int rowstart, rowend;
  int nnzloc,i0;
  bool symmetric_pattern=true;

	if (P>1){    
		rowstart = dist[rank];
		rowend = dist[rank+1]-1;
		nnzloc=irn[rowend+1]-irn[rowstart];

		int *ptrloc = new int[n_local+1];
		int *indloc = new int[nnzloc];  
		double *valloc = new double[nnzloc];

		i0 = irn[rowstart]; // total number of nnz values before me
		ptrloc[0] = 0;
#pragma omp for		
		for (int i=0; i<n_local; i++){
			ptrloc[i+1]=irn[rowstart+1+i] - i0;
		}

#pragma omp for
		for (int i=0; i<nnzloc; i++){
			indloc[i]=jcn[i+i0];	  
			valloc[i]=val[i+i0];
		}
		if (rank==0){
			std::cout<<"ptrloc: "<<ptrloc[0]<<" "<<ptrloc[n_local]<<std::endl;  
			std::cout<<"indloc: "<<indloc[0]<<" "<<indloc[nnzloc-1]<<std::endl;  
			std::cout<<"valloc: "<<valloc[0]<<" "<<valloc[nnzloc-1]<<std::endl;
		}
#ifdef NEWSPK	
		if (upd){
			std::cout<<"Updating matrix values"<<std::endl;
			spss->update_matrix_values(n_local, ptrloc, indloc, valloc, \
					dist, symmetric_pattern);
		}else
#endif			
		{
			spss->set_distributed_csr_matrix(n_local, ptrloc, indloc, valloc, dist);
		}
		delete ptrloc, indloc, valloc;
	}
	else
	{
#ifdef NEWSPK
		if (upd){
			std::cout<<"Updating matrix values"<<std::endl;
			spss->update_matrix_values(n_local, irn, jcn, val, dist, symmetric_pattern);
		}else
#endif
		{
			spss->set_distributed_csr_matrix(n_local, irn, jcn, val, dist);
		}
	}
	return;
}
//===============================================================================//
extern "C" void spk_reord(StrumpackSparseSolverMPIDist<double,int>** spss_,MPI_Fint* comm_) {
// reorder (distributed) matrix

  StrumpackSparseSolverMPIDist<double,int>* spss= *spss_;   

  MPI_Comm comm=MPI_Comm_f2c(*comm_);
  int rank,P;	
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &P);	
  std::chrono::steady_clock::time_point t0, t1;  

  // Reordering	
  t0 = std::chrono::steady_clock::now();    
  spss->reorder();
  t1 = std::chrono::steady_clock::now();
  if (!rank)
    std::cout<<"Time to reorder (s) = "<< std::chrono::duration_cast<
    std::chrono::microseconds>(t1 - t0).count()*1e-6 << std::endl;  
  
  return;
}
//===============================================================================//
extern "C" void spk_fact(StrumpackSparseSolverMPIDist<double,int>** spss_,MPI_Fint* comm_) {
// factorize (distributed) matrix

  StrumpackSparseSolverMPIDist<double,int>* spss= *spss_;   

  MPI_Comm comm=MPI_Comm_f2c(*comm_);
  int rank,P;	
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &P);	
  std::chrono::steady_clock::time_point t0, t1;  

  // Factorization	
  t0 = std::chrono::steady_clock::now();  
  spss->factor();	    
  t1 = std::chrono::steady_clock::now();
  if (!rank)
      std::cout<<"Time to factorize (s) = "<< std::chrono::duration_cast<
      std::chrono::microseconds>(t1 - t0).count()*1e-6 << std::endl;	

  return;
}
//===============================================================================//
extern "C" void spk_solve(int* n_, double** rhs_,
		StrumpackSparseSolverMPIDist<double,int>** spss_,MPI_Fint* comm_,int* phase) {

  int n=*n_;
  double* rhs=*rhs_;

  StrumpackSparseSolverMPIDist<double,int>* spss= *spss_;   

  MPI_Comm comm=MPI_Comm_f2c(*comm_);
  int thread_level,rank,P;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &P);	
  std::chrono::steady_clock::time_point t0, t1;  

  int* dist = distribute(n, P);
  int n_local = dist[rank+1]-dist[rank];

  // set local RHS    
	std::vector<double> b(n_local), x(n_local);

#pragma omp for
	for (int i=dist[rank]; i<dist[rank+1]; i++)
		b[i-dist[rank]]=rhs[i];    

	t0 = std::chrono::steady_clock::now();
	spss->solve(b.data(),x.data(),false);

// Gather the solution
	std::vector<double> x_glob(n), x_buf(n);
	x_glob.assign(n,0);
	x_buf.assign(n,0);

#pragma omp for
	for (int i=dist[rank]; i<dist[rank+1]; i++)
			x_buf[i]=x[i-dist[rank]];

	MPI_Allreduce(x_buf.data(), x_glob.data(), n, MPI_DOUBLE_PRECISION, MPI_SUM, comm);
		
	t1 = std::chrono::steady_clock::now();
	if (!rank){
		std::cout<<"Time to solve (s) = "<< std::chrono::duration_cast<
			std::chrono::microseconds>(t1 - t0).count()*1e-6 << std::endl;     
	}

#pragma omp for	
	for (int i=0;i<n;i++){
				(*rhs_)[i] = x_glob[i];
	}

	x.clear();
	b.clear();
	x_glob.clear();
	x_buf.clear(); 

	return;
}
//==========================================================================================//
extern "C" void spk_finalize(StrumpackSparseSolverMPIDist<double,int>** spss_,MPI_Fint* comm_) {
  StrumpackSparseSolverMPIDist<double,int>* spss= *spss_;
  MPI_Comm comm=MPI_Comm_f2c(*comm_);

	delete *spss_;
	scalapack::Cblacs_exit(1);
	return;
}  
//==========================================================================================//
extern "C" void convert2csr(int *indx_, int *n_, int *nnz_, int **irn, int **jcn, double **val)
{
  //int *rowptrE;
  int n =*n_, nnz=*nnz_, indx=*indx_;

  sparse_index_base_t    indexing;  
  sparse_matrix_t cooA;
  sparse_matrix_t csrA;  
  sparse_status_t stat;  

// create mkl coordinate sparse matrix
  if (indx==1){
    mkl_sparse_d_create_coo(&cooA, SPARSE_INDEX_BASE_ONE, n, n, nnz, *irn, *jcn, *val);  
  } else {
    mkl_sparse_d_create_coo(&cooA, SPARSE_INDEX_BASE_ZERO, n, n, nnz, *irn, *jcn, *val);
  }

// convert to csr format  
  mkl_sparse_convert_csr(cooA, SPARSE_OPERATION_NON_TRANSPOSE, &csrA);
  mkl_sparse_destroy(cooA);
  mkl_sparse_order(csrA); // important

// export csr values, rowptr (Begin and End counting) and colind
  MKL_INT *rowptrE, *rowptrB, *colind;
  double *values;
  mkl_sparse_d_export_csr(csrA, &indexing, &n, &n, &rowptrB, &rowptrE, &colind, &values);
  
  nnz = rowptrE[n-1] - indx; 
  if (nnz!=(*nnz_)) 
	  std::cout<<"New nnz: "<<nnz<<" Old nnz "<< *nnz_<<std::endl;    
  *nnz_ = nnz;

#pragma omp for  
  for (int i=0; i<n; i++){
	  (*irn)[i] = rowptrB[i] - indx;
  }
  
  (*irn)[n]=nnz;

#pragma omp for
  for (int i=0; i<nnz; i++){
	  (*jcn)[i] = colind[i] -indx;
	  (*val)[i] = values[i];
  }  

  mkl_sparse_destroy(csrA);

  return;
}
//=========================================================================================//
// distribute n rows among P ranks
int* distribute(int n, int P){
    int* dist;  
    int* nl;
    int n_local;

    dist = new int[P+1];
    nl = new int[P];      

    for (int i=0; i<P; i++){
  	  nl[i] = floor(n/P);
  	  if (i<n%P)
		  nl[i]+=1;
    }  
    dist[0]=0;
    for (int i=0; i<P; i++){
  	  dist[i+1] = dist[i] + nl[i];
    }
    return dist;
}
//=========================================================================================//
extern "C" void getmem(float *av, float *tot){
    long avpg = sysconf(_SC_AVPHYS_PAGES);
    long ppg = sysconf(_SC_PHYS_PAGES);
    long pgsize = sysconf(_SC_PAGESIZE);

    *tot = ppg*pgsize/1e9;
    *av = avpg*pgsize/1e9;

    return;
}

#endif

