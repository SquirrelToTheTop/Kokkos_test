#include <iostream>

// Hzlnt
#include "src/utils/array.hpp"

// Kokkos stuff
#include <Kokkos_Core.hpp>

// MPI
#include "mpi.h"

int main( int argc, char *argv[] ){

  std::cout << "\t> [Hzlnt] Hello !" << std::endl;

#ifdef WITH_MPI
  MPI_Init( &argc, &argv );
#endif

  Kokkos::initialize( argc, argv );
  {

    const size_t nvals = 885000000;
    DataArray1D<float> a( "test", nvals );

    // Initialize y vector.
    Kokkos::parallel_for( "a_init", nvals, KOKKOS_LAMBDA ( const int &idx ) {
      a[ idx ] = idx;
    });

    double sum = 0.0;
    Kokkos::parallel_reduce("kpr_sum", nvals, KOKKOS_LAMBDA ( const int& idx, double& lsum ) {
        lsum += a[ idx ] + 1 ;
    }, sum );


    std::cout << "\t> Sum of numbers from [0," << nvals << "] : " << sum << std::endl;
  }

  Kokkos::finalize();
  
#ifdef WITH_MPI
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;
}