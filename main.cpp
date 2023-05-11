#include <iostream>

// Hzlnt
#include "src/datamodel/collection.hpp"
#include "src/io/Reader_Hercule.hpp"
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

    Collection meshes = Collection(); 
    
    int timestepID = 12, contextID = 0;
    Reader_Hercule lecteur;

    lecteur.initializeReader();
    lecteur.open( "/home/squirrel/work/ramses_GrdCh/hercule_hdep" );
    uint32_t ndomains = lecteur.GetNumberOfDomains();

    // This cannot work because Kokkos requires "constness" ...;
    //
    // Kokkos::parallel_for( "read_amr", ndomains, KOKKOS_LAMBDA ( const int &idx ) {
    //   meshes.addItem( idx, lecteur.GetAMRData( timestepID, contextID, "Ramses3D" ) );
    // });

    for( uint32_t idom=0; idom<ndomains; ++idom ){
      meshes.addItem( idom, lecteur.GetAMRData( timestepID, contextID, "Ramses3D" ) );
    }

    const size_t nvals = 8850;
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
