#include <iostream>

// Hzlnt
#include "src/datamodel/collection.hpp"
#include "src/io/Reader_Hercule.hpp"
#include "src/utils/definition.hpp"

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
    
    int timestepID = 12;
    Reader_Hercule lecteur;

    lecteur.initializeReader();
    // lecteur.open( "/local/home/ls256408/work/ramses_GrdCh/blast2d/" );
    lecteur.open( "/home/squirrel/work/ramses_GrdCh/blast2d/" );

    uint32_t ndomains = lecteur.GetNumberOfDomains();

    // This cannot work because Kokkos requires "constness" and Hercule does probably not work
    // with openMP reads
    //
    // Kokkos::parallel_for( "read_amr", ndomains, KOKKOS_LAMBDA ( const int &idx ) {
    //   meshes.addItem( idx, lecteur.GetAMRData( timestepID, contextID, "Ramses3D" ) );
    // });

    // I/O -> doit se passer sur le Host
    for( uint32_t idom=0; idom<ndomains; ++idom ){
      HZL_TRACE("[Main] reading domain # " << idom );
      meshes.addItem( idom, lecteur.GetAMRData( timestepID, idom, "Ramses2D", 2 ) );
    }

    HZL_TRACE( "[Main] computing 'uncompress' ... " );
    meshes.uncompressAMRDesc();

    HZL_TRACE( "[Main] computing 'testKokkos' ... " );
    meshes.testKokkos();

    // const size_t nvals = 8850;
    // DataArray1D<float> a( "test", nvals );

    // // Initialize y vector: is the vector copied ?!
    // Kokkos::parallel_for( "a_init", nvals, KOKKOS_LAMBDA ( const int &idx ) {
    //   a[ idx ] = idx;
    // });

    // double sum = 0.0;
    // Kokkos::parallel_reduce("kpr_sum", nvals, KOKKOS_LAMBDA ( const int& idx, double& lsum ) {
    //   lsum += a[ idx ] + 1 ;
    // }, sum );


    // std::cout << "\t> Sum of numbers from [0," << nvals << "] : " << sum << std::endl;
  }

  Kokkos::finalize();
  
#ifdef WITH_MPI
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;
}
