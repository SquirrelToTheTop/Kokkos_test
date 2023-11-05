#include "../utils/hzlnt.hpp"
#include "Kokkos_UnorderedMap.hpp"

#include "lightAMR.hpp"
#include "mesh_utils.h"

#include <unordered_map>

class Collection;

// using Collection_container = Kokkos::UnorderedMap<uint32_t, LightAMR>;
using Collection_container = std::unordered_map<uint32_t, LightAMR>;

class Collection{

  public:

    Collection(){};
    Collection( const Collection & ) = delete;

    ~Collection(){};

    /**
     * @param key (in): identifiant unique pour l'élement LightAMR dans la collection (ID du domaine)
     * @param rv  (in): objet de type LightAMR
    */
    inline void addItem( uint32_t key, LightAMR && rv ) noexcept {
      HZL_TRACE( "[Collection]::addItem LightAMR with key : " << key );
      HZL_TRACE("[Collection]::addItem LightAMR Ncells # " << rv.getNumberOfCells() );
#ifdef USE_KOKKOS
      if( _collections.exists( key ) ){
        std::cerr << "\t>[Collection]::addItem key ('" << key << "') already exists !" << std::endl;
      }else{
        _collections.insert( key, rv );
        _keys.emplace_back( key );
      }
#endif

      if( _collections.find( key ) != _collections.end() ){
        std::cerr << "\t>[Collection]::addItem key ('" << key << "') already exists !" << std::endl;
      }else{
        _collections.emplace( key, std::move(rv) );
        _keys.emplace_back( key );
      }

    }

    inline uint32_t GetNumberOfItems() const { return static_cast<uint32_t>( _collections.size() ); }

    void uncompressAMRDesc() {
      for(auto & a : _collections )
        a.second.uncompress();
    }

    void testKokkos() {

      HZL_TRACE( "[Collection]::testKokkos " );
      // for(int i=0; i<_keys.size(); ++ i){
      Kokkos::parallel_for( "dump info", _keys.size(), KOKKOS_LAMBDA ( const int & i ) {
        
        HZL_TRACE("[Collection]::testKokkos Processing domain # " << i );

        const auto & lamr = _collections.at( _keys[ i ] );        
        WA_uc_ra_based( lamr, lamr.getNumbeOfLevels()+1 );

      });

      // }

      HZL_TRACE( "[Colleciton]::testKokkos Number of elements : " << _cells.size() );

    }


    void WA_uc_ra_based( const LightAMR & lamr, uint8_t maxLevel ) {

      int nchildren = lamr.getNumberOfChildPerLeaf();
      assert( nchildren == 8 || nchildren == 4); // warning dyablo

      // accessor to mesh data
      const auto &tree = lamr.getRefinementArrayConstPtr();
      const auto &mask = lamr.getOwnershipArrayConstPtr();
      const auto &ncpl = lamr.getNcplArrayConstPtr();

      auto ncplsum  = LightAMR_Utils::precomputeSumIndex( ncpl );
      auto nCCpl    = LightAMR_Utils::precomputeCoarsePerLevel( tree, ncpl );
      auto nCCplIdx = LightAMR_Utils::precomputeCoarseIndexPerLevel( tree, ncpl );

      auto firstChildArray = LightAMR_Utils::computeFirstChildIndex( tree, ncpl, nchildren );

      // counter of coarse cell per level
      // std::vector<uint64_t> nPcpl( ncpl.size(), 0 ); 

      // skip child if T-node; 0 ->  do not skip, 1 -> skip the cell
      std::vector<uint8_t> sctn( tree.size(), 0 ); 

      // could be openmp
      // omp_set_num_threads(1); // deactivate multithreading -> before add #pragma omp atomic on write in hashmap
      // #pragma omp parallel
      {
        // #pragma omp for
        for(uint32_t icell = 0; icell<tree.size(); ++icell){
          // HZL_LOG( "[WA_uc_ra_based]::Cell [" << icell );

          // skip children that were part of terminal node
          if( sctn[ icell ] == 1 ) continue;

          if( tree[ icell ] == 0 ){  // leaf cell
              
            if( mask[ icell ] == 0 ) { // belong to current domain and should not be skipped
              uint8_t clvl = LightAMR_Utils::findLevelFromIndex(ncplsum, icell);
              Logical_Pos_t ijk = LightAMR_Utils::getLogicalIJKFromLightAMRIndex( nCCpl, nCCplIdx, ncplsum, clvl, icell, nchildren );

              MeshUtils::Cell a = { ijk, icell, clvl, false };
              if( ! _cells.exists( a ) ) _cells.insert( a );
            }

          }else { // coarse cell

            // trois options: - on est pas au niveau max demandé et c'est pas un terminal node -> skip
            //                - on est pas au niveau max demandé et c'est un terminal node     -> add 
            //                - on est au niveau max demandé et elle nous appartient           -> add
            //                - on est au niveau max demandé et elle nous appartient pas       -> skip
            uint8_t clvl = LightAMR_Utils::findLevelFromIndex(ncplsum, icell);
            if( clvl > maxLevel ) break;

            // on a atteind le niveau max demandé donc si elle nous appartient on la considere feuille
            if ( clvl == maxLevel ){
            
              if( mask[ icell ] == 0 ) {
                Logical_Pos_t ijk = LightAMR_Utils::getLogicalIJKFromLightAMRIndex(nCCpl, nCCplIdx, ncplsum, clvl, icell, nchildren);
                MeshUtils::Cell a = { ijk, icell, clvl, false };
                if( ! _cells.exists( a ) ) _cells.insert( a );
              }
              
            }else{
                
              // uint64_t firstChildIdx = nchildren * nPcpl[ clvl ] + ncplsum[ clvl ];
              uint64_t firstChildIdx = firstChildArray[ icell ];

              // std::cerr << "\t> [WA_uc_ra_based] icell[" << icell << "] firstChildIdx @ " << firstChildIdx << std::endl;

              uint8_t sumRafLeaves = 0;
              uint8_t sumMaskLeaves = 0;
              for (size_t ileaf = 0; ileaf < nchildren; ++ileaf) {
                assert( firstChildIdx + ileaf < tree.size() &&  firstChildIdx + ileaf < mask.size() );
                sumRafLeaves += tree[ firstChildIdx + ileaf ];
                sumMaskLeaves += mask[ firstChildIdx + ileaf ];
              }

              if ( sumRafLeaves == 0 && sumMaskLeaves == 0 ){
                Logical_Pos_t ijk = LightAMR_Utils::getLogicalIJKFromLightAMRIndex( nCCpl, nCCplIdx, ncplsum, clvl, icell, nchildren );
                
                // std::cerr << "\t> [WA_uc_ra_based] icell[" << icell << "] T-node {" << ijk.i << ", " << ijk.j 
                //       << ", " << ijk.k << "} " << std::endl;
                MeshUtils::Cell a = { ijk, static_cast<uint32_t>( firstChildIdx ), clvl, true };
                if( ! _cells.exists( a ) ) _cells.insert( a );

                // set skip for children
                for( size_t ileaf=firstChildIdx; ileaf < firstChildIdx+nchildren; ++ileaf )
                  sctn[ ileaf ] = 1;

              }
              //else{
              //     std::cerr << "\t> [WA_uc_ra_based] icell[" << icell << "] skipping coarse " << std::endl;
              // }

            }

          }

        }
          
      }

    }
    
  private:

    Collection_container _collections;
    Kokkos::UnorderedMap<MeshUtils::Cell, uint32_t> _cells;

    std::vector<uint32_t> _keys;
  
};