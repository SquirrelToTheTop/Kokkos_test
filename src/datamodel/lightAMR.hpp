#pragma once

#include "../utils/array.hpp"

class LightAMR;

class LightAMR {

  public:

    LightAMR(){
      _isCompressed = false;
      _nbCells = 0;
    };

    LightAMR(const LightAMR &) { std::cout << "[lightAMR] Je fais une copie " << std::endl; };
    // LightAMR(LightAMR &&) = default;

    /* accesseur aux data pour l'AMR */
    const DataArray1D_host<uint8_t> & getRefinementArrayConstPtr() const { return _refinementArray; }
    const DataArray1D_host<uint8_t> & getOwnershipArrayConstPtr()  const { return _ownershipArray; }
    const DataArray1D_host<uint64_t> & getNcplArrayConstPtr()  const { return _ncpl; }

    DataArray1D_host<uint8_t> & getRefinementArrayPtr() { return _refinementArray; }
    DataArray1D_host<uint8_t> & getOwnershipArrayPtr()  { return _ownershipArray; }
    DataArray1D_host<uint64_t> & getNcplArrayPtr() { return _ncpl; }

    void setAMRCompressionState( bool compressed ) { _isCompressed = compressed; }
    bool getAMRCompressionState() const { return _isCompressed; }

    void setNumberOfCells( int64_t nbCells ) { assert( nbCells > 0 ); _nbCells = nbCells; }

  private:

    /* données AMR (modèle lightAMR, possible compression CPS52 de _refinementArray et _ownershipArray) */
    DataArray1D_host<uint8_t> _refinementArray, _ownershipArray;
    DataArray1D_host<uint64_t> _ncpl;
    
    // definis sur l'état des _refinementArray & _ownershipArray est compressé (CPS52)
    bool _isCompressed;

    // nombre de cellule (correspond à la taille de _refinementArray & _ownershipArray si décompressé )
    // toute cellule confondue (feuille ou noeud)
    uint64_t _nbCells;
};


  static std::vector<MeshUtils::Cell>
  WA_uc_ra_based( const LightAMR &lamr, uint8_t maxLevel ) {

    // vector mesh of current domain
    std::vector<MeshUtils::Cell> vmeshDom;
    vmeshDom.reserve( lamr.get_ncells()*0.875 );

    int nchildren = 1 << lamr.getDim();

    assert( lamr.getRefinementFactor() == 2 );
    assert( nchildren == 8 || nchildren == 4); // warning dyablo

    // accessor to mesh data
    auto &tree = lamr.getConstReferenceAMRtree();
    auto &mask = lamr.getConstReferenceMask();
    auto &ncpl = lamr.getConstReferenceNCellsPerLevel();

    std::vector<uint64_t> ncplsum = LightAMR_Utils::precomputeSumIndex(ncpl);
    std::vector<uint64_t> nCCpl = LightAMR_Utils::precomputeCoarsePerLevel(tree, ncpl);
    auto nCCplIdx = LightAMR_Utils::precomputeCoarseIndexPerLevel(tree, ncpl);

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

        // skip children that were part of terminal node
        if( sctn[ icell ] == 1 ) continue;

        if( tree[ icell ] == 0 ){  // leaf cell
            
          if( mask[ icell ] == 0 ) { // belong to current domain and should not be skipped
            uint8_t clvl = LightAMR_Utils::findLevelFromIndex(ncplsum, icell);
            Logical_Pos_t ijk = LightAMR_Utils::getLogicalIJKFromLightAMRIndex(nCCpl, nCCplIdx, ncplsum, clvl, icell, nchildren);

            MeshUtils::Cell a = { ijk, icell, clvl, false };
            vmeshDom.push_back( a );
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
              vmeshDom.push_back( a );
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
              Logical_Pos_t ijk = LightAMR_Utils::getLogicalIJKFromLightAMRIndex(nCCpl, nCCplIdx, ncplsum, clvl, icell, nchildren);
              
              // std::cerr << "\t> [WA_uc_ra_based] icell[" << icell << "] T-node {" << ijk.i << ", " << ijk.j 
              //       << ", " << ijk.k << "} " << std::endl;
              
              MeshUtils::Cell a = { ijk, static_cast<uint32_t>( firstChildIdx ), clvl, true };
              vmeshDom.push_back( a );

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

    // expect RVO
    return vmeshDom;
  }