#pragma once

#include "../utils/definition.hpp"

#include "mesh_utils.h"
#include "lightAMR_utils.h"

#include "encode.h"
#include "decode.h"
#include "decompression.h"

// class LightAMR;

class LightAMR {

  public:

    LightAMR( int dim ){
      _dim = dim;
      _isCompressed = false;
      _nbCells = 0;
      _nbChildPerLeaf = 8;
      _nbChildPerNode = 2;
    };

    LightAMR(const LightAMR &) { std::cout << "[lightAMR] Je fais une copie " << std::endl; };
    LightAMR(LightAMR &&) = default;

    /* accesseur aux data pour l'AMR */
    const DataArray1D_host<uint8_t>  & getRefinementArrayConstPtr() const { return _refinementArray; }
    const DataArray1D_host<uint8_t>  & getOwnershipArrayConstPtr()  const { return _ownershipArray; }
    const DataArray1D_host<uint64_t> & getNcplArrayConstPtr()       const { return _ncpl; }
    
    inline uint64_t getNumberOfCells()        const { return _nbCells; };
    inline int      getDimension()            const { return _dim; };
    inline int      getNumberOfChildPerNode() const { return _nbChildPerNode; };
    inline int      getNumberOfChildPerLeaf() const { return _nbChildPerLeaf; };
    inline int      getNumbeOfLevels()        const { return _ncpl.size(); };

    DataArray1D_host<uint8_t>  & getRefinementArrayPtr() { return _refinementArray; }
    DataArray1D_host<uint8_t>  & getOwnershipArrayPtr()  { return _ownershipArray; }
    DataArray1D_host<uint64_t> & getNcplArrayPtr()       { return _ncpl; }

    void setAMRCompressionState( bool compressed ) { _isCompressed = compressed; }
    bool getAMRCompressionState() const { return _isCompressed; }

    void setNumberOfCells( int64_t nbCells ) { assert( nbCells > 0 ); _nbCells = nbCells; }


    void uncompress() {

      int8_t * ptrt = reinterpret_cast<int8_t *>( _refinementArray.data() );
      int8_t * ptrm = reinterpret_cast<int8_t *>( _ownershipArray.data() );

      auto tmp_tree = uncompressAMR( ptrt, _refinementArray.size(), 0, _ncpl.size(), _nbCells );
      auto tmp_mask = uncompressAMR( ptrm, _ownershipArray.size(), 0, _ncpl.size(), _nbCells  );

      Kokkos::resize( _refinementArray, _nbCells );
      Kokkos::resize( _ownershipArray,  _nbCells );

      for(size_t ic=0; ic<_nbCells; ++ic ){
        _refinementArray[ ic ] = tmp_tree[ ic ];
        _ownershipArray [ ic ] = tmp_mask[ ic ];
      }

      // free( tmp_tree );
      // free( tmp_mask );

    }

  private:

    /* données AMR (modèle lightAMR, possible compression CPS52 de _refinementArray et _ownershipArray) */
    DataArray1D_host<uint8_t> _refinementArray, _ownershipArray;
    DataArray1D_host<uint64_t> _ncpl;
    
    // definis sur l'état des _refinementArray & _ownershipArray est compressé (CPS52)
    bool _isCompressed;

    // nombre de cellule (correspond à la taille de _refinementArray & _ownershipArray si décompressé )
    // toute cellule confondue (feuille ou noeud)
    uint64_t _nbCells;

    int _dim;
    int _nbChildPerNode;
    int _nbChildPerLeaf;
};