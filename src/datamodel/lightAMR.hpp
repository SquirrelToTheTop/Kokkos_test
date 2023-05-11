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