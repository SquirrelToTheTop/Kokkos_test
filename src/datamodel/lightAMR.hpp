#include "../utils/array.hpp"

class LightAMR;

class LightAMR {

  public:
    LightAMR(){};
    LightAMR( const LightAMR & ) = delete;
    LightAMR( LightAMR && ) = default;

    const DataArray1D_host<uint8_t> & getRefinementArrayConstPtr() const { return _refinementArray; }
    const DataArray1D_host<uint8_t> & getOwnershipArrayConstPtr()  const { return _ownershipArray; }
    const DataArray1D_host<uint64_t> & getNcplArrayConstPtr()  const { return _ncpl; }

    DataArray1D_host<uint8_t> & getRefinementArrayPtr() { return _refinementArray; }
    DataArray1D_host<uint8_t> & getOwnershipArrayPtr()  { return _ownershipArray; }
    DataArray1D_host<uint64_t> & getNcplArrayPtr() { return _ncpl; }

  private:

    DataArray1D_host<uint8_t> _refinementArray, _ownershipArray;
    DataArray1D_host<uint64_t> _ncpl;
};