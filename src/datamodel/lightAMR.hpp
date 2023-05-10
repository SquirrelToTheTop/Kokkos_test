#include "../utils/array.hpp"

class LightAMR;

class LightAMR {

  public:
    LightAMR(){};
    LightAMR(const LightAMR &) = delete;
    LightAMR(LightAMR &&) = default;

    const DataArray1D_host<uint8_t> & getRefinementArrayPtr() const { return _refinementArray; }
    const DataArray1D_host<uint8_t> & getOwnershipArrayPtr() const { return _ownershipArray; }

  private:
    DataArray1D_host<uint8_t> _refinementArray, _ownershipArray;
};