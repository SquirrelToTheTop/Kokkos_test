#pragma once

#include <string>

#include "hzlnt.hpp"

// array
template <typename type_t> using DataArray1D = Kokkos::View< type_t *, Layout, Device >;
template <typename type_t> using DataArray1D_host = Kokkos::View< type_t *, Layout, Host >;

// name type info
template<typename T> inline std::string getCustomTypeName(const T&){ return typeid(T).name(); }
template<> inline std::string getCustomTypeName<signed char>(const signed char&){ return "char"; }
template<> inline std::string getCustomTypeName<unsigned char>(const unsigned char&){ return "uchar"; }
template<> inline std::string getCustomTypeName<uint32_t>(const uint32_t&){ return "uint32_t"; }
template<> inline std::string getCustomTypeName<uint64_t>(const uint64_t&){ return "uint64_t"; }
template<> inline std::string getCustomTypeName<int32_t>(const int32_t&){ return "int32_t"; }
template<> inline std::string getCustomTypeName<int64_t>(const int64_t&){ return "int64_t"; }
template<> inline std::string getCustomTypeName<float>(const float&){ return "float"; }
template<> inline std::string getCustomTypeName<double>(const double&){ return "double"; }

template<typename T> inline std::string getCustomTypeName( const DataArray1D<T>& ref ){
  return getCustomTypeName( T() )+"["+std::to_string(ref.size())+"]";
}
