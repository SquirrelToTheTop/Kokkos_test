/* This file is part of the 'Hazelnut' package.
 *
 *  Author: L. Strafella
 */
#pragma once

#include <string>

// kokkos stuff
#include <Kokkos_Core.hpp>

// some definition for kokkos
using Device = Kokkos::DefaultExecutionSpace;
using Host   = Kokkos::DefaultHostExecutionSpace;
using Layout = Kokkos::LayoutRight;

#define HDF5_KEY_MAXLENGTH 128

#ifdef NDEBUG
#include <iostream>
// #define HZL_TRACE( msg ) do { std::cerr <<  "HZL::Trace[" << __FILE__ << "]::" << __LINE__ << ":" \
//                                           <<  msg << std::endl; } while (0)
#define HZL_TRACE( msg ) do { std::cerr <<  "HZL::Trace[line " << __LINE__ << "] :" \
                                          <<  msg << std::endl; } while (0)
#else
#define HZL_TRACE( msg )
#endif

#ifdef NDEBUG
#include <iostream>
#define HZL_LOG( msg ) do { std::cerr <<  "HZL::Log[line " << __LINE__ << "] :" \
                                          <<  msg << std::endl; } while (0)
#else
#define HZL_LOG( msg )
#endif

template<typename T> inline std::string getTypeName(){
  assert( "Unknown h5_type ! " && false);
  return "nuts";
};

template<> inline std::string getTypeName<std::string>(){ return "str"; }
template<> inline std::string getTypeName<char*>(){ return "char *"; }
template<> inline std::string getTypeName<char>(){ return "int8_t"; }
template<> inline std::string getTypeName<unsigned char>(){ return "uint8_t"; }
template<> inline std::string getTypeName<short>(){ return "int16_t"; }
template<> inline std::string getTypeName<int>(){ return "int32_t"; }
template<> inline std::string getTypeName<unsigned int>(){ return "uint32_t"; }
template<> inline std::string getTypeName<long int>(){ return "int64_t"; }
template<> inline std::string getTypeName<unsigned long int>(){ return "uint64_t"; }
template<> inline std::string getTypeName<float>(){ return "float32_t"; }
template<> inline std::string getTypeName<double>(){ return "float64_t"; }
template<> inline std::string getTypeName<bool>(){ return "bool"; }