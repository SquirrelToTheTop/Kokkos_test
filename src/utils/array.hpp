#pragma once

#include "hzlnt.hpp"

template <typename type_t> using DataArray1D = Kokkos::View< type_t *, Layout, Device >;

template <typename type_t> using DataArray1D_host = Kokkos::View< type_t *, Layout, Host >;