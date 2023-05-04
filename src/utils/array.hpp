#pragma once

#include "hzlnt.hpp"

template <typename type_t> using DataArray1D = Kokkos::View< type_t *, Layout, Device >;