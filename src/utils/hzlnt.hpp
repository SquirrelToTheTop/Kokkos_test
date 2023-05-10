#pragma once

// kokkos stuff
#include <Kokkos_Core.hpp>

// some definition for kokkos
using Device = Kokkos::DefaultExecutionSpace;
using Host   = Kokkos::DefaultHostExecutionSpace;
using Layout = Kokkos::LayoutRight;
