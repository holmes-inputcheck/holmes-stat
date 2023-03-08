cmake_minimum_required (VERSION 3.0)
project (ucb-holmes)

find_path(CMAKE_FOLDER NAMES cmake/emp-tool-config.cmake)
include(${CMAKE_FOLDER}/cmake/emp-base.cmake)

find_package(emp-ot REQUIRED)
include_directories(${EMP-OT_INCLUDE_DIRS})

find_package(emp-zk-holmes REQUIRED)
include_directories(${EMP-ZK-HOLMES_INCLUDE_DIRS})

find_package(GMP REQUIRED)
include_directories(${GMP_INCLUDE_DIR})
find_package(FLINT REQUIRED)
include_directories(${FLINT_INCLUDE_DIR}/flint)
find_package(MPFR REQUIRED)
include_directories(${MPFR_INCLUDE_DIR})


if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
   set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -ggdb -Wno-ignored-attributes -Xpreprocessor")
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Wno-ignored-attributes -Xpreprocessor")
else()
   set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -ggdb -Wno-ignored-attributes")
   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Wno-ignored-attributes")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/)

set(binaries_graph graph_t_test graph_z_test graph_chi_squared graph_f_test graph_t_test_dataset1 graph_z_test_dataset1 graph_f_test_dataset1 graph_chi_squared_dataset1)
foreach(bin ${binaries_graph})
   add_executable("${bin}" "graph/${bin}.cpp")
   target_link_libraries("${bin}" ${EMP-OT_LIBRARIES} ${EMP-ZK-HOLMES_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${FLINT_LIBRARIES})
endforeach()
