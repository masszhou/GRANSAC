CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(GRANSACSample)

IF(MSVC)
  ADD_DEFINITIONS(-DNOMINMAX)
  # Another MSVC quirk. http://stackoverflow.com/questions/6832666/lnk2019-when-including-asio-headers-solution-generated-with-cmake
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /MP") # For multi-threaded compilation on MSVC
ELSE()
  IF(APPLE)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++") # C++11 support on Apple clang
  ELSE()
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")# C++11 support on GCC
  ENDIF()
ENDIF()

# add opencv for visualization
set (OpenCV_DIR /usr/local/opencv/share/OpenCV)
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

# add eigen3 for least square
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

# add OpenMP for multi-threads support
find_package(OpenMP)
IF(OPENMP_FOUND)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  MESSAGE(STATUS "Found OpenMP")
ENDIF()

# add GRANSAC Header
include_directories(${CMAKE_SOURCE_DIR}/include)

# add executable
add_executable(line_example examples/line_fitting_example.cpp)
target_link_libraries(line_example ${OpenCV_LIBRARIES})

add_executable(quadratic_example examples/quadratic_fitting_example.cpp)
target_link_libraries(quadratic_example ${OpenCV_LIBRARIES})

add_executable(polynomial_example examples/polynomial_fitting_example.cpp)
target_link_libraries(polynomial_example ${OpenCV_LIBRARIES})

add_executable(lls_example examples/lls_fitting_example.cpp)
target_link_libraries(lls_example ${OpenCV_LIBRARIES})