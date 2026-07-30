cmake_minimum_required(VERSION 3.0.0)

project("PolyScope CGAL Feature")
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL")
#add_definitions(-DBOOST_NO_ZLIB)
# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()




#add_definitions(-DCGAL_ImageIO_USE_ZLIB)

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_STANDARD 20)
### Compiler options
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 ) # Emit a compile flags file to support completion engines 

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using Clang (linux or apple) or GCC
  message("Using clang/gcc compiler flags")
  SET(BASE_CXX_FLAGS "-std=c++20 -Wall -Wextra -g3")
  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces")
  SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message("Setting clang-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=5 -fcolor-diagnostics")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-fsanitize=address -fno-limit-debug-info")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
    message("Setting gcc-specific options")
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
  endif()


  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TRACE_INCLUDES}") # uncomment if you need to track down where something is getting included from
  SET(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -g3")
  SET(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
  SET(CMAKE_CXX_FLAGS_RELEASE        "-march=native -O3 -DNDEBUG")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message("Using Visual Studio compiler flags")
  set(BASE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /bigobj")
  set(BASE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") # parallel build
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")

  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions(-DNOMINMAX)
  add_definitions(-D_USE_MATH_DEFINES)

else()
  # unrecognized
  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()

if (MSVC)
  add_compile_options(/bigobj /MT)
else ()
  add_compile_options(-Wa,-mbig-obj)
endif ()
Set(Eigen  "C:/Users/mladeuil/Downloads/eigen-3.4.0/eigen-3.4.0/build")


# Find Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)



# Use Boost
#set(Boost_COMPILER "-vc141")
set(BOOST_ROOT "D:/library/boost_1_88_0")
set(Boost_INCLUDE_DIR "D:/library/boost_1_88_0")
set(Boost_LIBRARY_DIR "D:/library/boost_1_88_0/stage/lib")
find_package(Boost REQUIRED COMPONENTS system iostreams filesystem)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})


set(ZLIB_ROOT "D:/library/zlib-1.2.3-bin")
set(ZLIB_LIBRARY "D:/library/zlib-1.2.3-bin/lib/zlib.lib")
set(ZLIB_INCLUDE_DIR "D:/library/zlib-1.2.3-bin/include")


Set(DGtal_DIR  "D:/library/DGtal/build")
# Use DGtal
find_package(DGtal REQUIRED)
include_directories(${DGTAL_INCLUDE_DIRS})
link_directories(${DGTAL_LIBRARY_DIRS})



#find_package(Spectra 1.0.1 REQUIRED NO_MODULE)

# Add polyscope
message("\n\n == CMAKE recursively building Polyscope\n")
add_subdirectory("deps/polyscope")


# Create an executable
add_executable(
        CGAL-polyscope-example
        src/main.cpp
        )


Set(CGAL_DIR "D:/library/CGAL-6.1-beta1")
find_package(CGAL REQUIRED COMPONENTS ImageIO)
find_package(ZLIB REQUIRED)
set(CGAL_USE_ZLIB TRUE)
target_compile_definitions(CGAL-polyscope-example PUBLIC  CGAL_ImageIO_USE_ZLIB)

if(TARGET CGAL::CGAL_ImageIO)
if(CGAL_ImageIO_USE_ZLIB)
  message("OK ZLIB")
else ()
  message("PAS OK ZLIB")
endif()
  message("OK CGAL_ImageIO")
else()
message("PAS OK CGAL_ImageIO")
endif()


if(CGAL_USE_ZLIB)
  message("OK USE_ZLIB")
else()
  message("PAS OK USE_ZLIB")
endif()
# Include settings
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/deps/")
target_include_directories(CGAL-polyscope-example PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/libigl/include")
#target_include_directories(CGAL-polyscope-example PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/CGALHelper/include")
add_subdirectory("deps/geometry-central") # wherever you put it


# Link settings
target_compile_features(CGAL-polyscope-example PRIVATE cxx_std_17)
target_link_libraries(CGAL-polyscope-example ${DGTAL_LIBRARIES} CGAL::CGAL  CGAL::CGAL_ImageIO polyscope geometry-central ZLIB::ZLIB Eigen3::Eigen Boost::filesystem)

#target_link_libraries(CGAL-polyscope-example polyscope Spectra::Spectra)
