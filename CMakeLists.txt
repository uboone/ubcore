# ======================================================================
#  ubcore main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/lardataobj/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_package(cetmodules REQUIRED)
project(ubcore VERSION 09.79.00 LANGUAGES CXX)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

find_package( lardata REQUIRED EXPORT )
find_package( ubobj REQUIRED EXPORT )

# macros for dictionary and simple_plugin
include(ArtMake)
include(BuildPlugins)
include(CetTest)

# ADD SOURCE CODE SUBDIRECTORIES HERE
add_subdirectory(ubcore)

# Unit tests.
add_subdirectory(test)

# ups - table and config files

# packaging utility
cet_cmake_config()
