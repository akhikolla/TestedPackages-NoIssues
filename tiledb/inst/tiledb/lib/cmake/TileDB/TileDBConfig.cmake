#
# This file attempts to locate the TileDB library. If found, the following
# imported targets are created:
#   - TileDB::tiledb_shared
#   - TileDB::tiledb_static
# And the following variables are defined:
#   - TILEDB_FOUND
#   - TileDB_FOUND
#


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was Config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

include("${CMAKE_CURRENT_LIST_DIR}/TileDBTargets.cmake")
check_required_components("TileDB")

# Targets required for linking against tiledb_static:

# -- End required targets for static

# Define a convenience all-caps variable
if (NOT DEFINED TILEDB_FOUND)
  if (TARGET TileDB::tiledb_shared)
    set(TILEDB_FOUND TRUE)
  else()
    set(TILEDB_FOUND FALSE)
  endif()
endif()
