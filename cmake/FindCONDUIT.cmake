# - Try to find libconduit
# Once done this will define
#  CONDUIT_FOUND - System has libconduit
#  CONDUIT_INCLUDE_DIRS - The libconduit include directories
#  CONDUIT_LIBRARIES - The libraries needed to use libconduit

FIND_PATH(WITH_CONDUIT_PREFIX
    NAMES include/conduit/conduit.h
)

IF(ENABLE_MPI)
FIND_LIBRARY(CONDUIT_LIBRARIES
    NAMES conduit conduit_blueprint conduit_relay_mpi
    HINTS ${WITH_CONDUIT_PREFIX}/lib
)
ELSE(ENABLE_MPI)
FIND_LIBRARY(CONDUIT_LIBRARIES
    NAMES conduit conduit_blueprint conduit_relay
    HINTS ${WITH_CONDUIT_PREFIX}/lib
)
ENDIF(ENABLE_MPI)

FIND_PATH(CONDUIT_INCLUDE_DIRS
    NAMES conduit/conduit.h
    HINTS ${WITH_CONDUIT_PREFIX}/include
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CONDUIT DEFAULT_MSG
    CONDUIT_LIBRARIES
    CONDUIT_INCLUDE_DIRS
)

# Hide these vars from ccmake GUI
MARK_AS_ADVANCED(
	CONDUIT_LIBRARIES
	CONDUIT_INCLUDE_DIRS
)
