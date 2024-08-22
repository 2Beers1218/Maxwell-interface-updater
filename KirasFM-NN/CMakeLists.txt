# Set the name of the project and target:
SET( TARGET "KirasFM--minimal-example" )

# link the local source files
FILE ( GLOB_RECURSE TARGET_SRC  "source/*.cc" )
INCLUDE_DIRECTORIES ( include )

# add the executable
SET ( TARGET_SRC ${TARGET_SRC} main.cc )


# set Minimum
CMAKE_MINIMUM_REQUIRED( VERSION 2.8.12 )

# add deal.II
SET ( DEAL_II_DIR /home/tobias/deal.ii-candi/deal.II-v9.2.0/ )
FIND_PACKAGE( deal.II 9.2.0 QUIET
              HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
            )

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()
