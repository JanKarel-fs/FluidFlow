SET(NAME Sources)
SET(SOURCES initialisation.cpp inlet.cpp outlet.cpp setBoundaryConditions.cpp
setGrid.cpp settings.cpp slipWall.cpp noSlipWall.cpp symmetry.cpp timeStep.cpp
homogeneousNeumann.cpp inletJacobian.cpp outletJacobian.cpp
slipWallJacobian.cpp noSlipWallJacobian.cpp symmetryJacobian.cpp
homogeneousNeumannJacobian.cpp)

ADD_LIBRARY(${NAME} ${SOURCES})
TARGET_LINK_LIBRARIES(${NAME} Fvm)
