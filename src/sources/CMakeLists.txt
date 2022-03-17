SET(NAME Sources)
SET(SOURCES initialisation.cpp inlet.cpp outlet.cpp setBoundaryConditions.cpp
setGrid.cpp settings.cpp slipWall.cpp timeStep.cpp)

ADD_LIBRARY(${NAME} ${SOURCES})
TARGET_LINK_LIBRARIES(${NAME} Fvm)
