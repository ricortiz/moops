FIND_PACKAGE(VTK REQUIRED)
IF(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
    SET(APP_LIBS ${APP_LIBS} vtkIO)
ELSE(VTK_FOUND)
    MESSAGE(FATAL_ERROR "VTK library not found.")
ENDIF(VTK_FOUND)

# QT4_WRAP_CPP(gui_HEADERS_MOC gui.hpp)

add_library(gui gui.cpp vtkSurfaceSource.cpp)

# add_executable(qt_gui main.cpp ${gui_HEADERS_MOC})
# target_link_libraries(qt_gui gui QVTK ${QT_LIBRARIES} vtkRendering
#   vtkGraphics
#   vtkIO
#   vtkCommon)

