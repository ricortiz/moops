
INCLUDE(FindVTK)
IF(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
    SET(APP_LIBS ${APP_LIBS} QVTK vtkRendering vtkGraphics vtkIO vtkCommon)
ELSE(VTK_FOUND)
    MESSAGE(FATAL_ERROR "VTK library not found.")
ENDIF(VTK_FOUND)

SET(APP_LIBS ${APP_LIBS} gui ${QT_LIBRARIES})
QT4_WRAP_CPP(gui_HEADERS_MOC ../io/gui.hpp)
add_executable(heart_valveless_pump  heart_valveless_pump.cpp ${gui_HEADERS_MOC})

IF(USE_CUDA_FLUID_SOLVER)
    add_definitions(-DCUDA_FLUID_SOLVER)
    target_link_libraries(heart_valveless_pump ${APP_LIBS} cudaStokesSolver)
ELSE(USE_CUDA_FLUID_SOLVER)
    target_link_libraries(heart_valveless_pump ${APP_LIBS})
ENDIF(USE_CUDA_FLUID_SOLVER)
