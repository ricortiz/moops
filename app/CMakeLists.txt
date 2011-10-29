QT4_WRAP_CPP(HEADERS_MOC ../gui/gui.hpp)
ADD_EXECUTABLE(heart_valveless_pump heart_valveless_pump.cpp ${HEADERS_MOC})
TARGET_LINK_LIBRARIES(heart_valveless_pump ${APP_LIBS} gui)

