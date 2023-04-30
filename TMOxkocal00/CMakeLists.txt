CMAKE_MINIMUM_REQUIRED(VERSION 2.4.8)
SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_INCLUDE_CURRENT_DIR TRUE)

IF(UNIX AND NOT LINUX_SET)
 ADD_DEFINITIONS(-D LINUX)
ENDIF(UNIX AND NOT LINUX_SET)

INCLUDE_DIRECTORIES(
../tmolib/
.
)

LINK_LIBRARIES (
tmo
fftw3
)

SET(TMO_SOURCES
 TMOxkocal00.cpp
 TMOPlugin.cpp
 laplace.cpp
)


SET(TMOWARD_HEADERS
 TMOxkocal00.h
 TMOPlugin.h
 laplca.h
)


add_library( TMOxkocal00 SHARED ${TMO_SOURCES} )
add_custom_command( TARGET TMOxkocal00 POST_BUILD 
	COMMAND cp -f libTMOxkocal00.so ../TMOCmd/TMOxkocal00.tml
	COMMAND cp -f libTMOxkocal00.so ../TMOgui/TMOxkocal00.tml 
	COMMAND cp -f libTMOxkocal00.so ../TMOxkocal00.tml 
#	COMMAND rm -f libTMOxkocal00.so 
	)


