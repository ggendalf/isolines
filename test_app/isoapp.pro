QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17
QMAKE_CXXFLAGS += -bigobj -fp:precise

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

INCLUDEPATH += "../boost_lib" "../isobuild"
LIBS += "-L../boost_lib/stage/lib/"

SOURCES += \
    main.cpp \
	mainwindow.cpp \
	../isobuild/isolines.cpp


HEADERS += \
	mainwindow.h \
	../isobuild/matrix.hpp \
	../isobuild/isolines.h \
	../isobuild/simple.h

FORMS += \
	mainwindow.ui

TRANSLATIONS += \
	isoapp_ru_RU.ts

#CONFIG += lrelease
#CONFIG += embed_translations

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
