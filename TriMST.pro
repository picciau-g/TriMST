QT += core

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle

QMAKE_MAC_SDK = macosx10.11

QMAKE_CXXFLAGS_RELEASE += -fpermissive
QMAKE_CXXFLAGS_DEBUG += -fpermissive
QMAKE_CXXFLAGS += -std=c++0x
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3

SOURCES += main.cpp \
    LibTri/Vertex3D.cpp \
    LibTri/Vertex2D.cpp \
    LibTri/Triangle.cpp \
    LibTri/Timer.cpp \
    LibTri/Reader.cpp \
    LibTri/normals.cpp \
    LibTri/Edge.cpp \
    mstsegmenter.cpp

HEADERS += Vertex3D.h \
    LibTri/Vertex2D.h \
    LibTri/Triangle.h \
    LibTri/Timer.h \
    LibTri/Sorting.h \
    LibTri/Reader.h \
    LibTri/normals.h \
    LibTri/Mesh.h \
    LibTri/Edge.h \
    mstsegmenter.h

#include(deployment.pri)
#qtcAddDeployment()

