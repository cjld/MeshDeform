 TEMPLATE = app
 CONFIG += console

# TEMPLATE = lib
# CONFIG += static
# DEFINES += DBUILD_LIB

CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    Align.cpp\
    ARAPDeform.cpp\
    FeatureVector.cpp\
    LQSolver.cpp\
    SR_ARAP.cpp \
    FVAnalysis.cpp \
    DeformCaller.cpp \
    DeformInfo.cpp \
    deformsf.cpp

# matlab
LIBS += -LE:\matlab\R2013a\extern\lib\win64\microsoft \
    -llibmx -llibmex -llibmat -llibeng

INCLUDEPATH += E:\matlab\R2013a\extern\include


HEADERS += \
    DMEngine.h \
    Align.h\
    ARAPDeform.h\
    FeatureVector.h\
    LQSolver.h \
    FVAnalysis.h \
    DeformCaller.h \
    DeformInfo.h \
    deformsf.h

#eigen
INCLUDEPATH += E:\MyLib\eigen\eigen-eigen-6b38706d90a9

#openmesh
INCLUDEPATH += 'E:/MyLib/OpenMesh/OpenMesh 3.1/include'

LIBS += -L'E:/MyLib/OpenMesh/OpenMesh 3.1/lib'
Release:LIBS += -lOpenMeshCore -lOpenMeshTools
Debug:LIBS += -lOpenMeshCored -lOpenMeshToolsd

Debug:DEFINES += WRITE_ITER_MESH
Release:DEFINES += DRELEASE

DEFINES += _USE_MATH_DEFINES

# DEFINES += DUSE_OPENMP

QMAKE_CXXFLAGS += /openmp
