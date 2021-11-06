QT += quick

CONFIG += c++11

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

#helps with templated version but not much with virtual version
QMAKE_CXXFLAGS_DEBUG += -O9
QMAKE_CFLAGS_DEBUG += -O9

SOURCES += \
        MandelEvaluator.cpp \
        MandelImageCombiner.cpp \
        MandelMath.cpp \
        MandelModel.cpp \
        ShareableImageWrapper.cpp \
        double_double.cpp \
        main.cpp \
        multiprec.cpp

RESOURCES += qml.qrc

# Additional import path used to resolve QML modules in Qt Creator's code model
QML_IMPORT_PATH =

# Additional import path used to resolve QML modules just for Qt Quick Designer
QML_DESIGNER_IMPORT_PATH =

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

HEADERS += \
  MandelEvaluator.hpp \
  MandelImageCombiner.hpp \
  MandelMath.hpp \
  MandelModel.hpp \
  ShareableImageWrapper.hpp \
  double_double.hpp \
  multiprec.hpp
