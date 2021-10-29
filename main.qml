import QtQuick 2.12
import QtQuick.Window 2.12

import cz.seapraha.MandelImageCombiner 1.0;
import cz.seapraha.MandelModel 1.0;

Window {
    width: 640
    height: 480
    visible: true
    title: qsTr("Mandelbrot")

    Row {
        id: toprow;
        anchors.left: parent.left;
        anchors.right: parent.right;
        anchors.top: parent.top;
        Text {
            id: labelCX;
            width: 50;
            text: "fluff";
        }
        Text {
            id: labelCY;
            width: 50;
            text: "bluff";
        }
        Text {
            id: labelTimes;
            width: 600;
            text: "truff";
        }
    }

    MouseArea {
        anchors.left: parent.left;
        anchors.right: parent.right;
        anchors.top: toprow.bottom;
        anchors.bottom: parent.bottom;
        hoverEnabled: true;
        property bool pressed: false;
        property int pressx;
        property int pressy;
        onPositionChanged:
        {
            if (pressed)
            {
                mandelModel.drag(mouse.x-pressx, mouse.y-pressy);
                pressx=mouse.x;
                pressy=mouse.y;
            };
            //labelCX.text=mandelModel.pixelXtoRE(mouse.x);
            //labelCY.text=mandelModel.pixelYtoIM(mouse.y);
            labelCX.text=mouse.x;
            labelCY.text=mouse.y;
        }
        onPressed: {
            if (mouse.button==Qt.LeftButton)
            {
                pressed=true;
                pressx=mouse.x;
                pressy=mouse.y;
            };
        }
        onReleased: {
            pressed=false;
        }

        onWheel: {
            if (wheel.angleDelta.y>0)
              mandelModel.zoom(wheel.x, wheel.y, +1);
            else if (wheel.angleDelta.y<0)
              mandelModel.zoom(wheel.x, wheel.y, -1);
            wheel.accepted=true;
        }
        onHeightChanged: {
            imageCombiner.resetBgImage(width, height, 0xff000000);
            mandelModel.setImageSize(width, height);
            mandelModel.writeToImage(imageCombiner.getBaseImage());
            imageCombiner.update();
        }
        onWidthChanged: {
            imageCombiner.resetBgImage(width, height, 0xff000000);
            mandelModel.setImageSize(width, height);
            mandelModel.writeToImage(imageCombiner.getBaseImage());
            imageCombiner.update();
        }

        MandelImageCombiner {
            id: imageCombiner;
            Component.onCompleted: {
                resetBgImage(width, height, 0xff000000);
                mandelModel.setImageSize(width, height);
                //setBaseImage(mandelModel.getAsImage());
            }
            anchors.fill: parent;
        }
    }


    MandelModel {
        id: mandelModel
    }

    Timer {
        interval: 100;
        repeat: true;
        running: true;
        onTriggered: {
            mandelModel.writeToImage(imageCombiner.getBaseImage());
            imageCombiner.update();
            labelTimes.text=mandelModel.getTimes();
        }
    }
}
