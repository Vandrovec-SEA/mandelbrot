import QtQuick 2.12
import QtQuick.Controls 2.12
import QtQuick.Window 2.12

import cz.seapraha.MandelImageCombiner 1.0;
import cz.seapraha.MandelModel 1.0;
import cz.seapraha.LaguerreModel 1.0;

Window {
    width: 640
    height: 480
    visible: true
    title: qsTr("Mandelbrot")
    //Keys.forwardTo: mouseArea //can't attach Keys to Window

    Column {
        anchors.left: parent.left;
        //anchors.right: parent.right;
        anchors.top: parent.top;
        id: toprow;
        //width: 150
        Text {
            id: labelXY;
            text: "fluff";
        }
        Text {
            id: labelInfoGen;
            text: "bluff";
        }
        Text {
            id: labelInfoSpec;
            text: "truff";
        }
    }

    ButtonGroup { id: bgroupView }
    Column {
        anchors.right: comboColumn.left
        anchors.top: parent.top;
        id: rgroupView
        RadioButton {
            id: rbuttonViewMand
            checked: true
            text: "Mand"
            height: 13
            //implicitIndicatorHeight: 13
            ButtonGroup.group: bgroupView
        }
        RadioButton {
            id: rbuttonViewLagu
            text: "Lagu"
            height: 13
            //implicitIndicatorWidth: 13
            ButtonGroup.group: bgroupView
        }
    }

    Column {
        id: comboColumn
        anchors.top: parent.top;
        anchors.right: parent.right
        Row {
            anchors.right: parent.right
            Rectangle {
                id: busyIndicator
                width: 20
                height: 20
                radius: 10
                color: "blue"
            }
            ComboBox {
                textRole: "text"
                valueRole: "key"
                currentIndex: 1 //cls
                model: ListModel {
                    id: paintstyle_model
                    ListElement { text: "Kind"; key: 0 }//mandelModel.paintStyleKind }
                    ListElement { text: "Cls"; key: 1 }//mandelModel.paintStyleCls }
                    ListElement { text: "Exter"; key: 2 }//mandelModel.paintStyleExter  }
                    ListElement { text: "Inter"; key: 3 }//mandelModel.paintStyleInter  }
                    ListElement { text: "Near"; key: 4 }//mandelModel.paintStyleNear  }
                    ListElement { text: "d/dz"; key: 5 }//mandelModel.paintStyleFZ  }
                    ListElement { text: "d/dc"; key: 6 }//mandelModel.paintStyleFC  }
                }
                onActivated: mandelModel.selectedPaintStyle=paintstyle_model.get(currentIndex).key;
            }
        }
        Row {
            anchors.right: parent.right
            Text {
                id: busyIndicator2
                text: "0"
            }
            ComboBox {
                textRole: "text"
                valueRole: "key"
                currentIndex: 0 //double
                model: ListModel {
                    id: precision_model
                    ListElement { text: "Double"; key: 0 }
                    ListElement { text: "DDouble"; key: 1 }
                    //ListElement { text: "QDouble"; key: 2 }
                    ListElement { text: "Multi"; key: 3 }
                }
                onActivated:
                {
                    mandelModel.selectedPrecision=precision_model.get(currentIndex).key;
                    laguerreModel.selectedPrecision=precision_model.get(currentIndex).key;
                }
            }
        }
    }

    MouseArea {
        id: mouseArea
        anchors.left: parent.left;
        anchors.right: parent.right;
        anchors.top: toprow.bottom;
        anchors.bottom: parent.bottom;
        hoverEnabled: true;
        acceptedButtons: Qt.LeftButton | Qt.RightButton
        property int dragging: 0;
        property int drag_last_x;
        property int drag_last_y;
        Keys.onPressed: { //never triggers
            console.log(event);
            if (event.key===Qt.Key_E)
                mandelModel.paintOrbit(mandelImageCombiner.getOverlayImage(), mouse.x, mouse.y);
        }
        onPositionChanged:
        {
            if (dragging==1)
            {
                mandelModel.drag(mouse.x-drag_last_x, mouse.y-drag_last_y);
                drag_last_x=mouse.x;
                drag_last_y=mouse.y;
            };
            if (dragging==2)
            {
                laguerreModel.drag(mouse.x-drag_last_x, mouse.y-drag_last_y);
                drag_last_x=mouse.x;
                drag_last_y=mouse.y;
            };
            //TODO: else copy from    ((mandel.mousePt.c.re<>view.orbit.re) or (mandel.mousePt.c.im<>view.orbit.im)) then
            //labelCX.text=mandelModel.pixelXtoRE_str(mouse.x);
            //labelCY.text=mandelModel.pixelYtoIM_str(mouse.y);
            if (rbuttonViewMand.checked)
              mandelModel.paintOrbit(mandelImageCombiner.getOverlayImage(), mouse.x, mouse.y);
            if (rbuttonViewLagu.checked)
            {
              laguerreModel.paintOrbit(laguerreImageCombiner.getOverlayImage(), mouse.x, mouse.y);
              if (laguerreModel==null)
                console.log("laguModel missing in action");
            }
            //labelCX.text=mouse.x;
            //labelCY.text=mouse.y;
            laguerreModel.getTimes();
        }
        onPressed: {
            if (mouse.button==Qt.LeftButton)
            {
                if (rbuttonViewMand.checked)
                  dragging=1;
                if (rbuttonViewLagu.checked)
                  dragging=2;
                drag_last_x=mouse.x;
                drag_last_y=mouse.y;
                if (rbuttonViewMand.checked)
                  mandelModel.paintOrbit(mandelImageCombiner.getOverlayImage(), mouse.x, mouse.y);
                if (rbuttonViewLagu.checked)
                  laguerreModel.paintOrbit(laguerreImageCombiner.getOverlayImage(), mouse.x, mouse.y);
            };
        }
        onReleased: {
            dragging=0;
            if (mouse.button==Qt.RightButton)
            {
                mainPopupMenu.popup(mouse.x, mouse.y);
            }
        }

        onWheel: {
            if (wheel.angleDelta.y>0)
            {
              if (rbuttonViewMand.checked)
                mandelModel.zoom(wheel.x, wheel.y, +1);
              else if (rbuttonViewLagu.checked)
                laguerreModel.zoom(wheel.x, wheel.y, +1);
            }
            else if (wheel.angleDelta.y<0)
            {
              if (rbuttonViewMand.checked)
                mandelModel.zoom(wheel.x, wheel.y, -1);
              else if (rbuttonViewLagu.checked)
                laguerreModel.zoom(wheel.x, wheel.y, -1);
            }
            wheel.accepted=true;
        }
        onHeightChanged: {
            mandelImageCombiner.resetBgImage(width, height, 0xff000000);
            mandelModel.setImageSize(width, height);
            mandelModel.writeToImage(mandelImageCombiner.getBaseImage());
            mandelImageCombiner.update();
            laguerreImageCombiner.resetBgImage(width, height, 0xff000000);
            laguerreModel.setImageSize(width, height);
            laguerreModel.writeToImage(laguerreImageCombiner.getBaseImage());
            laguerreImageCombiner.update();
        }
        onWidthChanged: {
            mandelImageCombiner.resetBgImage(width, height, 0xff000000);
            mandelModel.setImageSize(width, height);
            mandelModel.writeToImage(mandelImageCombiner.getBaseImage());
            mandelImageCombiner.update();
            laguerreImageCombiner.resetBgImage(width, height, 0xff000000);
            laguerreModel.setImageSize(width, height);
            laguerreModel.writeToImage(laguerreImageCombiner.getBaseImage());
            laguerreImageCombiner.update();
        }

        MandelImageCombiner {
            id: mandelImageCombiner;
            visible: rbuttonViewMand.checked
            Component.onCompleted: {
                resetBgImage(width, height, 0xff000000);
                mandelModel.setImageSize(width, height);
                //setBaseImage(mandelModel.getAsImage());
            }
            anchors.fill: parent;
        }

        MandelImageCombiner {
            id: laguerreImageCombiner;
            visible: rbuttonViewLagu.checked
            Component.onCompleted: {
                resetBgImage(width, height, 0xff000000);
                laguerreModel.setImageSize(width, height);
                //setBaseImage(mandelModel.getAsImage());
            }
            anchors.fill: parent;
        }

        Menu {
            id: mainPopupMenu
            /*MenuItem {
                text: "test1"
                onTriggered: console.log("clicked test1");
            }*/
            MenuItem {
                text: "set Lagu"
                onTriggered: laguerreModel.setParams(mandelModel.viewInfo);
            }
            Menu {
                title: "Presets..."
                Repeater {
                    model: [
                        {viewRe: -0.5, viewIm: 0, viewZoom: 1/128, caption: 'Main view'},
                        {viewRe: -0.75, viewIm: 0, viewZoom: 1/134217728, caption: 'Zoom on the 1-2 joint'},
                        {viewRe: -0.125, viewIm: 0.649519052836944556, viewZoom: 1/268435456, caption: 'Zoom on the 1-3 joint'},
                        {viewRe: -1.401155524, viewIm: 0, viewZoom: 1/8589934592, caption: 'A 1536 atom on real axis'},
                        {viewRe: -0.7458232657523781425, viewIm: 0.1659563046792491885, viewZoom: 1/1152921504606846976.0, caption: 'An extremely small 980 atom'},
                        {viewRe: -0.16400146484, viewIm: 1.0978469849, viewZoom: 1/262144, caption: 'Detail of ray 2/9'},
                        {viewRe: -0.54418945312, viewIm: 0.65454101562, viewZoom: 1/2048, caption: 'Detail of head of 2/5 bulb'},//  julia.c.re:=-0.53540039062,  julia.c.im:=-0.66320800781,//}
                        {viewRe: -0.706765832036462260364, viewIm: 0.3540325897413834077558, viewZoom: 1/1125899906842624, caption: 'Just some zoom of a 408 atom'},
                        {viewRe: -1.262303108533615159097935, viewIm: 0.383593549747416773254917, viewZoom: 1/70368744177664, caption: 'just a zoom'},
                        {viewRe: -1.67440967369022279335495, viewIm: 4.716562790887584498e-5, viewZoom: 1/137438953472, caption: 'Embedded Julia set'},
                        {viewRe: -0.7498572, viewIm: -0.0150818, viewZoom: 1/8388608, caption: 'Deep pools'},
                        {viewRe: 0.2860166666955662541362, viewIm: 0.0115374014484412084383, viewZoom: 1/2305843009213693952, caption: 'Spiral'},
                        {viewRe: 0.28601560169483100000, viewIm: 0.01153748599423500000, viewZoom: 1/17592186044416, caption: 'Dense 4-spiral'},
                        {viewRe: -1.39473140503329, viewIm: 0.00450421402574, viewZoom: 1/17179869184, caption: 'Star'},
                        {viewRe: -1.39473141148726511, viewIm: 0.00450422164786152379, viewZoom: 1/1125899906842624, caption: 'Circled mandel'},
                        {viewRe: -1.94098287909277016167, viewIm: 0.00064812392023092638912, viewZoom: 1/2199023255552, caption: 'Minijulia'},
                        {viewRe: -0.74364388703715904315, viewIm: 0.13182590420531251074, viewZoom: 1/8796093022208, caption: 'Wikipedia'}
                    ]
                    MenuItem {
                        text: modelData.caption
                        onTriggered: mandelModel.setView_double(modelData.viewRe, modelData.viewIm, modelData.viewZoom);
                    }
                }
            }
        }
    }



    MandelModel {
        id: mandelModel
    }

    LaguerreModel {
        id: laguerreModel
    }

    Timer {
        interval: 100;
        repeat: true;
        running: true;
        onTriggered: {
            if (rbuttonViewMand.checked)
            {
              labelInfoSpec.text=mandelModel.writeToImage(mandelImageCombiner.getBaseImage());
              mandelImageCombiner.update();
              labelXY.text=mandelModel.getTextXY();
              labelInfoGen.text=mandelModel.getTextInfoGen();
              //labelInfoSpec.text=mandelModel.getTextInfoSpec();
              busyIndicator.color=(mandelModel.threadsWorking===0?"green":mandelModel.threadsWorking===mandelModel.threadsMax?"red":"yellow")
              busyIndicator2.text=mandelModel.threadsWorking
            };
            if (rbuttonViewLagu.checked)
            {
              laguerreModel.getTimes();
              laguerreModel.writeToImage(laguerreImageCombiner.getBaseImage());
              laguerreImageCombiner.update();
              labelXY.text=laguerreModel.getTextXY();
              labelInfoGen.text=laguerreModel.getTextInfoGen();
              labelInfoSpec.text=laguerreModel.getTextInfoSpec();
            };
        }
    }
}
