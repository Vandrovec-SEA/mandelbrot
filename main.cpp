#include <QGuiApplication>
#include <QQmlApplicationEngine>

#include "MandelImageCombiner.hpp"
#include "MandelModel.hpp"
#include "LaguerreModel.hpp"
#include "ShareableImageWrapper.hpp"

int main(int argc, char *argv[])
{
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
  QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
#endif

  QGuiApplication app(argc, argv);

  qmlRegisterType<MandelImageCombiner>("cz.seapraha.MandelImageCombiner", 1, 0, "MandelImageCombiner");
  qmlRegisterType<MandelModel>("cz.seapraha.MandelModel", 1, 0, "MandelModel");
  qmlRegisterType<LaguerreModel>("cz.seapraha.LaguerreModel", 1, 0, "LaguerreModel");
  qRegisterMetaType<ShareableImageWrapper>();
  qRegisterMetaType<ShareableViewInfo>();

  QQmlApplicationEngine engine;
  const QUrl url(QStringLiteral("qrc:/main.qml"));
  QObject::connect(&engine, &QQmlApplicationEngine::objectCreated,
                   &app, [url](QObject *obj, const QUrl &objUrl) {
    if (!obj && url == objUrl)
      QCoreApplication::exit(-1);
  }, Qt::QueuedConnection);
  engine.load(url);

  return app.exec();
}
