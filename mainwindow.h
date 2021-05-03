#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>
#include <QPainter>
#include <QImageWriter>

#include "trafficsystem.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow() override;

    void start();

    TrafficSystem traffic_system;

public slots:
    void update();
    void paintEvent(QPaintEvent*) override;

private:
    Ui::MainWindow *ui;

    QTimer* visualization_timer;
    double time;
};
#endif // MAINWINDOW_H
