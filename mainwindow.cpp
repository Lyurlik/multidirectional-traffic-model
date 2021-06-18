#include "mainwindow.h"
#include "ui_mainwindow.h"

double eta = 0.02;     // weighting parameter used for approximation of functions
double d0 = 70;       // Gaussian Kernel parameter. This valus specifies in which range does a car contribute to the density. d0 = 70 implies a range of 70 meters.

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    setGeometry(0, 0, 1580, 980); // size of visualization window

    // ------------- read network topology ------------ /////////////////////////////////
    traffic_system.network.loadIntersections("../ModelValidation/IntersectionTable.csv");
    traffic_system.network.loadRoads("../ModelValidation/RoadTable.csv");
    traffic_system.network.loadImportance("../ModelValidation/RoadFRC.csv");
    traffic_system.network.loadTurns("../ModelValidation/TurnTable.csv");

    // ------------- read real data & set Gaussian Kernel estimation parameter ------------- ///////////////////////
    traffic_system.Grenoble.loadData("../ModelValidation/Timestamp.csv", "../ModelValidation/Density.csv", "../ModelValidation/AllInflows.csv", "../ModelValidation/AllOutflows.csv");
    traffic_system.Grenoble.setDensityReconstructionParameter(d0);      // parameter d_0

    traffic_system.NSWE.processIntersections();   // find all the parameters like L, alpha, phi_max etc. in NSWE formulation for every intersection in a network
    traffic_system.NSWE.constructInterpolation(eta);   // contruct interpolation to define network and intersection parameters for all the cells (continuous plane)

    time = 3600 * 14;    // start simulation time, here it is 2pm (multiply by 3600 to have time in seconds
    traffic_system.init_simulation(time);
    traffic_system.setTimestep(0.1);     // time step for NSWE model
    traffic_system.NSWE.setInitialDensity(traffic_system.Grenoble.reconstructDensity(time));   // here we read data from Martin, divide roads by 10 and turn it into NSWE formulation
    traffic_system.Grenoble.setEMAtimewindow(100);   // use ema procedure to make the density more smooth in time. The larger the time window is, the smoother it looks

    // create timer for visualization
    visualization_timer = new QTimer(this);
    visualization_timer->setTimerType(Qt::PreciseTimer);
    connect(visualization_timer, SIGNAL(timeout()), this, SLOT(update()));    //connect visualization timer with update fct (concurrent threads)
    connect(this, SIGNAL(destroyed()), &traffic_system, SLOT(destroy_simulation()));
    start();

    // here call computation of the weighted SSIM index averaged over 9 zones (it is a unique value).
    std::cout << time << " " << traffic_system.NSWE.getSSIMDiff_mean_weighted(traffic_system.Grenoble.reconstructDensity(time)) << " " << traffic_system.NSWE.res_SSIM.transpose() << std::endl;

}

MainWindow::~MainWindow()
{
    traffic_system.destroy_simulation();
    delete ui;
}


void MainWindow::start()
{
    traffic_system.start_simulation(); // launches concurrent thread of traffic simulation
    visualization_timer->start(20);  // 20 milisec of real time (50 images per second)
}


// --------------- visual update (connected to timeout of visualization timer)
void MainWindow::update()
{
    time = traffic_system.getTime();     // how much time of simulation has passed
    QMainWindow::update();      // repaint whole window, update image (call paintEvent)
}


// ----------------------- here we specify how to draw the output of the code. These are usually two windows with density distributions:
 //  NSWE (left window) and real density (right window)
void MainWindow::paintEvent(QPaintEvent*)
{

    QImage main_window_image(geometry().size(), QImage::Format_RGB32);
    QRect window_rect(QPoint(0, 0), geometry().size());

    QPainter painter(&main_window_image);
    painter.fillRect(window_rect, QBrush(Qt::white));

    int width = geometry().width();
    int height = geometry().height();

    QString s = "beta = " + QString::number(eta) + ",   time = " + QString::number(time / 3600.0) + ", SSIM =" + QString::number(traffic_system.NSWE.getSSIMDiff_mean_weighted(traffic_system.Grenoble.reconstructDensity(time)));
    painter.setFont(QFont("Times", 12, QFont::DemiBold));
    painter.drawText(width / 2 - 30, 30, s);

    QRect left_rectangle(10, 50, (width - 110) / 2, height - 70);
    QRect right_rectangle((width - 110) / 2 + 30, 50, (width - 110) / 2, height - 70);

    traffic_system.NSWE.displayDensity(painter, left_rectangle, true);   // display NSWE result on the left
    double max_time = std::min(time, traffic_system.Grenoble.getMaximalTime());     // to prevent from crash by asking more data than it has
    traffic_system.Grenoble.displayDensity(max_time, painter, right_rectangle, true);  // display real data result on the right

    // compute weighted SSIM index to compare two density distributions every time step
    std::cout << time << " " << traffic_system.NSWE.getSSIMDiff_mean_weighted(traffic_system.Grenoble.reconstructDensity(time)) << " " << traffic_system.NSWE.res_SSIM.transpose() << std::endl;

    QRect colorbar_rectangle(width - 60, 50, 50, height - 70);
    traffic_system.displayColorbar(painter, colorbar_rectangle);
    QPainter painter_main_window(this);
    painter_main_window.drawImage(window_rect, main_window_image);

}
