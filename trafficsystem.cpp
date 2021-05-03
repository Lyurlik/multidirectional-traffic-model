#include "trafficsystem.h"
#include "utils.h"

TrafficSystem::TrafficSystem() :
    NSWE(&network),
    Grenoble(&network, &NSWE),
    started(false),
    aborted(false)
{

}

TrafficSystem::~TrafficSystem()
{
    destroy_simulation();
}


void TrafficSystem::setTimestep(double new_timestep)
{
    timestep = new_timestep;
}


void TrafficSystem::start_simulation()
{
    if (!started) {
        if (!isRunning()) {
            start(HighestPriority);
        }
        started = true;
    }
}


void TrafficSystem::init_simulation(double init_time)
{
    time = init_time;
}

double TrafficSystem::getTime()
{
    return time;
}

void TrafficSystem::destroy_simulation()
{
    qDebug() << "INFO Traffic System: STOPPING.";
    started = false;
    aborted = true;
    wait(1000);
    exit(0);
}


// What must be done in thread
void TrafficSystem::run()
{
    while(true) {

        while(started) {
            Grenoble.setFlow(NSWE.inflow, NSWE.outflow, time);   // for time-dependent inflows
            NSWE.update(timestep);
            time += timestep;
        }
        if (aborted) return;

    }
}


void TrafficSystem::displayColorbar(QPainter& painter, QRect rectangle)  // colorbar is the bar which varies from green for zero density to red for congestion
{
    QImage image(rectangle.size(), QImage::Format_RGB32);

    for (int j = 0; j < rectangle.height(); j++) {
        QColor color;
        double val = 1 - 1.0 * j / rectangle.height();

        color = getColorGreenRed(val);

        for (int i = 0; i < rectangle.width(); i++) {
            image.setPixelColor(i, j, color);
        }
    }

    painter.drawImage(rectangle, image);
}



