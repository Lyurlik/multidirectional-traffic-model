#ifndef TRAFFICSYSTEM_H
#define TRAFFICSYSTEM_H


#include <QMutex>
#include <QSize>
#include <QThread>
#include <QWaitCondition>

#include "grenobledata.h"
#include "nswemodel.h"
#include "urbannetwork.h"

class TrafficSystem : public QThread
{
    Q_OBJECT
public:
    TrafficSystem();
    ~TrafficSystem() override;

    void setTimestep(double new_timestep);
    double getTime();

    void displayColorbar(QPainter& painter, QRect rectangle);

    NSWEModel NSWE;
    GrenobleData Grenoble;
    UrbanNetwork network;

public slots:
    void start_simulation();
    void init_simulation(double init_time = 0);
    void destroy_simulation();

protected:
    void run() override;

private:

    double time;
    double timestep;

    bool started;
    bool aborted;
};

#endif // TRAFFICSYSTEM_H
