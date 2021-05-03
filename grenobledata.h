#ifndef GRENOBLEDATA_H
#define GRENOBLEDATA_H

#include <vector>
#include <QString>
#include <QDebug>
#include <QElapsedTimer>
#include <QPainter>
#include "urbannetwork.h"
#include "density.h"
#include "nswemodel.h"

class GrenobleData
{
public:
    GrenobleData(UrbanNetwork* network, NSWEModel* nswe);

    void loadData(QString filename_time, QString filename_density, QString filename_inflow, QString filename_outflow);

    void setDensityReconstructionParameter(double d0);
    double getDensityReconstructionParameter();

    void setEMAtimewindow(double timewindow);

    double getMaximalTime();
    Density4 reconstructDensity(double t);

    void displayDensity(double time, QPainter& painter, QRect rectangle, bool paint_network);

    void setFlow(std::vector<std::vector<Eigen::Vector4d> >& inflow, std::vector<std::vector<Eigen::Vector4d> >& outflow, double time);

    //------------------------

    UrbanNetwork* network;
    NSWEModel* nswe;

    Density4 ema_density;
    double ema_parameter;
    bool ema_initialized;

    double max_time;
    std::vector<double> scenario_timepoints;
    std::vector<std::vector<double> > scenario_densities;
    std::vector<std::vector<double> > scenario_velocities;
    std::vector<std::vector<double> > scenario_inflows;
    std::vector<std::vector<double> > scenario_outflows;

    double density_reconstruction_width;
    double maximal_density_ever;
};

#endif // GRENOBLEDATA_H
