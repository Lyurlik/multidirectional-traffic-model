#ifndef NSWEMODEL_H
#define NSWEMODEL_H

#include <QString>
#include <QDebug>
#include <QElapsedTimer>

#include "utils.h"
#include "urbannetwork.h"
#include "density.h"

#define CELL_NUMBER 60    // number of cells in the domain in one dimension, e.g. CELL_NUMBER 60 implies that there are 3600 cells in total
#define M_PI   3.14159265359    // constant pi
#define M_PI_2 1.57079632679    // constant pi/2

class NSWEModel
{
public:

    NSWEModel(UrbanNetwork* network);

    void processIntersections();    // recompute all intersection and road parameters in NSWE formulation
    void constructInterpolation(double strength = 0.05);   // approximate all parameters in every cell using Inverse Distance Weighting
    void update(double dt);              // update the density using Godunov scheme. This function contains the NSWE model.
    void displayDensity(QPainter& painter, QRect rectangle, bool paint_network);   // this function defines density visualization procedure
    Eigen::Vector4d interpolateDensity(double x, double y);
    double interpolateLength(double x, double y);
    double getSSIMDiff_mean_weighted(Density4 density);   // here the weighted SSIM idex is computed
    void setInitialDensity(Density4 density);    // needed to set the initial density in NSWE model equal to the reconstructed density

    UrbanNetwork* network;


    //-----data-per-intersection-------------------
    std::vector<Eigen::MatrixXd> P_in;
    std::vector<Eigen::MatrixXd> P_out;

    std::vector<Eigen::Matrix4d> alpha_intersec;
    std::vector<Eigen::Matrix4d> beta_intersec;

    std::vector<Eigen::Vector4d> sin_intersec;
    std::vector<Eigen::Vector4d> cos_intersec;
    std::vector<double> L_intersec;
    std::vector<bool> used_for_out_interpolation;

    std::vector<Eigen::Vector4d> vel_intersec;
    std::vector<Eigen::Vector4d> rho_max_intersec;

    double max_density_ever;
    double dx, dy;

    //-----data-per-grid-cell-------------------  (define for each cell: 4dim vectors for each space dimension or 4dim diagonal matrices)
    std::vector<std::vector<Eigen::Vector2d> > pos;

    std::vector<std::vector<Eigen::Matrix4d> > alpha;
    std::vector<std::vector<Eigen::Matrix4d> > beta;

    std::vector<std::vector<Eigen::Vector4d> > sin_grid;
    std::vector<std::vector<Eigen::Vector4d> > cos_grid;
    std::vector<std::vector<double> > L;

    std::vector<std::vector<Eigen::Vector4d> > vel;
    std::vector<std::vector<Eigen::Vector4d> > rho_max;

    std::vector<std::vector<Eigen::Vector4d> > inflow;
    std::vector<std::vector<Eigen::Vector4d> > outflow;

    std::vector<std::vector<Eigen::Vector4d> > Demand;
    std::vector<std::vector<Eigen::Vector4d> > Supply;
    std::vector<std::vector<Eigen::Vector4d> > rho;
    std::vector<std::vector<Eigen::Vector4d> > new_rho;

    Eigen::VectorXd res_SSIM;

    //-------------------------------------

    std::vector<std::vector<Eigen::Vector4d> > rho_ema;  // smoothed density
    double ema_parameter; // this visualization parameter can be tuned: the larger the parameter, the more oscillatory is the visualization
};

#endif // NSWEMODEL_H
