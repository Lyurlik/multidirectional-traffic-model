#ifndef DENSITY_H
#define DENSITY_H

#include <QDebug>
#include <QColor>
#include "Eigen/Eigen"
#include "urbannetwork.h"

#define CELL_NUMBER 60    // number of cells in the domain in one dimension, e.g. CELL_NUMBER 60 implies that there are 3600 cells in total

class Density
{
public:
    Density(int grid_size);

    void gaussian_blur(double standart_deviation, double network_scale = 1.0); // razmytie Gaussian Kernel (prakticheski odno i to zhe)
    void box_blur(int blur_size);

    void clear();
    void addValue(double x, double y, double value);
    void setDomain(UrbanNetwork* domain_network);

    double getNormalizedValue(double x, double y);

    UrbanNetwork* domain;

    int grid_size;
    Eigen::MatrixXd data;
};


struct Density4
{
    Density N, S, W, E;
    Density4(int grid_size) : N(grid_size), S(grid_size), W(grid_size), E(grid_size) {}

    void setDomain(UrbanNetwork* domain_network)
    {
        N.setDomain(domain_network);
        S.setDomain(domain_network);
        W.setDomain(domain_network);
        E.setDomain(domain_network);
    }
    void clear()
    {
        N.clear(); S.clear(); W.clear(); E.clear();
    }
    void gaussian_blur(double standart_deviation, double network_scale = 1.0)
    {
        N.gaussian_blur(standart_deviation, network_scale);
        S.gaussian_blur(standart_deviation, network_scale);
        W.gaussian_blur(standart_deviation, network_scale);
        E.gaussian_blur(standart_deviation, network_scale);
    }
};


//-------------------------------------------------
//CLASS DESCRIPTION ENDED
//-------------------------------------------------
//FUNCTIONS DEFINITIONS STARTED
//-------------------------------------------------




#endif // DENSITY_H
