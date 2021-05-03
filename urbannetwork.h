#ifndef URBANNETWORK_H
#define URBANNETWORK_H

#include <QString>
#include <QFile>
#include <QDebug>
#include <QPainter>
#include <vector>
#include <iostream>
#include "Eigen/Eigen"

class UrbanNetwork
{
public:
    UrbanNetwork();

    void loadIntersections(QString filename);
    void loadRoads(QString filename, double distance_between_cars = 6.0);
    void loadTurns(QString filename);

    void paint(QImage& image);


    double x_min;
    double y_min;
    double range_x;
    double range_y;

    std::vector<unsigned int> node_IDs;
    std::vector<bool> node_on_border;
    std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d> > nodes;

    std::vector<unsigned int> edge_IDs;
    std::vector<std::pair<unsigned int, unsigned int> > edge_list;
    std::vector<double> edge_max_velocity;
    std::vector<double> edge_max_density;
    std::vector<double> edge_max_flow;

    std::vector< std::vector<unsigned int> > incoming_edges;
    std::vector< std::vector<unsigned int> > outcoming_edges;

    std::vector< std::vector<unsigned int> > turn_table;

    Eigen::MatrixXd turning_ratios;
    double min_density;
    double min_velocity;
};

#endif // URBANNETWORK_H
