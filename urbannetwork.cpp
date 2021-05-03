#include "urbannetwork.h"
#include "utils.h"

UrbanNetwork::UrbanNetwork()
{
}

// ----------------- read intersections file: ------------------------
          // for the network we read intersections' coordinates, id and whether it is on border (on border means that external inflow or outflow is possible)
void UrbanNetwork::loadIntersections(QString filename)
{
    QFile* input_data = new QFile(filename);
    if (!input_data->open(QFile::ReadOnly | QFile::Text)) {
        qDebug() << "ERROR Network: Could not open network file!!!";
    }
    QTextStream* input_data_stream = new QTextStream(input_data);

    node_IDs.clear();
    node_on_border.clear();
    nodes.clear();

    //---------------------------
    //Read file
    //---------------------------

    input_data_stream->readLine();

    while(!input_data_stream->atEnd()) {
        QString line = input_data_stream->readLine();
        line.replace(',',' ');
        QTextStream ss(&line);

        Eigen::Vector2d v;
        unsigned int is_on_border, ID;
        ss >> v.x() >> v.y() >> ID >> is_on_border;

        node_on_border.push_back(is_on_border == 1);
        nodes.push_back(v);
        node_IDs.push_back(ID);
    }


    //---------------------------
    //Set domain properties
    //---------------------------

    x_min = nodes[0].x();
    y_min = nodes[0].y();
    double x_max = x_min;
    double y_max = y_min;
    for (unsigned int i = 1; i < nodes.size(); i++) {
        if (nodes[i].x() < x_min)
            x_min = nodes[i].x();
        if (nodes[i].y() < y_min)
            y_min = nodes[i].y();
        if (nodes[i].x() > x_max)
            x_max = nodes[i].x();
        if (nodes[i].y() > y_max)
            y_max = nodes[i].y();
    }
    range_x = x_max - x_min;
    range_y = y_max - y_min;

    qDebug() << range_x << range_y;
    qDebug() << "INFO Network: intersections file is loaded.";
}



// ----------------------------- read roads file --------------------------
void UrbanNetwork::loadRoads(QString filename, double distance_between_cars)
{
    QFile* input_data = new QFile(filename);
    if (!input_data->open(QFile::ReadOnly | QFile::Text)) {
        qDebug() << "ERROR Network: Could not open network file!!!";
    }
    QTextStream* input_data_stream = new QTextStream(input_data);

    edge_IDs.clear();
    edge_list.clear();
    edge_max_velocity.clear();
    min_density = 1e9;
    min_velocity = 1e9;

    //---------------------------
    //Read file
    //---------------------------

    input_data_stream->readLine();
    while(!input_data_stream->atEnd()) {
        QString line = input_data_stream->readLine();
        line.replace(',',' ');
        QTextStream ss(&line);

        unsigned int ID1, ID2, ID_road;
        double max_vel, length, x_data, y_data;
        int lanes;

        ss >> x_data >> y_data >> ID1 >> ID2 >> ID_road >> max_vel >> lanes >> length;  // id1 and id2 are intersections that road of id ID_road is connecting

        x_data *= 1000;
        y_data *= 1000;
        max_vel /= 3.6;

        unsigned int index_1 = 0;
        unsigned int index_2 = 0;
        for (unsigned int i = 0; i < node_IDs.size(); i++) {
            if (ID1 == node_IDs[i])
                index_1 = i;
            if (ID2 == node_IDs[i])
                index_2 = i;
        }

        edge_IDs.push_back(ID_road);
        edge_list.push_back(std::make_pair(index_1, index_2));
        edge_max_velocity.push_back(max_vel);
        edge_max_density.push_back(lanes / distance_between_cars);
        edge_max_flow.push_back(max_vel * lanes / distance_between_cars / 3);

        if (min_density > lanes / distance_between_cars) {
            min_density = lanes / distance_between_cars;
        }
        if (min_velocity > max_vel) {
            min_velocity = max_vel;
        }
    }

    //---------------------------
    //Build edges adjacency and default turning ratios lists
    //---------------------------
    //Default turning ratios: all outgoing roads have equal probability
    //---------------------------

    incoming_edges.clear();
    outcoming_edges.clear();
    incoming_edges.resize(nodes.size());
    outcoming_edges.resize(nodes.size());

    for (unsigned int i = 0; i < edge_list.size(); i++) {
        incoming_edges[edge_list[i].second].push_back(i);
        outcoming_edges[edge_list[i].first].push_back(i);
    }

    QFile* output_data_1 = new QFile("inflows_empty.txt");
    if (!output_data_1->open(QFile::WriteOnly | QFile::Text)) {
        qDebug() << "ERROR Network: Could not open network file!!!";
    }
    QTextStream* output_data_stream_1 = new QTextStream(output_data_1);

    for (unsigned int k = 0; k < nodes.size(); k++) {
        if (!node_on_border[k])
            continue;

        for (unsigned int j = 0; j < outcoming_edges[k].size(); j++) {
            *output_data_stream_1 << edge_IDs[outcoming_edges[k][j]] << ",0\n";
        }
    }
    output_data_1->close();

    qDebug() << "INFO Network: roads file is loaded.";
}


void UrbanNetwork::loadTurns(QString filename)
{
    QFile* input_data = new QFile(filename);
    if (!input_data->open(QFile::ReadOnly | QFile::Text)) {
        qDebug() << "ERROR Network: Could not open turns file!!!";
    }
    QTextStream* input_data_stream = new QTextStream(input_data);

    turn_table.clear();
    turn_table.resize(edge_list.size());

    //---------------------------
    //Read file
    //---------------------------

    turning_ratios.resize(edge_list.size(), edge_list.size());
    turning_ratios.setZero();

    input_data_stream->readLine();

    while(!input_data_stream->atEnd()) {
        QString line = input_data_stream->readLine();
        line.replace(',',' ');
        QTextStream ss(&line);

        unsigned int id1, id2, id_turn, id_node;
        double turn;
        ss >> id1 >> id2 >> id_turn >> id_node >> turn;     // id1 is incoming road, id2 is outgoing road and turn is the corresponding turning ratio between roads

        for (unsigned int i = 0; i < edge_IDs.size(); i++) {
            if (edge_IDs[i] != id1)
                continue;

            unsigned int k = edge_list[i].second;
            for (unsigned int j = 0; j < outcoming_edges[k].size(); j++) {
                if (edge_IDs[outcoming_edges[k][j]] != id2)
                    continue;

                turning_ratios(i, outcoming_edges[k][j]) = turn;   // turning ratio alpha

                break;
            }

            break;
        }
    }


    for (unsigned int k = 0; k < nodes.size(); k++) {

        //-----beta (supply coefficient)------
        for (unsigned int j = 0; j < outcoming_edges[k].size(); j++) {
            double sum_inflows = 0;
            for (unsigned int i = 0; i < incoming_edges[k].size(); i++) {
                sum_inflows += turning_ratios(incoming_edges[k][i], outcoming_edges[k][j]) * edge_max_flow[incoming_edges[k][i]];
            }
            if (sum_inflows < 1e-8)
                sum_inflows = 1;
            for (unsigned int i = 0; i < incoming_edges[k].size(); i++) {
                turning_ratios(outcoming_edges[k][j], incoming_edges[k][i]) = turning_ratios(incoming_edges[k][i], outcoming_edges[k][j])
                        * edge_max_flow[incoming_edges[k][i]] / sum_inflows;
            }
        }
    }
}



void UrbanNetwork::paint(QImage& image)
{
    QPainter painter(&image);

    QBrush brush(QColor(180, 180, 180, 200));
    painter.setBrush(brush);

    for (unsigned int i = 0; i < edge_list.size(); i++) {
        Eigen::Vector2d v1 = nodes[edge_list[i].first];
        Eigen::Vector2d v2 = nodes[edge_list[i].second];

        v1.x() = (v1.x() - x_min) * image.width() / range_x;
        v1.y() = image.height() - 1 - (v1.y() - y_min) * image.height() / range_y;
        v2.x() = (v2.x() - x_min) * image.width() / range_x;
        v2.y() = image.height() - 1 - (v2.y() - y_min) * image.height() / range_y;

        painter.setPen(QPen(brush, 8));
        painter.drawLine(v1.x(), v1.y(), v2.x(), v2.y());

        Eigen::Vector2d direction = (v2 - v1).normalized();
        Eigen::Vector2d mid_point_1 = v1 * 0.6 + v2 * 0.4;
        Eigen::Vector2d mid_point_2 = v1 * 0.4 + v2 * 0.6;

        double len = (mid_point_2 - mid_point_1).norm() / 2;
        Eigen::Vector2d mid_point_1_left = mid_point_1 + Eigen::Vector2d(-direction.y(), direction.x()) * len;
        Eigen::Vector2d mid_point_1_right = mid_point_1 - Eigen::Vector2d(-direction.y(), direction.x()) * len;

        QPolygonF polygon;
        polygon << QPointF(mid_point_2.x(), mid_point_2.y()) <<
                   QPoint(mid_point_1_left.x(), mid_point_1_left.y()) <<
                   QPoint(mid_point_1_right.x(), mid_point_1_right.y());
        painter.setPen(QPen(brush, 1));
        painter.drawPolygon(polygon);
    }
}

