#include "grenobledata.h"

GrenobleData::GrenobleData(UrbanNetwork* network, NSWEModel* nswe) :
    network(network),
    nswe(nswe),
    ema_density(CELL_NUMBER),
    ema_parameter(0),
    ema_initialized(false),
    density_reconstruction_width(-1),
    maximal_density_ever(1)
{

}


// -------------------------------------------- load data estimated from real sensors in Grenoble ------------------------------------------
void GrenobleData::loadData(QString filename_time, QString filename_density, QString filename_inflow, QString filename_outflow)
{
    QFile* input_data = new QFile(filename_time);
    if (!input_data->open(QFile::ReadOnly | QFile::Text)) {
        qDebug() << "ERROR GrenobleData: Could not open time file!!!";
    }

    QTextStream* input_data_stream = new QTextStream(input_data);

    long long first_timestamp = -1;
    while(!input_data_stream->atEnd()) {
        long long timestamp;
        (*input_data_stream) >> timestamp;    // time in seconds in unix time for which corresponding data will be read (here each minute)
        if (first_timestamp < 0)
            first_timestamp = timestamp;

        scenario_timepoints.push_back(timestamp - first_timestamp);
    }
    input_data->close();
    delete input_data_stream;
    delete input_data;

    max_time = scenario_timepoints.back();

    //--------------------------

    scenario_densities.resize(scenario_timepoints.size());
    scenario_inflows.resize(scenario_timepoints.size());
    scenario_outflows.resize(scenario_timepoints.size());

    for (int i = 0; i < scenario_timepoints.size(); i++) {
        scenario_densities[i].resize(network->edge_IDs.size());
        scenario_inflows[i].resize(network->edge_IDs.size());
        scenario_outflows[i].resize(network->edge_IDs.size());
    }

    //--------------------------

    input_data = new QFile(filename_density);
    if (!input_data->open(QFile::ReadOnly | QFile::Text)) {
        qDebug() << "ERROR GrenobleData: Could not open density file!!!";
    }

    input_data_stream = new QTextStream(input_data);

    while(!input_data_stream->atEnd()) {
        QString line = input_data_stream->readLine();
        line.replace(',',' ');
        QTextStream ss(&line);
        int road_id;
        ss >> road_id;

        int index = -1;
        for (int i = 0; i < network->edge_IDs.size(); i++) {
            if (network->edge_IDs[i] == road_id) {
                index = i;
                break;
            }
        }
        if (index == -1)
            continue;

        for (int i = 0; i < scenario_timepoints.size(); i++) {
            float density;
            ss >> density;                                   // density is given for each data time instant (here each minute) per road
            scenario_densities[i][index] = density;
        }
    }
    input_data->close();
    delete input_data_stream;
    delete input_data;

    //--------------------------

    input_data = new QFile(filename_inflow);
    if (!input_data->open(QFile::ReadOnly | QFile::Text)) {
        qDebug() << "ERROR GrenobleData: Could not open inflow file!!!";
    }

    input_data_stream = new QTextStream(input_data);

    while(!input_data_stream->atEnd()) {
        QString line = input_data_stream->readLine();
        line.replace(',',' ');
        QTextStream ss(&line);
        int road_id;
        ss >> road_id;

        int index = -1;
        for (int i = 0; i < network->edge_IDs.size(); i++) {
            if (network->edge_IDs[i] == road_id) {
                index = i;
                break;
            }
        }
        if (index == -1)
            continue;

        if (!network->node_on_border[network->edge_list[index].first])      // inflow can be set only for the node_on_border
            continue;

        for (int i = 0; i < scenario_timepoints.size(); i++) {
            float inflow;
            ss >> inflow;  // inflow is read for each time instant as veh/hour
            scenario_inflows[i][index] = inflow / 3600;    // we divide by 3600 to convert veh/hour into veh/sec
        }
    }

    //--------------------------

    input_data = new QFile(filename_outflow);
    if (!input_data->open(QFile::ReadOnly | QFile::Text)) {
        qDebug() << "ERROR GrenobleData: Could not open outflow file!!!";
    }

    input_data_stream = new QTextStream(input_data);

    while(!input_data_stream->atEnd()) {
        QString line = input_data_stream->readLine();
        line.replace(',',' ');
        QTextStream ss(&line);
        int road_id;
        ss >> road_id;

        int index = -1;
        for (int i = 0; i < network->edge_IDs.size(); i++) {
            if (network->edge_IDs[i] == road_id) {
                index = i;
                break;
            }
        }
        if (index == -1)
            continue;

        if (!network->node_on_border[network->edge_list[index].second])       // outflow can be set only for the node_on_border
            continue;

        for (int i = 0; i < scenario_timepoints.size(); i++) {
            float outflow;
            ss >> outflow;     // outflow is read for each time instant as veh/hour
            scenario_outflows[i][index] = outflow / 3600;  // we divide by 3600 to convert veh/hour into veh/sec
        }
    }
    input_data->close();
    delete input_data_stream;
    delete input_data;

    qDebug() << "INFO GrenobleData: data are loaded.";
}

void GrenobleData::setDensityReconstructionParameter(double d0)
{
    density_reconstruction_width = d0;
}
double GrenobleData::getDensityReconstructionParameter()
{
    return density_reconstruction_width;
}

void GrenobleData::setEMAtimewindow(double timewindow)
{
    ema_parameter = 2.0 / (1 + timewindow);
}

double GrenobleData::getMaximalTime()
{
    return max_time;
}


// ------------------------ real density given for every road is reconstructed for every cell in the domain -------------------------

Density4 GrenobleData::reconstructDensity(double t)
{
    unsigned int scenario_index = round(t / 60);

    Density4 estimated_density(CELL_NUMBER);
    estimated_density.setDomain(network);


    for (unsigned int i = 0; i < network->edge_IDs.size(); i++) {

        double density = scenario_densities[scenario_index][i];

        Eigen::Vector2d start = network->nodes[network->edge_list[i].first];
        Eigen::Vector2d end = network->nodes[network->edge_list[i].second];


        unsigned int index_of_outcoming_road = 0;
        for (unsigned int j = 0; j < network->outcoming_edges[network->edge_list[i].first].size(); j++)
        {
            if (network->outcoming_edges[network->edge_list[i].first][j] == i)
            {
                index_of_outcoming_road = j;
                break;
            }
        }
        Eigen::Vector4d p = nswe->P_out[network->edge_list[i].first].col(index_of_outcoming_road);

        int N = 10;    // separate each road in 10 parts and set densities on the boundary between this parts. This is needed for the Gaussian Kernel estimation for 2D
        density = density * (end - start).norm() / 1000.0 / N;
        for (int j = 1; j <= N; j++) {
            Eigen::Vector2d car_position = start + (end - start) * j / (N + 1);
            double est_L = nswe->interpolateLength(car_position.x(), car_position.y());

            // take matrix P for each road and define density for each direction
            estimated_density.N.addValue(car_position.x(), car_position.y(), density * p(0) / est_L);
            estimated_density.S.addValue(car_position.x(), car_position.y(), density * p(1) / est_L);
            estimated_density.W.addValue(car_position.x(), car_position.y(), density * p(2) / est_L);
            estimated_density.E.addValue(car_position.x(), car_position.y(), density * p(3) / est_L);
        }

    }
    //------Gaussian Kernel estimation: as we placed virtual cars with specified densities at some positions,
    //------we need to smooth the result such that each particular position has an influence on the cirle with radius specified by density_reconstruction_width
    estimated_density.gaussian_blur(density_reconstruction_width);

    if (ema_parameter > 1e-8) {
        if (!ema_initialized) {
            ema_density.setDomain(network);
            ema_density.N.data = estimated_density.N.data;
            ema_density.S.data = estimated_density.S.data;
            ema_density.W.data = estimated_density.W.data;
            ema_density.E.data = estimated_density.E.data;
        }
        else {
            ema_density.N.data = (1 - ema_parameter) * ema_density.N.data + ema_parameter * estimated_density.N.data;
            ema_density.S.data = (1 - ema_parameter) * ema_density.S.data + ema_parameter * estimated_density.S.data;
            ema_density.W.data = (1 - ema_parameter) * ema_density.W.data + ema_parameter * estimated_density.W.data;
            ema_density.E.data = (1 - ema_parameter) * ema_density.E.data + ema_parameter * estimated_density.E.data;
        }
        return ema_density;
    }

    return estimated_density;
}


//---------------displaying the real data density in the rectangle for GUI visualization
void GrenobleData::displayDensity(double time, QPainter& painter, QRect rectangle, bool paint_network)
{
    QImage image(rectangle.size(), QImage::Format_RGB32);

    Density4 density = reconstructDensity(time);

    for (int i = 0; i < rectangle.width(); i++) {
        for (int j = 0; j < rectangle.height(); j++) {
            QColor color(0, 0, 0);
            double x = 1.0 * i / rectangle.width();
            double y = 1 - 1.0 * (j + 1) / rectangle.height();

            Eigen::Vector4d val;     
            val(0) = density.N.getNormalizedValue(x, y);
            val(1) = density.S.getNormalizedValue(x, y);
            val(2) = density.W.getNormalizedValue(x, y);
            val(3) = density.E.getNormalizedValue(x, y);

            double sum_val = val.sum() / nswe->max_density_ever;
            sum_val *= 40;
            if (sum_val > 1) sum_val = 1;
            color = getColorGreenRed(sum_val);
            image.setPixelColor(i, j, color);
        }
    }

    if (paint_network) {
        network->paint(image);
    }

    // -------------START: here we paint dashed lines that separate Grenoble in 9 zones (needed for SSIM computation) ------------
    int x_size = 3;
    int y_size = 3;

    QPainter image_painter(&image);
    image_painter.setPen(QPen(Qt::black, 3, Qt::DashLine));
    double x_rect_size = rectangle.width() / x_size;
    double y_rect_size = rectangle.height() / y_size;

    for (int i = 1; i < x_size; i++) {
        image_painter.drawLine(i * x_rect_size, 0, i * x_rect_size, rectangle.height());
    }
    for (int j = 1; j < y_size; j++) {
        image_painter.drawLine(0, j * y_rect_size, rectangle.width(), j * y_rect_size);
    }
    image_painter.end();

    // -------------END: here we paint dashed lines that separate Grenoble in 9 zones (needed for SSIM computation) ------------

    painter.drawImage(rectangle, image);
}


//--------------------here arrays inflow and outflow are filled with real data inflows and outflows specified for every road for a given timestep
void GrenobleData::setFlow(std::vector<std::vector<Eigen::Vector4d> >& inflow, std::vector<std::vector<Eigen::Vector4d> >& outflow, double time)
{
    for (int i = 0; i < CELL_NUMBER; i++) {
        for (int j = 0; j < CELL_NUMBER; j++) {
            inflow[i][j].setZero();
            outflow[i][j].setZero();
        }
    }

    if (time > max_time)
        return;

    unsigned int scenario_index = round(time / 60);  // the index of minute

    for (unsigned int i = 0; i < network->edge_IDs.size(); i++) {
        double inflow_road = scenario_inflows[scenario_index][i];  // [timestamp][road]
        double outflow_road = scenario_outflows[scenario_index][i];  // [timestamp][road]


        Eigen::Vector2d start = network->nodes[network->edge_list[i].first];
        Eigen::Vector2d end = network->nodes[network->edge_list[i].second];

        unsigned int index_of_outcoming_road = 0;
        for (unsigned int j = 0; j < network->outcoming_edges[network->edge_list[i].first].size(); j++)
        {
            if (network->outcoming_edges[network->edge_list[i].first][j] == i)
            {
                index_of_outcoming_road = j;
                break;
            }
        }
        Eigen::Vector4d p = nswe->P_out[network->edge_list[i].first].col(index_of_outcoming_road);


        int int_x_start = round((start.x() - network->x_min) * CELL_NUMBER / network->range_x - 0.5);
        int int_y_start = round((start.y() - network->y_min) * CELL_NUMBER / network->range_y - 0.5);
        int int_x_end = round((end.x() - network->x_min) * CELL_NUMBER / network->range_x - 0.5);
        int int_y_end = round((end.y() - network->y_min) * CELL_NUMBER / network->range_y - 0.5);

        if (int_x_start < 1) int_x_start = 1;
        if (int_y_start < 1) int_y_start = 1;
        if (int_x_start >= CELL_NUMBER - 2) int_x_start = CELL_NUMBER - 3;
        if (int_y_start >= CELL_NUMBER - 2) int_y_start = CELL_NUMBER - 3;

        if (int_x_end < 1) int_x_end = 1;
        if (int_y_end < 1) int_y_end = 1;
        if (int_x_end >= CELL_NUMBER - 2) int_x_end = CELL_NUMBER - 3;
        if (int_y_end >= CELL_NUMBER - 2) int_y_end = CELL_NUMBER - 3;

        // We distribute in- and outflow among 4 cells uniformly
        for (int s = 0; s < 4; s++) {
            inflow[int_x_start][int_y_start](s) += inflow_road * p(s) / 4;
            inflow[int_x_start][int_y_start + 1](s) += inflow_road * p(s) / 4;
            inflow[int_x_start + 1][int_y_start](s) += inflow_road * p(s) / 4;
            inflow[int_x_start + 1][int_y_start + 1](s) += inflow_road * p(s) / 4;

            outflow[int_x_end][int_y_end](s) += outflow_road * p(s) / 4;
            outflow[int_x_end][int_y_end + 1](s) += outflow_road * p(s) / 4;
            outflow[int_x_end + 1][int_y_end](s) += outflow_road * p(s) / 4;
            outflow[int_x_end + 1][int_y_end + 1](s) += outflow_road * p(s) / 4;
        }
    }
}



