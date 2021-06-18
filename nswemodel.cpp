#include "nswemodel.h"

#define CRIT_TO_MAX 0.33333  // we set rho_c = 1/3 rho_max for the shape of triangular fundamental diagram.


// ------------------ calculate demand function as a function of density
Eigen::Vector4d demand(Eigen::Vector4d rho, Eigen::Vector4d rho_max, Eigen::Vector4d v_max)
{
    Eigen::Vector4d rho_crit = rho_max * CRIT_TO_MAX;
    Eigen::Vector4d result;
    for (int i = 0; i < 4; i++) {
        if (rho(i) < 0)
            result(i) = 0;
        else if (rho(i) < rho_crit(i))
            result(i) = v_max(i) * rho(i);
        else
            result(i) = v_max(i) * rho_crit(i);
    }
    return result;
}


// ---------------- calculate supply function as a function of density
Eigen::Vector4d supply(Eigen::Vector4d rho, Eigen::Vector4d rho_max, Eigen::Vector4d v_max)
{
    Eigen::Vector4d rho_crit = rho_max * CRIT_TO_MAX;

    Eigen::Vector4d result;
    for (int i = 0; i < 4; i++) {
        if (rho(i) < rho_crit(i))
            result(i) = v_max(i) * rho_crit(i);
        else if (rho(i) < rho_max(i)) {
            result(i) = v_max(i) * rho_crit(i) * (rho_max(i) - rho(i))  / (rho_max(i) - rho_crit(i));
        }
        else result(i) = 0;
    }
    return result;
}


NSWEModel::NSWEModel(UrbanNetwork* network) :
    network(network)
{
    ema_parameter = 0.0005;  // this parameter is to smooth the visualization
}


// ------------------------- calculate all network and intersection parameters in NSWE formulation for all intersections in the network -------------
void NSWEModel::processIntersections()
{
    P_in.resize(network->nodes.size());
    P_out.resize(network->nodes.size());
    alpha_intersec.resize(network->nodes.size());
    beta_intersec.resize(network->nodes.size());
    L_intersec.resize(network->nodes.size());
    sin_intersec.resize(network->nodes.size());
    cos_intersec.resize(network->nodes.size());
    vel_intersec.resize(network->nodes.size());
    rho_max_intersec.resize(network->nodes.size());
    used_for_out_interpolation.resize(network->nodes.size());

    // run over all intersections in a network
    for (unsigned int k = 0; k < network->nodes.size(); k++) {
        P_in[k].resize(4, network->incoming_edges[k].size());
        P_in[k].setZero();

         // process incoming roads of each intersection k (calculate Pin for incoming roads)
        for (unsigned int i = 0; i < network->incoming_edges[k].size(); i++) {
            Eigen::Vector2d dir = network->nodes[network->edge_list[network->incoming_edges[k][i]].second]
                            - network->nodes[network->edge_list[network->incoming_edges[k][i]].first];

            double ang = atan2(dir.y(), dir.x());  // theta in road formulation
            double cos_ang = cos(ang);
            double sin_ang = sin(ang);

            if (ang > M_PI_2) {     //NW
                P_in[k](0, i) = sin_ang / (sin_ang - cos_ang);
                P_in[k](2, i) = -cos_ang / (sin_ang - cos_ang);
            }
            else if (ang > 0) {     //NE
                P_in[k](0, i) = sin_ang / (sin_ang + cos_ang);
                P_in[k](3, i) = cos_ang / (sin_ang + cos_ang);
            }
            else if (ang > -M_PI_2) {     //SE
                P_in[k](1, i) = -sin_ang / (-sin_ang + cos_ang);
                P_in[k](3, i) = cos_ang / (-sin_ang + cos_ang);
            }
            else {     //SW
                P_in[k](1, i) = -sin_ang / (-sin_ang - cos_ang);
                P_in[k](2, i) = -cos_ang / (-sin_ang - cos_ang);
            }
        }


       // process outgoing roads of each intersection k (calculate Pout for outcoming roads)
        P_out[k].resize(4, network->outcoming_edges[k].size());
        P_out[k].setZero();

        for (unsigned int j = 0; j < network->outcoming_edges[k].size(); j++) {
            Eigen::Vector2d dir = network->nodes[network->edge_list[network->outcoming_edges[k][j]].second]
                            - network->nodes[network->edge_list[network->outcoming_edges[k][j]].first];
            double ang = atan2(dir.y(), dir.x());
            double cos_ang = cos(ang);
            double sin_ang = sin(ang);

            if (ang > M_PI_2) {     //NW
                P_out[k](0, j) = sin_ang / (sin_ang - cos_ang);
                P_out[k](2, j) = -cos_ang / (sin_ang - cos_ang);
            }
            else if (ang > 0) {     //NE
                P_out[k](0, j) = sin_ang / (sin_ang + cos_ang);
                P_out[k](3, j) = cos_ang / (sin_ang + cos_ang);
           }
            else if (ang > -M_PI_2) {     //SE
                P_out[k](1, j) = -sin_ang / (-sin_ang + cos_ang);
                P_out[k](3, j) = cos_ang / (-sin_ang + cos_ang);
            }
            else {     //SW
                P_out[k](1, j) = -sin_ang / (-sin_ang - cos_ang);
                P_out[k](2, j) = -cos_ang / (-sin_ang - cos_ang);
            }
        }


        // find cos(theta) and sin(theta) with line for each intersection
        cos_intersec[k].setZero();
        sin_intersec[k].setZero();
        Eigen::Vector4d max_flow = Eigen::Vector4d::Zero();

        for (unsigned int j = 0; j < network->outcoming_edges[k].size(); j++) {
             Eigen::Vector2d dir = network->nodes[network->edge_list[network->outcoming_edges[k][j]].second]
                                - network->nodes[network->edge_list[network->outcoming_edges[k][j]].first];
             double ang = atan2(dir.y(), dir.x());
             double cos_ang = cos(ang);  // is found for each road
             double sin_ang = sin(ang);
             for (int s1 = 0; s1 < 4; s1++) {
                  cos_intersec[k](s1) += P_out[k](s1, j) * cos_ang * network->edge_max_flow[network->outcoming_edges[k][j]];
                  sin_intersec[k](s1) += P_out[k](s1, j) * sin_ang * network->edge_max_flow[network->outcoming_edges[k][j]];
                  max_flow(s1) += P_out[k](s1, j) * network->edge_max_flow[network->outcoming_edges[k][j]];
            }
        }
        for (unsigned int i = 0; i < network->incoming_edges[k].size(); i++) {
             Eigen::Vector2d dir = network->nodes[network->edge_list[network->incoming_edges[k][i]].second]
                                - network->nodes[network->edge_list[network->incoming_edges[k][i]].first];
             double ang = atan2(dir.y(), dir.x());
             double cos_ang = cos(ang);
             double sin_ang = sin(ang);
             for (int s1 = 0; s1 < 4; s1++) {
                 cos_intersec[k](s1) += P_in[k](s1, i) * cos_ang * network->edge_max_flow[network->incoming_edges[k][i]];
                 sin_intersec[k](s1) += P_in[k](s1, i) * sin_ang * network->edge_max_flow[network->incoming_edges[k][i]];
                 max_flow(s1) += P_in[k](s1, i) * network->edge_max_flow[network->incoming_edges[k][i]];
             }
        }

        for (int s1 = 0; s1 < 4; s1++) {
            if (max_flow(s1) > 1e-7) {
                cos_intersec[k](s1) /= max_flow(s1);
                sin_intersec[k](s1) /= max_flow(s1);
                }
            else {
                cos_intersec[k](s1) = sin_intersec[k](s1) = 0;
            }
        }

        alpha_intersec[k].setZero();
        beta_intersec[k].setZero();

        // find alpha and beta for each intersection and each direction
        for (unsigned int i = 0; i < network->incoming_edges[k].size(); i++) {
            for (unsigned int j = 0; j < network->outcoming_edges[k].size(); j++) {
                double alpha_flow = network->turning_ratios(network->incoming_edges[k][i], network->outcoming_edges[k][j]) *
                        network->edge_max_flow[network->incoming_edges[k][i]];
                double beta_flow = network->turning_ratios(network->outcoming_edges[k][j], network->incoming_edges[k][i]) *
                        network->edge_max_flow[network->outcoming_edges[k][j]];

                for (int s1 = 0; s1 < 4; s1++) {
                    for (int s2 = 0; s2 < 4; s2++) {
                        alpha_intersec[k](s1, s2) += P_in[k](s1, i) * P_out[k](s2, j) * alpha_flow;
                        beta_intersec[k](s1, s2) += P_in[k](s1, i) * P_out[k](s2, j) * beta_flow;
                    }
                }

            }
        }

        for (int s1 = 0; s1 < 4; s1++) {
            double sum_alpha = 0;
            double sum_beta = 0;
            for (int s2 = 0; s2 < 4; s2++) {
                sum_alpha += alpha_intersec[k](s1, s2);
                sum_beta += beta_intersec[k](s2, s1);
            }
            if (sum_alpha < 1e-8) {
                sum_alpha = 1;
                alpha_intersec[k](s1, 0) = alpha_intersec[k](s1, 1) = alpha_intersec[k](s1, 2) = alpha_intersec[k](s1, 3) = 0.25;
            }
            if (sum_beta < 1e-8) {
                sum_beta = 1;
                beta_intersec[k](0, s1) = beta_intersec[k](1, s1) = beta_intersec[k](2, s1) = beta_intersec[k](3, s1) = 0.25;
            }
            for (int s2 = 0; s2 < 4; s2++) {
                alpha_intersec[k](s1, s2) /= sum_alpha;
                beta_intersec[k](s2, s1) /= sum_beta;
            }
        }


        // calculate L for each intersection
        L_intersec[k] = 0;
        double sum_rho = 0;
        for (unsigned int j = 0; j < network->outcoming_edges[k].size(); j++) {
            Eigen::Vector2d dir = network->nodes[network->edge_list[network->outcoming_edges[k][j]].second]
                            - network->nodes[network->edge_list[network->outcoming_edges[k][j]].first];

            L_intersec[k] += dir.norm() * network->edge_max_density[network->outcoming_edges[k][j]];
            sum_rho += network->edge_max_density[network->outcoming_edges[k][j]];

        }
        if (network->outcoming_edges[k].size() > 0)
            L_intersec[k] /= sum_rho;

        // calculate rho_max for each intersection and each direction
        rho_max_intersec[k].setZero();
        for (int s1 = 0; s1 < 4; s1++) {
            for (unsigned int j = 0; j < network->outcoming_edges[k].size(); j++) {
                rho_max_intersec[k](s1) += P_out[k](s1, j) *  network->edge_max_density[network->outcoming_edges[k][j]];
            }
            for (unsigned int i = 0; i < network->incoming_edges[k].size(); i++) {
                 rho_max_intersec[k](s1) += P_in[k](s1, i) *  network->edge_max_density[network->incoming_edges[k][i]];
            }
        }

        // calculate free-flow speed v for each intersection and each direction
        vel_intersec[k].setZero();
        for (int s1 = 0; s1 < 4; s1++) {
            for (unsigned int j = 0; j < network->outcoming_edges[k].size(); j++) {
                 vel_intersec[k](s1) += P_out[k](s1, j) *  network->edge_max_density[network->outcoming_edges[k][j]]
                       * network->edge_max_velocity[network->outcoming_edges[k][j]];
            }

            for (unsigned int i = 0; i < network->incoming_edges[k].size(); i++) {
                 vel_intersec[k](s1) += P_in[k](s1, i) * network->edge_max_density[network->incoming_edges[k][i]]
                        * network->edge_max_velocity[network->incoming_edges[k][i]];
            }

            if (rho_max_intersec[k](s1) > 1e-8)
                vel_intersec[k](s1) /= rho_max_intersec[k](s1);
            else
                vel_intersec[k](s1) = 0;
        }

        used_for_out_interpolation[k] = network->outcoming_edges[k].size() > 0;
    }  

    qDebug() << "NSWE: intersections processed.";
}



// ------------------------
void NSWEModel::constructInterpolation(double strength) // approximate all parameters to define them for each cell in a 2D plane

{
    pos.resize(CELL_NUMBER);
    alpha.resize(CELL_NUMBER);
    beta.resize(CELL_NUMBER);
    sin_grid.resize(CELL_NUMBER);
    cos_grid.resize(CELL_NUMBER);
    L.resize(CELL_NUMBER);
    vel.resize(CELL_NUMBER);
    rho_max.resize(CELL_NUMBER);
    inflow.resize(CELL_NUMBER);
    outflow.resize(CELL_NUMBER);
    Demand.resize(CELL_NUMBER);
    Supply.resize(CELL_NUMBER);
    rho.resize(CELL_NUMBER);
    new_rho.resize(CELL_NUMBER);
    rho_ema.resize(CELL_NUMBER);

    for (int i = 0; i < CELL_NUMBER; i++) {
        pos[i].resize(CELL_NUMBER);
        alpha[i].resize(CELL_NUMBER);
        beta[i].resize(CELL_NUMBER);
        sin_grid[i].resize(CELL_NUMBER);
        cos_grid[i].resize(CELL_NUMBER);
        L[i].resize(CELL_NUMBER);
        vel[i].resize(CELL_NUMBER);
        rho_max[i].resize(CELL_NUMBER);
        inflow[i].resize(CELL_NUMBER);
        outflow[i].resize(CELL_NUMBER);
        Demand[i].resize(CELL_NUMBER);
        Supply[i].resize(CELL_NUMBER);
        rho[i].resize(CELL_NUMBER);
        new_rho[i].resize(CELL_NUMBER);
        rho_ema[i].resize(CELL_NUMBER);
    }

    dx = network->range_x / CELL_NUMBER;
    dy = network->range_y / CELL_NUMBER;
    for (int i = 0; i < CELL_NUMBER; i++) {
        for (int j = 0; j < CELL_NUMBER; j++) {
            pos[i][j] = Eigen::Vector2d(network->x_min + (i + 0.5) * dx, network->y_min + (j + 0.5) * dy);
        }
    }


    max_density_ever = 0;
    for (int i = 0; i < CELL_NUMBER; i++) {
        for (int j = 0; j < CELL_NUMBER; j++) {

            alpha[i][j].setZero();
            beta[i][j].setZero();
            sin_grid[i][j].setZero();
            cos_grid[i][j].setZero();
            L[i][j] = 0;
            vel[i][j].setZero();
            rho_max[i][j].setZero();
            inflow[i][j].setZero();
            outflow[i][j].setZero();
            double weight_in_out = 0;
            double weight_out = 0;

            for (int k = 0; k < network->nodes.size(); k++) {
                double dist = (network->nodes[k] - pos[i][j]).norm();
                double kernel = exp(-dist * dist * strength * strength);
                dist = exp(-dist * strength) * network->node_weight[k];
                if (!network->node_on_border[k]) {
                    weight_in_out += dist;
                    alpha[i][j] += alpha_intersec[k] * dist;
                    beta[i][j] += beta_intersec[k] * dist;
                    sin_grid[i][j] += sin_intersec[k] * dist;
                    cos_grid[i][j] += cos_intersec[k] * dist;
                    rho_max[i][j] += rho_max_intersec[k] * dist;
                    vel[i][j] += vel_intersec[k] * dist;
                }
                if (used_for_out_interpolation[k]) {
                    weight_out += dist;
                    L[i][j] += L_intersec[k] * dist;
                }
            }


            alpha[i][j] /= weight_in_out;
            beta[i][j] /= weight_in_out;
            L[i][j] /= weight_out;
            sin_grid[i][j] /= weight_in_out;
            cos_grid[i][j] /= weight_in_out;
            rho_max[i][j] /= weight_in_out;
            vel[i][j] /= weight_in_out;

            if (rho_max[i][j].sum() > max_density_ever)
                max_density_ever = rho_max[i][j].sum();

            rho[i][j].setZero();
            new_rho[i][j] = rho[i][j];
            rho_ema[i][j] = rho[i][j];
        }
    }


    qDebug() << "NSWE: interpolation is ready.";
}



// ---------------------- density update is defined here (Godunov scheme) ---- implementation of the NSWE model
void NSWEModel::update(double dt)
{

    for (int i = 0; i < CELL_NUMBER; i++) {
        for (int j = 0; j < CELL_NUMBER; j++) {
            Demand[i][j] = demand(rho[i][j], rho_max[i][j], vel[i][j]);
            Supply[i][j] = supply(rho[i][j], rho_max[i][j], vel[i][j]);
        }
    }

    for (int i = 1; i < CELL_NUMBER - 1; i++) {
        for (int j = 1; j < CELL_NUMBER - 1; j++) {
            Eigen::Vector4d density_update = Eigen::Vector4d::Zero();

            Eigen::Matrix4d partial_flows;

            for (int s0 = 0; s0 < 4; s0++) {
                for (int s1 = 0; s1 < 4; s1++) {
                    partial_flows(s0, s1) = std::min(alpha[i][j](s0,s1) * Demand[i][j](s0), beta[i][j](s0,s1) * Supply[i][j](s1));
                }
            }

            // mixing term
            for (int s0 = 0; s0 < 4; s0++) {
                for (int s1 = 0; s1 < 4; s1++) {
                    density_update(s0) += partial_flows(s1, s0) / L[i][j];
                    density_update(s0) -= partial_flows(s0, s1) / L[i][j];
                }
            }

            // term of extended model (inflow + outflow external)
            for (int s1 = 0; s1 < 4; s1++) {
                density_update(s1) += 1/L[i][j] * (std::min(inflow[i][j](s1), Supply[i][j](s1))
                        - std::min(Demand[i][j](s1), outflow[i][j](s1)));
            }

            // transportation part
            for (int s1 = 0; s1 < 4; s1++) {

                double sin_prev = (sin_grid[i][j](s1) + sin_grid[i][j - 1](s1)) / 2;
                double sin_next = (sin_grid[i][j](s1) + sin_grid[i][j + 1](s1)) / 2;
                double cos_prev = (cos_grid[i][j](s1) + cos_grid[i - 1][j](s1)) / 2;
                double cos_next = (cos_grid[i][j](s1) + cos_grid[i + 1][j](s1)) / 2;

                // Godunov scheme, upwind scheme
                if (cos_prev > 0)
                    density_update(s1) += cos_prev * std::min(Demand[i - 1][j](s1), Supply[i][j](s1)) / dx;
                else
                    density_update(s1) += cos_prev * std::min(Demand[i][j](s1), Supply[i - 1][j](s1)) / dx;


                if (cos_next > 0)
                    density_update(s1) -= cos_next * std::min(Demand[i][j](s1), Supply[i + 1][j](s1)) / dx;
                else
                    density_update(s1) -= cos_next * std::min(Demand[i + 1][j](s1), Supply[i][j](s1)) / dx;


                if (sin_prev > 0)
                    density_update(s1) += sin_prev * std::min(Demand[i][j - 1](s1), Supply[i][j](s1)) / dy;
                else
                    density_update(s1) += sin_prev * std::min(Demand[i][j](s1), Supply[i][j - 1](s1)) / dy;


                if (sin_next > 0)
                    density_update(s1) -= sin_next * std::min(Demand[i][j](s1), Supply[i][j + 1](s1)) / dy;
                else
                    density_update(s1) -= sin_next * std::min(Demand[i][j + 1](s1), Supply[i][j](s1)) / dy;

            }

            new_rho[i][j] = rho[i][j] + dt * density_update;
        }
    }

    // smoothing density using EMA procedure
    for (int i = 1; i < CELL_NUMBER - 1; i++) {
        for (int j = 1; j < CELL_NUMBER - 1; j++) {
            rho[i][j] = new_rho[i][j];

            rho_ema[i][j] = (1 - ema_parameter) * rho_ema[i][j] + ema_parameter * rho[i][j];
        }
    }

}



//---------------displaying the NSWE density in the rectangle for GUI visualization
void NSWEModel::displayDensity(QPainter& painter, QRect rectangle, bool paint_network)
{
    QImage image(rectangle.size(), QImage::Format_RGB32);

    for (int i = 0; i < rectangle.width(); i++) {
        for (int j = 0; j < rectangle.height(); j++) {
            QColor color(0, 0, 0);
            double x = network->x_min + network->range_x * (i + 0.5) / rectangle.width();
            double y = network->y_min + network->range_y * (rectangle.height() - j - 0.5) / rectangle.height();


            Eigen::Vector4d density = interpolateDensity(x, y);

            double sum_val = density.sum() / max_density_ever;
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


//---------------interpolation of the NSWE density for any point between cells for the purpose of smooth visualization; biliniear interpolation
Eigen::Vector4d NSWEModel::interpolateDensity(double x, double y)
{
    x = (x - network->x_min) / dx - 0.5;
    y = (y - network->y_min) / dy - 0.5;

    if (x < 0) x = 0;
    if (y < 0) y = 0;
    if (x > CELL_NUMBER - 1 - 1e-8) x = CELL_NUMBER - 1 - 1e-8;
    if (y > CELL_NUMBER - 1 - 1e-8) y = CELL_NUMBER - 1 - 1e-8;

    int int_x = int(x);
    int int_y = int(y);
    double fl_x = x - int_x;
    double fl_y = y - int_y;


    Eigen::Vector4d val11;
    Eigen::Vector4d val12;
    Eigen::Vector4d val21;
    Eigen::Vector4d val22;

    val11 = rho_ema[int_x][int_y];
    val12 = rho_ema[int_x][int_y + 1];
    val21 = rho_ema[int_x + 1][int_y];
    val22 = rho_ema[int_x + 1][int_y + 1];

    return val11 * (1 - fl_x) * (1 - fl_y) + val12 * (1 - fl_x) * fl_y +
                    val21 * fl_x * (1 - fl_y) + val22 * fl_x * fl_y;
}


//---------------------this is used to get the value of L in NSWE model for real-data density reconstruction (see GrenobleData::reconstructDensity)
double NSWEModel::interpolateLength(double x, double y)
{
    x = (x - network->x_min) / dx - 0.5;
    y = (y - network->y_min) / dy - 0.5;

    if (x < 0) x = 0;
    if (y < 0) y = 0;
    if (x > CELL_NUMBER - 1 - 1e-8) x = CELL_NUMBER - 1 - 1e-8;
    if (y > CELL_NUMBER - 1 - 1e-8) y = CELL_NUMBER - 1 - 1e-8;

    int int_x = int(x);
    int int_y = int(y);
    double fl_x = x - int_x;
    double fl_y = y - int_y;


    double val11 = L[int_x][int_y];
    double val12 = L[int_x][int_y + 1];
    double val21 = L[int_x + 1][int_y];
    double val22 = L[int_x + 1][int_y + 1];

    return val11 * (1 - fl_x) * (1 - fl_y) + val12 * (1 - fl_x) * fl_y +
                    val21 * fl_x * (1 - fl_y) + val22 * fl_x * fl_y;
}


//------------------set density for initial time step (in mainwindow we set it equal to the reconstructed density from real data)
void NSWEModel::setInitialDensity(Density4 density)
{
    for (int i = 0; i < CELL_NUMBER; i++) {
        for (int j = 0; j < CELL_NUMBER; j++) {
            rho[i][j](0) = density.N.data(i, j);
            rho[i][j](1) = density.S.data(i, j);
            rho[i][j](2) = density.W.data(i, j);
            rho[i][j](3) = density.E.data(i, j);

            rho_ema[i][j] = rho[i][j];
        }
    }
}


// --------------- calculate the weighted SSIM index to compare two density distributions ---------------------

int N_BATCH = 20;  // corresponds to dividing the domain into 3x3 = 9 zones, since CELL_NUMBER/N_BATCH = 3
double C = 0.000000000001;   // small constant used to avoid division by zero

double NSWEModel::getSSIMDiff_mean_weighted(Density4 density)
{
    int num_batches = CELL_NUMBER / N_BATCH;

    Eigen::VectorXd weights(num_batches * num_batches);
    res_SSIM.resize(num_batches * num_batches);

    double res_ssim = 0;

    for (int q = 0; q < num_batches; q++) {
        for (int r = 0; r < num_batches; r++) {

            double res = 0;
            double mean_1 = 0;
            double mean_2 = 0;
            double var_1 = 0;
            double var_2 = 0;
            double cross_var = 0;

            for (int i = 0; i < N_BATCH; i++) {
                for (int j = 0; j < N_BATCH; j++) {
                    int i1 = q * N_BATCH + i;
                    int j1 = r * N_BATCH + j;
                    mean_1 += rho_ema[i1][j1].sum();
                    mean_2 += density.N.data(i1, j1) + density.S.data(i1, j1)
                                                    + density.W.data(i1, j1) + density.E.data(i1, j1);
                }
            }
            mean_1 /= N_BATCH * N_BATCH;
            mean_2 /= N_BATCH * N_BATCH;

            weights(q * num_batches + r) = mean_2 + 0.0001;   // we add a small value here to avoid division by zero

            for (int i = 0; i < N_BATCH; i++) {
                for (int j = 0; j < N_BATCH; j++) {
                    int i1 = q * N_BATCH + i;
                    int j1 = r * N_BATCH + j;
                    var_1 += (rho_ema[i1][j1].sum() - mean_1) * (rho_ema[i1][j1].sum() - mean_1);

                    var_2 += (density.N.data(i1, j1) + density.S.data(i1, j1)
                                                   + density.W.data(i1, j1) + density.E.data(i1, j1) - mean_2) *
                                                    (density.N.data(i1, j1) + density.S.data(i1, j1)
                                                     + density.W.data(i1, j1) + density.E.data(i1, j1) - mean_2);

                    cross_var += (rho_ema[i1][j1].sum() - mean_1) * (density.N.data(i1, j1)
                                   + density.S.data(i1, j1) + density.W.data(i1, j1) + density.E.data(i1, j1) - mean_2);
                }
            }
            var_1 /= N_BATCH * N_BATCH;
            var_2 /= N_BATCH * N_BATCH;
            cross_var /= N_BATCH * N_BATCH;

            var_1 = sqrt(var_1);
            var_2 = sqrt(var_2);

            res = (2 * mean_1 * mean_2 + C) * (2 * cross_var + C) / (mean_1 * mean_1 + mean_2 * mean_2 + C) / (var_1 * var_1 + var_2 * var_2 + C);
            res = (res + 1) / 2;
            res_ssim += res * weights(q * num_batches + r);
            res_SSIM(q * num_batches + r) = res;
        }
    }

    res_ssim /= weights.sum();
    return res_ssim;
}
