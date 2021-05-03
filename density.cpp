#include "density.h"
#include "utils.h"

Density::Density(int grid_size) :
    grid_size(grid_size),
    data(grid_size, grid_size)
{
    domain = nullptr;
    data.setZero();
}


//----------------gaussian distribution is approximated as 3 convolutions with uniform distribution; each of them is implemented in box_blur
void Density::gaussian_blur(double standart_deviation, double network_scale)
{
    int n_boxes = 3;

    standart_deviation *= grid_size / domain->range_x / network_scale;

    double wIdeal = sqrt((12 * standart_deviation * standart_deviation / n_boxes) + 1);  // Ideal averaging filter width

    int wl = static_cast<int>(floor(wIdeal));
    if(wl % 2 == 0)
        wl--;
    int wu = wl+2;

    double mIdeal = (12 * standart_deviation * standart_deviation - n_boxes * wl * wl -
                     4 * n_boxes * wl - 3 * n_boxes) / (-4 * wl - 4);
    int m = static_cast<int>(round(mIdeal));

    for(int i = 0; i < n_boxes; i++) {
        int blur_size = i < m ? wl : wu;
        this->box_blur(blur_size);
    }
}

void Density::box_blur(int blur_size)
{
    Eigen::MatrixXd accumulator(grid_size, grid_size);

    accumulator(0, 0) = data(0, 0);

    for (int i = 1; i < grid_size; i++) {
        accumulator(i, 0) = accumulator(i - 1, 0) + data(i, 0);
    }
    for (int j = 1; j < grid_size; j++) {
        accumulator(0, j) = accumulator(0, j - 1) + data(0, j);
    }
    for (int i = 1; i < grid_size; i++) {
        for (int j = 1; j < grid_size; j++) {
            accumulator(i, j) = accumulator(i - 1, j) + accumulator(i, j - 1) - accumulator(i - 1, j - 1) + data(i, j);
        }
    }

    int blur_r = (blur_size - 1) / 2;

    //00
    for (int i = 0; i < blur_r + 1; i++) {
        for (int j = 0; j < blur_r + 1; j++) {
            data(i, j) = accumulator(i + blur_r, j + blur_r);
        }
    }

    //01
    for (int i = 0; i < blur_r + 1; i++) {
        for (int j = blur_r + 1; j < grid_size - blur_r; j++) {
            data(i, j) = accumulator(i + blur_r, j + blur_r) - accumulator(i + blur_r, j - blur_r - 1);
        }
    }

    //02
    for (int i = 0; i < blur_r + 1; i++) {
        for (int j = grid_size - blur_r; j < grid_size; j++) {
            data(i, j) = accumulator(i + blur_r, grid_size - 1) - accumulator(i + blur_r, j - blur_r - 1);
        }
    }

    //10
    for (int i = blur_r + 1; i < grid_size - blur_r; i++) {
        for (int j = 0; j < blur_r + 1; j++) {
            data(i, j) = accumulator(i + blur_r, j + blur_r) - accumulator(i - blur_r - 1, j + blur_r);
        }
    }

    //11
    for (int i = blur_r + 1; i < grid_size - blur_r; i++) {
        for (int j = blur_r + 1; j < grid_size - blur_r; j++) {
            data(i, j) = accumulator(i + blur_r, j + blur_r) + accumulator(i - blur_r - 1, j - blur_r - 1)
                   - accumulator(i + blur_r, j - blur_r - 1) - accumulator(i - blur_r - 1, j + blur_r);
        }
    }

    //12
    for (int i = blur_r + 1; i < grid_size - blur_r; i++) {
        for (int j = grid_size - blur_r; j < grid_size; j++) {
            data(i, j) = accumulator(i + blur_r, grid_size - 1) + accumulator(i - blur_r - 1, j - blur_r - 1)
                   - accumulator(i + blur_r, j - blur_r - 1) - accumulator(i - blur_r - 1, grid_size - 1);
        }
    }

    //20
    for (int i = grid_size - blur_r; i < grid_size; i++) {
        for (int j = 0; j < blur_r + 1; j++) {
            data(i, j) = accumulator(grid_size - 1, j + blur_r) - accumulator(i - blur_r - 1, j + blur_r);
        }
    }

    //21
    for (int i = grid_size - blur_r; i < grid_size; i++) {
        for (int j = blur_r + 1; j < grid_size - blur_r; j++) {
            data(i, j) = accumulator(grid_size - 1, j + blur_r) + accumulator(i - blur_r - 1, j - blur_r - 1)
                   - accumulator(grid_size - 1, j - blur_r - 1) - accumulator(i - blur_r - 1, j + blur_r);
        }
    }

    //22
    for (int i = grid_size - blur_r; i < grid_size; i++) {
        for (int j = grid_size - blur_r; j < grid_size; j++) {
            data(i, j) = accumulator(grid_size - 1, grid_size - 1) + accumulator(i - blur_r - 1, j - blur_r - 1)
                   - accumulator(grid_size - 1, j - blur_r - 1) - accumulator(i - blur_r - 1, grid_size - 1);
        }
    }

    data /= blur_size * blur_size;
}

void Density::clear()
{
    data.setZero();
}

void Density::addValue(double x, double y, double value)
{
    if (x < domain->x_min || x > domain->x_min + domain->range_x ||
            y < domain->y_min || y > domain->y_min + domain->range_y) {
        qDebug() << "ERROR Density: value setting request is out of range!";
        exit(3);
    }

    double dx = domain->range_x / grid_size;
    double dy = domain->range_y / grid_size;

    x = (x - domain->x_min) / dx - 0.5;
    y = (y - domain->y_min) / dy - 0.5;

    if (x < 0) x = 0;
    if (y < 0) y = 0;
    if (x > grid_size - 1 - 1e-8) x = grid_size - 1 - 1e-8;
    if (y > grid_size - 1 - 1e-8) y = grid_size - 1 - 1e-8;


    int int_x = int(x);
    int int_y = int(y);
    double fl_x = x - int_x;
    double fl_y = y - int_y;

    data(int_x, int_y) += value * (1 - fl_x) * (1 - fl_y);
    data(int_x, int_y + 1) += value * (1 - fl_x) * fl_y;
    data(int_x + 1, int_y) += value * fl_x * (1 - fl_y);
    data(int_x + 1, int_y + 1) += value * fl_x * fl_y;
}



void Density::setDomain(UrbanNetwork* network_domain)
{
    domain = network_domain;
}


double Density::getNormalizedValue(double x, double y)
{
    x *= (grid_size - 1);
    y *= (grid_size - 1);
    int int_x = int(x);
    int int_y = int(y);
    double fl_x = x - int_x;
    double fl_y = y - int_y;

    double val11 = data(int_x, int_y);
    double val12 = data(int_x, int_y + 1);
    double val21 = data(int_x + 1, int_y);
    double val22 = data(int_x + 1, int_y + 1);

    return val11 * (1 - fl_x) * (1 - fl_y) + val12 * (1 - fl_x) * fl_y +
                    val21 * fl_x * (1 - fl_y) + val22 * fl_x * fl_y;

}

