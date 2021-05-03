#ifndef UTILS_H
#define UTILS_H

#include <QColor>
#include <QDebug>

inline QColor getColorGreenRed(double value)
{
    if (value < 0) value = 0;
    if (value > 1) value = 1;

    int r, g, b;

    if (value < 0.0625) {
        r = 128 * value / 0.0625;
        g = 128 + 64 * value / 0.0625;
        b = 102;
    }
    else if (value < 0.1875) {
        r = 128 + 127 * (value - 0.0625) / 0.125;
        g = 192 + 63 * (value - 0.0625) / 0.125;
        b = 102 - 51 * (value - 0.0625) / 0.125;
    }
    else if (value < 0.5625) {
        r = 255;
        g = 255 - 153 * (value - 0.1875) / 0.375;
        b = 51 - 51 * (value - 0.1875) / 0.375;
    }
    else if (value < 0.8125) {
        r = 255;
        g = 102 - 102 * (value - 0.5625) / 0.25;
        b = 0;
    }
    else {
        r = 255 - 76 * (value - 0.8125) / 0.1875;
        g = 0;
        b = 0;
    }

    return QColor(r, g, b);
}

#endif // UTILS_H
