#pragma once

#define _USE_MATH_DEFINES

#include <math.h>

float smoothingKernel(float radius, float dst) {
    if (dst < radius) {
        float volume = M_PI * std::pow(radius, 8) / 4;
        float value = radius * radius - dst * dst;
        return std::pow(value, 3) / volume;
    }
    return 0;
}