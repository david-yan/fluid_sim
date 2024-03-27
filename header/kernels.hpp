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

float smoothingKernelDerivative(float radius, float dst) {
    if (dst >= radius) return 0;

    float f = radius * radius - dst * dst;
    float scale = -24.f / (M_PI * std::pow(radius, 8));
    return scale * dst * f * f;
}

float spikyKernel(float radius, float dst) {
    if (dst < radius) {
        float volume = M_PI * std::pow(radius, 5) / 10;
        float value = radius - dst;
        return std::pow(value, 3) / volume;
    }
    return 0;
}

float spikyKernelDerivative(float radius, float dst) {
    if (dst >= radius) return 0;

    float f = radius - dst;
    float scale = -30.f / (M_PI * std::pow(radius, 5));
    return scale * dst * f * f;
}