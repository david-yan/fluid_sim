#pragma once

#define _USE_MATH_DEFINES

#include <math.h>

float smoothingKernel(float radius, float dst) {
    if (dst >= radius) return 0;

    float volume = M_PI * std::pow(radius, 8) / 4;
    float value = radius * radius - dst * dst;
    return std::pow(value, 3) / volume;
}

float smoothingKernelDerivative(float radius, float dst) {
    if (dst >= radius) return 0;

    float f = radius * radius - dst * dst;
    float scale = -24.f / (M_PI * std::pow(radius, 8));
    return scale * dst * f * f;
}

float spikyPow2Kernel(float radius, float dst) {
    if (dst >= radius) return 0;

    float volume = M_PI * std::pow(radius, 4) / 6;
    float value = radius - dst;
    return std::pow(value, 2) / volume;
}

float spikyPow2KernelDerivative(float radius, float dst) {
    if (dst >= radius) return 0;

    float f = radius - dst;
    float scale = -12.f / (M_PI * std::pow(radius, 4));
    return scale * f;
}

float spikyPow3Kernel(float radius, float dst) {
    if (dst >= radius) return 0;
    
    float volume = M_PI * std::pow(radius, 5) / 10;
    float value = radius - dst;
    return std::pow(value, 3) / volume;
}

float spikyPow3KernelDerivative(float radius, float dst) {
    if (dst >= radius) return 0;

    float f = radius - dst;
    float scale = -30.f / (M_PI * std::pow(radius, 5));
    return scale * f * f;
}