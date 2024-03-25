#define _USE_MATH_DEFINES

#include <math.h>

float smoothingKernel(float radius, float dst) {
    float volume = M_PI * std::pow(radius, 8) / 4;
    float value = std::max(0.f, static_cast<float>(radius * radius - dst * dst));
    return std::pow(value, 3) / volume;
}