#ifndef GRADIENT_H
#define GRADIENT_H

#include <vector>
#include <stdexcept>

template<typename T>
std::vector<T> calculateGradient(const std::vector<T>& v) {
    size_t n = v.size();
    if (n == 0) {
        throw std::invalid_argument("The input vector is empty.");
    }
    if (n == 1) {
        return std::vector<T>(1, 0); // Gradient of a single value is 0
    }

    std::vector<T> gradient(n);

    // Calculate gradient for the first element
    gradient[0] = v[1] - v[0];

    // Calculate gradient for elements in the middle
    for (size_t i = 1; i < n - 1; ++i) {
        gradient[i] = (v[i + 1] - v[i - 1]) / T(2);
    }

    // Calculate gradient for the last element
    gradient[n - 1] = v[n - 1] - v[n - 2];

    return gradient;
}

#endif