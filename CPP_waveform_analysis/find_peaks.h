// WaveformAnalysis.h

#ifndef FIND_PEAKS_H
#define FIND_PEAKS_H

#include <vector>
#include <stdexcept>
#include <iostream>
#include <functional>
#include <map>
#include <limits>
#include <cmath>
#include <tuple>

// Comparator functions
bool less(int a, int b);
bool greater(int a, int b);

// Modified boolrelextrema function to work with comparators
std::vector<bool> boolrelextrema(const std::vector<int>& data, std::function<bool(int, int)> comparator, int order = 1, const std::string& mode = "clip");

// argrelextrema function adapted for C++
std::vector<int> argrelextrema(const std::vector<int>& data, std::function<bool(int, int)> comparator, int order = 1, const std::string& mode = "clip");

// Wrapper functions for minima and maxima
std::vector<int> argrelmin(const std::vector<int>& data, int order = 1, const std::string& mode = "clip");
std::vector<int> argrelmax(const std::vector<int>& data, int order = 1, const std::string& mode = "clip");

// Function to calculate peak prominences
std::tuple<std::vector<double>, std::vector<int>, std::vector<int>> calculate_peak_prominences(
    const std::vector<double>& x, 
    const std::vector<int>& peaks, 
    int wlen = -1
);

// Function to calculate peak widths
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> calculate_peak_widths(
    const std::vector<double>& x,
    const std::vector<int>& peaks,
    double relHeight,
    const std::vector<double>& prominences,
    const std::vector<int>& leftBases,
    const std::vector<int>& rightBases
);

// Function to unpack condition arguments
std::pair<double, double> _unpack_condition_args(double interval, const std::vector<double>& x, const std::vector<int>& peaks);
std::pair<std::vector<double>, std::vector<double>> _unpack_condition_args(const std::pair<std::vector<double>, std::vector<double>>& interval, const std::vector<double>& x, const std::vector<int>& peaks);

// Function to select by property
std::vector<bool> _select_by_property(const std::vector<double>& peak_properties, double pmin, double pmax);

// Function to select peaks by threshold
std::tuple<std::vector<bool>, std::vector<double>, std::vector<double>> _select_by_peak_threshold(const std::vector<double>& x, const std::vector<int>& peaks, double tmin, double tmax);

// // Structure and function for finding local maxima
// struct LocalMaximaResult {
//     std::vector<int> midpoints;
//     std::vector<int> leftEdges;
//     std::vector<int> rightEdges;
// };

std::tuple<std::vector<int>, std::vector<int>, std::vector<int>> findLocalMaxima1D(const std::vector<double>& x);

// Main function to find peaks with various conditions
std::pair<std::vector<int>, std::map<std::string, std::vector<double>>> findPeaks(
    const std::vector<double>& x,
    double height = std::numeric_limits<double>::quiet_NaN(),
    double threshold = std::numeric_limits<double>::quiet_NaN(),
    int distance = 1,
    double prominence = std::numeric_limits<double>::quiet_NaN(),
    double width = std::numeric_limits<double>::quiet_NaN(),
    int wlen = -1,
    double relHeight = 0.5
);

#endif // FIND_PEAKS_H
