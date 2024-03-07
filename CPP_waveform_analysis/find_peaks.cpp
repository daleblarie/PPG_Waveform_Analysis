# include "find_peaks.h"

// Comparator functions
bool less(int a, int b) {
    return a < b;
}

// Function to compare two elements and return true if the first is greater than the second
bool greater(int a, int b) {
    return a > b;
}

// Modified boolrelextrema function to work with comparators
std::vector<bool> boolrelextrema(const std::vector<int>& data, std::function<bool(int, int)> comparator, int order, const std::string& mode) {
    if (order < 1) {
        throw std::invalid_argument("Order must be an int >= 1");
    }

    int dataSize = data.size();
    std::vector<bool> extrema(dataSize, true); // Initialize all to true

    for (int i = 0; i < dataSize; ++i) {
        for (int shift = 1; shift <= order; ++shift) {
            int plusIndex = i + shift;
            int minusIndex = i - shift;

            if (mode == "clip") {
                plusIndex = std::min(plusIndex, dataSize - 1);
                minusIndex = std::max(minusIndex, 0);
            } else if (mode == "wrap") {
                plusIndex = (plusIndex + dataSize) % dataSize;
                minusIndex = (minusIndex + dataSize) % dataSize;
            } else {
                throw std::invalid_argument("Invalid mode. Use 'clip' or 'wrap'.");
            }

            extrema[i] = extrema[i] && comparator(data[i], data[plusIndex]) && comparator(data[i], data[minusIndex]);
            if (!extrema[i])
                break;
        }
    }

    return extrema;
}

// argrelextrema function adapted for C++
std::vector<int> argrelextrema(const std::vector<int>& data, std::function<bool(int, int)> comparator, int order, const std::string& mode) {
    auto extremaBool = boolrelextrema(data, comparator, order, mode);
    std::vector<int> extremaIndices;
    for (size_t i = 0; i < extremaBool.size(); ++i) {
        if (extremaBool[i]) {
            extremaIndices.push_back(i);
        }
    }
    return extremaIndices;
}

// Wrapper functions for minima and maxima
std::vector<int> argrelmin(const std::vector<int>& data, int order, const std::string& mode) {
    return argrelextrema(data, less, order, mode);
}

std::vector<int> argrelmax(const std::vector<int>& data, int order, const std::string& mode) {
    return argrelextrema(data, greater, order, mode);
}


std::tuple<std::vector<double>, std::vector<int>, std::vector<int>> calculate_peak_prominences(
    const std::vector<double>& x, 
    const std::vector<int>& peaks, 
    int wlen
) {
    if (x.empty() || peaks.empty()) {
        // std::cout<<"Input vectors x or peaks are empty."<<std::endl;
        std::vector<double> prominences;
        std::vector<int> left_bases;
        std::vector<int> right_bases;
        return {prominences, left_bases, right_bases};
    }

    std::vector<double> prominences(peaks.size(), 0.0);
    std::vector<int> left_bases(peaks.size(), 0);
    std::vector<int> right_bases(peaks.size(), 0);

    int x_size = static_cast<int>(x.size());

    for (size_t peak_nr = 0; peak_nr < peaks.size(); ++peak_nr) {
        int peak = peaks[peak_nr];
        if (peak < 0 || peak >= x_size) {
            throw std::invalid_argument("Peak index out of range.");
        }

        int i_min = 0;
        int i_max = x_size - 1;

        if (wlen >= 2) {
            // Adjust window around the evaluated peak (within bounds)
            i_min = std::max(peak - wlen / 2, i_min);
            i_max = std::min(peak + wlen / 2, i_max);
        }

        // Find the left base
        int i = peak;
        double left_min = x[peak];
        left_bases[peak_nr] = i;
        while (i_min <= i && x[i] <= x[peak]) {
            if (x[i] < left_min) {
                left_min = x[i];
                left_bases[peak_nr] = i;
            }
            i--;
        }

        // Find the right base
        i = peak;
        double right_min = x[peak];
        right_bases[peak_nr] = i;
        while (i <= i_max && x[i] <= x[peak]) {
            if (x[i] < right_min) {
                right_min = x[i];
                right_bases[peak_nr] = i;
            }
            i++;
        }

        prominences[peak_nr] = x[peak] - std::max(left_min, right_min);
    }

    return {prominences, left_bases, right_bases};
}


std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> calculate_peak_widths(
    const std::vector<double>& x,
    const std::vector<int>& peaks,
    double relHeight,
    const std::vector<double>& prominences,
    const std::vector<double>& leftBase,
    const std::vector<double>& rightBase) {
    // std::cout<<"calculate_peak_widths 1"<<std::endl;
    std::vector<int> leftBases = std::vector<int>(leftBase.begin(), leftBase.end());
    std::vector<int> rightBases = std::vector<int>(rightBase.begin(), rightBase.end());
    if (relHeight < 0) {
        throw std::invalid_argument("`relHeight` must be greater or equal to 0.0");
    }
    // std::cout<<"calculate_peak_widths 2"<<std::endl;
    // std::cout<<"peaks.size()"<<peaks.size()<<std::endl;
    // std::cout<<"prominences.size()"<<prominences.size()<<std::endl;
    // std::cout<<"leftBases.size()"<<leftBases.size()<<std::endl;
    // std::cout<<"rightBases.size()"<<rightBases.size()<<std::endl;
    if (!(peaks.size() == prominences.size() && peaks.size() == leftBases.size() && peaks.size() == rightBases.size())) {
        std::cout<<"INPUTS ARE NOT ALL THE SAME SIZE"<<std::endl;
        throw std::invalid_argument("All input arrays must have the same size");
    }
    // std::cout<<"calculate_peak_widths 3"<<std::endl;

    std::vector<double> widths(peaks.size(), 0.0);
    std::vector<double> widthHeights(peaks.size(), 0.0);
    std::vector<double> leftIps(peaks.size(), 0.0);
    std::vector<double> rightIps(peaks.size(), 0.0);
    // std::cout<<"calculate_peak_widths 4"<<std::endl;

    for (size_t p = 0; p < peaks.size(); ++p) {
        int i_min = leftBases[p];
        int i_max = rightBases[p];
        int peak = peaks[p];
        if (!(0 <= i_min && i_min <= peak && peak <= i_max && i_max < static_cast<int>(x.size()))) {
            throw std::invalid_argument("Prominence data is invalid for peak at index " + std::to_string(peak));
        }
        double height = widthHeights[p] = x[peak] - prominences[p] * relHeight;

        // Find intersection point on left side
        int i = peak;
        while (i_min < i && height < x[i]) {
            --i;
        }
        double left_ip = static_cast<double>(i);
        if (x[i] < height) {
            left_ip += (height - x[i]) / (x[i + 1] - x[i]);
        }

        // Find intersection point on right side
        i = peak;
        while (i < i_max && height < x[i]) {
            ++i;
        }
        double right_ip = static_cast<double>(i);
        if (x[i] < height) {
            right_ip -= (height - x[i]) / (x[i - 1] - x[i]);
        }

        widths[p] = right_ip - left_ip;
        leftIps[p] = left_ip;
        rightIps[p] = right_ip;
    }
    // std::cout<<"calculate_peak_widths 5"<<std::endl;

    // Note: The original Python code includes a warning for peaks with a width of 0.
    // In C++, you might handle this with logging or other mechanisms as appropriate.

    return {widths, widthHeights, leftIps, rightIps};
}

std::pair<double, double> _unpack_condition_args(double interval, const std::vector<double>& x, const std::vector<int>& peaks) {
    // Single value case: imin is interval, imax is not used
    return {interval, std::numeric_limits<double>::quiet_NaN()};
}

std::pair<std::vector<double>, std::vector<double>> _unpack_condition_args(const std::pair<std::vector<double>, std::vector<double>>& interval, const std::vector<double>& x, const std::vector<int>& peaks) {
    if (interval.first.size() != x.size() || interval.second.size() != x.size()) {
        throw std::invalid_argument("Array size of interval borders must match x");
    }
    
    std::vector<double> imin(peaks.size()), imax(peaks.size());
    for (size_t i = 0; i < peaks.size(); ++i) {
        imin[i] = interval.first[peaks[i]];
        imax[i] = interval.second[peaks[i]];
    }
    
    return {imin, imax};
}



template <typename T>
std::vector<int> argsort(const std::vector<T>& v) {
    // Initialize original indices
    std::vector<int> idx(v.size());
    for(size_t i = 0; i < idx.size(); ++i) {
        idx[i] = i;
    }

    // Sort indices based on comparing values in v
    std::sort(idx.begin(), idx.end(),
        [&v](int i1, int i2) {return v[i1] < v[i2];});

    // std::cout<<"ARGSORT TIME"<<std::endl;
    // for (auto & el: idx){
    //     std::cout<<el<<" ";
    // }

    return idx;
}


std::vector<bool> selectByPeakDistance(const std::vector<int>& peaks,
                                       const std::vector<double>& priority,
                                       double distance) {
    size_t peaksSize = peaks.size();
    std::vector<bool> keep(peaksSize, true);
    std::vector<int> priorityToPosition = argsort(priority);
    
    
    // Convert distance to an integer value, assuming distance is in units of array indices
    int distanceInt = static_cast<int>(std::ceil(distance));
    
    for (int i = peaksSize - 1; i >= 0; --i) {
        int j = priorityToPosition[i];
        
        if (!keep[j]) continue; // Skip if already marked for exclusion
        
        int k = j - 1;
        while (k >= 0 && peaks[j] - peaks[k] < distanceInt) {
            keep[k] = false;
            --k;
        }
        
        k = j + 1;
        while (k < static_cast<int>(peaksSize) && peaks[k] - peaks[j] < distanceInt) {
            keep[k] = false;
            ++k;
        }
    }
    
    return keep;
}

std::vector<bool> _select_by_property(const std::vector<double>& peak_properties, double pmin, double pmax) {
    std::vector<bool> keep(peak_properties.size(), true);
    for (size_t i = 0; i < peak_properties.size(); ++i) {
        if (pmin != std::numeric_limits<double>::quiet_NaN()) {
            keep[i] = keep[i] && (pmin <= peak_properties[i]);
        }
        if (pmax != std::numeric_limits<double>::quiet_NaN()) {
            keep[i] = keep[i] && (peak_properties[i] <= pmax);
        }
    }
    return keep;
}

std::tuple<std::vector<bool>, std::vector<double>, std::vector<double>> _select_by_peak_threshold(const std::vector<double>& x, const std::vector<int>& peaks, double tmin, double tmax) {
    std::vector<bool> keep(peaks.size(), true);
    std::vector<double> left_thresholds(peaks.size()), right_thresholds(peaks.size());
    
    for (size_t i = 0; i < peaks.size(); ++i) {
        int peak = peaks[i];
        left_thresholds[i] = x[peak] - x[peak - 1];
        right_thresholds[i] = x[peak] - x[peak + 1];
        
        double min_threshold = std::min(left_thresholds[i], right_thresholds[i]);
        double max_threshold = std::max(left_thresholds[i], right_thresholds[i]);
        
        if (tmin != std::numeric_limits<double>::quiet_NaN()) {
            keep[i] = keep[i] && (tmin <= min_threshold);
        }
        if (tmax != std::numeric_limits<double>::quiet_NaN()) {
            keep[i] = keep[i] && (max_threshold <= tmax);
        }
    }
    
    return {keep, left_thresholds, right_thresholds};
}



std::tuple<std::vector<int>, std::vector<int>, std::vector<int>> findLocalMaxima1D(const std::vector<double>& x) {
    std::vector<int> midpoints;
    std::vector<int> leftEdges;
    std::vector<int> rightEdges;
    // std::cout<<"findLocalMaxima1D 1"<<std::endl;
    int i = 1;  // Start from the second element
    int iMax = x.size() - 1;  // The last element can't be a maximum
    // std::cout<<"findLocalMaxima1D 2"<<std::endl;

    while (i < iMax) {
        if (x[i - 1] < x[i]) {
            int iAhead = i + 1;
            // Find next sample that is not equal to x[i]
            while (iAhead < iMax && x[iAhead] == x[i]) {
                iAhead++;
            }
            // A maximum is found if the next nonequal sample is smaller than x[i]
            if (x[iAhead] < x[i]) {
                leftEdges.push_back(i);
                rightEdges.push_back(iAhead - 1);
                midpoints.push_back((leftEdges.back() + rightEdges.back()) / 2);
                // Skip the samples that can't be a maximum
                i = iAhead;
            }
        }
        i++;
    }
    // std::cout<<"findLocalMaxima1D 4"<<std::endl;
    // for (const auto& el : leftEdges) {
        // std::cout << el << " ";
    // }
    // std::cout<<std::endl;
    // std::cout<<"findLocalMaxima1D 5"<<std::endl;
    // for (const auto& el : rightEdges) {
        // std::cout << el << " ";
    // }
    // std::cout<<std::endl;

    // std::cout<<"findLocalMaxima1D 6"<<std::endl;

    // for (const auto& el : midpoints) {
    //     std::cout << el << " ";
    // }
    // std::cout<<std::endl;
    // std::cout<<"findLocalMaxima1D 7"<<std::endl;

    return {midpoints, leftEdges, rightEdges};
}



// Assuming all the provided C++ functions and structures are defined above this point

std::tuple<std::vector<int>, std::map<std::string, std::vector<double>>> filterPeaksByProperty(const std::vector<int>& peaks, const std::vector<double>& properties, double minValue, double maxValue, std::map<std::string, std::vector<double>>& all_properties_for_all_points) {
    std::map<std::string, std::vector<double>> filtered_properties;
    std::vector<int> filteredPeaks;
    for (size_t i = 0; i < peaks.size(); ++i) {
        if ((std::isnan(minValue) || properties[i] >= minValue) && (std::isnan(maxValue) || properties[i] <= maxValue)) {
            filteredPeaks.push_back(peaks[i]);
            // map iterator created 
            // iterator pointing to start of map 
            std::map<std::string, std::vector<double>>::iterator it = all_properties_for_all_points.begin(); 
        
            // Iterating over the map using Iterator till map end. 
            while (it != all_properties_for_all_points.end()) { 
                // Accessing the key 
                std::string word = it->first; 
                // std::cout<<"Getting properties and working on "<<word<<std::endl;
                // Accessing the value 
                std::vector<double> vector = it->second; 

                filtered_properties[word].push_back(vector[i]);
                // std::cout << word << " size is " <<  filtered_properties[word].size() << std::endl; 
                // iterator incremented to point next item 
                it++; 
            } 
            // std::cout<<"on to the next"<<std::endl;
        }
    }
    return {filteredPeaks, filtered_properties};
}

std::pair<std::vector<int>, std::map<std::string, std::vector<double>>> findPeaks(
    const std::vector<double>& x,
    double height,
    double threshold,
    int distance,
    double prominence,
    double width,
    int wlen,
    double relHeight
) {

    std::tuple<std::vector<int>, std::map<std::string, std::vector<double>>> filtered_data;


    if (distance < 1) {
        throw std::invalid_argument("`distance` must be greater or equal to 1");
    }

    auto [discovered_peaks, left_edges, right_edges]= findLocalMaxima1D(x);

    // std::cout<<"FOUND LOCAL MAXMIMA"<<std::endl;
    if (discovered_peaks.size() == 0){
        // std::cout<<"No peaks found"<<std::endl;
    } else{
        // std::cout<<"first found peak is"<<discovered_peaks[0]<<std::endl;
    }

    // Evaluate distance condition if required
    if (distance >= 1) {
        std::vector<double> peakPriorities(discovered_peaks.size());
        for (size_t i = 0; i < discovered_peaks.size(); ++i) {
            peakPriorities[i] = x[discovered_peaks[i]]; // Example: priority based on peak height
            // std::cout<<x[discovered_peaks[i]]<<std::endl;
        }
        auto keep = selectByPeakDistance(discovered_peaks, peakPriorities, static_cast<double>(distance));
        
        // Filter peaks based on the 'keep' mask
        std::vector<int> filteredPeaks;
        for (size_t i = 0; i < keep.size(); ++i) {
            if (keep[i]) filteredPeaks.push_back(discovered_peaks[i]);
        }
        discovered_peaks = filteredPeaks;
    }

    // std::cout<<"FOUND discovered_peaks are"<<discovered_peaks[0]<<std::endl;
    std::map<std::string, std::vector<double>> properties;
    // Evaluate conditions based on the provided optional arguments
    if (height) {
        // std::cout<<"FILTERING WITH HEIGHT"<<std::endl;
        if (discovered_peaks.size() == 0){
            // std::cout<<"pre filtering height No peaks found"<<std::endl;
        } else{
            // std::cout<<"pre filtering height first found peak is"<<discovered_peaks[0]<<std::endl;
        }
        std::vector<double> peakHeights(discovered_peaks.size());
        for (size_t i = 0; i < discovered_peaks.size(); ++i) {
            peakHeights[i] = x[discovered_peaks[i]];
        }
        properties["peak_heights"] = peakHeights;
        filtered_data = filterPeaksByProperty(discovered_peaks, peakHeights, height, std::numeric_limits<double>::max(), properties);
        discovered_peaks = std::get<0>(filtered_data);
        properties = std::get<1>(filtered_data);
        if (discovered_peaks.size() == 0){
            // std::cout<<"post filtering height No peaks found after filtering by height"<<std::endl;
        } else{
            // std::cout<<"post filtering height first found peak is"<<discovered_peaks[0]<<std::endl;
        }
    }

    // Additional conditions such as threshold, prominence, width can be evaluated similarly
    // using _select_by_property, _select_by_peak_threshold, and other helper functions
    // For simplicity, these steps are omitted here

    if (prominence || width) {
        // std::cout<<"Calculating WITH prominence"<<std::endl;
        // if (discovered_peaks.size() == 0){
            // std::cout<<" pre Calculating WITH prominence No peaks found after filtering by height"<<std::endl;
        // } else{
            // std::cout<<"pre Calculating WITH prominence first found peak is"<<discovered_peaks[0]<<std::endl;
        // }
        // Calculate prominence if required
        auto [prominences, leftBases, rightBases] = calculate_peak_prominences(x, discovered_peaks, wlen);

        // std::cout<<"Prominences";
        // for (auto& el : prominences){
        //     std::cout<<el<<" ";
        // }
        // std::cout<<std::endl;
        // std::cout<<"leftBases";
        // for (auto& el : leftBases){
        //     std::cout<<el<<" ";
        // }
        // std::cout<<std::endl;


        if (prominence) {
            // std::cout<<"filtereing WITH prominence"<<std::endl;
            properties["prominences"] = prominences;
            properties["left_bases"] = std::vector<double>(leftBases.begin(), leftBases.end());
            properties["right_bases"] = std::vector<double>(rightBases.begin(), rightBases.end());
            // for (auto& el : properties["left_bases"]){
            //     std::cout<<el<<" ";
            // }
            // std::cout<<std::endl<<std::endl;
            filtered_data = filterPeaksByProperty(discovered_peaks, prominences, prominence, std::numeric_limits<double>::max(), properties);
            discovered_peaks = std::get<0>(filtered_data);
            properties = std::get<1>(filtered_data);




            // std::map<std::string, std::vector<double>>::iterator it = properties.begin(); 
        
            // // Iterating over the map using Iterator till map end. 
            // while (it != properties.end()) { 
            //     // Accessing the key 
            //     std::string word = it->first; 
            //     // std::cout<<"Getting properties and working on "<<word<<std::endl;
            //     // Accessing the value 
            //     std::cout << word << " size is " <<  properties[word].size() << std::endl; 
            //     // iterator incremented to point next item 
            //     it++; 
            // } 











            // std::cout<<"leftBases";
            // for (auto& el : properties["left_bases"]){
            //     std::cout<<el<<" ";
            // }
            // std::cout<<std::endl<<std::endl;
            // if (discovered_peaks.size() == 0){
            //     std::cout<<" post filtering WITH prominence No peaks found after filtering by height"<<std::endl;
            // } else{
            //     std::cout<<" post filtering WITH prominence first found peak is"<<discovered_peaks[0]<<std::endl;
            // }
        }

        // Calculate widths if required
        if (width) {
            // std::cout<<"filtereing WITH width"<<std::endl;
            // std::cout<<"filtereing WITH width 0"<<std::endl;

            auto [widths, widthHeights, leftIps, rightIps] = calculate_peak_widths(x, discovered_peaks, relHeight, properties["prominences"], properties["left_bases"], properties["right_bases"]);
            // std::cout<<"filtereing WITH width 1"<<std::endl;
            // if (discovered_peaks.size() == 0){
            //     std::cout<<" pre filtering WITH width No peaks found after filtering by height"<<std::endl;
            // } else{
            //     std::cout<<"pre filtering WITH width first found peak is"<<discovered_peaks[0]<<std::endl;
            // }
            properties["widths"] = widths;
            // std::cout<<"filtereing WITH width 3"<<std::endl;
            properties["width_heights"] = widthHeights;
            // std::cout<<"filtereing WITH width 4"<<std::endl;
            properties["left_ips"] = leftIps;
            // std::cout<<"filtereing WITH width 5"<<std::endl;
            properties["right_ips"] = rightIps;


            filtered_data = filterPeaksByProperty(discovered_peaks, widths, width, std::numeric_limits<double>::max(), properties);
            discovered_peaks = std::get<0>(filtered_data);
            properties = std::get<1>(filtered_data);
            // std::cout<<"filtereing WITH width 2"<<std::endl;
            // if (discovered_peaks.size() == 0){
            //     std::cout<<" post filtering WITH width No peaks found after filtering by height"<<std::endl;
            // } else{
            //     std::cout<<"post filtering WITH width first found peak is"<<discovered_peaks[0]<<std::endl;
            // }
        }
    }
    // std::cout<<"RETURNING FROM MAIN"<<std::endl;

    // for (auto &el : discovered_peaks){
    //     std::cout<<el<<" ";
    // }
    // std::cout<<std::endl;
    return {discovered_peaks, properties};
}