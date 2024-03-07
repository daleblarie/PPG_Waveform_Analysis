// waveform_analysis.cpp
#include <vector>
#include <stdexcept>
#include <iostream>
#include <functional>
#include <map>
#include <limits>
#include <cmath>
#include <tuple>

#include "find_peaks.h"
#include "gradient.h"


#include <chrono>
using namespace std::chrono;






extern "C" {
    // extern double _calculate_PI(const std::vector<double>& pleth_signal);
    double* gradient(double* input_arr, double* return_arr, int input_array_length){
        std::vector<double> input_arr_vector(input_arr, input_arr + input_array_length);
        auto signal_gradient = calculateGradient(input_arr_vector);
        int i =0;
        for (auto val: signal_gradient){
            return_arr[i] = val;
            i++;
        }
        
        return return_arr;
    }

    int32_t* filter_diastolic_peaks(int32_t* systolic_peaks, int systolic_peaks_length, int32_t* potential_diastolic_peaks, int potential_diastolic_peaks_length, int32_t* return_array) {
           // getting rid of false peaks from the gradient peaks to get true diastolic peaks
        std::vector<int> diastolic_peaks;
        size_t dias_index = 0;
        
        for (size_t sys_index = 0; sys_index < systolic_peaks_length; ++sys_index) {
            int sys_peak = systolic_peaks[sys_index];
            while (dias_index < potential_diastolic_peaks_length - 1 &&
                potential_diastolic_peaks[dias_index] < sys_peak + 10) {
                ++dias_index;
            }
            
            if (dias_index < potential_diastolic_peaks_length) {
                diastolic_peaks.push_back(potential_diastolic_peaks[dias_index]);
            }
        }

        for (int i=0; i<systolic_peaks_length; i++){
            return_array[i] = diastolic_peaks[i];
        }

        // std::cout<<"Time for filter diastolic "<<duration_cast<microseconds>(high_resolution_clock::now() - start).count()<<std::endl;
        return return_array;
    }


    int32_t* find_dicrotic_notches(double* pleth_signal, int32_t* systolic_peaks, int systolic_peaks_length, int32_t* diastolic_peaks, int diastolic_peaks_length, int32_t* return_array){
        // Calculate dicrotic notch
        // std::cout<<"diastolic_peaks_length"<<std::endl;
        // std::cout<<diastolic_peaks_length<<std::endl;
        // std::cout<<systolic_peaks_length<<std::endl;
        std::vector<int> dicrotic_notch;
        for (size_t i = 0; i < systolic_peaks_length; ++i) {
            // std::cout<<"1 "<<i<<std::endl;
            int sys_peak = systolic_peaks[i];
            int dias_peak = diastolic_peaks[i];
            // std::cout<<"2 "<<i<<std::endl;

            int n_points = dias_peak - sys_peak;
            // check to make sure that our matching isnt messed 
            if (n_points <=0 ){ continue;}
            // std::cout<<"2.5 "<<n_points<<std::endl;
            auto subtraction_line = std::vector<double>(n_points);
            // std::cout<<"3 "<<i<<std::endl;

            // Manually calculate the subtraction line
            double slope = (pleth_signal[dias_peak] - pleth_signal[sys_peak]) / n_points;
            for (int j = 0; j < n_points; ++j) {
                subtraction_line[j] = pleth_signal[sys_peak] + slope * j;
            }
            // std::cout<<"4 "<<i<<std::endl;

            auto sys_to_dias_window = std::vector<double>(pleth_signal + sys_peak, pleth_signal + dias_peak);
            auto min_element_iter = std::min_element(sys_to_dias_window.begin(), sys_to_dias_window.end());
            int notch = std::distance(sys_to_dias_window.begin(), min_element_iter) + sys_peak;
            // std::cout<<"5 "<<i<<std::endl;

            dicrotic_notch.push_back(notch);
        }
        // std::cout<<"dicrotic_notch filled"<<std::endl;

        for(int i=0; i< dicrotic_notch.size(); i++){
            return_array[i] = dicrotic_notch[i];
        }
        return return_array;
    }


    double* calculate_SI_from_peaks(int32_t* systolic_peaks, int32_t* diastolic_peaks, int num_peaks, double* return_array, int return_array_length,   double patient_height_meters, double sampling_freq) {
        auto start = high_resolution_clock::now();

        // CALCULATE SI VALUE
        // Initial conditions
        int prev_peak = 0;
        double SI_value = 0.0;

        for (size_t i = 0; i < num_peaks; ++i) {
            // Set the SI value for the range from prev_peak to the current diastolic peak
            for (int j =prev_peak; j<diastolic_peaks[i]; j++){
                return_array[j] = SI_value;
            }

            // Calculate the new SI value
            SI_value = patient_height_meters / ((diastolic_peaks[i] - systolic_peaks[i]) / sampling_freq);

            // Update prev_peak for the next iteration
            prev_peak = diastolic_peaks[i];
        }
        for (int j =prev_peak; j<return_array_length; j++){
            return_array[j] = SI_value;
        }

        // std::cout<<"Time for calculate SI values "<<duration_cast<microseconds>(high_resolution_clock::now() - start).count()<<std::endl;
        start = high_resolution_clock::now();
        return return_array;
    }

    double* calculate_RI_from_peaks(double* pleth_signal, int32_t* systolic_peaks, int32_t* diastolic_peaks, int num_peaks, double* return_array, int return_array_length) {
        // Calculate RI value
        int prev_peak = 0;
        double RI_value = 0.0;

        for (size_t i = 0; i < num_peaks; ++i) {
            for (int j = prev_peak; j < diastolic_peaks[i]; ++j) {
                return_array[j] = RI_value;
            }
            
            RI_value = 100 * (pleth_signal[diastolic_peaks[i]] / pleth_signal[systolic_peaks[i]]);
            prev_peak = diastolic_peaks[i];
        }

        for (int j = prev_peak; j < return_array_length; ++j) {
            return_array[j] = RI_value;
        }

        return return_array;
    }


    double* calculate_AI_from_peaks(double* pleth_signal, int32_t* systolic_peaks, int32_t* diastolic_peaks, int num_peaks, double* return_array, int return_array_length) {

        auto start = high_resolution_clock::now();

            // Calculate AI value
        int prev_peak = 0;
        double AI_value = 0.0;

        for (size_t i = 0; i < num_peaks; ++i) {
            for (int j = prev_peak; j < diastolic_peaks[i]; ++j) {
                return_array[j] = AI_value;
            }
            
            AI_value = 100 * ((pleth_signal[systolic_peaks[i]] - pleth_signal[diastolic_peaks[i]]) / pleth_signal[systolic_peaks[i]]);
            prev_peak = diastolic_peaks[i];
        }

        for (int j = prev_peak; j < return_array_length; ++j) {
            return_array[j] = AI_value;
        }
        // std::cout<<"Time for calculate AI values "<<duration_cast<microseconds>(high_resolution_clock::now() - start).count()<<std::endl;
        start = high_resolution_clock::now();
        return return_array;
    }


    double* calculate_IPA_from_peaks(double* pleth_signal, int32_t* systolic_peaks, int num_systolic_peaks, int32_t* dicrotic_notch, double* return_array, int return_array_length){


        std::vector<double> pleth_signal_vector(pleth_signal, pleth_signal + return_array_length);

        std::vector<int> sys_peak_vector(systolic_peaks, systolic_peaks + num_systolic_peaks);

        std::vector<int> notch_vector(dicrotic_notch, dicrotic_notch + num_systolic_peaks);

        double ipa_value = -1;
        int prev_trough = 0;
        // std::cout<<std::endl;
        for (int i = 0; i < sys_peak_vector.size() - 1; i++) {
            int sys_peak = sys_peak_vector[i];
            int next_sys_peak = sys_peak_vector[i+1];
            int notch = dicrotic_notch[i];


            int trough = std::distance(pleth_signal_vector.begin() + sys_peak,
                std::min_element(pleth_signal_vector.begin() + sys_peak, pleth_signal_vector.begin() + next_sys_peak)) + sys_peak;

            for (int j = prev_trough; j < trough; j++) {
                if (j > return_array_length){std::cout<<"WE ARE GOING PAS OUR ARRAY BOUNDS"<<std::endl;}
                return_array[j] = ipa_value;
            }
   
            double A1 = 0.0;
            for (int idx = prev_trough; idx < notch; ++idx) {
                A1 += pleth_signal_vector[idx];
            }
            double A2 = 0.0;
            for (int idx = notch; idx < trough; ++idx) {
                A2 += pleth_signal_vector[idx];
            }
      
            if (A1 == 0) {
                ipa_value = -1;
            } else {
                ipa_value = A2 / A1;
            }
            prev_trough = trough;
            // signal_troughs[i] = trough;
        }

        if (prev_trough >= return_array_length){std::cout<<"PREV TROUGH IS PAST INDEX"<<std::endl;};
        for (int j = prev_trough; j < return_array_length; j++) {
            if (j > return_array_length){std::cout<<"WE ARE GOING PAS OUR ARRAY BOUNDS"<<std::endl;}
            return_array[j] = ipa_value;
        }

        return return_array;
    }


    double* calculate_SI(double* pleth_signal, double* return_array, int array_length, double patient_height_meters, double sampling_freq){
        auto start = high_resolution_clock::now();
        // #TODO 
        // convert input pleth signal and return array to be vectors isntead, see if thats possible in ctypes and numpy
// 
        // converting to vector
        // https://stackoverflow.com/questions/8777603/what-is-the-simplest-way-to-convert-array-to-vector
        std::vector<double> pleth_signal_vector(pleth_signal, pleth_signal + array_length);

        // TIMING
        // std::cout<<"allocate pleth_signal_vector "<<duration_cast<microseconds>(high_resolution_clock::now() - start).count()<<std::endl;
        start = high_resolution_clock::now();
        // getting systolic peaks
        auto [systolic_peaks, systolic_properties] = findPeaks(pleth_signal_vector,
                    std::numeric_limits<double>::quiet_NaN(),   //height
                    std::numeric_limits<double>::quiet_NaN(),   //threshold
                    10,                                         // distance
                    0.15,                                       //prominence
                    std::numeric_limits<double>::quiet_NaN(),   //width
                    std::numeric_limits<int>::quiet_NaN(),      //wlen
                    0.5);  //relHeight
        // TIMING
        // std::cout<<"Time for systolic peaks "<<duration_cast<microseconds>(high_resolution_clock::now() - start).count()<<std::endl;
        start = high_resolution_clock::now();
        // getting potential diastolic peaks
        auto pleth_gradient = calculateGradient(pleth_signal_vector);
        // std::cout<<"Time for get pleth gradient "<<duration_cast<microseconds>(high_resolution_clock::now() - start).count()<<std::endl;
        start = high_resolution_clock::now();
        auto [potential_diastolic_peaks, potential_diastolic_properties] = findPeaks(pleth_gradient,
                    std::numeric_limits<double>::quiet_NaN(),   //height
                    std::numeric_limits<double>::quiet_NaN(),   //threshold
                    25,                                         // distance
                    std::numeric_limits<double>::quiet_NaN(),                                       //prominence
                    std::numeric_limits<double>::quiet_NaN(),   //width
                    std::numeric_limits<int>::quiet_NaN(),      //wlen
                    0.5);  //relHeight
        // std::cout<<"Time for potential diastolic peaks "<<duration_cast<microseconds>(high_resolution_clock::now() - start).count()<<std::endl;
        start = high_resolution_clock::now();


        // getting rid of false peaks from the gradient peaks to get true diastolic peaks
        std::vector<int> diastolic_peaks;
        size_t dias_index = 0;
        
        for (size_t sys_index = 0; sys_index < systolic_peaks.size(); ++sys_index) {
            int sys_peak = systolic_peaks[sys_index];
            while (dias_index < potential_diastolic_peaks.size() - 1 &&
                potential_diastolic_peaks[dias_index] < sys_peak + 10) {
                ++dias_index;
            }
            
            if (dias_index < potential_diastolic_peaks.size()) {
                diastolic_peaks.push_back(potential_diastolic_peaks[dias_index]);
            }
        }
        // std::cout<<"Time for filter diastolic "<<duration_cast<microseconds>(high_resolution_clock::now() - start).count()<<std::endl;
        start = high_resolution_clock::now();
        // CALCULATE SI VALUE
        // Initial conditions
        int prev_peak = 0;
        double SI_value = 0.0;

        for (size_t i = 0; i < systolic_peaks.size(); ++i) {
            // Set the SI value for the range from prev_peak to the current diastolic peak
            for (int j =prev_peak; j<diastolic_peaks[i]; j++){
                return_array[j] = SI_value;
            }

            // Calculate the new SI value
            SI_value = patient_height_meters / ((diastolic_peaks[i] - systolic_peaks[i]) / sampling_freq);

            // Update prev_peak for the next iteration
            prev_peak = diastolic_peaks[i];
        }
        // std::cout<<"Time for calculate SI values "<<duration_cast<microseconds>(high_resolution_clock::now() - start).count()<<std::endl;
        start = high_resolution_clock::now();
        // Set the remaining elements of the return_array to the last SI_value
        for (int j =prev_peak; j<array_length; j++){
            return_array[j] = SI_value;
        }
        
        // std::cout<<"Time for fill rest of array "<<duration_cast<microseconds>(high_resolution_clock::now() - start).count()<<std::endl;
        start = high_resolution_clock::now();
        // // RETURNING GRADIENT PEAKS
        // int j=0;
        // for (const auto& ind : diastolic_peaks){
        //     return_array[j] = ind;
        //     j++;
        //     // std::cout<<ind<<std::endl;
        // }
        // // std::cout<<"THESE ARE ALL THE INDECES"<<std::endl;

        // for (int i = j; i < array_length; i++){
        //     return_array[i] = -1;
        // }

        // // RETURNING POTENTIAL GRADIENT PEAKS
        // int j=0;
        // for (const auto& ind : potential_diastolic_peaks){
        //     return_array[j] = ind;
        //     j++;
        //     // std::cout<<ind<<std::endl;
        // }
        // // std::cout<<"THESE ARE ALL THE INDECES"<<std::endl;

        // for (int i = j; i < array_length; i++){
        //     return_array[i] = -1;
        // }
        
        // RETURNING GRADIENT
        // int j = 0;
        // for (const auto& ind : pleth_gradient){
        //     return_array[j] = ind;
        //     j++;
        //     // std::cout<<ind<<std::endl;
        // }
        // for (int i = j; i < array_length; i++){
        //     return_array[i] = -1;
        // }

        return return_array;
    };
    // Assuming calculateGradient and findPeaks functions are similar to calculate_SI example

    double* calculate_RI(double* pleth_signal, double* return_array, int array_length, double sampling_freq) {
        // Convert input array to vector
        std::vector<double> pleth_signal_vector(pleth_signal, pleth_signal + array_length);

                // getting systolic peaks
        auto [systolic_peaks, systolic_properties] = findPeaks(pleth_signal_vector,
                    std::numeric_limits<double>::quiet_NaN(),   //height
                    std::numeric_limits<double>::quiet_NaN(),   //threshold
                    10,                                         // distance
                    0.15,                                       //prominence
                    std::numeric_limits<double>::quiet_NaN(),   //width
                    std::numeric_limits<int>::quiet_NaN(),      //wlen
                    0.5);  //relHeight

        // getting potential diastolic peaks
        auto pleth_gradient = calculateGradient(pleth_signal_vector);
        auto [potential_diastolic_peaks, potential_diastolic_properties] = findPeaks(pleth_gradient,
                    std::numeric_limits<double>::quiet_NaN(),   //height
                    std::numeric_limits<double>::quiet_NaN(),   //threshold
                    25,                                         // distance
                    std::numeric_limits<double>::quiet_NaN(),   //prominence
                    std::numeric_limits<double>::quiet_NaN(),   //width
                    std::numeric_limits<int>::quiet_NaN(),      //wlen
                    0.5);  //relHeight



        std::vector<int> diastolic_peaks;
        size_t dias_index = 0;

        for (size_t sys_index = 0; sys_index < systolic_peaks.size(); ++sys_index) {
            int sys_peak = systolic_peaks[sys_index];
            while (dias_index < potential_diastolic_peaks.size() - 1 && potential_diastolic_peaks[dias_index] < sys_peak + 10) {
                ++dias_index;
            }
            
            if (dias_index < potential_diastolic_peaks.size()) {
                diastolic_peaks.push_back(potential_diastolic_peaks[dias_index]);
            }
        }

        // Calculate RI value
        int prev_peak = 0;
        double RI_value = 0.0;

        for (size_t i = 0; i < systolic_peaks.size(); ++i) {
            for (int j = prev_peak; j < diastolic_peaks[i]; ++j) {
                return_array[j] = RI_value;
            }
            
            RI_value = 100 * (pleth_signal_vector[diastolic_peaks[i]] / pleth_signal_vector[systolic_peaks[i]]);
            prev_peak = diastolic_peaks[i];
        }

        for (int j = prev_peak; j < array_length; ++j) {
            return_array[j] = RI_value;
        }

        return return_array;
    }

    double* calculate_AI(double* pleth_signal, double* return_array, int array_length) {
        // Convert input array to vector
        std::vector<double> pleth_signal_vector(pleth_signal, pleth_signal + array_length);

        // getting systolic peaks
        auto [systolic_peaks, systolic_properties] = findPeaks(pleth_signal_vector,
                    std::numeric_limits<double>::quiet_NaN(),   //height
                    std::numeric_limits<double>::quiet_NaN(),   //threshold
                    10,                                         // distance
                    0.15,                                       //prominence
                    std::numeric_limits<double>::quiet_NaN(),   //width
                    std::numeric_limits<int>::quiet_NaN(),      //wlen
                    0.5);  //relHeight

        // getting potential diastolic peaks
        auto pleth_gradient = calculateGradient(pleth_signal_vector);
        auto [potential_diastolic_peaks, potential_diastolic_properties] = findPeaks(pleth_gradient,
                    std::numeric_limits<double>::quiet_NaN(),   //height
                    std::numeric_limits<double>::quiet_NaN(),   //threshold
                    25,                                         // distance
                    std::numeric_limits<double>::quiet_NaN(),   //prominence
                    std::numeric_limits<double>::quiet_NaN(),   //width
                    std::numeric_limits<int>::quiet_NaN(),      //wlen
                    0.5);  //relHeight

        std::vector<int> diastolic_peaks;
        size_t dias_index = 0;

        for (size_t sys_index = 0; sys_index < systolic_peaks.size(); ++sys_index) {
            int sys_peak = systolic_peaks[sys_index];
            while (dias_index < potential_diastolic_peaks.size() - 1 && potential_diastolic_peaks[dias_index] < sys_peak + 10) {
                ++dias_index;
            }
            
            if (dias_index < potential_diastolic_peaks.size()) {
                diastolic_peaks.push_back(potential_diastolic_peaks[dias_index]);
            }
        }

        // Calculate AI value
        int prev_peak = 0;
        double AI_value = 0.0;

        for (size_t i = 0; i < systolic_peaks.size(); ++i) {
            for (int j = prev_peak; j < diastolic_peaks[i]; ++j) {
                return_array[j] = AI_value;
            }
            
            AI_value = 100 * ((pleth_signal_vector[systolic_peaks[i]] - pleth_signal_vector[diastolic_peaks[i]]) / pleth_signal_vector[systolic_peaks[i]]);
            prev_peak = diastolic_peaks[i];
        }

        for (int j = prev_peak; j < array_length; ++j) {
            return_array[j] = AI_value;
        }

        return return_array;
    }

    double* calculate_IPA(double* pleth_signal, double* return_array, int array_length, double sampling_freq) {
        // Convert input array to vector
        std::vector<double> pleth_signal_vector(pleth_signal, pleth_signal + array_length);

        // getting systolic peaks
        auto [systolic_peaks, systolic_properties] = findPeaks(pleth_signal_vector,
                    std::numeric_limits<double>::quiet_NaN(),   //height
                    std::numeric_limits<double>::quiet_NaN(),   //threshold
                    10,                                         // distance
                    0.15,                                       //prominence
                    std::numeric_limits<double>::quiet_NaN(),   //width
                    std::numeric_limits<int>::quiet_NaN(),      //wlen
                    0.5);  //relHeight

        // getting potential diastolic peaks
        auto pleth_gradient = calculateGradient(pleth_signal_vector);
        auto [potential_diastolic_peaks, potential_diastolic_properties] = findPeaks(pleth_gradient,
                    std::numeric_limits<double>::quiet_NaN(),   //height
                    std::numeric_limits<double>::quiet_NaN(),   //threshold
                    25,                                         // distance
                    std::numeric_limits<double>::quiet_NaN(),   //prominence
                    std::numeric_limits<double>::quiet_NaN(),   //width
                    std::numeric_limits<int>::quiet_NaN(),      //wlen
                    0.5);  //relHeight
        // std::cout<<"Got gradients"<<std::endl;
        std::vector<int> diastolic_peaks;
        size_t dias_index = 0;

        // Filter false peaks from gradient peaks to get true diastolic peaks
        for (const auto& sys_peak : systolic_peaks) {
            while (dias_index < potential_diastolic_peaks.size() - 1 && potential_diastolic_peaks[dias_index] < sys_peak + 10) {
                ++dias_index;
            }
            
            if (dias_index < potential_diastolic_peaks.size()) {
                diastolic_peaks.push_back(potential_diastolic_peaks[dias_index]);
            }
        }
        // std::cout<<"Filtered False Peaks"<<std::endl;


        // Calculate dicrotic notch
        std::vector<int> dicrotic_notch;
        for (size_t i = 0; i < systolic_peaks.size(); ++i) {
            // std::cout<<"1 "<<i<<std::endl;
            int sys_peak = systolic_peaks[i];
            int dias_peak = diastolic_peaks[i];
            // std::cout<<"2 "<<i<<std::endl;

            int n_points = dias_peak - sys_peak;
            // check to make sure that our matching isnt messed 
            if (n_points <=0 ){ continue;}
            // std::cout<<"2.5 "<<n_points<<std::endl;
            auto subtraction_line = std::vector<double>(n_points);
            // std::cout<<"3 "<<i<<std::endl;

            // Manually calculate the subtraction line
            double slope = (pleth_signal[dias_peak] - pleth_signal[sys_peak]) / n_points;
            for (int j = 0; j < n_points; ++j) {
                subtraction_line[j] = pleth_signal[sys_peak] + slope * j;
            }
            // std::cout<<"4 "<<i<<std::endl;

            auto sys_to_dias_window = std::vector<double>(pleth_signal + sys_peak, pleth_signal + dias_peak);
            auto min_element_iter = std::min_element(sys_to_dias_window.begin(), sys_to_dias_window.end());
            int notch = std::distance(sys_to_dias_window.begin(), min_element_iter) + sys_peak;
            // std::cout<<"5 "<<i<<std::endl;

            dicrotic_notch.push_back(notch);
        }
        // std::cout<<"Got dicrotic notches"<<std::endl;

        // std::vector<double> IPA_array(array_length, 0);
        // std::vector<int> signal_troughs(systolic_peaks.size(), 0);
        double ipa_value = -1;
        int prev_trough = 0;
        // std::cout<<"Got dicrotic 1notches"<<std::endl;
        // std::cout<<systolic_peaks.size()<<std::endl;
        int sys_peak = systolic_peaks[0];
        if (systolic_peaks.size() > 0){
            // Only do this is we have detected peaks
            for (size_t i = 0; i < systolic_peaks.size() - 1; i++) {
                // std::cout<<"i "<<i<<" "<<systolic_peaks.size() - 1<<std::endl;
                // std::cout<<"TEST 1.1"<<std::endl;
                sys_peak = systolic_peaks[i];
                // std::cout<<"TEST 1.2"<<std::endl;
                int next_sys_peak = systolic_peaks[i+1];
                // std::cout<<"TEST 1.3"<<std::endl;
                int notch = dicrotic_notch[i];

                // std::cout<<"TEST 1"<<std::endl;

                int trough = std::distance(pleth_signal_vector.begin() + sys_peak,
                    std::min_element(pleth_signal_vector.begin() + sys_peak, pleth_signal_vector.begin() + next_sys_peak)) + sys_peak;
                // std::cout<<"TEST 2"<<std::endl;

                for (int j = prev_trough; j < trough; j++) {
                    if (j > array_length){std::cout<<"WE ARE GOING PAS OUR ARRAY BOUNDS"<<std::endl;}
                    return_array[j] = ipa_value;
                }
                // std::cout<<"TEST 3"<<std::endl;
    
                double A1 = 0.0;
                for (int idx = prev_trough; idx < notch; ++idx) {
                    A1 += pleth_signal_vector[idx];
                }
                double A2 = 0.0;
                for (int idx = notch; idx < trough; ++idx) {
                    A2 += pleth_signal_vector[idx];
                }
                // std::cout<<"TEST 4"<<std::endl;
        
                if (A1 == 0) {
                    ipa_value = -1;
                } else {
                    ipa_value = A2 / A1;
                }
                prev_trough = trough;
                // signal_troughs[i] = trough;
            }
            // std::cout<<"Got FINISHED FILLING"<<std::endl;
        }

        if (prev_trough >= array_length){std::cout<<"PREV TROUGH IS PAST INDEX"<<std::endl;};
        for (int j = prev_trough; j < array_length; j++) {
            if (j > array_length){std::cout<<"WE ARE GOING PAS OUR ARRAY BOUNDS"<<std::endl;}
            return_array[j] = ipa_value;
        }

        return return_array;
    }

    double _calculate_PI(const std::vector<double>& pleth_signal) {
        if (pleth_signal.empty()) {
            std::cerr << "Error: Pleth signal is empty." << std::endl;
            return -1;
        }

        // Calculate minimum value (DC component) of the signal
        double min_value = *std::min_element(pleth_signal.begin(), pleth_signal.end());
        
        // Calculate sum of DC component
        double DC_sum = min_value * pleth_signal.size();
        if (DC_sum == 0) {
            // std::cerr << "Error: DC sum is zero, division by zero is not allowed." << std::endl;
            return -1;
        }
        
        // Calculate sum of AC components manually
        double AC_sum = 0.0;
        for (const double& current : pleth_signal) {
            AC_sum += (current - min_value);
        }
        
        // Calculate Perfusion Index (PI) value
        double PI_value = 100.0 * AC_sum / DC_sum;

        return PI_value;
    }

    double calculate_PI(double* pleth_signal, int array_length){
        std::vector<double> pleth_signal_vector(pleth_signal, pleth_signal + array_length);
        return (_calculate_PI(pleth_signal_vector));
    }

    void calculate_PVI(double* pleth_signal, double* return_array, int array_length, int window, int step) {
        // Convert input array to vector
        std::vector<double> pleth_signal_vector(pleth_signal, pleth_signal + array_length);
        double PI_min = std::numeric_limits<double>::infinity();
        double PI_max = 0;
        double current_PI;
        double pvi_value;

        int PVI_index = window;
        int start_index = 0;

        while (start_index + window < array_length) {
            std::vector<double> window_signal(pleth_signal_vector.begin() + start_index, pleth_signal_vector.begin() + start_index + window);
            current_PI = _calculate_PI(window_signal); // Ensure you have a suitable implementation for _calculate_PI that matches this context
            if (current_PI > PI_max) {
                PI_max = current_PI;
            }
            if (current_PI < PI_min && current_PI >= 0) { // Ensure PI_min doesn't consider error code (-1) from _calculate_PI
                PI_min = current_PI;
            }

            if  (current_PI == -1){ // the case where calculate_pi returns -1 (aka there is a 0 in the pleth wave)
                pvi_value = -1;
            } else{
                pvi_value = (PI_max == 0) ? -1 : 100 * (PI_max - PI_min) / PI_max;

            }
            for (int i = PVI_index; i < std::min(PVI_index + step, array_length); ++i) {
                return_array[i] = pvi_value;
            }
            // std::cout<<start_index<<", "<<PVI_index<<std::endl;;

            start_index += step;
            PVI_index += step;
        }
        // Fill the remaining part of the return_array with the last calculated PVI value
        std::fill(return_array + PVI_index-step, return_array + array_length, pvi_value);

        // This function doesn't return anything as it directly modifies the passed in return_array
    }


}



int main() {
    // Example usage
    std::vector<double> x = {1, 3, 2, 4, 1, 5, 6, 4};
    auto [peaks, properties] = findPeaks(x, 3, std::numeric_limits<double>::quiet_NaN(), 1, 1, std::numeric_limits<double>::quiet_NaN(), -1, 0.5);

    std::cout << "Peaks at indices: ";
    for (const auto& peak : peaks) {
        std::cout << peak << " ";
    }
    std::cout << std::endl;

    // Optionally print properties
    if (!properties["peak_heights"].empty()) {
        std::cout << "Peak heights: ";
        for (const auto& height : properties["peak_heights"]) {
            std::cout << height << " ";
        }
        std::cout << std::endl;
    }

    // Continue for other properties as needed...

    return 0;
}
