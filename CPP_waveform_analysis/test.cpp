// #include <chrono>
// using namespace std::chrono;


// extern "C" {
//     double* calculate_SI(double* pleth_signal, double* return_array, int array_length, double patient_height_meters, double sampling_freq){
//         auto start = high_resolution_clock::now();
//         // #TODO 
//         // convert input pleth signal and return array to be vectors isntead, see if thats possible in ctypes and numpy
// // 
//         // converting to vector
//         // https://stackoverflow.com/questions/8777603/what-is-the-simplest-way-to-convert-array-to-vector
//         std::vector<double> pleth_signal_vector(pleth_signal, pleth_signal + array_length);

//         // TIMING
//         std::cout<<"allocate pleth_signal_vector "<<duration_cast<microseconds>(high_resolution_clock::now() - start).count()<<std::endl;
//         start = high_resolution_clock::now();
//         // getting systolic peaks
//         auto [systolic_peaks, systolic_properties] = findPeaks(pleth_signal_vector,
//                     std::numeric_limits<double>::quiet_NaN(),   //height
//                     std::numeric_limits<double>::quiet_NaN(),   //threshold
//                     10,                                         // distance
//                     0.15,                                       //prominence
//                     std::numeric_limits<double>::quiet_NaN(),   //width
//                     std::numeric_limits<int>::quiet_NaN(),      //wlen
//                     0.5);  //relHeight
//         // TIMING
//         std::cout<<"Time for systolic peaks "<<duration_cast<microseconds>(high_resolution_clock::now() - start).count()<<std::endl;
//         start = high_resolution_clock::now();
//         // getting potential diastolic peaks
//         auto pleth_gradient = calculateGradient(pleth_signal_vector);
//         std::cout<<"Time for get pleth gradient "<<duration_cast<microseconds>(high_resolution_clock::now() - start).count()<<std::endl;
//         start = high_resolution_clock::now();
//         auto [potential_diastolic_peaks, potential_diastolic_properties] = findPeaks(pleth_gradient,
//                     std::numeric_limits<double>::quiet_NaN(),   //height
//                     std::numeric_limits<double>::quiet_NaN(),   //threshold
//                     25,                                         // distance
//                     std::numeric_limits<double>::quiet_NaN(),                                       //prominence
//                     std::numeric_limits<double>::quiet_NaN(),   //width
//                     std::numeric_limits<int>::quiet_NaN(),      //wlen
//                     0.5);  //relHeight
//         std::cout<<"Time for potential diastolic peaks "<<duration_cast<microseconds>(high_resolution_clock::now() - start).count()<<std::endl;
//         start = high_resolution_clock::now();


//         // getting rid of false peaks from the gradient peaks to get true diastolic peaks
//         std::vector<int> diastolic_peaks;
//         size_t dias_index = 0;
        
//         for (size_t sys_index = 0; sys_index < systolic_peaks.size(); ++sys_index) {
//             int sys_peak = systolic_peaks[sys_index];
//             while (dias_index < potential_diastolic_peaks.size() - 1 &&
//                 potential_diastolic_peaks[dias_index] < sys_peak + 10) {
//                 ++dias_index;
//             }
            
//             if (dias_index < potential_diastolic_peaks.size()) {
//                 diastolic_peaks.push_back(potential_diastolic_peaks[dias_index]);
//             }
//         }
//         std::cout<<"Time for filter diastolic "<<duration_cast<microseconds>(high_resolution_clock::now() - start).count()<<std::endl;
//         start = high_resolution_clock::now();
//         // CALCULATE SI VALUE
//         // Initial conditions
//         int prev_peak = 0;
//         double SI_value = 0.0;

//         for (size_t i = 0; i < systolic_peaks.size(); ++i) {
//             // Set the SI value for the range from prev_peak to the current diastolic peak
//             for (int j =prev_peak; j<diastolic_peaks[i]; j++){
//                 return_array[j] = SI_value;
//             }

//             // Calculate the new SI value
//             SI_value = patient_height_meters / ((diastolic_peaks[i] - systolic_peaks[i]) / sampling_freq);

//             // Update prev_peak for the next iteration
//             prev_peak = diastolic_peaks[i];
//         }
//         std::cout<<"Time for calculate SI values "<<duration_cast<microseconds>(high_resolution_clock::now() - start).count()<<std::endl;
//         start = high_resolution_clock::now();
//         // Set the remaining elements of the return_array to the last SI_value
//         for (int j =prev_peak; j<array_length; j++){
//             return_array[j] = SI_value;
//         }
        
//         std::cout<<"Time for fill rest of array "<<duration_cast<microseconds>(high_resolution_clock::now() - start).count()<<std::endl;
//         start = high_resolution_clock::now();

//         return return_array;
//     };



//     double* calculate_SI_from_peaks(double* systolic_peaks, double* diastolic_peaks, double* return_array, int return_array_length, int num_peaks,  double patient_height_meters, double sampling_freq) {
//         // CALCULATE SI VALUE
//         // Initial conditions
//         int prev_peak = 0;
//         double SI_value = 0.0;

//         for (size_t i = 0; i < num_peaks; ++i) {
//             // Set the SI value for the range from prev_peak to the current diastolic peak
//             for (int j =prev_peak; j<diastolic_peaks[i]; j++){
//                 return_array[j] = SI_value;
//             }

//             // Calculate the new SI value
//             SI_value = patient_height_meters / ((diastolic_peaks[i] - systolic_peaks[i]) / sampling_freq);

//             // Update prev_peak for the next iteration
//             prev_peak = diastolic_peaks[i];
//         }

//         for (int j =prev_peak; j<return_array_length; j++){
//             return_array[j] = SI_value;
//         }

//         return return_array;
//     }
// }



