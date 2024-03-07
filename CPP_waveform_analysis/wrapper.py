
import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
import pandas as pd
import time

import scipy
import scipy.signal

import matplotlib.pyplot as plt

# compile command
# g++ -fPIC -Wall -shared -o Waveform_analysis/CPP_waveform_analysis/waveform_calculations.so Waveform_analysis/CPP_waveform_analysis/*.cpp -std=c++11 -O3

lib = ctypes.CDLL('Waveform_Analysis/CPP_waveform_analysis/waveform_calculations.so')

class PPG_METRICS_CTYPES:
    def __init__(self, pure_cpp=False):
        self.pure_cpp=pure_cpp  # determines which implementation to use for ctypes
        self.c_double_p = ctypes.POINTER(ctypes.c_double)
        self.c_int_p = ctypes.POINTER(ctypes.c_int)
        self.patient_height_meters_type = ctypes.c_double
        self.sampling_rate_type = ctypes.c_double
        self.array_length_type = ctypes.c_int
        self.PI_restype = ctypes.c_double
        self.window_type = ctypes.c_int
        self.step_type = ctypes.c_int
        self.num_peaks_type = ctypes.c_int
        
        lib.gradient.argtypes = (self.c_double_p, self.c_double_p, self.array_length_type)
        lib.gradient.restype = self.c_double_p

        lib.filter_diastolic_peaks.argtypes = (
            self.c_int_p, 
            self.array_length_type,
            self.c_int_p,
            self.array_length_type,
            self.c_int_p)
        lib.filter_diastolic_peaks.restype = self.c_int_p

        lib.find_dicrotic_notches.argtypes = (
            self.c_double_p, 
            self.c_int_p,
            self.num_peaks_type,
            self.c_int_p,
            self.num_peaks_type,
            self.c_int_p)
        lib.find_dicrotic_notches.restype = self.c_int_p


        lib.calculate_SI_from_peaks.argtypes = (
                self.c_int_p, 
                self.c_int_p,
                self.num_peaks_type,
                self.c_double_p,
                self.array_length_type,
                self.patient_height_meters_type,
                self.sampling_rate_type)
        lib.calculate_SI_from_peaks.restype = self.c_double_p

        lib.calculate_RI_from_peaks.argtypes = (
            self.c_double_p, 
            self.c_int_p, 
            self.c_int_p, 
            self.num_peaks_type,
            self.c_double_p,
            self.array_length_type,)
        lib.calculate_RI_from_peaks.restype = self.c_double_p

        lib.calculate_AI_from_peaks.argtypes = (
            self.c_double_p, 
            self.c_int_p, 
            self.c_int_p, 
            self.num_peaks_type,
            self.c_double_p,
            self.array_length_type,)
        lib.calculate_AI_from_peaks.restype = self.c_double_p

        lib.calculate_IPA_from_peaks.argtypes = (
            self.c_double_p, 
            self.c_int_p, 
            self.num_peaks_type,
            self.c_int_p, 
            self.c_double_p,
            self.array_length_type,)
        lib.calculate_IPA_from_peaks.restype = self.c_double_p

        lib.calculate_SI.argtypes = (self.c_double_p, self.c_double_p, self.array_length_type, self.patient_height_meters_type, self.sampling_rate_type)
        lib.calculate_SI.restype = self.c_double_p

        lib.calculate_RI.argtypes = (self.c_double_p, self.c_double_p, self.array_length_type, self.sampling_rate_type)
        lib.calculate_RI.restype = self.c_double_p

        lib.calculate_AI.argtypes = (self.c_double_p, self.c_double_p, self.array_length_type)
        lib.calculate_AI.restype = (self.c_double_p)

        lib.calculate_IPA.argtypes = (self.c_double_p, self.c_double_p, self.array_length_type, self.sampling_rate_type)
        lib.calculate_IPA.restype = self.c_double_p

        lib.calculate_PI.argtypes = (self.c_double_p, self.array_length_type)
        lib.calculate_PI.restype = self.PI_restype

        lib.calculate_PVI.argtypes = (self.c_double_p, self.c_double_p, self.array_length_type, self.window_type, self.step_type)
        lib.calculate_PVI.restype = self.c_double_p


    """function aliases"""
    def calculate_SI(self, pleth_array, patient_height_meters=1, sampling_rate = 62.4725):
        if self.pure_cpp:
            return self.test_SI(pleth_array, patient_height_meters, sampling_rate)
        else:
            return self.test_SI_from_peaks(pleth_array, patient_height_meters, sampling_rate)
    
    def calculate_RI(self, pleth_array):
        if self.pure_cpp:
            return self.test_RI(pleth_array)
        else:
            return self.test_RI_from_peaks(pleth_array)
        
    def calculate_AI(self, pleth_array):
        if self.pure_cpp:
            return self.test_AI(pleth_array)
        else:
            return self.test_AI_from_peaks(pleth_array)
            
    def calculate_IPA(self, pleth_array):
        if self.pure_cpp:
            return self.test_IPA(pleth_array)
        else:
            return self.test_IPA_from_peaks(pleth_array)
                
    def calculate_PI(self, pleth_array):
        return self.test_PI(pleth_array)
    
    def calculate_PVI(self, pleth_array):
        return self.test_PVI(pleth_array)
    
    """faster functions that take a kind of hybrid approach"""
    def test_SI_from_peaks(self, pleth_array, patient_height_meters=1, sampling_rate = 62.4725):
        # getting it as right type
        pleth_array = pleth_array.astype(np.float64)
        start = time.perf_counter()
        systolic_peaks = scipy.signal.find_peaks(pleth_array, distance=10, prominence=0.15)[0].astype(np.int32)
        # print("Time for systolic peaks {:.0f}".format((time.perf_counter() - start)* 1_000_000))
        start = time.perf_counter()
        if systolic_peaks.shape[0] == 0:
            # print("nopeaks detected in test_SI_from_peaks")
            return np.ones(shape = pleth_array.shape).astype(np.float64) * -1


        pleth_array_pointer = pleth_array.ctypes.data_as(self.c_double_p)
        grad = np.zeros(shape=pleth_array.shape).astype(np.float64)
        grad_pointer = grad.ctypes.data_as(self.c_double_p)
        lib.gradient(pleth_array_pointer, grad_pointer, pleth_array.shape[0])
        
        
        # print("Time for get pleth gradient {:.0f}".format((time.perf_counter() - start)* 1_000_000))
        start = time.perf_counter()
        potential_diastolic_peaks = scipy.signal.find_peaks(grad, distance=25)[0].astype(np.int32)
        # print("Time for potential diastolic peaks {:.0f}".format((time.perf_counter() - start)* 1_000_000))
        start = time.perf_counter()



        # getting rid of false peaks from the gradient peaks
        systolic_peaks_pointer = systolic_peaks.ctypes.data_as(self.c_int_p)
        potential_diastolic_peaks_pointer = potential_diastolic_peaks.ctypes.data_as(self.c_int_p)
        diastolic_peaks = np.ones(shape = systolic_peaks.shape).astype(np.int32) * -1
        diastolic_peaks_pointer = diastolic_peaks.ctypes.data_as(self.c_int_p)


        lib.filter_diastolic_peaks(systolic_peaks_pointer, systolic_peaks.shape[0], potential_diastolic_peaks_pointer, potential_diastolic_peaks.shape[0], diastolic_peaks_pointer)

        # print("Time for filter diastolic {:.0f}".format((time.perf_counter() - start)* 1_000_000))
        start = time.perf_counter()

        return_array = np.ones(shape = pleth_array.shape).astype(np.float64) * -1
        return_array_pointer = return_array.ctypes.data_as(self.c_double_p)

        lib.calculate_SI_from_peaks(systolic_peaks_pointer, diastolic_peaks_pointer, systolic_peaks.shape[0],return_array_pointer, return_array.shape[0], patient_height_meters, sampling_rate)

        return return_array

    def test_RI_from_peaks(self, pleth_array):
        # getting it as right type
        pleth_array = pleth_array.astype(np.float64)
        start = time.perf_counter()
        systolic_peaks = scipy.signal.find_peaks(pleth_array, distance=10, prominence=0.15)[0].astype(np.int32)
        # print("Time for systolic peaks {:.0f}".format((time.perf_counter() - start)* 1_000_000))
        start = time.perf_counter()

        if systolic_peaks.shape[0] == 0:
            # print("nopeaks detected in test_RI_from_peaks")
            return np.ones(shape = pleth_array.shape).astype(np.float64) * -1


        pleth_array_pointer = pleth_array.ctypes.data_as(self.c_double_p)
        grad = np.zeros(shape=pleth_array.shape).astype(np.float64)
        grad_pointer = grad.ctypes.data_as(self.c_double_p)
        lib.gradient(pleth_array_pointer, grad_pointer, pleth_array.shape[0])
        
        
        # print("Time for get pleth gradient {:.0f}".format((time.perf_counter() - start)* 1_000_000))
        start = time.perf_counter()
        potential_diastolic_peaks = scipy.signal.find_peaks(grad, distance=25)[0].astype(np.int32)
        # print("Time for potential diastolic peaks {:.0f}".format((time.perf_counter() - start)* 1_000_000))
        start = time.perf_counter()



        # getting rid of false peaks from the gradient peaks
        systolic_peaks_pointer = systolic_peaks.ctypes.data_as(self.c_int_p)
        potential_diastolic_peaks_pointer = potential_diastolic_peaks.ctypes.data_as(self.c_int_p)
        diastolic_peaks = np.ones(shape = systolic_peaks.shape).astype(np.int32) * -1
        diastolic_peaks_pointer = diastolic_peaks.ctypes.data_as(self.c_int_p)


        lib.filter_diastolic_peaks(systolic_peaks_pointer, systolic_peaks.shape[0], potential_diastolic_peaks_pointer, potential_diastolic_peaks.shape[0], diastolic_peaks_pointer)

        # print("Time for filter diastolic {:.0f}".format((time.perf_counter() - start)* 1_000_000))
        start = time.perf_counter()

        # CALCULATING AI NOW
        return_array = np.ones(shape = pleth_array.shape).astype(np.float64) * -1
        return_array_pointer = return_array.ctypes.data_as(self.c_double_p)

        lib.calculate_RI_from_peaks(pleth_array_pointer, systolic_peaks_pointer, diastolic_peaks_pointer, systolic_peaks.shape[0], return_array_pointer, return_array.shape[0])

        return return_array

    def test_AI_from_peaks(self, pleth_array):
        start = time.perf_counter()
        systolic_peaks = scipy.signal.find_peaks(pleth_array, distance=10, prominence=0.15)[0]
        # getting it as right type
        pleth_array = pleth_array.astype(np.float64)
        start = time.perf_counter()
        systolic_peaks = scipy.signal.find_peaks(pleth_array, distance=10, prominence=0.15)[0].astype(np.int32)
        # print("Time for systolic peaks {:.0f}".format((time.perf_counter() - start)* 1_000_000))
        start = time.perf_counter()

        if systolic_peaks.shape[0] == 0:
            # print("nopeaks detected in test_AI_from_peaks")
            return np.ones(shape = pleth_array.shape).astype(np.float64) * -1
        
        pleth_array_pointer = pleth_array.ctypes.data_as(self.c_double_p)
        grad = np.zeros(shape=pleth_array.shape).astype(np.float64)
        grad_pointer = grad.ctypes.data_as(self.c_double_p)
        lib.gradient(pleth_array_pointer, grad_pointer, pleth_array.shape[0])
        
        
        # print("Time for get pleth gradient {:.0f}".format((time.perf_counter() - start)* 1_000_000))
        start = time.perf_counter()
        potential_diastolic_peaks = scipy.signal.find_peaks(grad, distance=25)[0].astype(np.int32)
        # print("Time for potential diastolic peaks {:.0f}".format((time.perf_counter() - start)* 1_000_000))
        start = time.perf_counter()



        # getting rid of false peaks from the gradient peaks
        systolic_peaks_pointer = systolic_peaks.ctypes.data_as(self.c_int_p)
        potential_diastolic_peaks_pointer = potential_diastolic_peaks.ctypes.data_as(self.c_int_p)
        diastolic_peaks = np.ones(shape = systolic_peaks.shape).astype(np.int32) * -1
        diastolic_peaks_pointer = diastolic_peaks.ctypes.data_as(self.c_int_p)


        lib.filter_diastolic_peaks(systolic_peaks_pointer, systolic_peaks.shape[0], potential_diastolic_peaks_pointer, potential_diastolic_peaks.shape[0], diastolic_peaks_pointer)

        # print("Time for filter diastolic {:.0f}".format((time.perf_counter() - start)* 1_000_000))
        start = time.perf_counter()

        # CALCULATING AI NOW
        return_array = np.ones(shape = pleth_array.shape).astype(np.float64) * -1
        return_array_pointer = return_array.ctypes.data_as(self.c_double_p)


        lib.calculate_AI_from_peaks(pleth_array_pointer, systolic_peaks_pointer, diastolic_peaks_pointer, systolic_peaks.shape[0], return_array_pointer, return_array.shape[0])

        return return_array

    def test_IPA_from_peaks(self, pleth_array):
        # getting it as right type
        pleth_array = pleth_array.astype(np.float64)
        start = time.perf_counter()
        systolic_peaks = scipy.signal.find_peaks(pleth_array, distance=10, prominence=0.15)[0].astype(np.int32)
        # print("Time for systolic peaks {:.0f}".format((time.perf_counter() - start)* 1_000_000))
        start = time.perf_counter()

        if systolic_peaks.shape[0] == 0:
            # print("nopeaks detected in test_IPA_from_peaks")
            return np.ones(shape = pleth_array.shape).astype(np.float64) * -1

        pleth_array_pointer = pleth_array.ctypes.data_as(self.c_double_p)
        grad = np.zeros(shape=pleth_array.shape).astype(np.float64)
        grad_pointer = grad.ctypes.data_as(self.c_double_p)
        lib.gradient(pleth_array_pointer, grad_pointer, pleth_array.shape[0])
        
        
        # print("Time for get pleth gradient {:.0f}".format((time.perf_counter() - start)* 1_000_000))
        start = time.perf_counter()
        potential_diastolic_peaks = scipy.signal.find_peaks(grad, distance=25)[0].astype(np.int32)
        # print("Time for potential diastolic peaks {:.0f}".format((time.perf_counter() - start)* 1_000_000))
        start = time.perf_counter()



        # getting rid of false peaks from the gradient peaks
        systolic_peaks_pointer = systolic_peaks.ctypes.data_as(self.c_int_p)
        potential_diastolic_peaks_pointer = potential_diastolic_peaks.ctypes.data_as(self.c_int_p)
        diastolic_peaks = np.ones(shape = systolic_peaks.shape).astype(np.int32) * -1
        diastolic_peaks_pointer = diastolic_peaks.ctypes.data_as(self.c_int_p)


        lib.filter_diastolic_peaks(systolic_peaks_pointer, systolic_peaks.shape[0], potential_diastolic_peaks_pointer, potential_diastolic_peaks.shape[0], diastolic_peaks_pointer)

        # print("Time for filter diastolic {:.0f}".format((time.perf_counter() - start)* 1_000_000))
        start = time.perf_counter()

        # print(systolic_peaks.shape[0], diastolic_peaks.shape[0])
        dicrotic_notches = np.ones(shape = systolic_peaks.shape).astype(np.int32) * -1
        dicrotic_notches_pointer = dicrotic_notches.ctypes.data_as(self.c_int_p)
        # print("Time for dicrotic notches {:.0f}".format((time.perf_counter() - start)* 1_000_000))
        start = time.perf_counter()

        systolic_peaks_pointer = systolic_peaks.ctypes.data_as(self.c_int_p)
        # print(systolic_peaks)
        # print(systolic_peaks.shape)
        lib.find_dicrotic_notches(pleth_array_pointer, systolic_peaks_pointer, systolic_peaks.shape[0], diastolic_peaks_pointer, diastolic_peaks.shape[0], dicrotic_notches_pointer)
        # print("Time for getting notches {:.0f}".format((time.perf_counter() - start)* 1_000_000))
        start = time.perf_counter()

        # CALCULATE IPA NOW
        return_array = np.zeros(shape=pleth_array.shape[0]).astype(np.float64)
        return_array_pointer = return_array.ctypes.data_as(self.c_double_p)
        

        # print("getting_IPA")
        lib.calculate_IPA_from_peaks(pleth_array_pointer, systolic_peaks_pointer, systolic_peaks.shape[0], dicrotic_notches_pointer, return_array_pointer, return_array.shape[0])

        return return_array



    """Full C++ implementations"""

    def test_SI(self, pleth_array, patient_height_meters=1, sampling_rate = 62.4725):
        pleth_array_pointer = pleth_array.ctypes.data_as(self.c_double_p)
        return_array = np.ones(shape = pleth_array.shape) * -1
        return_array_pointer = return_array.ctypes.data_as(self.c_double_p)

        lib.calculate_SI(pleth_array_pointer, return_array_pointer, pleth_array.shape[0], patient_height_meters, sampling_rate)

        return return_array

    def test_RI(self, pleth_array, sampling_rate=62.4725):
        pleth_array_pointer = pleth_array.ctypes.data_as(self.c_double_p)
        return_array = np.ones(shape=pleth_array.shape) * -1
        return_array_pointer = return_array.ctypes.data_as(self.c_double_p)


        lib.calculate_RI(pleth_array_pointer, return_array_pointer, pleth_array.shape[0], sampling_rate)
        return return_array

    def test_AI(self, pleth_array, sampling_rate=62.4725):
        pleth_array_pointer = pleth_array.ctypes.data_as(self.c_double_p)
        return_array = np.ones(shape=pleth_array.shape) * -1
        return_array_pointer = return_array.ctypes.data_as(self.c_double_p)


        lib.calculate_AI(pleth_array_pointer, return_array_pointer, pleth_array.shape[0])
        return return_array

    def test_IPA(self, pleth_array, sampling_rate=62.4725):
        pleth_array_pointer = pleth_array.ctypes.data_as(self.c_double_p)
        return_array = np.ones(shape=pleth_array.shape) * -1
        return_array_pointer = return_array.ctypes.data_as(self.c_double_p)

        
        lib.calculate_IPA(pleth_array_pointer, return_array_pointer, return_array.shape[0], sampling_rate)
        return return_array

    def test_PI(self, pleth_array):
        pleth_array_pointer = pleth_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

        pi_value = lib.calculate_PI(pleth_array_pointer, pleth_array.shape[0])
        return pi_value

    def test_PVI(self, pleth_array, window=500, step=10):
        pleth_array_pointer = pleth_array.ctypes.data_as(self.c_double_p)
        return_array = np.ones(shape=pleth_array.shape) * -1
        return_array_pointer = return_array.ctypes.data_as(self.c_double_p)


        # Call the C function
        lib.calculate_PVI(pleth_array_pointer, return_array_pointer, pleth_array.shape[0], window, step)

        return return_array


if __name__ == "__main__":
    import wfdb
    # hadm id 24239751
    # | subject_id | hadm_id  | stay_id  | caregiver_id | charttime           | storetime           | itemid | value | valuenum | valueuom | warning |
    # |   10079700 | 24239751 | 35381116 |        62802 | 2115-10-11 07:22:00 | 2115-10-11 14:15:00 | 226707 | 74    |       74 | Inch     |       0 |
    test_height_meters = 74 * 0.0254 #convert inches to meters
    signal, info = wfdb.rdsamp("./Waveform_Analysis/p100/p10079700/85594648/85594648", sampfrom=80000)
    sampling_freq = info["fs"]
    print(sampling_freq)
    window_time = 150 #seconds
    window_size = int(window_time*sampling_freq)


    pleth_index = info["sig_name"].index("Pleth")
    test_pleth_signal = signal[:window_size,pleth_index]

    ctypes_obj = PPG_METRICS_CTYPES()

    import time
    overall_start = time.perf_counter()
    SI_value = ctypes_obj.test_SI_from_peaks(test_pleth_signal, test_height_meters, sampling_freq)
    print("Total Time for partial Ctypes SI calculation {:.0f}".format((time.perf_counter() - overall_start)* 1_000_000))
    print()

    overall_start = time.perf_counter()
    RI_value = ctypes_obj.test_RI_from_peaks(test_pleth_signal)
    print("Total Time for partial Ctypes RI calculation {:.0f}".format((time.perf_counter() - overall_start)* 1_000_000))
    print()

    overall_start = time.perf_counter()
    AI_value = ctypes_obj.test_AI_from_peaks(test_pleth_signal)
    print("Total Time for partial Ctypes AI calculation {:.0f}".format((time.perf_counter() - overall_start)* 1_000_000))
    print()

    overall_start = time.perf_counter()
    IPA_value = ctypes_obj.test_IPA_from_peaks(test_pleth_signal)
    print("Total Time for partial Ctypes IPA calculation {:.0f}".format((time.perf_counter() - overall_start)* 1_000_000))
    print()


    overall_start = time.perf_counter()
    SI_value = ctypes_obj.test_SI(test_pleth_signal, test_height_meters, sampling_freq)
    print("Total Time for Full Ctypes SI calculation {:.0f}".format((time.perf_counter() - overall_start)* 1_000_000))
    print()
    start = time.perf_counter()
    PLOTTING = False
    if PLOTTING:
        fig, ax = plt.subplots(3,1, sharex=True)

        time_x = np.arange(len(test_pleth_signal))/sampling_freq
        ax[0].set_title("Pleth Signal")
        ax[0].plot(time_x, test_pleth_signal)
        ax[0].set_xlabel('Time')
        ax[1].set_title("Grad of Pleth signal")
        ax[1].plot(time_x, np.gradient(test_pleth_signal))
        ax[1].scatter(SI_value/sampling_freq, np.gradient(test_pleth_signal)[SI_value.astype(int)], c="r")

        ax[1].set_xlabel('Time')
        ax[2].plot(time_x, SI_value)
        ax[2].set_xlabel("Time")
        ax[2].set_title("Stiffness Index")
        plt.tight_layout()
        
        print("GOT SI")


    RI_value = ctypes_obj.test_RI(test_pleth_signal.astype(np.double), sampling_freq)

    if PLOTTING:
        fig, ax = plt.subplots(3,1, sharex=True)
        ax[0].set_title("Pleth Signal")
        ax[0].plot(time_x, test_pleth_signal)
        ax[1].set_title("Grad of Pleth signal")
        ax[1].plot(time_x, np.gradient(test_pleth_signal))
        ax[2].plot(time_x, RI_value)
        ax[2].set_title("Reflection Index")
        plt.tight_layout()
        print("GOT RI")


        AI_value = ctypes_obj.test_AI(test_pleth_signal.astype(np.double), sampling_freq)

    if PLOTTING:
        fig, ax = plt.subplots(3,1, sharex=True)
        ax[0].set_title("Pleth Signal")
        ax[0].plot(time_x, test_pleth_signal)
        ax[1].set_title("Grad of Pleth signal")
        ax[1].plot(time_x, np.gradient(test_pleth_signal))
        ax[2].plot(time_x, AI_value)
        ax[2].set_title("Augmentation Index")
        plt.tight_layout()
        print("GOT AI")


    IPA_value = ctypes_obj.test_IPA(test_pleth_signal.astype(np.double), sampling_freq)

    if PLOTTING:
        fig, ax = plt.subplots(3,1, sharex=True)
        ax[0].set_title("Pleth Signal")
        ax[0].plot(time_x, test_pleth_signal)
        ax[1].set_title("Grad of Pleth signal")
        ax[1].plot(time_x, np.gradient(test_pleth_signal))
        ax[2].plot(time_x, IPA_value)
        ax[2].set_title("Inflection point area values")
        plt.tight_layout()
        print("GOT IPA")

    PVI_value = ctypes_obj.test_PVI(test_pleth_signal.astype(np.double))

    if PLOTTING:
        fig, ax = plt.subplots(1)
        ax.set_title("PVI")
        ax.plot(time_x, PVI_value)
        plt.tight_layout()
        print("GOT PVI")

    if PLOTTING:
        plt.show()


    # test_data = pd.read_csv("Waveform_Analysis/analysis_speed_test_signals.csv")[1:]
    # print("data_loaded")

    # num_repeats = 10
    # test_height_meters = 1.75
    # sampling_freq = 62.4725

    # # Test for SI calculation
    # start = time.time()
    # for i in range(num_repeats):
    #     for series_name, series in test_data.items():
    #         SI_value = ctypes_obj.test_SI(series.to_numpy().astype(np.double), test_height_meters, sampling_freq)
    # print("Total time for SI calculation {} python: {:.2f}".format("Ctypes ", time.time() - start))

    # # Test for RI calculation
    # start = time.time()
    # for i in range(num_repeats):
    #     for series_name, series in test_data.items():
    #         RI_value = ctypes_obj.test_RI(series.to_numpy().astype(np.double), sampling_freq)
    # print("Total time for RI calculation {} python: {:.2f}".format("Ctypes ", time.time() - start))

    # # Test for AI calculation
    # start = time.time()
    # for i in range(num_repeats):
    #     for series_name, series in test_data.items():
    #         AI_value = ctypes_obj.test_AI(series.to_numpy().astype(np.double), sampling_freq)
    # print("Total time for AI calculation {} python: {:.2f}".format("Ctypes ", time.time() - start))

    # # Test for IPA calculation
    # start = time.time()
    # for i in range(num_repeats):
    #     for series_name, series in test_data.items():
    #         IPA_value = ctypes_obj.test_IPA(series.to_numpy().astype(np.double), sampling_freq)
    # print("Total time for IPA calculation {} python: {:.2f}".format("Ctypes ", time.time() - start))

    # # Test for PI calculation
    # start = time.time()
    # for i in range(num_repeats):
    #     for series_name, series in test_data.items():
    #         PI_value = ctypes_obj.test_PI(series.to_numpy().astype(np.double))
    # print("Total time for PI calculation {} python: {:.2f}".format("Ctypes ", time.time() - start))

    # # Test for PVI calculation
    # start = time.time()
    # for i in range(num_repeats):
    #     for series_name, series in test_data.items():
    #         PVI_value = ctypes_obj.test_PVI(series.to_numpy().astype(np.double))
    # print("Total time for PVI calculation {} python: {:.2f}".format("Ctypes ", time.time() - start))

