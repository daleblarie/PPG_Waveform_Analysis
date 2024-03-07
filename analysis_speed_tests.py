import scipy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import time
import sys

from ppg_metrics_numba_fullclass import *
from CPP_waveform_analysis.wrapper import *

from ppg_related_metrics import PPG_METRICS


# numpy forward fill from 
# https://stackoverflow.com/questions/41190852/most-efficient-way-to-forward-fill-nan-values-in-numpy-array
def numpy_fill(arr):
    '''Solution provided by Divakar.'''
    mask = np.isnan(arr)
    idx = np.where(~mask,np.arange(mask.shape[1]),0)
    np.maximum.accumulate(idx,axis=1, out=idx)
    out = arr[np.arange(idx.shape[0])[:,None], idx]
    return out

if __name__=="__main__":
    # load test data
    try:
        output_df = pd.read_csv("Waveform_Analysis/analysis_speed_test_signals.csv")
    except FileNotFoundError:
        all_files = os.listdir("Waveform_Analysis/analysis_speed_test_signals/")
        output_df = None
        for filename in all_files:
            print(filename)
            data = np.genfromtxt("Waveform_Analysis/analysis_speed_test_signals/" + filename)
            # data = np.genfromtxt("Waveform_Analysis/analysis_speed_test_signals/p011232-2162-11-12-10-04.csv")
            # print(data.shape)
            # print(np.count_nonzero(~np.isnan(data)))
            # if more than 3/4 nan, just drop it tbh
            if np.count_nonzero(np.isnan(data)) > (data.shape[0] *3 /4):
                continue
            if np.isnan(data[0]):
                data[0] = 0
            data = pd.DataFrame(data)
            data = data.ffill()
            # data = data[10000:20000].reset_index(drop=True)

            if output_df is None:
                output_df = data
            else:
                output_df.insert(output_df.shape[1], str(output_df.shape[1]), data)

        output_df.to_csv("Waveform_Analysis/analysis_speed_test_signals.csv", index_col=0)
        
    print(output_df)
    # plt.figure()
    # plt.plot(data)
    # plt.show()

    # TESTING RAW PYTHON HERE
    ppg_functions = PPG_METRICS(False)

    ppg_functions_numba = PPG_METRICS_NUMBA()

    ppg_functions_ctypes = PPG_METRICS_CTYPES()
    ppg_functions_ctypes_pure_cpp = PPG_METRICS_CTYPES(True)

    def run_speed_test(ppg_function_object, input_df, ppg_function_object_type="raw", num_repeats=10):

        start = time.time()
        for i in range(num_repeats):
            for series_name, series in input_df.items():
                SI_value = ppg_function_object.calculate_SI(series.to_numpy())
                # print("Average SI Value: {:.2f}".format(np.mean(SI_value)))
        print("Total time for SI calculation {} python: {:.2f}".format(ppg_function_object_type, time.time() - start))

        start = time.time()
        for i in range(num_repeats):
            for series_name, series in input_df.items():
                RI_value = ppg_function_object.calculate_RI(series.to_numpy())
                # print("Average RI Value: {:.2f}".format(np.mean(RI_value)))
        print("Total time for RI calculation {} python: {:.2f}".format(ppg_function_object_type, time.time() - start))

        start = time.time()
        for i in range(num_repeats):
            for series_name, series in input_df.items():
                AI_value = ppg_function_object.calculate_AI(series.to_numpy())
                # print("Average AI Value: {:.2f}".format(np.mean(AI_value)))
        print("Total time for AI calculation {} python: {:.2f}".format(ppg_function_object_type, time.time() - start))

        start = time.time()
        for i in range(num_repeats):
            for series_name, series in input_df.items():
                IPA_value = ppg_function_object.calculate_IPA(series.to_numpy())
                # print("Average IPA Value: {:.2f}".format(np.mean(IPA_value)))
        print("Total time for IPA calculation {} python: {:.2f}".format(ppg_function_object_type, time.time() - start))

        start = time.time()
        for i in range(num_repeats):
            for series_name, series in input_df.items():
                PI_value = ppg_function_object.calculate_PI(series.to_numpy())
                # print("Average PI Value: {:.2f}".format(np.mean(PI_value)))
        print("Total time for PI calculation {} python: {:.2f}".format(ppg_function_object_type, time.time() - start))

        start = time.time()
        for i in range(num_repeats):
            for series_name, series in input_df.items():
                PVI_value = ppg_function_object.calculate_PVI(series.to_numpy())
                # print("Average PVI Value: {:.2f}".format(np.mean(PVI_value)))
        print("Total time for PVI calculation {} python: {:.2f}".format(ppg_function_object_type, time.time() - start))

        # start = time.time()
        # for i in range(num_repeats):
        #     for series_name, series in input_df.items():
        #         PSQI_value = ppg_function_object.calculate_PSQI(series.to_numpy())
        #         # print("Average PSQI Value: {:.2f}".format(np.mean(PSQI_value)))
        # print("Total time for PSQI calculation {} python: {:.2f}".format(ppg_function_object_type, time.time() - start))

        # start = time.time()
        # for i in range(num_repeats):
        #     for series_name, series in input_df.items():
        #         ESQI_value = ppg_function_object.calculate_ESQI(series.to_numpy())
        #         # print("Average ESQI Value: {:.2f}".format(np.mean(ESQI_value)))
        # print("Total time for ESQI calculation {} python: {:.2f}".format(ppg_function_object_type, time.time() - start))

        # start = time.time()
        # for i in range(num_repeats):
        #     for series_name, series in input_df.items():
        #         filtered_breath_rate = ppg_function_object.calculate_resp_rate(series.to_numpy())
        #         # print("RESP RATE: {:.2f}".format(filtered_breath_rate))
        # print("Total time for RESP RATE calculation {} python: {:.2f}".format(ppg_function_object_type, time.time() - start))

        # start = time.time()
        # for i in range(num_repeats):
        #     for series_name, series in input_df.items():
        #         PSD_value = ppg_function_object.calculate_PSD(series.to_numpy())
        #         # print("PSD Values: {}".format(PSD_value))
        # print("Total time for PSD calculation {} python: {:.2f}".format(ppg_function_object_type, time.time() - start))

        # start = time.time()
        # for i in range(num_repeats):
        #     for series_name, series in input_df.items():
        #         Entropy_value = ppg_function_object.calculate_entropy(series.to_numpy())
        #         # print("Entropy Values: {}".format(Entropy_value))
        # print("Total time for Entropy calculation {} python: {:.2f}".format(ppg_function_object_type, time.time() - start))

        plt.show()
        print("Total Test Samples Number: {:d}\n and Total Num Repeats{:d}".format(input_df.shape[1], num_repeats))


    def run_raw_speedtest(output_df, num_repeats=10, type_name="NUMBA"):
        if type_name=="NUMBA WITH SOME NESTING":
            from ppg_metrics_numba_nested import calculate_SI, calculate_RI, calculate_AI, calculate_IPA, calculate_PI, calculate_PVI,calculate_PSQI, calculate_ESQI, calculate_resp_rate, calculate_PSD, calculate_entropy

        if type_name=="NUMBA NO NEST":
            from ppg_metrics_numba_nonesting import calculate_SI, calculate_RI, calculate_AI, calculate_IPA, calculate_PI, calculate_PVI,calculate_PSQI, calculate_ESQI, calculate_resp_rate, calculate_PSD, calculate_entropy


        start = time.time()
        for i in range(num_repeats):
            for series_name, series in output_df.items():
                SI_value = calculate_SI(series.to_numpy())
                # print("Average SI Value: {:.2f}".format(np.mean(SI_value)))
        print("Total time for SI calculation {} python: {:.2f}".format(type_name, time.time() - start))

        start = time.time()
        for i in range(num_repeats):
            for series_name, series in output_df.items():
                RI_value = calculate_RI(series.to_numpy())
                # print("Average RI Value: {:.2f}".format(np.mean(RI_value)))
        print("Total time for RI calculation {} python: {:.2f}".format(type_name, time.time() - start))

        start = time.time()
        for i in range(num_repeats):
            for series_name, series in output_df.items():
                AI_value = calculate_AI(series.to_numpy())
                # print("Average AI Value: {:.2f}".format(np.mean(AI_value)))
        print("Total time for AI calculation {} python: {:.2f}".format(type_name, time.time() - start))

        start = time.time()
        for i in range(num_repeats):
            for series_name, series in output_df.items():
                IPA_value = calculate_IPA(series.to_numpy())
                # print("Average IPA Value: {:.2f}".format(np.mean(IPA_value)))
        print("Total time for IPA calculation {} python: {:.2f}".format(type_name, time.time() - start))

        start = time.time()
        for i in range(num_repeats):
            for series_name, series in output_df.items():
                PI_value = calculate_PI(series.to_numpy())
                # print("Average PI Value: {:.2f}".format(np.mean(PI_value)))
        print("Total time for PI calculation {} python: {:.2f}".format(type_name, time.time() - start))

        start = time.time()
        for i in range(num_repeats):
            for series_name, series in output_df.items():
                PVI_value = calculate_PVI(series.to_numpy())
                # print("Average PVI Value: {:.2f}".format(np.mean(PVI_value)))
        print("Total time for PVI calculation {} python: {:.2f}".format(type_name, time.time() - start))

        # start = time.time()
        # for i in range(num_repeats):
        #     for series_name, series in output_df.items():
        #         PSQI_value = calculate_PSQI(series.to_numpy())
        #         # print("Average PSQI Value: {:.2f}".format(np.mean(PSQI_value)))
        # print("Total time for PSQI calculation {} python: {:.2f}".format(type_name, time.time() - start))

        # start = time.time()
        # for i in range(num_repeats):
        #     for series_name, series in output_df.items():
        #         ESQI_value = calculate_ESQI(series.to_numpy())
        #         # print("Average ESQI Value: {:.2f}".format(np.mean(ESQI_value)))
        # print("Total time for ESQI calculation {} python: {:.2f}".format(type_name, time.time() - start))

        # start = time.time()
        # for i in range(num_repeats):
        #     for series_name, series in output_df.items():
        #         filtered_breath_rate = calculate_resp_rate(series.to_numpy())
        #         # print("RESP RATE: {:.2f}".format(filtered_breath_rate))
        # print("Total time for RESP RATE calculation {} python: {:.2f}".format(type_name, time.time() - start))

        # start = time.time()
        # for i in range(num_repeats):
        #     for series_name, series in output_df.items():
        #         PSD_value = calculate_PSD(series.to_numpy())
        #         # print("PSD Values: {}".format(PSD_value))
        # print("Total time for PSD calculation {} python: {:.2f}".format(type_name, time.time() - start))

        # start = time.time()
        # for i in range(num_repeats):
        #     for series_name, series in output_df.items():
        #         Entropy_value = calculate_entropy(series.to_numpy())
        #         # print("Entropy Values: {}".format(Entropy_value))
        # print("Total time for Entropy calculation {} python: {:.2f}".format(type_name, time.time() - start))

        plt.show()
        print("Total Test Samples Number: {:d}\n and Total Num Repeats{:d}".format(output_df.shape[1], num_repeats))

    print("\n\n\n--------------------")
    print("BEGINNING SPEED TEST\n\n\n")
    run_speed_test(ppg_functions_ctypes, input_df=output_df, ppg_function_object_type="Ctypes")
    print("--------------------\n\n\n")
    run_speed_test(ppg_functions_ctypes_pure_cpp, input_df=output_df, ppg_function_object_type="Ctypes_Pure_Cpp")
    print("--------------------\n\n\n")
    run_raw_speedtest(output_df, type_name="NUMBA WITH SOME NESTING")
    print("\n\n\n--------------------\n\n\n")
    run_raw_speedtest(output_df, type_name="NUMBA NO NEST")
    print("\n\n\n--------------------\n\n\n")
    run_speed_test(ppg_functions, output_df, "Raw")
    print("\n\n\n--------------------\n\n\n")
    run_speed_test(ppg_functions_numba, input_df=output_df, ppg_function_object_type="Numba Class")
