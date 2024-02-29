import scipy
import matplotlib.pyplot as plt
import numpy as np

class PPG_METRICS:
    # default_sampling_freq = 62.4725
    def __init__(self, PLOTTING=False):
        self.PLOTTING = PLOTTING
        default_sampling_freq = 62.4725
    # class to contain all the calculation functions
        
    # Stiffness index
    def calculate_SI(self, pleth_signal, patient_height_meters=1, sampling_freq=62.4725):
        systolic_peaks = scipy.signal.find_peaks(pleth_signal, distance=10, prominence=0.15)[0]
        potential_diastolic_peaks = scipy.signal.find_peaks(np.gradient(pleth_signal), distance=25)[0]
        SI_array = np.zeros(shape=len(pleth_signal))

        # getting rid of false peaks from the gradient peaks
        diastolic_peaks = []
        dias_index = 0
        for sys_index in range(len(systolic_peaks)):
            sys_peak = systolic_peaks[sys_index]
            while potential_diastolic_peaks[dias_index] < sys_peak + 10 and dias_index < potential_diastolic_peaks.shape[0]-1:
                dias_index += 1
            diastolic_peaks.append(potential_diastolic_peaks[dias_index])

        # calculating SI over time with a 0 order hold
        prev_peak = 0
        SI_value = 0
        i=0
        for i in range(len(systolic_peaks)):
            SI_array[prev_peak:diastolic_peaks[i]] = SI_value
            SI_value = patient_height_meters / ((diastolic_peaks[i] - systolic_peaks[i])/sampling_freq)
            prev_peak = diastolic_peaks[i]
        SI_array[prev_peak:] = SI_value

        if self.PLOTTING:
            fig, ax = plt.subplots(3,1, sharex=True)

            time = np.arange(len(pleth_signal))/sampling_freq
            ax[0].set_title("Pleth Signal")
            ax[0].plot(time, pleth_signal)
            ax[0].set_xlabel('Time')
            ax[0].scatter([peak/sampling_freq for peak in systolic_peaks], pleth_signal[systolic_peaks], c="r")
            ax[0].scatter([peak/sampling_freq for peak in diastolic_peaks], pleth_signal[diastolic_peaks], c="b")
            ax[1].set_title("Grad of Pleth signal")
            ax[1].plot(time, np.gradient(pleth_signal))
            ax[1].set_xlabel('Time')
            ax[1].scatter([peak/sampling_freq for peak in systolic_peaks], np.gradient(pleth_signal)[systolic_peaks], c="r")
            ax[1].scatter([peak/sampling_freq for peak in diastolic_peaks], np.gradient(pleth_signal)[diastolic_peaks], c="b")
            ax[2].plot(time, SI_array)
            ax[2].set_xlabel("Time")
            ax[2].set_title("Stiffness Index")
            plt.tight_layout()

        return SI_array

    def calculate_RI(self, pleth_signal):
        systolic_peaks = scipy.signal.find_peaks(pleth_signal, distance=10, prominence=0.15)[0]
        potential_diastolic_peaks = scipy.signal.find_peaks(np.gradient(pleth_signal), distance=25)[0]
        RI_array = np.zeros(shape=len(pleth_signal))

        # getting rid of false peaks from the gradient peaks
        diastolic_peaks = []
        dias_index = 0
        for sys_index in range(len(systolic_peaks)):
            sys_peak = systolic_peaks[sys_index]
            while potential_diastolic_peaks[dias_index] < sys_peak + 10 and dias_index < potential_diastolic_peaks.shape[0]-1:
                dias_index += 1
            diastolic_peaks.append(potential_diastolic_peaks[dias_index])

        # calculating RI over time with a 0 order hold
        prev_peak = 0
        RI_value = 0
        i=0
        for i in range(len(systolic_peaks)):
            RI_array[prev_peak:diastolic_peaks[i]] = RI_value
            RI_value = 100 * (pleth_signal[diastolic_peaks[i]] / pleth_signal[systolic_peaks[i]])
            prev_peak = diastolic_peaks[i]
        RI_array[prev_peak:] = RI_value

        if self.PLOTTING:
            fig, ax = plt.subplots(3,1, sharex=True)
            ax[0].set_title("Pleth Signal")
            ax[0].plot(pleth_signal)
            ax[0].scatter(systolic_peaks, pleth_signal[systolic_peaks], c="r")
            ax[0].scatter(diastolic_peaks, pleth_signal[diastolic_peaks], c="b")
            ax[1].set_title("Grad of Pleth signal")
            ax[1].plot(np.gradient(pleth_signal))
            ax[1].scatter(systolic_peaks, np.gradient(pleth_signal)[systolic_peaks], c="r")
            ax[1].scatter(diastolic_peaks, np.gradient(pleth_signal)[diastolic_peaks], c="b")
            ax[2].plot(RI_array)
            ax[2].set_title("Reflection Index")
            plt.tight_layout()

        return RI_array

    def calculate_AI(self, pleth_signal):
        systolic_peaks = scipy.signal.find_peaks(pleth_signal, distance=10, prominence=0.15)[0]
        potential_diastolic_peaks = scipy.signal.find_peaks(np.gradient(pleth_signal), distance=25)[0]
        AI_array = np.zeros(shape=len(pleth_signal))

        # getting rid of false peaks from the gradient peaks
        diastolic_peaks = []
        dias_index = 0
        for sys_index in range(len(systolic_peaks)):
            sys_peak = systolic_peaks[sys_index]
            while potential_diastolic_peaks[dias_index] < sys_peak + 10 and dias_index < potential_diastolic_peaks.shape[0]-1:
                dias_index += 1
            diastolic_peaks.append(potential_diastolic_peaks[dias_index])

        # calculating AI over time with a 0 order hold
        prev_peak = 0
        AI_value = 0
        i=0
        for i in range(len(systolic_peaks)):
            AI_array[prev_peak:diastolic_peaks[i]] = AI_value
            AI_value = 100 * ((pleth_signal[systolic_peaks[i]] - pleth_signal[diastolic_peaks[i]]) / pleth_signal[systolic_peaks[i]])
            prev_peak = diastolic_peaks[i]
        AI_array[prev_peak:] = AI_value

        if self.PLOTTING:
            fig, ax = plt.subplots(3,1, sharex=True)
            ax[0].set_title("Pleth Signal")
            ax[0].plot(pleth_signal)
            ax[0].scatter(systolic_peaks, pleth_signal[systolic_peaks], c="r")
            ax[0].scatter(diastolic_peaks, pleth_signal[diastolic_peaks], c="b")
            ax[1].set_title("Grad of Pleth signal")
            ax[1].plot(np.gradient(pleth_signal))
            ax[1].scatter(systolic_peaks, np.gradient(pleth_signal)[systolic_peaks], c="r")
            ax[1].scatter(diastolic_peaks, np.gradient(pleth_signal)[diastolic_peaks], c="b")
            ax[2].plot(AI_array)
            ax[2].set_title("Augmentation Index")
            plt.tight_layout()

        return AI_array

    def calculate_IPA(self, pleth_signal):
        systolic_peaks = scipy.signal.find_peaks(pleth_signal, distance=10, prominence=0.15)[0]
        potential_diastolic_peaks = scipy.signal.find_peaks(np.gradient(pleth_signal), distance=25)[0]
        # getting rid of false peaks from the gradient peaks
        diastolic_peaks = []
        dias_index = 0
        for sys_index in range(len(systolic_peaks)):
            sys_peak = systolic_peaks[sys_index]
            while potential_diastolic_peaks[dias_index] < sys_peak + 10 and dias_index < potential_diastolic_peaks.shape[0]-1:
                dias_index += 1
            diastolic_peaks.append(potential_diastolic_peaks[dias_index])
        diastolic_peaks = np.array(diastolic_peaks)
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7309072/
        # The dicrotic notch is an essential feature of the PPG signal. Figure 8 describes the algorithm to detect the dicrotic notch. To do so, a line was drawn from the systolic peak to the diastolic peak. The minimum of the subtraction of the straight line from the signal is the dicrotic notch. However, to make it more robust, the fix index was used, which calculates the local minima within a given window (in this case 50 ms) around a given point. Reliable detection of the dicrotic notch in various situations is shown in Figure 9.

        dicrotic_notch = []
        for i in range(systolic_peaks.shape[0]):
            sys_peak, dias_peak = systolic_peaks[i], diastolic_peaks[i]
            subtraction_line = np.linspace(pleth_signal[sys_peak], pleth_signal[dias_peak], num=(dias_peak-sys_peak))
            sys_to_dias_window = pleth_signal[sys_peak:dias_peak]
            notch = np.argmin(sys_to_dias_window-subtraction_line) + sys_peak

            dicrotic_notch.append(notch)


        IPA_array = np.zeros(shape=len(pleth_signal))
        signal_troughs = np.zeros(shape=systolic_peaks.shape).astype(int)
        ipa_value = 0
        prev_trough = -1
        for i in range(len(dicrotic_notch)-1):
            sys_peak, notch = systolic_peaks[i], dicrotic_notch[i]
            trough = np.argmin(pleth_signal[sys_peak:systolic_peaks[i+1]]) + sys_peak
            IPA_array[prev_trough:trough] = ipa_value
            
            A1 = np.sum(pleth_signal[prev_trough:notch])
            if A1 ==0:
                ipa_value = -1
            else:
                A2 = np.sum(pleth_signal[notch:trough])
                ipa_value = A2/A1

            prev_trough = trough
            signal_troughs[i] = trough
        IPA_array[prev_trough:] = ipa_value


        if self.PLOTTING:
            fig, ax = plt.subplots(3,1, sharex=True)
            ax[0].set_title("Pleth Signal")
            ax[0].plot(pleth_signal)
            ax[0].scatter(signal_troughs, pleth_signal[signal_troughs], c="r")
            ax[0].scatter(dicrotic_notch, pleth_signal[dicrotic_notch], c="b")
            ax[1].set_title("Pleth Signal gradient")
            ax[1].plot(np.gradient(pleth_signal))
            ax[1].scatter(signal_troughs, np.gradient(pleth_signal)[signal_troughs], c="r")
            ax[1].scatter(dicrotic_notch, np.gradient(pleth_signal)[dicrotic_notch], c="b")
            ax[2].plot(IPA_array)
            ax[2].set_title("IPA Values")
            plt.tight_layout()

        return IPA_array

    def calculate_PI(self, pleth_signal):   #calculate perfusion index
        min_value = min(pleth_signal)
        DC_sum = min_value*len(pleth_signal)
        if DC_sum == 0:
            return -1
        AC_sum = np.sum(pleth_signal-min_value)
        PI_value = 100 * AC_sum/DC_sum



        # fig, ax = plt.subplots(1)

        # ax.plot(pleth_signal)
        # ax.fill_between(np.linspace(0, len(pleth_signal), num=len(pleth_signal)), min_value, pleth_signal, color="b", alpha=0.4)
        # ax.hlines(min_value,0, len(pleth_signal))
        # ax.fill_between(np.linspace(0, len(pleth_signal), num=len(pleth_signal)), 0, min_value, color="r", alpha=0.4)

        # ax.set_ylim([0,1])

        return PI_value

    def calculate_PVI(self, pleth_signal, window=500, step=10, PI_min=None, PI_max=None):
        PVI = np.zeros(len(pleth_signal))
        PVI_index = window
        start_index = 0
        if PI_min is None:
            PI_min = np.inf
        if PI_max is None:
            PI_max = 0
        while start_index + window < len(pleth_signal):
            window_signal = pleth_signal[start_index:start_index + window]
            current_PI = self.calculate_PI(window_signal)
            if current_PI > PI_max:
                PI_max = current_PI
            if current_PI < PI_min:
                PI_min = current_PI

            if PI_max==0:
                PVI[PVI_index:PVI_index+step] = -1
            else:
                PVI[PVI_index:PVI_index+step] = 100 * (PI_max - PI_min)/PI_max
            start_index += step
            PVI_index += step

        if self.PLOTTING:
            fig, ax = plt.subplots(1)
            ax.set_title("PVI")
            ax.plot(PVI)
            plt.tight_layout()

        return PVI, PI_min, PI_max

    def calculate_PSQI(self, pleth_signal, sampling_freq=62.4725):
        from scipy.signal import butter, lfilter, freqz
        def moving_average(data, n=3):
            window = np.ones(n)/n
            rolling_mean = np.convolve(data, window, mode="valid")
            return rolling_mean
        
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5597264/
        # Perfusion (PSQI): As previously mentioned, this is the gold standard for assessing PPG signal quality [14,15,16,17]. The perfusion index is the ratio of the pulsatile blood flow to the nonpulsatile or static blood in peripheral tissue. In other words, it is the difference of the amount of light absorbed through the pulse of when light is transmitted through the finger, which can be defined as follows:
        # PSQI=[(ymax−ymin)/|x¯|]×100
        # (2)
        # where PSQI is the perfusion index, x¯ is the statistical mean of the x signal (raw PPG signal), and y is the filtered PPG signal.
        
        # getting window size for moving mean filter
        # either window of 3 samples, or 0.25 seconds
        window_size = np.max((3, int(0.25 * sampling_freq)))
        filtered_data = moving_average(pleth_signal, window_size)

        y_max = max(filtered_data)
        y_min = min(filtered_data)
        mean = np.mean(pleth_signal)
        PSQI = ((y_max - y_min) / mean) * 100
        return PSQI

    def calculate_ESQI(self, pleth_signal):
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5597264/
        # Entropy (ESQI): Recently, Selvaraj et al. [19] found that entropy is a good indicator for PPG signal quality. Entropy quantifies how much the probability density function (PDF) of the signal differs from a uniform distribution and thus provides a quantitative measure of the uncertainty present in the signal [24], which is defined [25] as:
        # ESQI=−∑n=1Nx[n]^2loge(x[n]^2)
        # where x signal is the raw PPG signal and N is the number of data points.
        
        return (-1 * np.sum(np.square(pleth_signal) * np.log(np.square(pleth_signal))))

    def calculate_resp_rate(self, pleth_signal, sampling_freq=62.4725):
        from scipy.signal import butter, lfilter, filtfilt, welch, find_peaks

        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9056464/
        # Respiratory Rate Estimation The peak detection algorithm may also detect some false peaks [44, 45]. The proposed algorithm ensures that the algorithm detects only true peaks by removing the false peaks. The peak is considered to be a false peak if the R–R interval between two adjacent R-peaks is less than 30% of the mean R–R interval of the analysed signal, as mentioned in [22] and any false peak is dropped. The R–R intervals are calculated again after dropping any false peaks. The new values of the RR interval are then used to calculate different time-series measurements. These measurements include heart rate in beats per minute (BPM), RR-intervals or inter-beat interval (IBI) and estimated respiratory rate in breaths per min. The respiratory rate is calculated using Welch’s method [46, 47]. Welch’s method divides the inter-beat intervals (signal) into overlapping segments and computes a modified periodogram for each segment. Then an average of all periodograms along with an array of frequencies are returned at the output. The respiratory rate is the maximum frequency within the frequency band (Hz) that is, where PSD is maximum. To calculate the final respiratory rate (RespR) the frequency is multiplied by 60 to convert the respiratory frequency band from Hz to breaths per minute. Mathematically it is represented as,
        # Respirationrate=f×argmaxiP(i)×60

        def butter_bandpass_filter(data, lowcut, highcut, sampling_freq, order=3):
            nyq = 0.5 * sampling_freq
            low = lowcut / nyq
            high = highcut / nyq
            b, a = butter(order, [low, high], btype='band', analog=False)
            y = lfilter(b, a, data)
            return y
        
        def peak_enhancement(data, low=0, high=1024):
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9056464/
            # Peak Enhancement To enhance the signal to noise ratio and to get better detection of the peaks, the algorithm performs peak enhancement using the peak enhancement algorithm. This is a crucial step as all the further calculations will be dependent on the detected peak. Thus, accurate detection of peaks will play a significant role in this algorithm. This function makes the higher peaks more dominant while suppressing the smaller peaks, for better detection of the peaks. This function scales the signal to the specified lower and upper range. The formula is as in Eq. (1):
            # 𝑃𝑒𝑛ℎ=𝑟𝑏∗((𝑥−𝑚𝑖𝑛(𝑥))/𝑟𝑎𝑛𝑔𝑒(𝑥))+𝑙_𝑙𝑡
            # where 𝑃𝑒𝑛ℎ is peak enhancement, 𝑥 is the raw signal, 𝑟𝑏 is the range of upper limit and 𝑙_𝑙𝑡 (lower limits), given by the user (by default it is set at 1024 and 0, respectively). The 𝑟𝑎𝑛𝑔𝑒(𝑥) is a valued range calculated by subtracting the maximum value of the analysed signal from the minimum value.

            ## idk why they feel the need to call it peak enhancement, theyre just normalizing to a range, its nothing special, and doesnt actually do anything
            enhanced_data = (high-low) * ((data-min(data))/max(data) - min(data)) + low
            return enhanced_data

        resp_freq_low = 0.1 #hz (6 breaths/min)
        resp_freq_high = 0.4 #hz (24 breaths/min)

        filtered_signal = butter_bandpass_filter(pleth_signal, resp_freq_low, resp_freq_high, sampling_freq)
        peak_enhanced_signal = peak_enhancement(filtered_signal)


        # find peaks
        RR_peaks = find_peaks(filtered_signal)[0]
        RR_peaks_enhanced = find_peaks(peak_enhanced_signal)[0]

        # remove false peaks
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9056464/
        # The peak is considered to be a false peak if the R–R interval between two adjacent R-peaks is less than 30% of the mean R–R interval of the analysed signal, as mentioned in [22] and any false peak is dropped. The R–R intervals are calculated again after dropping any false peaks.
        average_RR_distance = np.mean(np.diff(RR_peaks))
        average_RR_distance_enhanced = np.mean(np.diff(RR_peaks))


        # re-find true peaks
        RR_peaks = find_peaks(filtered_signal, distance=average_RR_distance)[0]
        RR_peaks_enhanced = find_peaks(peak_enhanced_signal, distance=average_RR_distance_enhanced)[0]


        # get breath rate for entire signal using welches method. Get average breath length (peak to peak RR interval), and use that as num samples per segment
        # get RR windows
        RR_length = np.mean(np.diff(RR_peaks))
        RR_length_enhanced = np.mean(np.diff(RR_peaks_enhanced))
        
        if np.isnan(RR_length) or np.isnan(RR_length_enhanced):
            print("No Pleth peak detected. Returning Nan")
            if self.PLOTTING:
                fig, ax = plt.subplots(1,1, sharex="col")
                time = np.arange(len(pleth_signal))/sampling_freq
                ax.plot(time, pleth_signal)
                ax.set_xlabel('Time')
                ax.set_ylabel('Signal')
                ax.set_title('Raw Signal')  
                plt.tight_layout()
            return np.nan, np.nan

        # get power spectra for signals, using average rr window
        filtered_signal_welch_freq, filtered_signal_welch_psd = welch(filtered_signal, fs=sampling_freq, nperseg=RR_length)
        peak_enhanced_welch_freq, peak_enhanced_welch_psd = welch(peak_enhanced_signal, fs=sampling_freq, nperseg=RR_length_enhanced)

        max_power_arg = np.argmax(filtered_signal_welch_psd)
        max_power_freq = filtered_signal_welch_freq[max_power_arg]
        filtered_breath_rate = max_power_freq * 60 #hz. Gives max power frequency in terms of cycles per minute

        max_power_arg = np.argmax(peak_enhanced_welch_psd)
        max_power_freq = peak_enhanced_welch_freq[max_power_arg]
        enhanced_breath_rate = max_power_freq * 60 #hz. Gives max power frequency in terms of cycles per minute
        

        if self.PLOTTING:
            def calculate_fft(signal, sampling_freq):
                N = len(signal)  # Number of points in the signal
                fft_output = np.fft.fft(signal)  # Take FFT
                freqs = np.fft.fftfreq(N, 1/sampling_freq)  # Frequency bins
                
                # Only take the first half of the spectrum
                positive_freqs = freqs[:N//2]
                positive_fft_output = np.abs(fft_output[:N//2])  # Take the magnitude for plotting
                return positive_freqs, positive_fft_output
            
            raw_signal_freq, raw_signal_fft = calculate_fft(pleth_signal, sampling_freq)
            filtered_signal_freq, filtered_signal_fft = calculate_fft(filtered_signal, sampling_freq)
            peak_enhanced_signal_freq, peak_enhanced_signal_fft = calculate_fft(peak_enhanced_signal, sampling_freq)

            max_display_freq = 10
            max_freq_index = np.max(np.argwhere(max_display_freq > filtered_signal_freq))
            resp_freq_high_display = 2 * resp_freq_high
            resp_freq_high_display_index = np.max(np.argwhere(resp_freq_high_display > filtered_signal_freq))

            fig, ax = plt.subplots(3,2)
            time = np.arange(len(pleth_signal))/sampling_freq

            ax[0,0].plot(time, pleth_signal)
            ax[0,0].set_xlabel('Time')
            ax[0,0].set_ylabel('Signal')
            ax[0,0].set_title('Raw Signal')   

            ax[0,1].set_title("Raw Signal FFT")
            ax[0,1].plot(raw_signal_freq[1:max_freq_index], raw_signal_fft[1:max_freq_index]) #ignoring DC Component
            ax[0,1].set_xlabel('Frequency (Hz)')
            ax[0,1].set_ylabel('Magnitude')
            ax[0,1].set_title('Positive Frequency Spectrum')   

            ax[1,0].plot(time, filtered_signal)
            ax[1,0].set_xlabel('Time')
            ax[1,0].set_ylabel('Signal')
            ax[1,0].set_title('Filtered Signal')   
            ax[1,0].scatter(RR_peaks/sampling_freq, filtered_signal[RR_peaks], marker="x", c="r")

            ax[1,1].set_title("Filtered Signal FFT")
            ax[1,1].plot(filtered_signal_freq[1:resp_freq_high_display_index], filtered_signal_fft[1:resp_freq_high_display_index]) #ignoring DC Component
            ax[1,1].scatter(filtered_signal_freq[np.argmax(filtered_signal_fft[1:])], np.max(filtered_signal_fft[1:]), marker="x", c="r")
            ax[1,1].set_xlabel('Frequency (Hz)')
            ax[1,1].set_ylabel('Magnitude')
            ax[1,1].set_title('Positive Frequency Spectrum')  


            ax[2,0].plot(time, peak_enhanced_signal)
            ax[2,0].set_xlabel('Time')
            ax[2,0].set_ylabel('Signal')
            ax[2,0].set_title('Peak Enhanced Signal')   
            ax[2,0].scatter(RR_peaks_enhanced/sampling_freq, peak_enhanced_signal[RR_peaks_enhanced], marker="x", c="r")

            ax[2,1].set_title("Peak Enhanced FFT")
            ax[2,1].plot(peak_enhanced_signal_freq[1:resp_freq_high_display_index], peak_enhanced_signal_fft[1:resp_freq_high_display_index]) #ignoring DC Component
            ax[2,1].scatter(peak_enhanced_signal_freq[np.argmax(peak_enhanced_signal_fft[1:])], np.max(peak_enhanced_signal_fft[1:]), marker="x", c="r")
            ax[2,1].set_xlabel('Frequency (Hz)')
            ax[2,1].set_ylabel('Magnitude')
            ax[2,1].set_title('Positive Frequency Spectrum')   
            plt.tight_layout()

        return filtered_breath_rate 

    def calculate_PSD(self, pleth_signal, sampling_freq=62.4725):    # calculating power spectral density
        resp_freq_low = 0.1 #hz (6 breaths/min)
        resp_freq_high = 0.4 #hz (24 breaths/min)
        hr_freq_low = 0.5 # hz (30 beats/min)
        hr_freq_high = 3 #hz (180 beats/min)
        bands = {
            "resp": (resp_freq_low, resp_freq_high),
            "hr": (hr_freq_low, hr_freq_high)
        }

        

        def bandpower(x, fs, fmin, fmax):
            # https://stackoverflow.com/questions/44547669/python-numpy-equivalent-of-bandpower-from-matlab
            f, Pxx = scipy.signal.periodogram(x, fs=fs)
            ind_min = np.argmax(f > fmin) - 1
            ind_max = np.argmax(f > fmax) - 1
            return scipy.trapz(Pxx[ind_min: ind_max], f[ind_min: ind_max])

        power_bands = {}
        for band, (f_low, f_high) in bands.items():
            power = bandpower(pleth_signal, sampling_freq, f_low, f_high)
            power_bands[band] = power

        if self.PLOTTING:
            # calculating FFT
            N = len(pleth_signal)  # Number of points in the signal
            fft_output = np.fft.fft(pleth_signal)  # Take FFT
            freqs = np.fft.fftfreq(N, 1/sampling_freq)  # Frequency bins
            # Only take the first half of the spectrum
            positive_freqs = freqs[:N//2]
            positive_fft_output = np.abs(fft_output[:N//2])  # Take the magnitude for plotting

            # taking PSD
            frequency_components, PSD = scipy.signal.periodogram(pleth_signal, sampling_freq, scaling="density")

            fig, ax = plt.subplots(2)
            ax[0].set_title("FFT")
            ax[0].plot(positive_freqs[1:], positive_fft_output[1:]) #ignoring DC Component
            ax[0].set_xlabel('Frequency (Hz)')
            ax[0].set_ylabel('Magnitude')
            ax[0].set_title('Positive Frequency Spectrum')



            ax[1].set_title("PSD")
            ax[1].semilogy(frequency_components[1:], PSD[1:]) #ignoring DC Component
            plt.tight_layout()

        # for label in power_bands.keys():
            # print(label, power_bands[label])
        
        return power_bands

    def calculate_entropy(self, pleth_signal):
        #calculate probability distribution
        _, value_counts = np.unique(pleth_signal, return_counts=True)
        probabilities = value_counts / len(pleth_signal)

        # Calculate Shannon entropy
        entropy = -np.sum(probabilities * np.log2(probabilities))
        return entropy


if __name__ =="__main__":
    import wfdb
    ppg_metrics = PPG_METRICS()
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
    
    # stiffnes_index = ppg_metrics.calculate_SI(test_pleth_signal, test_height_meters, sampling_freq)
    # reflection_index = ppg_metrics.calculate_RI(test_pleth_signal)
    # augmentation_index = ppg_metrics.calculate_AI(test_pleth_signal)
    # inflection_point_area = ppg_metrics.calculate_IPA(test_pleth_signal)
    # prefusion_index = ppg_metrics.calculate_PI(test_pleth_signal)
    # prefusion_index = ppg_metrics.calculate_PVI(test_pleth_signal)
    # power_spectral_density = ppg_metrics.calculate_PSD(test_pleth_signal, sampling_freq)
    resp_rate = ppg_metrics.calculate_resp_rate(test_pleth_signal, sampling_freq)
    plt.show()
