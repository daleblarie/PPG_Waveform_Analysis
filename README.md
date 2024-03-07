## About the Different Implementations

Some common clinically useful metrics extracted from a PPG waveform.
There are a few versions here. One version is implemented in base python/numpy/scipy and is able to display plots. This is `Raw Python`. It is a class that contains the relevant calculation functions, and all functions are implemented.

There are 3 versions that use numba for JIT compliation to increase execution speed. They differ in how the code is structured. 

`Numba Class Python` is a class implemented in python, with parts of its methods @JIT compiled.

`Numba With some Nesting` is a series of functions that exist outside a class, with parts of the functions @JIT compiled.. These are imported and called on their own instead of being managed by an object.

`Numba No Nesting` is a series of functions that exist outside a class, with @JIT compiled functions, and no function declaration nesting. These are imported and called on their own instead of being managed by an object.

Finally there are two Ctypes implementations. Because of the effort and increased complexity required, only some of the functions are implemented. 
One is a `pure C++` implementation, where all functionality is replicated in C++.

The other is a `hybrid python/C++` implementation. Scipy functions and some other things are implemented in python because they are already optimized and very fast, while other parts are implemented in C++ for speed gains.

Overall, `Numba No Nesting` is the recommended implementation to use because it is overall the fastest, and most stable.
Use it by importing like so `from ppg_metrics_numba_nonesting import calculate_SI, calculate_RI, calculate_AI, calculate_IPA, calculate_PI, calculate_PVI,calculate_PSQI, calculate_ESQI, calculate_resp_rate, calculate_PSD, calculate_entropy`

# Speed Tests

These were speed tested ten times total on a data set of 90, one hour long PPG waveforms sampled at 62.5hz taken from the Mimic 3 waveform database (a total of 900 signals each with ~225,000 points)

## SPEEDUP RESULTS
| |Raw Python |Numba Class python |Numba With some Nesting |Numba No Nesting |Pure Ctypes |Ctypes/Python Hybrid 
| --- |	--- |	--- |	--- |	--- |	--- |	--- |
|SI, RI, AI, IPA, PI, PVI Calculation Total Time (seconds)|	706.67|	1606.48|	1000.55|	61.88|	420.7|	65.92|
|SI, RI, AI, IPA, PI, PVI SPEEDUP|	1|	0.4398872068|	0.7062815452|	11.42000646|	1.679748039|	10.72011529|
|Total Time(seconds)|	818.74|	1924.32	|1175.59	|157.98 | --| --|		
|Total Time (minutes)|	13.65|	32.07|	19.59|	2.63| --| --|		
|OVERALL SPEEDUP|	1|	0.4254697763|	0.696450293|	5.182554754		| -- | --|

![Speed Comparison of Different Implementations of PPG Calculations](https://github.com/daleblarie/PPG_Waveform_Analysis/assets/33942693/9369d4f4-2ae7-4981-a12c-01fa43c9913c)


## RAW TIME VALUES (seconds)		


| |Raw Python |Numba Class python |Numba With some Nesting |Numba No Nesting |Pure Ctypes |Ctypes/Python Hybrid 
| --- |	--- |	--- |	--- |	--- |	--- |	--- |
|SI CALCULATION|	15.1	|211.59|	194.27|	9.91|	100.4|	9.34|
|RI CALCULATION|	13.11	|209.94|	191.74|	9.92|	100.21|	9.36|
|AI CALCULATION|	13.31	|215.4|	198.13|	9.93|	99.35|	16.15|
|IPA CALCULATION|	59.83|	468.47|	394.73|	10.9|	99.69|	10|
|PI CALCULATION|	8.66|	122.85|	1.12|	0.45|	0.32|	0.33|
|PVI CALCULATION|	596.66|	378.23|	20.56|	20.77|	20.73|	20.74|
|PSQI CALCULATION|	20.16|	95.92|	84.44|	3.93| -- | --|		
|ESQI CALCULATION|	2.11	|115.17|	1.5	|1.51| --| --|		
|RESP RATE CALCULATION|	46.56|	58.17|	46.66|	47.57|		
|PSD CALCULATION|	37.94|	43.35|	37.25|	37.91| --| --|		
|Entropy CALCULATION|	5.3	|5.23	|5.19	|5.18|	

