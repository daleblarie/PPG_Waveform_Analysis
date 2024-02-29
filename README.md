Some common clinically useful metrics extracted from a PPG waveform.
One version is implemented in base python/numpy/scipy and is able to display plots. The other is setup to run with numba for JIT compliation to increase execution speed.



## RAW TIME VALUES (seconds)		


| -- | Raw Python | Numba Implementation |
|---|---|---|
| SI CALCULATION | 17.15	| 10.25 |
| RI CALCULATION | 15.03 | 10.15 |
| AI CALCULATION | 15.4	| 10.23 |

| IPA CALCULATION | 33.96 | 21.8 |
| PI CALCULATION | 8.71 |	0.45 |
| PVI CALCULATION	| 600.15 | 20.71 |
| PSQI CALCULATION |20.16 | 3.93 |
| ESQI CALCULATION | 2.11 | 1.51 |
| RESP RATE CALCULATION | 46.56 | 47.57 |
| PSD CALCULATION | 37.94 |	37.91 |
| Entropy CALCULATION | 5.3	| 5.18 |
--
SPEEDUP
Total Time(seconds)  |712.67|	79.03
Total Time (minutes)	|11.88	|1.32
OVERALL SPEEDUP	|1|	9.017714792

![Comparison of Implementation Version Speed](https://github.com/daleblarie/PPG_Waveform_Analysis/assets/33942693/d40c1902-43e2-46e7-aead-5dc9058afa7d)
