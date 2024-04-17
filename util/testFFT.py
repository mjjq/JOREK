#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy import interpolate

def high_freuqncy_pass(signal, cutoff_freq, sampling_rate):
    # Step 1: Apply FFT to the signal
    fft_result = np.fft.fft(signal)
    
    # Step 2: Frequency domain representation
    # Generate frequency axis
    n = len(signal)
    freqs = np.fft.fftfreq(n, d=1/sampling_rate)
    
    # Identify the indices corresponding to frequencies
    low_freq_indices = np.where(np.abs(freqs) < cutoff_freq)[0]
    
    # Step 3: Remove or attenuate low-frequency components
    # Set low-frequency components to zero
    fft_result[low_freq_indices] = 0
    
    # Step 4: Apply inverse FFT
    filtered_signal = np.fft.ifft(fft_result)
    
    return filtered_signal.real  # Return real part of the filtered signal

def low_freuqncy_pass(signal, cutoff_freq, sampling_rate):
    # Step 1: Apply FFT to the signal
    fft_result = np.fft.fft(signal)
    
    # Step 2: Frequency domain representation
    # Generate frequency axis
    n = len(signal)
    freqs = np.fft.fftfreq(n, d=1/sampling_rate)
    
    # Identify the indices corresponding to frequencies
    low_freq_indices = np.where(np.abs(freqs) > cutoff_freq)[0]
    
    # Step 3: Remove or attenuate low-frequency components
    # Set low-frequency components to zero
    fft_result[low_freq_indices] = 0
    
    # Step 4: Apply inverse FFT
    filtered_signal = np.fft.ifft(fft_result)
    
    return filtered_signal.real  # Return real part of the filtered signal

def window_fft(timeInput, signalInput, window_size):
    data_len = len(signalInput)
    if data_len != len(timeInput):
        print("Error, lengths mismatch, len(time) = ", len(timeInput), ", len(signal) = ", len(signalInput))
        exit(-1)
    elif data_len == 0:
        print("Error, zero length, len(time/signal) = ", len(timeInput))
        exit(-1)
    elif window_size > data_len or window_size == 0:
        print("Error, window is zero or too large, window_size = ", window_size, ", len(time/signal) = ", len(timeInput))
        exit(-1)

    norm_factor = window_size
    time = np.linspace(np.min(timeInput), np.max(timeInput), data_len)
    f_interp = interpolate.interp1d(timeInput, signalInput)
    signal = f_interp(time)

    id_mid  = int((window_size+1)/2-1)
    freq_1d=np.fft.fftfreq(window_size, time[1] - time[0])[:id_mid + 1]
    time_1d = time[:-(window_size-1)] #time[window_size-1:]
    spec_2d = []
    for i in tqdm(range(data_len-window_size+1), desc="Processing window FFT"):
        sample=signal[i:i+window_size]
        spec=np.fft.fft(sample)/norm_factor
        spec=np.abs(spec[:id_mid + 1])
        spec[1:-1] = 2.e0 * spec[1:-1]
        spec_2d.append(spec)

    freq_2d, times_2d = np.meshgrid(freq_1d, time_1d)
    spec_2d = np.array(spec_2d)
    return time, signal, times_2d, freq_2d, spec_2d

if __name__ == "__main__":
    n_len = 1000
    time = np.linspace(0, 10, n_len)  # Example time points from 0 to 10 with 1000 samples
    f1 = 2.0
    f2 = 5.0
    signal = np.sin(2 * np.pi * f1 * time) + 0.5 * np.sin(2 * np.pi * f2 * time)

    time_new, signal_new, times_2d, freq_2d, spec_2d = window_fft(time, signal, window_size = 100)
    
    # Plot the signal
    plt.figure(figsize=(10, 7))  # Adjust the figure size as needed
    plt.subplot(2, 1, 1)
    plt.plot(time, signal, label='Time Domain Signal (original)')
    plt.plot(time_new, signal_new, '.', label='Time Domain Signal (interpolated)')
    plt.xlim([min(time), max(time)])
    plt.legend()
    plt.ylabel('Signal Amplitude')
    plt.grid(True)
    plt.subplot(2, 1, 2)
    plt.contourf(times_2d, freq_2d, spec_2d, cmap='viridis')
    plt.xlim([min(time), max(time)])
    plt.ylim(0,10)
    plt.xlabel('Time')
    plt.ylabel('Frequency')
    plt.colorbar(orientation='horizontal', shrink=0.5)
    plt.tight_layout()
    plt.show()    

