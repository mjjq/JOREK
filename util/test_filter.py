import numpy as np

def filter_low_frequency(signal, cutoff_freq, sampling_rate):
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

# Example usage
# Generate a sample signal (e.g., sine wave)
sampling_rate = 1000  # Hz
duration = 1  # Second
frequency1 = 20  # Hz
frequency2 = 40  # Hz
t = np.linspace(0, duration, sampling_rate * duration)
signal = np.sin(2 * np.pi * frequency1 * t) + np.sin(2 * np.pi * frequency2 * t) 

# Apply the filter
cutoff_freq = 30  # Hz, adjust this value as needed
filtered_signal = filter_low_frequency(signal, cutoff_freq, sampling_rate)

# Plotting
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plt.plot(t, signal, label='Original Signal')
plt.plot(t, filtered_signal, label='Filtered Signal')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Original and Filtered Signal')
plt.legend()
plt.grid(True)
plt.show()
