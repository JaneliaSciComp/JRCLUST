function outdata = low_pass_resample(indata, cutoff_freq, sample_rate)
% LOW_PASS_RESAMPLE is an example of how to use the Matlab filters and
% resampling to modify the raw Grapevine data which is sampled at 30kHz
% and minimally filtered.  This should only be run on raw data.
% usage: LOW_PASS_RESAMPLE(X, CUTOFF_FREQ [Hz], SAMPLE_RATE[Hz])

% This parameter requires that this function only be run on raw data.
RAW_SAMPLE_RATE = 30000; % Hz
% Calculate the Nyquist freqency for raw data
Nq = RAW_SAMPLE_RATE / 2;
% create filter parameters for a 4th order low pass butterworth filter.
[z, p] = butter(4, cutoff_freq/Nq, 'low');
% apply this filter to the raw data.
% filtfilt applies the filter twice (one forward and back) and results
% in no phase difference.
filtered_data = filtfilt(z, p, indata);
% you may also apply filter instead.
% filtered_data = filter(z, p, indata);
% resample the result using the new sampling rate.
outdata = resample(filtered_data, sample_rate, RAW_SAMPLE_RATE);
