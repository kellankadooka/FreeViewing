%BUTTERWORTH FILTER ON POSITIONAL DATA
SF = 30; % sampling frequency
NF = SF/2; % Nyquist frequency
CF = 9; % Cut-off frequency
% initalize normalized cut-off frequency Wn with a value between 0 and 1
Wn = CF/NF; % == 9Hz/30Hz = 0.3
% run butter
[b,a] = butter(2, Wn, 'low'); % 2nd order, low-pass filter
%filter data using filtfilt (zero-phase digital filtering)
var_filtered = filtfilt(b,a,var_median);

plot(var_median)
hold on
plot(var_filtered)