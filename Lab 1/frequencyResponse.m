clear
clc
%% set parameter of ButterWorth Filter
% n = 2; % order of the filter
fc = 1.2*10e3; % cutoff frequency
f0 = 10e3; % fundamental frequency
f_o = 3.5*10e3; % output frequency, 23dB below 

% List of the parameter for ButterWorth Filter from Wikipedia
% n = 1: Bn = (s+1)
% n = 2: Bn = s^{2}+1.414s+1
% n = 3: Bn = (s+1)(s^{2}+s+1)
% n = 4: Bn = (s^{2}+0.7654s+1)(s^{2}+1.8478s+1)
% n = 5: Bn = (s+1)(s^{2}+0.6180s+1)(s^{2}+1.6180s+1)

% for denominator
a = [1 1.414 1];
 
% for nominator
b = 1;


%% Build the plot of ButterworthFilter
step = linspace(0,8);
h = freqs(b,a,step); % a for poles, b for zeros, step for normalized freq range from 0 to 8 rad/s
mag = 20*log10(abs(h)); % convert magnitude to dB

% plot frequency response
figure(1)

% frequency response

plt = plot(fc*step,mag);
xlim([0 60000]);

grid on
title('Frequency Response (Magnitude)');
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
xline(fc,'b',{'Cutoff frequency'}); % 3dB frequency (at cutoff frequency)
xline(f0,'r',{'Fundamental frequency'}); % fundamental frequency
xline(f_o,'g',{'Output frequency'}); % 23dB below


%% Verification of the constraints
% obtain the interpolation points
cq = interp1(fc*step, mag, [fc f0 f_o]);

message = ['Attenuation at cutoff frequency: ', num2str(cq(1)), ' dB'];
disp(message);

message = ['Attenuation at fundamental frequency: ', num2str(cq(2)), ' dB'];
disp(message);

message = ['Attenuation at output harmonic: ', num2str(cq(3)), ' dB'];
disp(message);

if (abs(cq(2)) < 2)
    flag = 'Yes';
else
    flag = 'No';
end


message = ['The loss of the fundamental frequency due to filtering is less than 2 dB: ', flag];
disp(message);

message = ['Output of harmonics is ', num2str(-(cq(3)-cq(2))), ' dB below fundamental frequency'];
disp(message);

if (abs(cq(3)-cq(2)) > 13.5)
    flag = 'Yes';
else
    flag = 'No';
end

message = ['Does the design meet the required attentuation from filtering (at least 13.5 (23-9.5)dB): ', flag];
disp(message);

















