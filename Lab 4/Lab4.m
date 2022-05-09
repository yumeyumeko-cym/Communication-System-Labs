clear 
format long e
tend = 10;
tbeg = -10;
N=100000;
tstep = (tend-tbeg)/N;
sampling_rate = 1/tstep;

tt = tbeg:tstep:tend-tstep;

%% Experiment 1
yt1 = load('lab4_num_expt1');
lag1 = [100 200 500];

for i = 1:3
    figure(i)
    tiledlayout(1,2);
    
    maxlag = lag1(i);
    %Autocorrelation of yt
    Ry  = xcorr(yt1.yt,yt1.yt,maxlag);
    %tau vector
    tau_vec = -(maxlag*tstep):tstep:maxlag*tstep;
    %Abs. PSD corresponding to yt
    Sy = abs(fftshift(fft(fftshift(Ry))));
    %define the frequency vector corresponding to tau_vec
    Ntau = length(tau_vec);
    %Nyquist sampling rate
    fmax = sampling_rate/2; 
    fmin = -fmax;
    fstep = (fmax-fmin)/Ntau;
    %Frequency window
    freq = fmin:fstep:fmax-fstep;
    
    nexttile;
    plot(tau_vec, Ry);
    title("Autocorrelation");
    subtitle("max lag = " + maxlag);
    xlabel("\tau (s)");
    ylabel("R_y");
    
    nexttile;
    plot(freq, Sy);
    title("Absolute PSD");
    subtitle("max lag = " + maxlag);
    xlabel("Frequency (Hz)");
    ylabel("|S_y|");
 
end

%% Experiment 2
yt2 = load('lab4_num_expt2');
lag2 = [100 200 20000];

for i = 4:6
    figure(i);
    tiledlayout(1,2);

    maxlag = lag2(i-3);
    %Autocorrelation of yt
    Ry  = xcorr(yt2.yt,yt2.yt,maxlag);
    %tau vector
    tau_vec = -(maxlag*tstep):tstep:maxlag*tstep;
    %Abs. PSD corresponding to yt
    Sy = abs(fftshift(fft(fftshift(Ry))));
    %define the frequency vector corresponding to tau_vec
    Ntau = length(tau_vec);
    %Nyquist sampling rate
    fmax = sampling_rate/2; 
    fmin = -fmax;
    fstep = (fmax-fmin)/Ntau;
    %Frequency window
    freq = fmin:fstep:fmax-fstep;

    nexttile;
    plot(tau_vec, Ry);
    title("Autocorrelation R_y");
    subtitle("max lag = " + maxlag);
    xlabel("\tau (s)");
    ylabel("R_y");

    nexttile;
    plot(freq, Sy);
    title("Absolute PSD");
    subtitle("max lag = " + maxlag);
    xlabel("Frequency (Hz)");
    ylabel("|S_y|");


end

% Signal in time domain
figure(7);
plot(tt, yt2.yt);

xlim([-100*tstep 100*tstep]);
title("Signal yt(n): Time Domain");
subtitle("$y(t)$", 'interpreter', 'latex')
xlabel("Time (s)");
ylabel("y(t)");

% Signal in frequency domain
figure(8);

Nyt = length(yt2.yt);
% Nyquist sampling rate
fmax = sampling_rate/2; 
fmin = -fmax;
fstep = (fmax-fmin)/Nyt;
% Frequency window
freq = fmin:fstep:fmax-fstep;

plot(freq, abs(fftshift(fft(fftshift(yt2.yt)))));
title("Signal yt(n): Magnitude Spectrum");
xlabel("Frequency (Hz)");
ylabel("Magnitude");

%% Experiment 3
value3 = load('lab4_num_expt3');
xt3 = value3.xt;
yt3 = value3.yt;

figure(9);
tiledlayout(1,2);

nexttile;
plot(tt, xt3);
title("Signal x: Time Domain");
subtitle("$x(t)$", 'interpreter', 'latex');
xlabel("Time (s)");
ylabel("x(t)");

nexttile;
plot(tt, yt3);
title("Signal y: Time Domain");
subtitle("$y(t)$", 'interpreter', 'latex');
xlabel("Time (s)");
ylabel("y(t)");


maxlag = 20000;
R_xy = xcorr(yt3,xt3,maxlag);
tau_vec = -(maxlag*tstep):tstep:maxlag*tstep;
%Abs. PSD corresponding to yt
S_xy = abs(fftshift(fft(fftshift(R_xy))));
%define the frequency vector corresponding to tau_vec
Ntau = length(tau_vec);
%Nyquist sampling rate
fmax = sampling_rate/2; 
fmin = -fmax;
fstep = (fmax-fmin)/Ntau;
%Frequency window
freq = fmin:fstep:fmax-fstep;

fig = figure(10);
plot(tau_vec, R_xy);
title("Autocorrelation R_{xy}");
subtitle("max lag = " + maxlag)
xlabel("\tau (s)");
ylabel("R_{xy}");

































