clear
clc

%% Input Signal
figure(1);
layout_input = tiledlayout(3,1);
title(layout_input, 'Input Plot - YIMING CHEN & RUIYI DENG', 'FontSize', 16);


% Input - Square Wave Signal
f0 = 10000; % asked period is .1 ms
T0 = 1/f0; % period
t_s = 0.001*T0;
samples = 6*T0/t_s + 1;
range_s = -3*T0:t_s:3*T0; % x-axis

% generate the sqaure wave
input_square = square(range_s*2*pi*f0, 50);

% plot the square wave
nexttile([1 1]);
plot(range_s,input_square,'LineWidth',2);
ax = gca;
set(ax,'Fontsize',12)
xlim([0 range_s(length(range_s))]);

title('Input - Time Domain')
xlabel('Time (s)'); ylabel('Voltage (V)');
%% apply Fourier Transform on the Square Wave
N=100; % number of harmonics
nvec = -N:N;
c_in = zeros(size(nvec));
for n = nvec
    m = n+N+1;
    if (mod(n,2))
      c_in(m) = sinc(n/2); % fourier transformation of the sqaure wave
    else
      c_in(m) = 0.0;
    end
end
% parameter
f = nvec*f0; % frequency vector
mag = abs(c_in);
phase = angle(c_in);

%% plot of magnitude
nexttile;
stem(f,mag,'LineWidth',2);


axis([-8*f0 8*f0 0 max(mag)])

ax = gca;
set(ax,'Fontsize',12)

title('Magnitude Spectrum')
xlabel('Frequency (Hz)'); ylabel('Amplitude');
%% plot of frequency
nexttile;
stem(f,phase,'LineWidth',2);

axis([-8*f0 8*f0 -pi pi])

Ha = gca;
set(Ha,'Fontsize',12)

title('Phase Spectrum')
xlabel('Frequency (Hz)'); ylabel('Phase');





%% Output Signal
figure(2);
layout_input = tiledlayout(3,1);
title(layout_input, 'Onput Plot - YIMING CHEN & RUIYI DENG', 'FontSize', 16);

%% Implementation of the 2nd-order ButterWorth Filter
fc = 1.2*10e3; % cutoff frequency found from the previous section
Hf = 1 ./ (1+1.414*(1i*f/fc) + (1i*f/fc).^2); % transfer function
c_out = c_in .* Hf; % Fourier coefficients of the filter output

%% Fourier Coeff
A = zeros(2*N+1,ceil(samples));
for n = nvec
    m=n+N+1;
    A(m,:) = c_out(m) .* exp(1i*2*pi*n*f0*range_s);
end
gp_out = sum(A);

%% Plot the output signal 
nexttile([1 1]);
plot(range_s,input_square,'r',range_s,real(gp_out),'g','LineWidth',2);

ax = gca;
set(ax,'Fontsize',12)
title('Filter Input and Output - Time Domain')
xlabel('Time (s)'); ylabel('Voltage (V)');
set(ax,'Fontsize',12)
legend('Input','Output')

%% plot of magnitude
nexttile;
stem(f,abs(c_in),'r','LineWidth',2);
hold on
stem(f,abs(c_out),'g','LineWidth',2);
hold off
axis([-8*f0 8*f0 0 max(abs(c_in))])

ax = gca;
set(ax,'Fontsize',12)

title('Magnitude Spectrum')
xlabel('Frequency (Hz)'); ylabel('Amplitude');

%% plot of frequency
phase_o = angle(c_out);
nexttile;
stem(f,phase,'r','LineWidth',2);
hold on
stem(f,phase_o,'g','LineWidth',2);
hold off
axis([-8*f0 8*f0 -pi pi])

ax = gca;
set(ax,'Fontsize',12)

title('Phase Spectrum')
xlabel('Frequency (Hz)'); ylabel('Phase');



%% end



















