clear
hold off
format long e

% time samples
N = 2^16; %No. of FFT samples
sampling_rate = 40e4; %unit Hz
tstep = 1/sampling_rate;
tmax = N*tstep/2;
tmin = -tmax;
tt = tmin:tstep:tmax-tstep;

%freq samples
fmax = sampling_rate/2; 
fmin = -fmax;
fstep = (fmax-fmin)/N;
freq = fmin:fstep:fmax-fstep;

%% Modulation

%carrier
fc = 20e3;
Ac = 1;
ct = Ac*cos(2*pi*fc*tt);

Tm = 0.0005;
mt = -2*sinc(tt/Tm); % message signal

%% Q1: Plotting the message signal

tiledlayout(1,2);

% time domain
nexttile;
plot(tt, mt, 'LineWidth', 2);
tim_dom_ax = gca;
set(tim_dom_ax);
xlabel('Time (s)');
ylabel('m(t) (V)');
title('Time Domain of m(t)');
axis([-2e-3 2e-3 min(mt) max(mt)]);

% frequency domain
Mf1 = fft(fftshift(mt));
Mf = fftshift(Mf1);
abs_Mf = abs(Mf);

nexttile;
plot(freq, abs_Mf, 'LineWidth', 2);
freq_dom_ax = gca;
set(freq_dom_ax);
xlabel('Frequency (Hz)');
ylabel('|M(f)|');
title('Frequency Domain of m(t)');
axis ([-5e3 5e3 0 max(abs(Mf))]);


% Measure highest frequency component from positive side of spectrum
[max_val, max_freq_index] = max(abs_Mf((length(abs_Mf)/2+1):end));
highest_freq = freq(max_freq_index + length(freq)/2);

fprintf('highest freqeuncy = %f\n', highest_freq)




%% Q2: Mod and Demod
mt_max = max(abs(mt));

ka = 0.5/mt_max; % 50% modulation
% AM modulation
st = (1+ka*mt).*ct;

% Plot the modulated signal
modulated_sig(st, mt, ka, tt, freq, 2, ka);


%% Q2 (i) - (iii)

% Plotting output of envelope detector and output of DC removal for time constant RC = 1/fc
RC = 1/fc;
yt = envelope_detect(st, RC, tt);
output_fig(yt, mt, ka, tt, 3, 0.5, "1/fc");

% Plotting output of envelope detector and output of DC removal for time constant RC = 10*Tm
RC = 10*Tm;
yt = envelope_detect(st, RC, tt);
output_fig(yt, mt, ka, tt, 4, 0.5, "10*Tm");

% Plotting output of envelope detector and output of DC removal for time
% constant RC = 1.5*Tm
RC = 1.5*Tm;
yt = envelope_detect(st, RC, tt);
output_fig(yt, mt, ka, tt, 5, 0.5, "1.5*Tm");




%% Q3: 200% Modulation
%For 200% modulation
ka=2/mt_max;

%AM signal
st = (1+ka*mt).*ct;

modulated_sig(st, mt, ka, tt, freq, 6, 2);

RC = 0.5*(Tm + (1/fc));
yt = envelope_detect(st, RC, tt);
output_fig(yt, mt, ka, tt, 7, 2, "0.5(T_m + 1/f_c)");



%% Function Definition
% Function definitions in a script must appear at the end of the file.
% The following are the functions used in lab 2.
function modulated_sig(signal, message, ka, tt, freq_vector, fig_num, percent_modulation)
    figure(fig_num)
    % Plotting the modulated signal
    tiledlayout(2,1);

    % time domain plotting s(t)
    nexttile;
    plot(tt, signal, 'LineWidth', 2);
    
    % plotting envelope
    hold on
    plot(tt, ka*message+1,'Color', 'r', 'LineStyle', '--', 'LineWidth', 2); % +1 for DC
    plot(tt, -ka*message-1,'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    hold off
    
    
    legend('Modulated Signal', 'Envelope');
    tim_dom_ax = gca;
    set(tim_dom_ax);
    xlabel('Time (s)');
    ylabel('Modulated Signal s(t) (V)');
    
    title_name = "Modulated Signal in Time Domain (" + percent_modulation*100 + "% Modulation)";
    title(title_name);
    
    axis([-2e-3 2e-3 min(signal) max(signal)]);

    % frequency domain
    Sf1 = fft(fftshift(signal));
    Sf = fftshift(Sf1);

    nexttile;
    plot(freq_vector, abs(Sf), 'LineWidth', 2);
    freq_dom_ax = gca;
    set(freq_dom_ax);
    xlabel('Frequency (Hz)');
    ylabel('|M(f)|');
    title_name = "Magnitude Spectrum of the Modulated Signal (" + percent_modulation*100 + "% Modulation)";
    title(title_name);
    axis([-25e3 25e3 0 max(abs(Sf))]);

end


function yt = envelope_detect(st, RC, tt)
    yt = st;
    n = 1;
    for t = tt
        if (n > 1)
            if (yt(n-1) > st(n))
                yt0 = yt(n-1);
                % time when C starts discharging
                tc = tt(n-1);
                yt(n) = yt0*exp(-(t-tc)/RC);
            end
        end
        n = n+1;
    end

end

function output_fig(signal, original, ka, time_vector, fig_num, percent_modulation, RC_name) 
    % Parameter lpf is used to control the Low-Pass Filter
    figure(fig_num);
    
    % To specify the figure for q(iii)

    tlayout = tiledlayout(1,2);

    
    title_name = "Output signals for " + percent_modulation*100 + "% Modulation with RC = " + RC_name;
    title(tlayout, title_name);

    % output of envelope detector
    nexttile;
    plot(time_vector, signal,'LineWidth',2);
    envelope_det_ax = gca;
    set(envelope_det_ax);
    xlabel('Time (s)');
    ylabel('y(t) (V)');
    title('After the envelope detector');
    axis([-2e-3 2e-3 0 max(signal)]);

    % dc removal and ka scale removal
    yt1 = (signal - 1) / ka;

    nexttile;
    plot(time_vector, yt1,'g', time_vector, original,'k','LineWidth',2);
    legend('after DC removal','message signal');
    output_signal_ax = gca;
    set(output_signal_ax);
    xlabel('Time (s)');
    ylabel('y1(t) (V)');
    title('After the DC removal');
    axis([-2e-3 2e-3 min(original) max(original)]);
      
   
end

































