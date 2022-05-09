clear
hold off
format long e
N = 4096; %No. of FFT samples
sampling_rate = 1000.0e3; %unit Hz
tstep = 1/sampling_rate;
tmax = N*tstep/2;

tmin = -tmax;
tt = tmin:tstep:tmax-tstep;
fmax = sampling_rate/2; 
fmin = -fmax;
fstep = (fmax-fmin)/N;
freq = fmin:fstep:fmax-fstep;

fc=64e3;
Ac = 1;
ct=Ac*cos(2*pi*fc*tt);
fm = 4e3;
mt = cos(2*pi*fm*tt);
st = mt.*ct;


% Plot of resulting signal when s(t) is again multiplied by cos(2pifct)
%st = mt.*ct.*cos(2*pi*fc*tt);

figure(1)
Hp1 = plot(tt,ct);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('Carrier c(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Carrier : Time domain');
axis([-1e-5 1e-5 -1.1 1.1])
pause(1)

figure(2)
Hp1 = plot(tt,mt);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('message  m(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('message signal : Time domain');
axis([-0.001 0.001 -1.1 1.1])
pause(1)
figure(3)
Hp1 = plot(tt,st);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('s(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('DSB-SC modulated wave : Time domain');
axis([-10/fc 10/fc min(st) max(st)])
pause(1)
Mf = fftshift(fft(fftshift(mt)))/(2*fmax);
figure(4)
%The amplitude of the spectrum is different from the Fourier transform
%amplitude due to discretization of discrete Fourier transform
Hp1=plot(freq,abs(Mf));
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Frequency (Hz) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('|M(f)|');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Spectrum of the message signal');
axis ([-15e3 15e3 0 max(abs(Mf))])
pause(1)


Sf = fftshift(fft(fftshift(st)))/(2*fmax);
figure(5)
Hp1=plot(freq,abs(Sf));
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Frequency (Hz) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('|S(f)|');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Spectrum of the DSB-SC wave S(f)');
axis ([-150e3 150e3 0 max(abs(Sf))])
pause(1)

%DSB-SC demodulation

%Local oscillator at the receiver perfectly synchronized
thet=0;
lo = cos(2*pi*fc*tt + thet); 
st1 = st .* lo;
figure(6)
Hp1=plot(tt,st1);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('  s hat(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('signal after remodulation at Rx, s hat(t) : Time domain');
axis([-10/fc 10/fc min(st1) max(st1)])
pause(1)
Sf1 = fftshift(fft(fftshift(st1)))/(2*fmax);
figure(7)
Hp1=plot(freq,abs(Sf1));
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Frequency (Hz) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('|S hat(f)|');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Spectrum S hat(f)');
axis ([-50e3 50e3 0 max(abs(Sf1))])
pause(1)
%Low pass filtering
f_cutoff = 30e3;
%ideal low pass filter
n=1;
for f = freq
    if abs(f) < f_cutoff
        Hf(n) = 1;
    else
        Hf(n) = 0;
    end
n=n+1;
end
Mf1 = Sf1 .* Hf;
mt1 = 2*fmax*fftshift(ifft(fftshift(Mf1)));
figure(8)
Hp1=plot(tt,mt1,'r',tt,mt*0.5,'g.');
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('m hat(t)(Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Output of low pass filter,  m hat(t)  : Time domain');
axis([-0.001 0.001 min(mt*0.5) max(mt*0.5)])
legend('LPF output', 'message sig');

% DSBSC signal with carrier power 50% of power of two SBs
figure(9)
st2 = (sqrt((1/2)/2) + mt).*ct;
Hp2 = plot(tt,st2);
set(Hp2,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('s(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('DSB-SC modulated wave : Time domain(carrier=0.5)');
axis([-10/fc 10/fc min(st2) max(st2)])
pause(1)
Sf2 = fftshift(fft(fftshift(st2)))/(2*fmax);
figure(10)
Hp2=plot(freq,abs(Sf2));
set(Hp2,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Frequency (Hz) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('|S(f)|');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Spectrum of the DSB-SC wave S(f)');
axis ([-150e3 150e3 0 max(abs(Sf2))])
pause(1)

% DSBSC signal with carrier power 300% of power of two SBs
figure(11)
st3 = (sqrt(3/2) + mt).*ct;
Hp3 = plot(tt,st3);
set(Hp3,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('s(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('DSB-SC modulated wave : Time domain(carrier=sqrt(3/2)');
axis([-10/fc 10/fc min(st3) max(st3)])
pause(1)
Sf3 = fftshift(fft(fftshift(st3)))/(2*fmax);
figure(12)
Hp3=plot(freq,abs(Sf3));
set(Hp3,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Frequency (Hz) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('|S(f)|');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Spectrum of the DSB-SC wave S(f)');
axis ([-150e3 150e3 0 max(abs(Sf3))])
pause(1)

% DSBSC signal with carrier power 500% (more than three times) of power of two SBs
figure(13)
st4 = (sqrt(5/2) + mt).*ct;
Hp4 = plot(tt,st4);
set(Hp4,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('s(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('DSB-SC modulated wave : Time domain(carrier=sqrt(5/2))');
axis([-10/fc 10/fc min(st4) max(st4)])
pause(1)
Sf4 = fftshift(fft(fftshift(st4)))/(2*fmax);
figure(14)
Hp4=plot(freq,abs(Sf4));
set(Hp4,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Frequency (Hz) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('|S(f)|');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Spectrum of the DSB-SC wave S(f)');
axis ([-150e3 150e3 0 max(abs(Sf4))])
pause(1)