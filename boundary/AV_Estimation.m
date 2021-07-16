function [IQ, V, OMEGA] = AV_Estimation(dz, WaveLength, InitSignal, InitTime, RefSignal, RefTime, TranSignal, TranTime, RecDist, Freq, PlotFlag, Vp)
% function estimates the attenuation (inverse quality factor)
% and phase velocity in considered layer
% using recorded traces before and after this layer

% Input variables:
% dz - grid step
% WaveLength - size to define domain sizes - �?СПРАВ�?ТЬ
% InitSignal - initial signal
% RefSignal - reflected signal (to estimate transmission losses)
% TranSignal - transmitted signal
% InitTime, RefTime, TranTime - corresponding time for each signal
% Rec1 - coordinate of first receiver line
% Rec2 - coordinate of second receiver line
% Freq - central frequency of the initial signal
% PlotFlag - 1, if you want to plot fourier transform results etc
% Vp - background velocity

% Output variables:
% OMEGA - frequencies array
% IQ - Inverse quality factor estimation
% V - Phase velocity estimation

%% Prepare variables
w0 = Freq*2*pi;
dw = w0/500;
w = [w0/2:dw:w0*2]; % frequency range, where we recover the attenuation and velocity
dt = InitTime(2) - InitTime(1); % time step

%% Plot the signals
if (PlotFlag == 1)
    figure;
    plot(InitTime, InitSignal, 'b');
    hold on;
    plot(TranTime, TranSignal, 'r');
    hold on;
    plot(RefTime, RefSignal, 'k');
    grid on;
    title('Signals in time domain');
    xlabel('Time, s');
    legend('Initial signal', 'Transmitted signal', 'Reflected signal');
end

%% correct the time to avoid background pieces between the receivers
% here we shift times for receiver 2 on time, needed for wave to propagate
% through homogeneous part of the media between receivers
% also we shift reflected signal in time to correctly estimate reflection coefficient

HomogDist = 3*WaveLength*dz;
Delta_t = HomogDist/real(Vp);
RecDist = RecDist - HomogDist;
[M,AmpInd1] = max(abs(InitSignal));
[M,AmpInd2] = max(abs(TranSignal));
EffVel = RecDist/(TranTime(AmpInd2)-Delta_t-InitTime(AmpInd1));

%% Plot the signals
% if (PlotFlag == 1)
%     figure;
%     plot(InitTime, InitSignal, 'b');
%     hold on;
%     plot(TranTime, TranSignal, 'r');
%     hold on;
%     plot(RefTime, RefSignal, 'k');
%     grid on;
%     title('Signals in time domain (Background-corrected)');
%     xlabel('Time, s');
%     legend('Initial signal', 'Transmitted signal', 'Reflected signal');
% end

%% Apply Fourier transform
InitFT(1:length(w)) = 0;
RefFT(1:length(w)) = 0;
TranFT(1:length(w)) = 0;
for k = 1:length(w)
    InitFT(k) = InitSignal.'*exp(1i*w(k)*InitTime)*dt; % Fourier image for initial signal
    TranFT(k) = TranSignal.'*exp(1i*w(k)*TranTime)*dt; % Fourier image for transmitted signal
    RefFT(k) = RefSignal.'*exp(1i*w(k)*RefTime)*dt;    % Fourier image for reflected signal
end

%% Plot Fourier images
if (PlotFlag == 1)
    figure;
    plot(w/2/pi,real(InitFT),'b');
    hold on;
    plot(w/2/pi,real(TranFT),'r');
    hold on;
    plot(w/2/pi,real(RefFT),'k');
    grid on;
    title('Real part of Fourier images');
    xlabel('Frequency, Hz');
    legend('Initial signal', 'Transmitted signal', 'Reflected signal');    
end
if (PlotFlag == 1)
    figure;
    plot(w/2/pi,imag(InitFT),'b');
    hold on;
    plot(w/2/pi,imag(TranFT),'r');
    hold on;
    plot(w/2/pi,imag(RefFT),'k');
    grid on;
    title('Imaginary part of Fourier images');
    xlabel('Frequency, Hz');
    legend('Initial signal', 'Transmitted signal', 'Reflected signal');    
end
if (PlotFlag == 1)
    figure;
    plot(w/2/pi,abs(InitFT),'b');
    hold on;
    plot(w/2/pi,abs(TranFT),'r');
    hold on;
    plot(w/2/pi,abs(RefFT),'k');
    grid on;
    title('Absolute value of Fourier images');
    xlabel('Frequency, Hz');
    legend('Initial signal', 'Transmitted signal', 'Reflected signal');    
end

%% for homog test - complex background velocity
BackgroundVelocity = real(Vp);
dist1 = 2*WaveLength*dz;
dist2 = WaveLength*dz;
TransmissionCorrection = exp(-1i*(dist1+dist2)*w.*BackgroundVelocity./(abs(BackgroundVelocity).^2));
TranFT = TranFT.*TransmissionCorrection;

%% for homog test - Elastic Reflection coefficient
% ReflectionCorrection = exp(-1i*2*dist1*w.*BackgroundVelocity./(abs(BackgroundVelocity).^2));
% RefFT = RefFT.*ReflectionCorrection;
% RefCoef = RefFT.*conj(InitFT)./(abs(InitFT).^2); % here we divide complex by real
% InitFT = InitFT.*(1-RefCoef.^2);

%% Plot Fourier images
if (PlotFlag == 1)
    figure;
    plot(w/2/pi,real(InitFT),'b');
    hold on;
    plot(w/2/pi,real(TranFT),'r');
    hold on;
    plot(w/2/pi,real(RefFT),'k');
    grid on;
    title('Real part of Fourier images');
    xlabel('Frequency, Hz');
    legend('Initial signal', 'Transmitted signal', 'Reflected signal');    
end
if (PlotFlag == 1)
    figure;
    plot(w/2/pi,imag(InitFT),'b');
    hold on;
    plot(w/2/pi,imag(TranFT),'r');
    hold on;
    plot(w/2/pi,imag(RefFT),'k');
    grid on;
    title('Imaginary part of Fourier images');
    xlabel('Frequency, Hz');
    legend('Initial signal', 'Transmitted signal', 'Reflected signal');    
end
if (PlotFlag == 1)
    figure;
    plot(w/2/pi,abs(InitFT),'b');
    hold on;
    plot(w/2/pi,abs(TranFT),'r');
    hold on;
    plot(w/2/pi,abs(RefFT),'k');
    grid on;
    title('Absolute value of Fourier images');
    xlabel('Frequency, Hz');
    legend('Initial signal', 'Transmitted signal', 'Reflected signal');    
end

% % plot frequency-dependent reflection coefficient
% if (PlotFlag == 1)
%     figure;
%     plot(w,real(RefCoef),'b');
%     hold on;
%     plot(w,imag(RefCoef),'r');
%     hold on;
%     plot(w,abs(RefCoef),'g');
%     title('Reflection Coefficient');
% end

% % plot corrected initial signal
% for k = 1:length(InitTime)
%     InitSignalCorrected(k) = InitFT*exp(-1i*w.'*InitTime(k))*dw; % Fourier image for initial signal
% end

%% Compute the phase to correct the estimation

% compute the cross-correltation
CC = InitFT.*conj(TranFT)./(InitFT.*conj(InitFT));

if (PlotFlag == 1)
    figure;
    plot(w/2/pi, real(CC), 'b', w/2/pi, imag(CC), 'r');
    grid on;
    title('Cross-correlation in frequency domain');
    legend('Real part', 'Imaginary part');
end

% compute the phase
%% here i tried to recover the algorithm from the VichMet paper
% Phase = unwrap(angle(TranFT./InitFT));
Phase = -unwrap(angle(CC),[],2);

PhasePreEst = w/EffVel*RecDist;
CorFlag = 0;

% plot phases comparison before correction
if (PlotFlag == 1)
    figure;
    plot(w/2/pi, Phase, 'b', w/2/pi, PhasePreEst, 'k');
    hold on;
end

while (CorFlag == 0)
    bb = norm(Phase - PhasePreEst);
    bbp = norm(Phase - PhasePreEst + 2*pi);
    bbm = norm(Phase - PhasePreEst - 2*pi);
    if (bb <= bbp)&&(bb <= bbm)
        CorFlag = 1;
    end
    if (bbp < bb)&&(bbp < bbm)
        Phase = Phase + 2*pi;
    end
    if (bbm < bb)&&(bbm < bbp)
        Phase = Phase - 2*pi;
    end
end

% plot phase after correction
if (PlotFlag == 1)
    plot(w/2/pi, Phase, 'r');
    grid on;
    title('Phase before and after correction');
    legend('Before', 'Pre-estimate', 'After');
end

%% obtain estimations
wmat = diag(w);
Re_s = RecDist./Phase*wmat; 
Im_s = (log(abs(InitFT)./abs(TranFT))/RecDist)*inv(wmat); % obtain imaginary part of slowness (some kind of inverted)
V = 1./(1./Re_s-1i*Im_s); % obtain complex phase velocity

%%
% plot velocity
if (PlotFlag == 1)
    figure;
    plot(w/2/pi, real(V));
    grid on;
    title('Real part of the velocity');
    xlabel('frequency, Hz');
    ylabel('\Re V, m\s');
    
    figure;
    plot(w/2/pi, imag(V));
    grid on;
    title('Imaginary part of the velocity');
    xlabel('frequency, Hz');
    ylabel('\Im V, m\s');
end

VV = V.*V;
V = real(V);
IQ = imag(VV)./real(VV);
OMEGA = w/2/pi;

% plot the inversed quality factor
if (PlotFlag == 1)
    figure;
    plot(OMEGA, IQ, 'b');
    grid on;
    title('Inversed quality factor');
    xlabel('Frequency, Hz');
    ylabel('Q^{-1}');
end

end