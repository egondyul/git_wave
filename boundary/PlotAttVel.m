 %function PlotAttVel
% function estimates inversed quality factor and phase velocity
% for the set of central frequencies and plot figures with these
% estimations

% clear all;
% below we consider our .mat-data format:
% out - <NFreq x NCase struct>
% out(j,i).freq - central frequency of the initial signal
% out(j,i).src - source line coordinate (in meters)
% out(j,i).rec1 - first receiver line coordinate (in meters)
% out(j,i).rec2 - second receiver line coordinate (in meters)
% out(j,i).LayerWidth - length of the considered layer of interest (homogeneous, layered, fractured etc.)
% out(j,i).time - time vector for both traces from receiver lines 1 and 2
% out(j,i).trace1 - trace recorded at the first receiver line
% out(j,i).trace2 - trace recorded at the second receiver line
% out(j,i).Vp - phase velocity in the background
% out(j,i).interface - start coordinate of the layer of interest
% out(j,i).WaveLength - size needed to define z-length of the domain
% out(j,i).dz - grid step in z-direction
MatName='1.mat';
load(MatName);

NFreq = size(out,1); % number of considered central frequencies
NCase = size(out,2); % number of considered cases (connectivity, permeability, etc.)

% create figures for attenuation and velocity
FVel = figure;
hold on;
FiQ = figure;
hold on;

% set colors to plot different cases
ccM=[0 0 0; 0 0 1; 0 1 0; 1 0 0; 1 1 0; 1 0 1; 0 1 1; 0.5 0.5 0.5; 1 0.5 0.5; 0.5 1 0.5];

% set logarithmic scale for figures
figure(FVel);
set(gca,'XScale','Log');
grid on;
figure(FiQ);
set(gca,'XScale','Log');
grid on;
%% Main cycle
for j=1:NFreq
    j  % print number of frequency
    FigTrace(j) = figure; % create plot for the set of traces for one central frequency
    hold on;
    grid on;
    for i=1:NCase
        i % print number of case for considered frequency
        cc = ccM(i,:); % set color for considered data to plot
        clear Time Trace1 Trace2;
        Time = out(j,i).time;
        Dt = Time(3)-Time(2); % time step - random neighbouring time values
        Trace1 = out(j,i).trace1;
        Trace2 = out(j,i).trace2;
        Nu0 = out(j,i).freq;
        RecDist = out(j,i).rec2 - out(j,i).rec1; % the distance between receivers (in meters)
        
        % plot traces for one central frequency
        % shift every next trace from the previous one at 1.5*amplitude
        figure(FigTrace(j));
        plot(Trace1+max(abs(Trace1))*i*1.5,Time,'b')  
        plot(Trace2+max(abs(Trace1))*i*1.5,Time,'r','Linewidth',2);
        TicksTrace(i,1)=max(abs(Trace1))*i*1.5; % remember ticks for traces plot
    
        %% Distinguish the signals - initial, transmitted and reflected
%        [M,MaxInd1] = max(abs(Trace1)); % indeces where maximal amplitude is reached for 1st and 2nd traces
        [M,MaxInd2] = max(abs(Trace2));
%         TravelTime = Time(MaxInd2) - Time(MaxInd1); % time spent by the wave to travel from receiver line 1 to receiver line 2
%         EffVel = RecDist/TravelTime; % Effective velocity of the wave between receiver lines

     % distinguish the initial signal
         TravelTime=abs(out(j,i).rec1-out(j,i).src)/out(j,i).Vp+3/out(j,i).freq;
        [M,IndTravelTime]=min(abs(Time-TravelTime));
         Alph=2; %2 period that grab the region
         InitSignal=Trace1(IndTravelTime-floor(Alph/Nu0/Dt):min(IndTravelTime + floor(Alph/Nu0/Dt), length(Trace1)));
         InitTime = Time(IndTravelTime - floor(Alph/Nu0/Dt):min(IndTravelTime + floor(Alph/Nu0/Dt), length(Trace1)));
%
        % distinguish the initial signal
%         Alph = 2; % number of periods, which we want to leave from the amplitude (minimal - 0.5)
%         InitSignal = Trace1(MaxInd1 - floor(Alph/Nu0/Dt):min(MaxInd1 + floor(Alph/Nu0/Dt), length(Trace1)));
%         InitTime = Time(MaxInd1 - floor(Alph/Nu0/Dt):min(MaxInd1 + floor(Alph/Nu0/Dt), length(Trace1)));
% %         

        % distinguish the transmitted signal
%         TranlTime=abs(out(j,i).rec2-out(j,i).src)/out(j,i).Vp+3/out(j,i).freq;
%         [M,IndTranTime]=min(abs(Time-TranlTime));
        Alph = 2;
        TranSignal = Trace2(MaxInd2 - floor(Alph/Nu0/Dt):min(MaxInd2 + floor(Alph/Nu0/Dt), length(Trace2)));
        TranTime = Time(MaxInd2 - floor(Alph/Nu0/Dt):min(MaxInd2 + floor(Alph/Nu0/Dt), length(Trace2)));      
        
        % distinguish the reflected signal
        Alph = 2;
%         src = out(j,i).src; % coordinate of the source line (in meters) - for new mat-files
%         RefToRecTime = 3/Nu0 + (out(j,i).interface - 2*src + out(j,i).rec1)/out(j,i).Vp;
        RefToRecTime = 3/Nu0 + ((out(j,i).interface - out(j,i).src) + (out(j,i).interface - out(j,i).rec1))/out(j,i).Vp;
        [M,IndRefTime]=min(abs(Time-RefToRecTime));
        RefSignal=Trace1(IndRefTime-floor(Alph/Nu0/Dt):min(IndRefTime + floor(Alph/Nu0/Dt), length(Trace1)));
        RefTime = Time(IndRefTime - floor(Alph/Nu0/Dt):min(IndRefTime + floor(Alph/Nu0/Dt), length(Trace1)));
%         [MM,MMax] = abs(Time - RefToRecTime);
%         [M,MaxIndRef] = max(abs(Trace1(min(MaxInd1 + floor(Alph/Nu0/Dt), length(Trace1)):length(Trace1))));
%         MaxIndRef = MaxIndRef + length(1:min(MaxInd1 + floor(Alph/Nu0/Dt), length(Trace1)));
%         RefSignal = Trace1(MaxIndRef - floor(Alph/Nu0/Dt):min(MaxIndRef + floor(Alph/Nu0/Dt),length(Trace1)));
%         RefTime = Time(MaxIndRef - floor(Alph/Nu0/Dt):min(MaxIndRef + floor(Alph/Nu0/Dt),length(Trace1)));        
                
        %% Estimation of the velocity and attenuation
        % main function for signal deconvolution
        PlotFlag = 1;
        [RES(j,i).IQ, RES(j,i).V, RES(j,i).OMEGA] = AV_Estimation(out(j,i).dz, out(j,i).WaveLength, InitSignal, InitTime, RefSignal, RefTime, TranSignal, TranTime, RecDist, Nu0, PlotFlag, out(j,i).Vp);
        
        % here we choose the frequency from our array closest to the central one
        [M,ClosestFreq] = min(abs(RES(j,i).OMEGA-Nu0));
        Frequency(j) = RES(j,i).OMEGA(ClosestFreq);
        Velocity(j,i) = RES(j,i).V(ClosestFreq);
        InverseQualityFactor(j,i) = RES(j,i).IQ(ClosestFreq);
        
        % plot the velocities and attenuation in frequency ranges of estimation
        figure(FVel);
        plot(RES(j,i).OMEGA,RES(j,i).V,'Color',cc);
        figure(FiQ);
        plot(RES(j,i).OMEGA,RES(j,i).IQ,'Color',cc);
    end
    
    figure(FigTrace(j));
    title(['\nu_0 = ' num2str(Nu0)], 'FontSize', 14);
    xlabel('Case','FontSize',14);
    ylabel('Time, s','FontSize',14);
    set(gca,'XTick',TicksTrace)
    set(gca,'XTickLabel',char('1', '2', '3', '4', '5', '6', '7', '8'));
    clear TicksTrace;
    set(gca,'FOntSize',14);
end

% plot the overall estimations
ResFigV = figure;
set(gca,'XScale','Log');
grid on;
hold on;

ResFigIQ = figure;
set(gca,'XScale','Log');
grid on;
hold on;

% plot the estimations
for i=1:NCase
    cc=ccM(i,:);
    figure(ResFigV);
    plot(Frequency, Velocity(:,i)','v','Color',cc);
    legend;
    figure(ResFigIQ);
    plot(Frequency, InverseQualityFactor(:,i)','v','Color',cc);
    legend;
end

% plot splines
SplineFrequency = (min(Frequency):(max(Frequency)-min(Frequency))/1000:max(Frequency));
for i=1:NCase
    cc=ccM(i,:);
    figure(ResFigV);
    SplineVelocity = spline(Frequency,Velocity(:,i)',SplineFrequency);
    plot(SplineFrequency,SplineVelocity,'Color',cc,'Linewidth',2);
    figure(ResFigIQ);
    SplineIQF = spline(Frequency,InverseQualityFactor(:,i)',SplineFrequency);
    plot(SplineFrequency,SplineIQF,'Color',cc,'Linewidth',2);
end

figure(ResFigV);
set(gca,'FontSize',14);
xlabel('Frequency, Hz');
ylabel('Velocity, m/s');

figure(ResFigIQ);
set(gca,'FontSize',14);
xlabel('Frequency, Hz');
ylabel('Q^{-1}');