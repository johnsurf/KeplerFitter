clear all; 
close all; 
clc;
restoredefaultpath;
CurrentPath = pwd;
Utilities = fullfile(CurrentPath,'Utilities')
%genpath(Utilities)
addpath(Utilities)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iUnits = 0
units = Units(iUnits);
% For Canonical Units
TU = units.TU;
DU = units.DU;
VU = units.VU;
AU = units.AU;
mu = units.mu;
twopi = units.twopi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%mu     =  1.0;         % Canonical Graviational Parameter = GM
%Rearth =  1.0;
%J2     =  1082.6269e-6;   %J2 pertubative constant
%J3     =      -2.30e-6;
%J4     =      -2.12e-6;
%J5     =      -0.20e-6;
%J6     =       1.00e-6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input / Output
%Input_filepath = 'InputDirectory';
filename       = 'TwoGoldenEyes'
%filename       = 'ThreeGoldenEyes'
%filename       = 'FourGoldenEyes'
in_dir         = pwd
out_dir        = fullfile(in_dir,'los_files')
%fid            = fopen(fullfile(Input_filepath, filename),'r');
if ~isfolder(out_dir)
    mkdir(out_dir);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Master Clock Start Time: 
%tEpoch          = 200/TU;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time   = [];
iFirst = 0;
for tSeconds = 600:900
    time = [time; tSeconds];    
end
TimeSize  = size(time);
timeSorted = sortrows(time);
if iFirst == 0
    tFirst   = timeSorted(1);
end
iMid           = floor(numel(time)/2);
TFit           = time(iMid)
%TFit            = median(time); 
tEpoch          = TFit/TU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iKepler = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a Sensor Orbit called Golden Eye1
[GoldenEye1, Kepler] = GenOrbit(tEpoch, units);
iKepler = [iKepler; Kepler];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a Sensor Orbit called Golden Eye 2
[GoldenEye2, Kepler] = GenOrbit(tEpoch, units);
iKepler = [iKepler; Kepler];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a Sensor Orbit called Golden Eye 3
%[GoldenEye3, Kepler] = GenOrbit(tEpoch, units);
%iKepler    = [iKepler; Kepler];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a Sensor Orbit called Golden Eye 4
%[GoldenEye4, Kepler] = GenOrbit(tEpoch, units);
%iKepler    = [iKepler; Kepler];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a Target called "The Thief"
[Thief, KeplerThief] = GenOrbit(tEpoch, units);
iKepler = [iKepler; KeplerThief];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KeplerThief

% Let us Digitize the Theif's Orbit's LOS vector with respect to the GoldenEye Sensors  

icnt = 0;
imax = 0;
imin = 0;

time = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ISensor    = 1:10;   % Allows for Sensor Satellite Numbering. 
track_data = [];
iFirst     = 0;

%for tSeconds = 600:PeriodGoldenEye1/4000.0:900
%for tSeconds = 0:PeriodGoldenEye1/4000.0:400
for tSeconds = 600:900
    
    %timeArray = tEpoch*TU + 300*rand(3,1);
    %timeArray = 300*rand(3,1);
    %time      = timeArray(1); 
    
    time = [time, tSeconds];
    icnt = icnt + 1;
    ISat = 0;
    
    ISat      = ISat + 1;
    timeMeas = tSeconds;
    %timeMeas = timeArray(1);
    %timeMeas = tSeconds + 2.0*(0.5 - rand(1))*50.;
    
    [row] = NineColumn(GoldenEye1, Thief, timeMeas/TU, ISensor(ISat)); 
    % %[rSensor, vSensor] = extrapolate(GoldenEye1, timeMeas/TU);
    % [rSensor, vSensor] = extrapolate(GoldenEye1, timeMeas/TU);
    % [rThief, vThief]   = extrapolate(Thief,   timeMeas/TU);
    % los   = rThief - rSensor;
    % range = norm(los);
    % los   = los/range;
    % row   = [ISensor(ISat), timeMeas, rSensor*DU, los, range*DU];
    track_data = [track_data; row];
    
    ISat     = ISat     + 1;
    timeMeas = tSeconds + 5;
    %timeMeas = timeArray(2);
    %timeMeas = tSeconds + 2.0*(0.5 - rand(1))*50.;
    
    [row] = NineColumn(GoldenEye2, Thief, timeMeas/TU, ISensor(ISat)); 
    track_data = [track_data; row];
    
%     ISat     = ISat     + 1;
%     timeMeas = tSeconds + 10;
%     %timeMeas = timeArray(3);
%     %timeMeas = tSeconds + 2.0*(0.5 - rand(1))*50.;
%     
%     [row] = NineColumn(GoldenEye3, Thief, timeMeas/TU, ISensor(ISat)); 
%     track_data = [track_data; row];
% 
%     ISat     = ISat     + 1;
%     timeMeas = tSeconds + 10;
%     %timeMeas = timeArray(3);
%     %timeMeas = tSeconds + 2.0*(0.5 - rand(1))*50.;
%     
%     [row] = NineColumn(GoldenEye4, Thief, timeMeas/TU, ISensor(ISat));
%     track_data = [track_data; row];
    
end

TimeArray = size(time);
timeSorted = sortrows(track_data,2);
if iFirst == 0
    tFirst   = timeSorted(1,2);
    iFirst = -1;
end
track_data      = timeSorted;
%track_data(:,2) = track_data(:,2) - tFirst;
%iMid           = floor(numel(track_data(:,2))/2);
%TFit           = tEpoch*TU
%TFit           = track_data(iMid,2)
%TFit            = median(track_data(:,2))

saveName = sprintf('%s_los_data.txt',strrep(filename,'.txt',''));
save(fullfile(out_dir,saveName),'track_data','-ascii','-double');
%save(fullfile(out_dir,'track_data.txt'),'track_data','-ascii','-double');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check Out Polynomial Fits   Sensor Orbits
state_vector_sat = cell(9,1);
kepler_fit_sat   = cell(6,1);
state_vector_los = cell(9,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%plot3(track_data(:,3), track_data(:,4),track_data(:,4))
Satellites = unique(track_data(:,1));
ChiSqSat      = [];
epsilon = 0.0000000001*ones(6,1);
Observed = [zeros(3,1)];
losM     = 0.000001*ones(3,1);
InvCovUV    = eye(2)*10^(8); 
for iSat=1:numel(Satellites)
    Index      = find(track_data(:,1) == Satellites(iSat));
    TrackOrbit = track_data(Index,:);
    plot3(TrackOrbit(:,3),TrackOrbit(:,4),TrackOrbit(:,5),'linewidth',2);
    hold on
    rows       = size(TrackOrbit,1);
    Iorder     = 4;
    %Iorder     = 1;
    Object     = 1; % for Sensor
    iType      = 1;  % iType is hard-wired at the moment
    %TFit       = median(TrackOrbit(:,2))
    %TFit       = 0.;
    [state_vector_fit, fitSat, Iorder, ChiSq] = PolyFit(Object, TrackOrbit,TFit,Iorder,iType);
    plot3(state_vector_fit(1,:),state_vector_fit(2,:),state_vector_fit(3,:),'linewidth',2);
    hold on
    Diff = state_vector_fit(1:3,:) - TrackOrbit(:,3:5)';
    ChiSqSat = [ChiSqSat, sum(sum(Diff.*Diff))];
    state_vector_sat{iSat} = PolyExtrapToFitTime(fitSat, TFit);
    time  = TFit/TU;
    Rpos =  state_vector_sat{iSat}(1:3)/DU;
    Rdot =  state_vector_sat{iSat}(4:6)/VU;
    [e a Inclination omega Omega Mp Period Jacobian] = kepler(time, Rpos, Rdot, mu);
    Kepler  = [e a Inclination omega Omega Mp]
    KeplerSat  = LagrangePlanetary(Kepler, units);
    Sensor   = [Rpos+epsilon(1:3); zeros(6,1)];
    [Rinit,Vinit,ParamList,Jacobian,Hessian,GravityCan] = KeplerSat.OrbitDerivatives(time, Sensor, losM, InvCovUV);    
    
    deltaPar     = zeros(6,1);
    ChiSqArray   = [];
    TangentArray = [];
    InvCovR = eye(3);
    deltaPar   = zeros(6,1);
    ChiSqArray = [];
    InvCovR = eye(3);
    
    Mind  = 6;
    oparm = zeros(Mind,1);
    
    for NR = 1:10
        
        nfree    = 0;
        oparm(1) = true;
        oparm(2) = true;
        oparm(3) = true;
        oparm(4) = true;
        oparm(5) = true;
        oparm(6) = true;
        
        IndexRows = find(oparm == 1);
        IndexCols = [IndexRows; Mind+1];
        
        for ipar=1:Mind
            if oparm(ipar) == true
                nfree = nfree+1;
            end
        end
        
        ChiSq  = 0.0;
        betaNR  = zeros(6,1);
        hessNR  = zeros(6);
        derivChiSq   = zeros(3,1);
        hessianChiSq = zeros(3);
        KeplerSat  = LagrangePlanetary(Kepler, units);
        rSatKepler = [];
        for ii =1:rows
            time   = TrackOrbit(ii,2)/TU;
            Rdata  = [TrackOrbit(ii,3:5)']/DU;
            Sensor = [Rdata+epsilon(1:3); zeros(6,1)];
            %Sensor = [Rdata; zeros(6,1)];
            % Measurement Angles:
            %losM    = state_vector_losArray(1:3,ii);
            losM     = TrackOrbit(ii,6:8)';
            losM     = losM/norm(losM);   % los
            rangeM   = norm(losM);
            thetaM   = acos(losM(3)/rangeM);
            phiM     = atan2(losM(2), losM(1));
            Observed =[thetaM; phiM];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % [Rinit,Vinit,ParamList,Jacobian,Hessian,GravityCan] = KeplerSat.OrbitDerivatives(time, Sensor, losM, InvCovUV);
            % deltaR   = Rinit - Rdata;
            % rSatKepler = [rSatKepler, Rinit];
            % %deltaR   = Rinit - Rdata + Jacobian(1:3,2:7)'*delZ;
            % ChiSq    = ChiSq  + 0.5*deltaR'*InvCovR*deltaR;
            % betaNR   = betaNR + Jacobian(1:3,2:7)'*InvCovR*deltaR;
            % hessPart = zeros(6);
            % for k = 1:3
            %     hessPart = hessPart + Hessian{k}(2:7,2:7)*deltaR(k);
            % end
            % %hessNR  = hessNR + Jacobian(1:3,2:7)'*InvCovR*Jacobian(1:3,2:7);
            % hessNR  = hessNR + hessPart + Jacobian(1:3,2:7)'*InvCovR*Jacobian(1:3,2:7);
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % [Rinit,Vinit,ParamList] = OrbitAtTime(KeplerSat, time, Sensor, losM);
            % deltaR   = Rinit - Rdata;
            %
            % Hessian = cell(6, 6);
            % Jacobian = [];
            % %LineCheck    = 50;                % Rx
            % [H_I, Jacob, Hess] = WorkOrder(KeplerSat, 50);
            % Hessian{1}  = Hess;
            % Jacobian      = [Jacobian; Jacob];   % row 1
            % %LineCheck    = 51;                % Ry
            % [H_I, Jacob, Hess] = WorkOrder(KeplerSat, 51);
            % Hessian{2}  = Hess;
            % Jacobian      = [Jacobian; Jacob];   % row 1
            % %LineCheck    = 52;                % Rz
            % [H_I, Jacob, Hess] = WorkOrder(KeplerSat, 52);
            % Hessian{3}  = Hess;
            % Jacobian      = [Jacobian; Jacob];   % row 1
            %
            % rSatKepler = [rSatKepler, Rinit];
            % %deltaR   = Rinit - Rdata + Jacobian(1:3,2:7)'*delZ;
            % ChiSq    = ChiSq  + 0.5*deltaR'*InvCovR*deltaR;
            % betaNR   = betaNR + Jacobian(1:3,2:7)'*InvCovR*deltaR;
            % hessPart = zeros(6);
            % for k = 1:3
            %     hessPart = hessPart + Hessian{k}(2:7,2:7)*deltaR(k);
            % end
            % %hessNR  = hessNR + Jacobian(1:3,2:7)'*InvCovR*Jacobian(1:3,2:7);
            % hessNR  = hessNR + hessPart + Jacobian(1:3,2:7)'*InvCovR*Jacobian(1:3,2:7);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [Rinit,Vinit,ParamList] = OrbitAtTime(KeplerSat, time, Sensor, losM, InvCovUV);
            deltaR   = Rinit - Rdata;
            % Approximate Gravity Perturbations: 
            %deltaR    = Rinit - Rdata + Jacobian(1:3,2:7)'*delZ;
            %norm(deltaR) 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [H_I, Jacob, Hess] = WorkOrder(KeplerSat, 105);            
            betaNR     = betaNR + Jacob(2:7)';
            hessNR     = hessNR + Hess(2:7,2:7);
            %H_I
            % rSatKepler = [rSatKepler, Rinit];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ChiSq      = ChiSq  + 0.5*deltaR'*InvCovR*deltaR;
            %rSatKepler = [rSatKepler, Rinit];
        end
        ChiSq
        ChiSqArray  =  [ChiSqArray, log(ChiSq)];
        %deltaPar    = -hessNR\betaNR
        
        AugmentedMatrix = [hessNR,-betaNR];
        AugmentedMatrix = AugmentedMatrix(IndexRows, IndexCols);
        hitzero   = 0;
        opzero = 0.0000000000001;
        %disp("Entering NRMin ");
        pinc = NRMin(AugmentedMatrix, nfree, opzero, hitzero);
        if(hitzero == 1)
            disp("  hit operational zero it NRMin  ")
        end
        ii = 1;
        stepSize = 0.0;
        for i=1:Mind
            if oparm(i)
                delp(i)   = pinc(ii);
                ii = ii + 1;
                stepSize = stepSize + delp(i)*delp(i);
            end
        end
        
        TangentMag   = log(norm(betaNR))
        TangentArray = [TangentArray, TangentMag];
        %Kepler      =  Kepler + deltaPar'
        Kepler      =  Kepler + delp;
    end
    disp([' GoldenEye ', num2str(iSat),' Kepler      = ', num2str(Kepler)])
    disp([' Generated Kepler Params = ', num2str(iKepler(iSat,:))])
    kepler_fit_sat{iSat}  = LagrangePlanetary(Kepler, units);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check Out Polynomial Fits   Line of Sight Polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%plot3(track_data(:,3), track_data(:,4),track_data(:,4))
Satellites = unique(track_data(:,1));
ChiSqSat      = [];
for iSat=1:numel(Satellites)
    Index      = find(track_data(:,1) == Satellites(iSat));
    TrackOrbit = track_data(Index,:);
    plot3(TrackOrbit(:,6),TrackOrbit(:,7),TrackOrbit(:,8),'linewidth',2);
    hold on
    rows       = size(TrackOrbit,1);
    Iorder     = 3;
    %Iorder     = 1;
    Object     = 2; % for los
    iType      = 1;  % iType is hard-wired at the moment
    %TFit       = median(TrackOrbit(:,2))
    %TFit       = 0.;
    [state_vector_fit, fitLos, Iorder, ChiSq] = PolyFit(Object, TrackOrbit,TFit,Iorder,iType);
    plot3(state_vector_fit(1,:),state_vector_fit(2,:),state_vector_fit(3,:),'linewidth',2);
    hold on
    Diff = state_vector_fit(1:3,:) - TrackOrbit(:,6:8)';
    ChiSqSat = [ChiSqSat, sum(sum(Diff.*Diff))];
    state_vector_los{iSat} = PolyExtrapToFitTime(fitLos, TFit);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ChiSqMin    = 10^100;
KeplerAve   = [];
KeplerChiSq = [];
for iSat=1:numel(Satellites)
    ChiSqSatMin = 10^100;
    % time        = TFit/TU;
    % Rpos        = state_vector_sat{iSat}(1:3)/DU;
    % Rdot        = state_vector_sat{iSat}(4:6)/VU;
    % % comment out 506 to see what happens. !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % [e a Inclination omega Omega Mp Period Jacobian] = kepler(time, Rpos, Rdot, mu);
    % KeplerSensor = [e a Inclination omega Omega Mp]; % Canonical Units
    % KeplerSat  = LagrangePlanetary(KeplerSensor, units);
    KeplerSat  = kepler_fit_sat{iSat};
    
    [StateVectorCell, RootsLaplace] = LaplaceInECI(state_vector_sat{iSat}, state_vector_los{iSat});
    for icell = 1:size(StateVectorCell,2)        
        StateVector = StateVectorCell{icell};
        time        = TFit/TU;
        Rpos        = StateVector(1:3)/DU;
        Rdot        = StateVector(4:6)/VU;
        ChiSq       = 0.0;
        % comment out 506 to see what happens. !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        [e a Inclination omega Omega Mp Period Jacobian] = kepler(time, Rpos, Rdot, mu);
        KeplerTest = [e a Inclination omega Omega Mp]; % Canonical Units
        Imaginary  = numel(find(imag(KeplerTest)));
        
        if Imaginary == 0
            KeplerObj  = LagrangePlanetary(KeplerTest, units);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Index = find(track_data(:,1) == iSat);
            track_Sat = track_data(Index,:);
            rows       = size(track_Sat,1);
            for icell = 1:size(StateVectorCell,2)
                StateVector = StateVectorCell{icell};
                time        = TFit/TU;
                Rpos        = StateVector(1:3)/DU;
                Rdot        = StateVector(4:6)/VU;
                ChiSq       = 0.0; 
                % comment out 506 to see what happens. !!!!!!!!!!!!!!!!!!!!!!!!!!!!
                [e a Inclination omega Omega Mp Period Jacobian] = kepler(time, Rpos, Rdot, mu);
                KeplerTest = [e a Inclination omega Omega Mp]; % Canonical Units
                
                Imaginary  = numel(find(imag(KeplerTest)));
                
                if Imaginary == 0
                    KeplerObj  = LagrangePlanetary(KeplerTest, units);
                    for ii =1:rows
                        tRecord = track_Sat(ii,2)/TU;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %losM    = state_vector_losArray(1:3,ii);
                        losM     = track_Sat(ii,6:8)';
                        losM     = losM/norm(losM);   % los
                        rangeM   = norm(losM);
                        thetaM   = acos(losM(3)/rangeM);
                        phiM     = atan2(losM(2), losM(1));
                        Observed = [thetaM; phiM];
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % %Sat Position at measured position
                        Rsensor  = [track_Sat(ii,3:5)']/DU;
                        Sensor   = [Rsensor; zeros(6,1)];
                        
                        % %Sat Position from Polynomial at TFit
                        % Sensor = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU]
                        
                        % %Sat Position from Orbit Fit to Kepler elements:
                        %[Rsensor, Vsensor, ParamList, Jacobian, Hessian, GravityCan] = KeplerSat.OrbitDerivatives(tRecord, Sensor, losM);
                        [Rsensor, Vsensor, ParamList] = OrbitAtTime(KeplerSat, tRecord, Sensor, losM, InvCovUV);
                        Sensor   = [Rsensor; Vsensor; zeros(3,1)];
                        %[Rinit,Vinit,ParamList,Jacobian,Hessian,GravityCan] = KeplerObj.OrbitDerivatives(tRecord, Sensor, losM);
                        [Rinit,   Vinit,   ParamList] = OrbitAtTime(KeplerObj, tRecord, Sensor, losM, InvCovUV);
                        los       = [ParamList(7); ParamList(8); ParamList(9)];
                        los       = los/norm(los);
                        AngleDiff = abs(1 - dot(los,losM));
                        ChiSq     = ChiSq + AngleDiff;
                    end
                    if ChiSq < ChiSqSatMin
                        ChiSqSatMin      = ChiSq;
                        KeplerSatBest    = KeplerTest;
                        StateVectorAve   = StateVector;
                    end
                    if ChiSq < ChiSqMin
                        ChiSqMin         = ChiSq;
                        KeplerBest       = KeplerTest;
                        StateVectorECI   = StateVector;
                    end
                end
            end
        end
    end
    if ChiSqSatMin < 10^(20)
        KeplerChiSq = [KeplerChiSq; ChiSqSatMin];
        KeplerAve   = [KeplerAve;   KeplerSatBest];
    end
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %      Begin Finite Differences Check
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %losM    = state_vector_losArray(1:3,ii);
% losM     = track_Sat(iMid,6:8)';
% losM     = losM/norm(losM);   % los
% rangeM   = norm(losM);
% thetaM   = acos(losM(3)/rangeM);
% phiM     = atan2(losM(2), losM(1));
% Observed = [thetaM; phiM];
% % Sat Position at measured position
% Rsensor  = [track_Sat(iMid,3:5)']/DU;
% Sensor   = [Rsensor; zeros(6,1)];
% KeplerObj  = LagrangePlanetary(KeplerBest, units);
% [Rextrap,Vextrap,ParamList] = OrbitAtTime(KeplerObj, tEpoch, Sensor, losM);
% for LineCheck = 1:131
%     [Jacob, JacobFinite, Hess, HessFinite] = DigitalJustice(KeplerObj, LineCheck, Sensor, losM);
%     Jacob
%     JacobFinite
%     Hess
%     HessFinite
%     Sum = sum(sum(Hess - HessFinite));
%     disp([' Line = ', num2str(LineCheck),'  sum = ', num2str(Sum)])
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%KeplerAve = sum(KeplerAve)/size(KeplerAve,1)
KeplerChiSq
KeplerAve
NumChiSq = size(KeplerChiSq,1);
KeplerCombined = zeros(1, size(KeplerAve,2))
InvCov         = 0.0;
for i = 1:NumChiSq
    CovInv = 1.0/KeplerChiSq(i);
    InvCov  = InvCov + CovInv;
    KeplerCombined = KeplerCombined + CovInv*KeplerAve(i,:);    
end
KeplerCombined = (1.0/InvCov)*KeplerCombined
% KeplerAveECI   = LagrangePlanetary(KeplerCombined, units);
% [rInitial, vInitial] = extrapolate(KeplerAveECI, TFit/TU);
% StateVectorECI     = [rInitial'*DU; vInitial'*VU; zeros(3,1)]
disp([' KeplerBest  = ', num2str(KeplerBest)])
disp([' KeplerThief = ',num2str(KeplerThief)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[KeplerFit1, CovFit1, TFit, SV7Fit1] = LeastSquaresAllSatellites(track_data, TFit, 1, 0, StateVectorECI, InvCovUV, units)
%[KeplerFit1, CovFit1, TFit, SV7Fit1] = LeastSquares(track_data, tFirst, 1, 0, zeros(7,1) )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[KeplerFit1, CovFit1, TFit, SV7Fit1] = LeastSquaresAllSatellites(track_data, TFit, 1, 0, StateVectorECI )
%[KeplerFit1, CovFit1, TFit, SV7Fit1] = LeastSquares(track_data, tFirst, 1, 0, zeros(7,1) )
