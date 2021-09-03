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
twopi = units.twopi;
% For Canonical Units
TU = units.TU;
DU = units.DU;
VU = units.VU;
AU = units.AU;
mu = units.mu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input / Output
%Input_filepath = 'InputDirectory';
%filename       = 'BatchLOSHistories'
%filename       = 'TwoGoldenEyes'
%filename       = 'ThreeGoldenEyes'
filename       = 'FourGoldenEyes'
%filename       = 'FiveGoldenEyes'
in_dir         = pwd
data_dir        = fullfile(in_dir,'los_files');
%data_dir        = fullfile(in_dir,'Data_Favorites');
%fid            = fopen(fullfile(Input_filepath, filename),'r');
InputName      = sprintf('%s_los_data.txt',strrep(filename,'.txt',''));
DataInputFile  = fullfile(data_dir, InputName);
track_data     = load(DataInputFile);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Master Clock Start Time:
%tEpoch          = 200/TU;
timeSorted = sortrows(track_data,2);
track_data      = timeSorted;
time   = track_data(:,2);
iFirst = 0;
if iFirst == 0
    tFirst   = track_data(1,2);
    iFirst = -1;
end
TimeSize  = size(time);
timeSorted = sortrows(time);
%track_data(:,2) = track_data(:,2) - tFirst;
%iFit           = floor(numel(track_data(:,2))/2);
%TFit           = tEpoch*TU
%TFit           = track_data(iFit,2)
%TFit            = median(track_data(:,2)

iFit           = floor(numel(time)/2);
%TFit          = median(time);
%tEpoch         = TFit/TU
%iFit           = numel(time);
TBeg           = time(1)
TMid           = time(iFit)
TEnd           = timeSorted(end);
ThiefPos       = track_data(end,10:12);
ThiefVel       = track_data(end,13:15);

% 1st Attempt
% TFit           = TMid
% %TFit           = TEnd
% TExtrap        = TEnd;
% TPolyFit       = TMid;
% %TPolyFit      = (TMid+TEnd)/2.0;
% %TPolyFit      = (TMid+2.0*TEnd)/3.0;

% % 2nd Attempt
% %TFit           = TMid
% TFit           = TEnd
% TExtrap        = TEnd;
% %TPolyFit       = TMid;
% TPolyFit      = (TMid+TEnd)/2.0;
% %TPolyFit      = (TMid+2.0*TEnd)/3.0;

% 3rd Attempt
%TFit           = TMid
TFit           = TEnd
TExtrap        = (TMid+TEnd)/2.0;
%TPolyFit       = TMid;
TPolyFit       = TMid
%TPolyFit      = (TMid+2.0*TEnd)/3.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check Out Polynomial Fits   Sensor Orbits
state_vector_sat = cell(10,1);
kepler_fit_sat   = cell(6,1);
state_vector_los = cell(10,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%plot3(track_data(:,3), track_data(:,4),track_data(:,4))
Satellites = unique(track_data(:,1));
ChiSqSat      = [];
epsilon = 0.001*ones(6,1);
%Observed = [zeros(3,1)];
losM     = ones(3,1);
losM     = losM/norm(losM);
InvCovUV    = eye(2)*10^(8);
%InvCovUV    = eye(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iKepler = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iSat=1:numel(Satellites)
    Index      = find(track_data(:,1) == Satellites(iSat));
    TrackOrbit = track_data(Index,:);
    plot3(TrackOrbit(:,3),TrackOrbit(:,4),TrackOrbit(:,5),'linewidth',2);
    hold on
    rows       = size(TrackOrbit,1);
    Iorder     = 5;
    %Iorder     = 1;
    Object     = 1; % for Sensor
    iType      = 1;  % iType is hard-wired at the moment
    %TPolyFit       = median(TrackOrbit(:,2))
    %TPolyFit       = 0.;
    [state_vector_fit, fitSat, Iorder, ChiSq] = PolyFit(Object, TrackOrbit, TPolyFit, Iorder, iType);
    plot3(state_vector_fit(1,:),state_vector_fit(2,:),state_vector_fit(3,:),'linewidth',2);
    hold on
    Diff = state_vector_fit(1:3,:) - TrackOrbit(:,3:5)';
    ChiSqSat = [ChiSqSat, sum(sum(Diff.*Diff))];
    state_vector_sat{iSat} = PolyExtrapToFitTime(fitSat, TPolyFit);
    time  = TPolyFit/TU;
    Rpos =  state_vector_sat{iSat}(1:3)/DU;
    Rdot =  state_vector_sat{iSat}(4:6)/VU;
    [e a Inclination omega Omega Mp Period Jacobian] = kepler(time, Rpos, Rdot, mu)
    KeplerObject = KeplerFromECI(time, Rpos, Rdot, units);
    InMotion = KeplerObject.Extrapolator(time);
    KeplerObject.JacECI_2_Kepler*KeplerObject.JacKepler_2_ECI
    
    % InMotion = KeplerObject.Extrapolator(2*time);
    % KeplerObject.JacECI_2_Kepler*KeplerObject.JacKepler_2_ECI
    % InMotion = KeplerObject.Extrapolator(3*time);
    % KeplerObject.JacECI_2_Kepler*KeplerObject.JacKepler_2_ECI
    % InMotion = KeplerObject.Extrapolator(4*time);
    % KeplerObject.JacECI_2_Kepler*KeplerObject.JacKepler_2_ECI
    % InMotion = KeplerObject.Extrapolator(5*time);
    % KeplerObject.JacECI_2_Kepler*KeplerObject.JacKepler_2_ECI
    % InMotion = KeplerObject.Extrapolator(6*time);
    % KeplerObject.JacECI_2_Kepler*KeplerObject.JacKepler_2_ECI
    % InMotion = KeplerObject.Extrapolator(time);
    % KeplerObject.JacECI_2_Kepler*KeplerObject.JacKepler_2_ECI
    
    Kepler  = [e a Inclination omega Omega Mp]
    XKepler = [Kepler'; time];
    [XEquin, XECI, JacKepler_2_Equinoctial] =  Kepler_2_All(XKepler, units);
    
    KeplerSat  = LagrangePlanetary(Kepler, units);
    iKepler = [iKepler; Kepler];
    
    %Sensor   =  [Rpos+epsilon(1:3); zeros(6,1)];
    Rdata  = [TrackOrbit(1,3:5)']/DU;
    Sensor = [Rdata+1000.0*epsilon(1:3); zeros(6,1)];
    % Measurement Angles:
    %losM    = state_vector_losArray(1:3,ii);
    losM     = TrackOrbit(1,6:8)';
    losM     = losM/norm(losM);   % los
    
    [Rinit,Vinit,ParamList,JacKepler_2_ECI,Hessian,GravityCan] = KeplerSat.OrbitDerivatives(time, Sensor, losM, InvCovUV);
    [Rinit,Vinit,ParamList] = OrbitAtTime(KeplerSat, time, Sensor, losM, InvCovUV);
    % Check that the Inverse Transformation should produce the inverse Jacobian matrix.
    %The following matrix produce should be the 6x6 Identity matrix:
    KeplerObject.JacECI_2_Kepler*JacKepler_2_ECI(1:6,2:7)
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Begin Finite Differences Check
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [Rextrap,Vextrap,ParamList] = OrbitAtTime(KeplerSat, tEpoch, Sensor, losM, InvCovUV);
    % HessDifference = 0.0;
    % for LineCheck = 1:175
    %     [Jacob, JacobFinite, Hess, HessFinite] = DigitalJustice(KeplerSat, LineCheck, Sensor, losM);
    %     %Jacob
    %     %JacobFinite
    %     %Hess
    %     %HessFinite
    %     Sum = sum(sum(abs(Hess - HessFinite)));
    %     HessDifference = HessDifference + Sum;
    %     disp([' Line = ', num2str(LineCheck),'  sum = ', num2str(Sum)])
    % end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % disp([' Total Differences Hessian = ', num2str(HessDifference)])
    %
    % HessDifference    = 0.0;
    % ForwardDifference = 0.0;
    % for LineCheck = 1:115
    %     [Jacob, JacobFinite, Hess, HessFinite, ForDerivs, DerivsDataList] = DigitalJustice(KeplerObject, LineCheck);
    %     SumForward = sum(sum(ForDerivs - DerivsDataList));
    %     ForwardDifference = ForwardDifference + abs(SumForward);
    %     %Jacob
    %     %JacobFinite
    %     SumJac = sum(Jacob - JacobFinite);
    %
    %     %Hess
    %     %HessFinite
    %     Sum = sum(sum(abs(Hess - HessFinite)));
    %     HessDifference = HessDifference + Sum;
    %     disp([' Line = ', num2str(LineCheck),'  Difference Jacobian = ',num2str(SumJac),'  Difference Hessian = ', num2str(Sum)])
    % end
    % disp([' Total Differences Forward / Hessian = ',num2str(ForwardDifference),'  ',num2str(HessDifference)])
    
    % Kepler  = [Rpos(1), Rpos(2), Rpos(3), Rdot(1), Rdot(2), Rdot(3)];
    % Kdelta  = KeplerOrbit(time, Rpos, Rdot);
    % epsilon = 0.000000001;
    %
    % for I = 1:90
    %
    %         H0 = Kdelta.DataList(I);
    %         Jacob = Kdelta.JacobianAll(I,:);
    %
    %         % Do time Increment Separately
    %         KeplerDelta = epsilon*eye(6);
    %         Kdelta = KeplerOrbit(time+epsilon, Kepler(1:3), Kepler(4:6));
    %         HT     = Kdelta.DataList(I);
    %         JacobT = Kdelta.JacobianAll(I,:);
    %
    %         Kepler1     = Kepler + KeplerDelta(1,:);
    %         Kdelta1 = KeplerOrbit(time, Kepler1(1:3), Kepler1(4:6));
    %         H1     = Kdelta1.DataList(I);
    %         Jacob1 = Kdelta1.JacobianAll(I,:);
    %
    %         Kepler2     = Kepler + KeplerDelta(2,:);
    %         Kdelta2 = KeplerOrbit(time, Kepler2(1:3), Kepler2(4:6));
    %         H2     = Kdelta2.DataList(I);
    %         Jacob2 = Kdelta2.JacobianAll(I,:);
    %
    %         Kepler3     = Kepler + KeplerDelta(3,:);
    %         Kdelta3 = KeplerOrbit(time, Kepler3(1:3), Kepler3(4:6));
    %         H3     = Kdelta3.DataList(I);
    %         Jacob3 = Kdelta3.JacobianAll(I,:);
    %
    %         Kepler4     = Kepler + KeplerDelta(4,:);
    %         Kdelta4 = KeplerOrbit(time, Kepler4(1:3), Kepler4(4:6));
    %         H4     = Kdelta4.DataList(I);
    %         Jacob4 = Kdelta4.JacobianAll(I,:);
    %
    %         Kepler5     = Kepler + KeplerDelta(5,:);
    %         Kdelta5 = KeplerOrbit(time, Kepler5(1:3), Kepler5(4:6));
    %         H5     = Kdelta5.DataList(I);
    %         Jacob5 = Kdelta5.JacobianAll(I,:);
    %
    %         Kepler6     = Kepler + KeplerDelta(6,:);
    %         Kdelta6 = KeplerOrbit(time, Kepler6(1:3), Kepler6(4:6));
    %         H6     = Kdelta6.DataList(I);
    %         Jacob6 = Kdelta6.JacobianAll(I,:);
    %
    %         %HT-H0
    %         I
    %         Jacob
    %         JacobFinite = [H1-H0, H2-H0, H3-H0, H4-H0, H5-H0, H6-H0]/epsilon
    %
    % end
    %
    
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
        %rSatKepler = [];
        for ii =1:rows
            time   = TrackOrbit(ii,2)/TU;
            Rdata  = [TrackOrbit(ii,3:5)']/DU;
            Sensor = [Rdata+epsilon(1:3); zeros(6,1)];
            % Measurement Angles:
            %losM    = state_vector_losArray(1:3,ii);
            losM     = TrackOrbit(ii,6:8)';
            losM     = losM/norm(losM);   % los
            rangeM   = norm(losM);
            thetaM   = acos(losM(3)/rangeM);
            phiM     = atan2(losM(2), losM(1));
            %Observed =[thetaM; phiM];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % [Rinit,Vinit,ParamList,Jacobian,Hessian,GravityCan] = KeplerSat.OrbitDerivatives(time, Sensor, losM);
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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            % deltaR    = Rinit - Rdata + Jacobian(1:3,2:7)'*delZ;
            % norm(deltaR)
            [H_I, Jacob, Hess] = WorkOrder(KeplerSat, 105);
            betaNR     = betaNR + Jacob(2:7)';
            hessNR     = hessNR + Hess(2:7,2:7);
            % H_I
            % rSatKepler = [rSatKepler, Rinit];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %rSatKepler = [rSatKepler, Rinit];
            ChiSq      = ChiSq  + 0.5*deltaR'*InvCovR*deltaR;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        %ChiSq
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if delp(1) > 0
            scale = (1 - Kepler(1))/delp(1);
        else
            scale = (Kepler(1) - 1)/delp(1);
        end
        if scale > 1
            scale = 1;
        end
        ChiSqMin = 10^100;
        dScale = scale*0.1;
        KeplerBest = Kepler;
        fractionBest = 1.0;
        KeplerTest = Kepler;
        for isearch = 0:10
            fraction = isearch*dScale;
            KeplerTest = Kepler + fraction*delp;
            KeplerSat1  = LagrangePlanetary(KeplerTest, units);
            ChiSqTest      = 0.0;
            for ii =1:rows
                time   = TrackOrbit(ii,2)/TU;
                Rdata  = [TrackOrbit(ii,3:5)']/DU;
                Sensor = [Rdata+epsilon(1:3); zeros(6,1)];
                % Measurement Angles:
                %losM    = state_vector_losArray(1:3,ii);
                losM     = TrackOrbit(ii,6:8)';
                losM     = losM/norm(losM);   % los
                rangeM   = norm(losM);
                thetaM   = acos(losM(3)/rangeM);
                phiM     = atan2(losM(2), losM(1));
                %Observed =[thetaM; phiM];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [Rinit,Vinit,ParamList] = OrbitAtTime(KeplerSat1, time, Sensor, losM, InvCovUV);
                deltaR   = Rinit - Rdata;
                ChiSqTest      = ChiSqTest  + 0.5*deltaR'*InvCovR*deltaR;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            if ChiSqTest < ChiSqMin
                ChiSqMin = ChiSqTest;
                fractionBest = fraction;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %ChiSq
        TangentMag   = log(norm(betaNR))
        TangentArray = [TangentArray, TangentMag];
        %Kepler      =  Kepler + deltaPar'
        %if NR == 1
        fractionBest = 1;
        %end
        Kepler      =  Kepler + fractionBest*delp;
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
    %TPolyFit       = median(TrackOrbit(:,2))
    %TPolyFit       = 0.;
    [state_vector_fit, fitLos, Iorder, ChiSq] = PolyFit(Object, TrackOrbit, TPolyFit,Iorder,iType);
    plot3(state_vector_fit(1,:),state_vector_fit(2,:),state_vector_fit(3,:),'linewidth',2);
    hold on
    Diff = state_vector_fit(1:3,:) - TrackOrbit(:,6:8)';
    ChiSqSat = [ChiSqSat, sum(sum(Diff.*Diff))];
    state_vector_los{iSat} = PolyExtrapToFitTime(fitLos, TPolyFit);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Examing the roots coming out of Laplace's Orbit Determination Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Search for the Statevector that predicts the LOS vector making the smallest
% angle with the direction of the Unit Vector LOS (measured).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ChiSqMin     = 10^100;
KeplerAve    = [];
KeplerChiSq  = [];
SVBest       = [];
SVBuild      = [];
SVestimated1 = [];
SVestimated2 = [];
StateVector  = [];
track_Sat    = track_data;
rows         = size(track_Sat,1);
Iorder       = 5;

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
    
    [StateVectorCell,RootsLaplace] = LaplaceInECI(state_vector_sat{iSat}, state_vector_los{iSat});
    for icell = 1:size(StateVectorCell,2)
        StateVector = StateVectorCell{icell};
        time        = StateVector(10)/TU;
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
            %Index = find(track_data(:,1) == Satellites(iSat));
            %track_Sat = track_data(Index,:);
            for icell = 1:size(StateVectorCell,2)
                StateVector = StateVectorCell{icell};
                time        = StateVector(10)/TU;
                Rpos        = StateVector(1:3)/DU;
                Rdot        = StateVector(4:6)/VU;
                ChiSq       = 0.0;
                % comment out 506 to see what happens. !!!!!!!!!!!!!!!!!!!!!!!!!!!!
                [e a Inclination omega Omega Mp Period Jacobian] = kepler(time, Rpos, Rdot, mu);
                KeplerTest = [e a Inclination omega Omega Mp]; % Classical Orbital Elements
                
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
                        %Observed = [thetaM; phiM];
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % %Sat Position at measured position
                        Rsensor  = [track_Sat(ii,3:5)']/DU;
                        Sensor   = [Rsensor; zeros(6,1)];
                        
                        % %Sat Position from Polynomial at TFit
                        % Sensor = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU]
                        
                        % %Sat Position from Orbit Fit to Kepler elements:
                        %[Rsensor, Vsensor, ParamList, Jacobian, Hessian, GravityCan] = KeplerSat.OrbitDerivatives(tRecord, Sensor, losM);
                        %[Rsensor, Vsensor, ParamList] = OrbitAtTime(KeplerSat, tRecord, Sensor, losM, InvCovUV);
                        %Sensor   = [Rsensor; Vsensor; zeros(3,1)];
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
                        KeplerBest       = KeplerTest
                        SVBest           = StateVector
                        [SVestimated, CovSV1, fit, CovFit, ChiSqFit] = PolyKeplerFit(SVBest, TExtrap, TPolyFit, track_data, Iorder, iType, units, InvCovUV);
                        InitialState     = InitialStateFit(Iorder, TPolyFit, fit, CovFit, SVestimated);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  JRS To Here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Roots       = [];
StateVectorSats = cell(10);
RTarget  = cell(3);
for iSat=1:numel(Satellites)
    [StateVectorCell,RootsLaplace] = LaplaceInECI(state_vector_sat{iSat}, state_vector_los{iSat});
    Roots = [Roots, RootsLaplace];
    for icell = 1:size(StateVectorCell,2)
        StateVector = StateVectorCell{icell};
        %time        = TPolyFit/TU;
        %Rpos        = StateVector(1:3)/DU;
        %Rdot        = StateVector(4:6)/VU;
        StateVectorSats{iSat} = [StateVectorSats{iSat}, StateVector];
        RTarget{iSat}        = [RTarget{iSat}, StateVector(1:3)];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the Average StateVector which agrees with each Sensor the closes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ichecks = 0;
if numel(Satellites) == 2
    ResSqMin = 10^(50);
    for iSat1 = 1:size(StateVectorSats{1}, 2)
        SV1 = StateVectorSats{1}(:,iSat1);
        for iSat2 = 1:size(StateVectorSats{2}, 2)
            SV2 = StateVectorSats{2}(:,iSat2);
            Ichecks= Ichecks + 1;
            SVAVE = (SV1 + SV2)/2.0;
            DiffSq1 = sum((SV1-SVAVE).*(SV1-SVAVE));
            DiffSq2 = sum((SV2-SVAVE).*(SV2-SVAVE));
            ResidSQ = DiffSq1 + DiffSq2;
            if ResidSQ < ResSqMin;
                ResSqMin = ResidSQ;
                SVBuild = SVAVE;
            end
        end
    end
end
if numel(Satellites) == 3
    ResSqMin = 10^(50);
    for iSat1 = 1:size(StateVectorSats{1}, 2)
        SV1 = StateVectorSats{1}(:,iSat1);
        for iSat2 = 1:size(StateVectorSats{2}, 2)
            SV2 = StateVectorSats{2}(:,iSat2);
            for iSat3 = 1:size(StateVectorSats{3}, 2)
                SV3 = StateVectorSats{3}(:,iSat3);
                Ichecks= Ichecks + 1;
                SVAVE = (SV1 + SV2 + SV3)/3.0;
                DiffSq1 = sum((SV1-SVAVE).*(SV1-SVAVE));
                DiffSq2 = sum((SV2-SVAVE).*(SV2-SVAVE));
                DiffSq3 = sum((SV3-SVAVE).*(SV3-SVAVE));
                ResidSQ = DiffSq1 + DiffSq2 + DiffSq3;
                if ResidSQ < ResSqMin
                    ResSqMin = ResidSQ;
                    SVBuild = SVAVE;
                end
            end
        end
    end
end
if numel(Satellites) == 4
    ResSqMin = 10^(50);
    for iSat1 = 1:size(StateVectorSats{1}, 2)
        SV1 = StateVectorSats{1}(:,iSat1);
        for iSat2 = 1:size(StateVectorSats{2}, 2)
            SV2 = StateVectorSats{2}(:,iSat2);
            for iSat3 = 1:size(StateVectorSats{3}, 2)
                SV3 = StateVectorSats{3}(:,iSat3);
                for iSat4 = 1:size(StateVectorSats{4}, 2)
                    Ichecks= Ichecks + 1;
                    SV4 = StateVectorSats{4}(:,iSat4);
                    SVAVE = (SV1 + SV2 + SV3 + SV4)/4.0;
                    DiffSq1 = sum((SV1-SVAVE).*(SV1-SVAVE));
                    DiffSq2 = sum((SV2-SVAVE).*(SV2-SVAVE));
                    DiffSq3 = sum((SV3-SVAVE).*(SV3-SVAVE));
                    DiffSq4 = sum((SV4-SVAVE).*(SV4-SVAVE));
                    ResidSQ = DiffSq1 + DiffSq2 + DiffSq3 + DiffSq4;
                    if ResidSQ < ResSqMin
                        ResSqMin = ResidSQ;
                        SVBuild = SVAVE;
                    end
                end
            end
        end
    end
end

if numel(SVBuild) == 10
    time        = SVBuild(10)/TU;
    Rpos        = SVBuild(1:3)/DU;
    Rdot        = SVBuild(4:6)/VU;
    ChiSq       = 0.0;
    % comment out 506 to see what happens. !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    [e a Inclination omega Omega Mp Period Jacobian] = kepler(time, Rpos, Rdot, mu);
    KeplerTest = [e a Inclination omega Omega Mp]; % Classical Orbital Elements
    
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
            %Observed = [thetaM; phiM];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %Sat Position at measured position
            Rsensor  = [track_Sat(ii,3:5)']/DU;
            Sensor   = [Rsensor; zeros(6,1)];
            
            % %Sat Position from Polynomial at TFit
            % Sensor = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU]
            
            % %Sat Position from Orbit Fit to Kepler elements:
            %[Rsensor, Vsensor, ParamList, Jacobian, Hessian, GravityCan] = KeplerSat.OrbitDerivatives(tRecord, Sensor, losM);
            %[Rsensor, Vsensor, ParamList] = OrbitAtTime(KeplerSat, tRecord, Sensor, losM, InvCovUV);
            %Sensor   = [Rsensor; Vsensor; zeros(3,1)];
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
            KeplerBest       = KeplerTest
            SVBest           = SVBuild
            [SVestimated, CovSV1, fit, CovFit, ChiSqFit] = PolyKeplerFit(SVBest, TExtrap, TPolyFit, track_data, Iorder, iType, units, InvCovUV);
            InitialState     = InitialStateFit(Iorder, TPolyFit, fit, CovFit, SVestimated);
        end
    end
end

%SVBest = [SVBest; TPolyFit];
%Ichecks
%StateVectorECI = SVBest;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Iorder = 5;
iType  = 1;
Object = 1;

if numel(SVBest) == 10
[SVestimated1, CovSV1, fit, CovFit, ChiSqFit] = PolyKeplerFit(SVBest, TExtrap, TPolyFit, track_data, Iorder, iType, units, InvCovUV);
ChiSq
ThiefPos' - SVestimated1(1:3)
ThiefVel' - SVestimated1(4:6)
time      = SVestimated1(10)/TU;
Rpos      = SVestimated1(1:3)/DU;
Rdot      = SVestimated1(4:6)/VU;
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
            %Observed = [thetaM; phiM];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %Sat Position at measured position
            Rsensor  = [track_Sat(ii,3:5)']/DU;
            Sensor   = [Rsensor; zeros(6,1)];
            
            % %Sat Position from Polynomial at TFit
            % Sensor = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU]
            
            % %Sat Position from Orbit Fit to Kepler elements:
            %[Rsensor, Vsensor, ParamList, Jacobian, Hessian, GravityCan] = KeplerSat.OrbitDerivatives(tRecord, Sensor, losM);
            %[Rsensor, Vsensor, ParamList] = OrbitAtTime(KeplerSat, tRecord, Sensor, losM, InvCovUV);
            %Sensor   = [Rsensor; Vsensor; zeros(3,1)];
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
            KeplerBest       = KeplerTest
            SVBest           = SVestimated1
            InitialState     = InitialStateFit(Iorder, TPolyFit, fit, CovFit, SVBest);
        end
    end

end

Iorder = 5;
iType  = 1;
Object = 1;

[Iorder, fit, CovFit, ChiSqFit] = PolyCorrFit(Object, track_data, TPolyFit, Iorder, iType, units, InvCovUV);
[state_vector, CovSV] = PolyCorrExtrapToFitTime(Iorder, TExtrap, TPolyFit, fit, CovFit, units);
if numel(state_vector) == 10
ChiSq
ThiefPos' - state_vector(1:3)
ThiefVel' - state_vector(4:6)

time      = state_vector(10)/TU;
Rpos      = state_vector(1:3)/DU;
Rdot      = state_vector(4:6)/VU;
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
            %Observed = [thetaM; phiM];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %Sat Position at measured position
            Rsensor  = [track_Sat(ii,3:5)']/DU;
            Sensor   = [Rsensor; zeros(6,1)];
            
            % %Sat Position from Polynomial at TFit
            % Sensor = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU]
            
            % %Sat Position from Orbit Fit to Kepler elements:
            %[Rsensor, Vsensor, ParamList, Jacobian, Hessian, GravityCan] = KeplerSat.OrbitDerivatives(tRecord, Sensor, losM);
            %[Rsensor, Vsensor, ParamList] = OrbitAtTime(KeplerSat, tRecord, Sensor, losM, InvCovUV);
            %Sensor   = [Rsensor; Vsensor; zeros(3,1)];
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
            StateVectorAve   = state_vector;
        end
        if ChiSq < ChiSqMin
            ChiSqMin         = ChiSq;
            KeplerBest       = KeplerTest
            SVBest           = state_vector
            InitialState     = InitialStateFit(Iorder, TPolyFit, fit, CovFit, SVBest);
        end
    end
end

%InitialState = InitialStateFit(Iorder, TPolyFit, TEnd, TExtrap, fit, CovFit);
Iorder = 5;
iType  = 1;
Object = 1;

if numel(SVBest) == 10
[SVestimated2, CovSV1, fit, CovFit, ChiSqFit] = PolyKeplerFit(SVBest, TExtrap, TPolyFit, track_data, Iorder, iType, units, InvCovUV);
ChiSqFit
ThiefPos' - SVestimated2(1:3)
ThiefVel' - SVestimated2(4:6)

time      = SVestimated2(10)/TU;
Rpos      = SVestimated2(1:3)/DU;
Rdot      = SVestimated2(4:6)/VU;
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
            %Observed = [thetaM; phiM];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %Sat Position at measured position
            Rsensor  = [track_Sat(ii,3:5)']/DU;
            Sensor   = [Rsensor; zeros(6,1)];
            
            % %Sat Position from Polynomial at TFit
            % Sensor = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU]
            
            % %Sat Position from Orbit Fit to Kepler elements:
            %[Rsensor, Vsensor, ParamList, Jacobian, Hessian, GravityCan] = KeplerSat.OrbitDerivatives(tRecord, Sensor, losM);
            %[Rsensor, Vsensor, ParamList] = OrbitAtTime(KeplerSat, tRecord, Sensor, losM, InvCovUV);
            %Sensor   = [Rsensor; Vsensor; zeros(3,1)];
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
            StateVectorAve   = state_vector;
        end
        if ChiSq < ChiSqMin
            ChiSqMin         = ChiSq;
            KeplerBest       = KeplerTest
            SVBest           = SVestimated2
            InitialState     = InitialStateFit(Iorder, TPolyFit, fit, CovFit, SVBest);
        end
    end
end

% function  InitialState = InitialStateFit(Iorder, TPolyFit, TEnd, TExtrap, fit, CovFit)
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     InitialState.Iorder    = Iorder;
%     InitialState.TPolyFit  = TPolyFit;
%     InitialState.fit       = fit;
%     InitialState.CovFit    = CovFit;
%
%     InitialState.TExtrap   = TExtrap;
%     InitialState.TEnd      = TEnd;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end

% Esimated Accuracy of the Inital Guess
% SizeInitialState = [];
% rows = size(track_data,1);
% for i = 1:rows
%     Time = track_data(i,2);
%     [state_vector, CovSV] = PolyCorrExtrapToFitTime(Iorder, Time, TFit, fit, CovFit, units);
%     sizeCov = sqrt(trace(CovSV(1:3,1:3)));
%     SizeInitialState = [SizeInitialState, sizeCov];
% end
% figure
% plot(track_data(:,2), SizeInitialState);
% title('Square_root of the trace of the position-space Initial State Covariance Matrix');
% ylim([0 300])
% ylabel('Meters');
% xlabel('Time (seconds)');
Nbins          = 1400;
tBins          = 1400/Nbins;
Time           = 0.0 - tBins/2.0;
SizeInitialStatePos = [];
SizeInitialStateVel = [];
TimeArray = [];
for I = 1:Nbins
    Time       = Time + tBins;
    TimeArray = [TimeArray, Time];
    [state_vector, CovSV] = PolyCorrExtrapToFitTime(Iorder, Time, TPolyFit, fit, CovFit, units);
    sizeCov = sqrt(trace(CovSV(1:3,1:3)));
    SizeInitialStatePos = [SizeInitialStatePos, sizeCov];
    sizeCov = sqrt(trace(CovSV(4:6,4:6)));
    SizeInitialStateVel = [SizeInitialStateVel, sizeCov];
end
figure
plot(TimeArray, SizeInitialStatePos);
title('Square_root of the trace of the position-space Initial State Covariance Matrix');
ylim([0 1000])
ylabel('Meters');
xlim([0 1400])
xlabel('Time (seconds)');

figure
plot(TimeArray, SizeInitialStateVel);
title('Square_root of the trace of the Velocity-space Initial State Covariance Matrix');
ylim([0 1000])
ylabel('Meters per Second');
xlim([0 1400])
xlabel('Time (seconds)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ClosestRoots = [];
% RSqMin = 10^50
% Pairs = [];
% for iSat = 1:numel(Satellites)-1
%     for jSat = 1:numel(Satellites)
%         if iSat ~= jSat
%             for SViIndex = 1:size(RTarget{iSat},2)
%                 for SVjIndex = 1:size(RTarget{jSat},2)
%                     diff = RTarget{iSat}(:,SViIndex) - RTarget{jSat}(:,SVjIndex)
%                     distSq = sum(diff.*diff)
%                     if distSq < RSqMin
%                         Pair = [jSat, SVjIndex]
%                         RSqMin = distSq
%                         StateVectorBest = StateVectorSats{jSat}(:,SVjIndex)
%                     end
%                 end
%             end
%         end
%     end
% end
% iSat = Pair(1)
% SViIndex = Pair(2)
% BuildPairs = Pair;
% RSqMin = 10^50
% for jSat = 1:numel(Satellites)
%     if iSat ~= jSat
%         for SVjIndex = 1:size(RTarget{jSat},2)
%             diff = RTarget{iSat}(:,SViIndex) - RTarget{jSat}(:,SVjIndex)
%             distSq = sum(diff.*diff)
%             if distSq < RSqMin
%                 Pair = [jSat, SVjIndex]
%                 RSqMin = distSq
%                 StateVectorBest = StateVectorSats{jSat}(:,SVjIndex)
%             end
%         end
%     end
% end
%
% iSat = Pair(1)
% SViIndex = Pair(2)
% BuildPairs = [BuildPairs; Pair]
% RSqMin = 10^50
% for jSat = 1:numel(Satellites)
%     if numel(find(find(BuildPairs(:,1))==jSat)) == 0
%         for SVjIndex = 1:size(RTarget{jSat},2)
%             diff = RTarget{iSat}(:,SViIndex) - RTarget{jSat}(:,SVjIndex)
%             distSq = sum(diff.*diff)
%             if distSq < RSqMin
%                 Pair = [jSat, SVjIndex]
%                 RSqMin = distSq
%                 StateVectorBest = StateVectorSats{jSat}(:,SVjIndex)
%             end
%         end
%     end
% end
% BuildPairs = [BuildPairs; Pair]
%
% StateVectorECI = StateVectorBest;
% format long
% RootArray = reshape(Roots,[8*numel(Satellites),1])
% RootArray = RootArray(find(imag(RootArray) == 0))
% RootArray = RootArray/abs(max(RootArray))
% format short
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Search for the Statevector that predicts the LOS vector making the smallest
% angle with the direction of the Unit Vector LOS (measured).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for iSat=1:numel(Satellites)
%     ChiSqSatMin = 10^100;
%     % time        = TFit/TU;
%     % Rpos        = state_vector_sat{iSat}(1:3)/DU;
%     % Rdot        = state_vector_sat{iSat}(4:6)/VU;
%     % % comment out 506 to see what happens. !!!!!!!!!!!!!!!!!!!!!!!!!!!!
%     % [e a Inclination omega Omega Mp Period Jacobian] = kepler(time, Rpos, Rdot, mu);
%     % KeplerSensor = [e a Inclination omega Omega Mp]; % Canonical Units
%     % KeplerSat  = LagrangePlanetary(KeplerSensor, units);
%     KeplerSat  = kepler_fit_sat{iSat};
%
%     [StateVectorCell,RootsLaplace] = LaplaceInECI(state_vector_sat{iSat}, state_vector_los{iSat});
%     for icell = 1:size(StateVectorCell,2)
%         StateVector = StateVectorCell{icell};
%         time        = TFit/TU;
%         Rpos        = StateVector(1:3)/DU;
%         Rdot        = StateVector(4:6)/VU;
%         ChiSq       = 0.0;
%         % comment out 506 to see what happens. !!!!!!!!!!!!!!!!!!!!!!!!!!!!
%         [e a Inclination omega Omega Mp Period Jacobian] = kepler(time, Rpos, Rdot, mu);
%         KeplerTest = [e a Inclination omega Omega Mp]; % Canonical Units
%         Imaginary  = numel(find(imag(KeplerTest)));
%
%         if Imaginary == 0
%             KeplerObj  = LagrangePlanetary(KeplerTest, units);
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Index = find(track_data(:,1) == iSat);
%             track_Sat = track_data(Index,:);
%             rows       = size(track_Sat,1);
%             for icell = 1:size(StateVectorCell,2)
%                 StateVector = StateVectorCell{icell};
%                 time        = TFit/TU;
%                 Rpos        = StateVector(1:3)/DU;
%                 Rdot        = StateVector(4:6)/VU;
%                 ChiSq       = 0.0;
%                 % comment out 506 to see what happens. !!!!!!!!!!!!!!!!!!!!!!!!!!!!
%                 [e a Inclination omega Omega Mp Period Jacobian] = kepler(time, Rpos, Rdot, mu);
%                 KeplerTest = [e a Inclination omega Omega Mp]; % Classical Orbital Elements
%
%                 Imaginary  = numel(find(imag(KeplerTest)));
%
%                 if Imaginary == 0
%                     KeplerObj  = LagrangePlanetary(KeplerTest, units);
%                     for ii =1:rows
%                         tRecord = track_Sat(ii,2)/TU;
%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %losM    = state_vector_losArray(1:3,ii);
%                         losM     = track_Sat(ii,6:8)';
%                         losM     = losM/norm(losM);   % los
%                         rangeM   = norm(losM);
%                         thetaM   = acos(losM(3)/rangeM);
%                         phiM     = atan2(losM(2), losM(1));
%                         Observed = [thetaM; phiM];
%                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         % %Sat Position at measured position
%                         Rsensor  = [track_Sat(ii,3:5)']/DU;
%                         Sensor   = [Rsensor; zeros(6,1)];
%
%                         % %Sat Position from Polynomial at TFit
%                         % Sensor = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU]
%
%                         % %Sat Position from Orbit Fit to Kepler elements:
%                         %[Rsensor, Vsensor, ParamList, Jacobian, Hessian, GravityCan] = KeplerSat.OrbitDerivatives(tRecord, Sensor, losM);
%                         [Rsensor, Vsensor, ParamList] = OrbitAtTime(KeplerSat, tRecord, Sensor, losM, InvCovUV);
%                         Sensor   = [Rsensor; Vsensor; zeros(3,1)];
%                         %[Rinit,Vinit,ParamList,Jacobian,Hessian,GravityCan] = KeplerObj.OrbitDerivatives(tRecord, Sensor, losM);
%                         [Rinit,   Vinit,   ParamList] = OrbitAtTime(KeplerObj, tRecord, Sensor, losM, InvCovUV);
%                         los       = [ParamList(7); ParamList(8); ParamList(9)];
%                         los       = los/norm(los);
%                         AngleDiff = abs(1 - dot(los,losM));
%                         ChiSq     = ChiSq + AngleDiff;
%                     end
%                     if ChiSq < ChiSqSatMin
%                         ChiSqSatMin      = ChiSq;
%                         KeplerSatBest    = KeplerTest;
%                         StateVectorAve   = StateVector;
%                     end
%                     if ChiSq < ChiSqMin
%                         ChiSqMin         = ChiSq;
%                         KeplerBest       = KeplerTest
%                         StateVectorECI   = StateVector
%                     end
%                 end
%             end
%         end
%     end
%     if ChiSqSatMin < 10^(20)
%         KeplerChiSq = [KeplerChiSq; ChiSqSatMin];
%         KeplerAve   = [KeplerAve;   KeplerSatBest];
%     end
%
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Begin Finite Differences Check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %losM    = state_vector_losArray(1:3,ii);
% losM     = track_data(iFit,6:8)';
% losM     = losM/norm(losM);   % los
% rangeM   = norm(losM);
% thetaM   = acos(losM(3)/rangeM);
% phiM     = atan2(losM(2), losM(1));
% Observed = [thetaM; phiM];
% % Sat Position at measured position
% Rsensor  = [track_data(iFit,3:5)']/DU;
% Sensor   = [Rsensor; zeros(6,1)];
% KeplerObj  = LagrangePlanetary(KeplerBest, units);
% [Rextrap,Vextrap,ParamList] = OrbitAtTime(KeplerObj, tEpoch, Sensor, losM);
% for LineCheck = 130:133
%     [Jacob, JacobFinite, Hess, HessFinite] = DigitalJustice(KeplerObj, LineCheck, Sensor, losM);
%     Jacob
%     JacobFinite
%     Hess
%     HessFinite
%     Sum = sum(sum(Hess - HessFinite));
%     disp([' Line = ', num2str(LineCheck),'  sum = ', num2str(Sum)])
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %KeplerAve = sum(KeplerAve)/size(KeplerAve,1)
% KeplerChiSq
% KeplerAve
% NumChiSq = size(KeplerChiSq,1);
% KeplerCombined = zeros(1, size(KeplerAve,2))
% InvCov         = 0.0;
% for i = 1:NumChiSq
%     CovInv = 1.0/KeplerChiSq(i);
%     InvCov  = InvCov + CovInv;
%     KeplerCombined = KeplerCombined + CovInv*KeplerAve(i,:);
% end
% KeplerCombined = (1.0/InvCov)*KeplerCombined
% % KeplerAveECI   = LagrangePlanetary(KeplerCombined, units);
% % [rInitial, vInitial] = extrapolate(KeplerAveECI, TFit/TU);
% % StateVectorECI     = [rInitial'*DU; vInitial'*VU; zeros(3,1)]
% disp([' KeplerBest  = ', num2str(KeplerBest)])
% %disp([' KeplerThief = ',num2str(KeplerThief)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try to Integrate the Equations of Motion for the Classical Elements
tBeg = TBeg/TU
tFit = TFit/TU
tEnd = TEnd/TU
Y         = zeros(6,1);
KLagrange = LagrangePlanetary(KeplerBest, units);
[Rinit,   Vinit,   ParamList] = OrbitAtTime(KLagrange, tFit, Sensor, losM, InvCovUV);
options   = odeset('RelTol', 1e-10, 'AbsTol', 1e-13);
dJdt      = perturbations_odefun(KLagrange,tFit,Y,units);
Y         = zeros(6,1);
% Propagate from fit Reference time tFit back to tBeg
[Tout, Yout] = ode45(@(t,y) perturbations_odefun(KLagrange, t, y, units), [tFit, tBeg], Y, options);
%Yout(end,:)
Perturbed = [Tout, Yout];
Perturbed = sortrows(Perturbed,1);
Y = Perturbed(end,2:7);
% Propagate from fit Reference time tFit back to tEnd
if tFit < tEnd
    [Tout, Yout] = ode45(@(t,y) perturbations_odefun(KLagrange, t, y, units), [tFit, tEnd], Y, options);
    %Yout(end,:)
    Perturbed(end,:) = [];
    Perturbed = [Perturbed; Tout, Yout];
    Perturbed = sortrows(Perturbed,1);
end
%KeplerInterpolated = interp1(Perturbed(:,1),Perturbed(:,2:7),timeArray);
figure
plot(Perturbed(:,1), Perturbed(:,2:7))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All these State Vectors have been Extrapolated to the latest data point
%StateVectorECI = state_vector
%StateVectorECI
%StateVectorECI = SVestimated1
SVBest
ChiSqMin
%[KeplerFit1, CovFit1, TFit, SV7Fit1] = LeastSquares(track_data, tFirst, 1, 0, zeros(7,1) )
[KeplerFit1, CovFit1, TFit, SV7Fit1] = LeastSquaresAllSatellites(track_data, TFit, 1, 0, InvCovUV, units, InitialState)