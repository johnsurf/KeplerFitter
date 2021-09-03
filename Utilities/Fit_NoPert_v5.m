clc, clear all
restoredefaultpath
file_path  = mfilename('fullpath')
in_dir  = fileparts(file_path)
out_dir = strcat(in_dir,'\Converted')
orbitparams = wgs84Constants;
TU = orbitparams.TU;
DU = orbitparams.DU;
VU = orbitparams.VU;
AU = orbitparams.AU;
%mu = orbitparams.mu;
mu = orbitparams.Grav;
twopi  = 2.0*pi;         % 2*pi
Set = 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Input_filepath = 'Input Directory'
filename       = 'Input File.txt'
out_dir = fullfile(in_dir,'Converted')
fid = fopen(fullfile(Input_filepath, filename),'r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[message] = null_los_message_v1();
trackID   = 1

if ~isfolder(out_dir)
    mkdir(out_dir);
end
%%
count     = 0;
SVInitial = zeros(9,1);
iFirst    = 0; 
while ~feof(fid)
    tline = fgetl(fid);
    line = textscan(tline,'%s','Delimiter',' ','MultipleDelimsAsOne',1);
    newline = line{1,1}';
    % If header row, start new block read
    if contains(tline,'Track ID') | contains(tline,'Sequence')
        % if isempty(message(1).trackID)
        %     count = 1;
        % end
        % % % Header Row Data Fields
        % message(count).trackID   = str2double(newline{1,4});
        % %message(count).blockID   = str2double(newline{1,6});
        % %message(count).availTime = str2double(newline{1,9});
        if count ~= 0
            disp(' analyze block of LOS data ')
            % Convert messages into Our format
            rows       = size(message,2);
            %track_data = zeros(rows,121);
            track_data = zeros(rows,8);
            for ii = 1:rows
                time                 = message(ii).timeUTC;
                if Set == 1
                    state                = message(ii).satVec
                elseif Set == 2
                    % Adding in conversion from km to m if necessary!!!
                    state                = message(ii).satVec';
                    state                = 1e3*state;
                    los                  = message(ii).losVec';
                    % no point scaling the Cov9x9 matrix since it is "truth" and we're
                    % using unit matrices for the covariances.
                end
                % Convert ECEF to ECI at some arbitrary time:
                [W T] = Earth_Rotation_Matrix(time);
                T     = T';
                %new_state         = state;
                %new_los           = los;
                new_state           = T*state;
                new_los             = T*los;
                new_los             = new_los/norm(new_los);
                %C.data(ii,2:4)       = new_9state(1:3);
                %C.data(ii,5:7)       = new_9state(4:6);
                %C.data(ii,8:10)      = new_9state(7:9);
                %Cov9x9_ECI           = ECR_to_ECI_PVA_cov(Cov9x9,time);
                
                %track_data(ii,1)     = trackID;
                track_data(ii,1)     = message(ii).satNum;
                track_data(ii,2)     = time;
                % Adding in conversion from km to m
                track_data(ii,3:5)   = new_state(1:3);
                track_data(ii,6:8)   = new_los(1:3);
                %track_data(ii,9:11)  = new_9state(7:9);
                %track_data(ii,12:92) = reshape(Cov9x9_ECI,[1 81]);
                       
            end
            timeSorted = sortrows(track_data,2)
            if iFirst == 0
                tFirst   = timeSorted(1,2)
                iFirst = -1;
            end
            track_data(:,2) = track_data(:,2) - tFirst
            
            %[Zfit, tFit] = LeastSquares(track_data, tFirst, 402, 0 )
            
            disp(' reset for next block ')
            count = 0
        end
        
    continue
    else
        count = count + 1;
        % Message Block, row 2
        message(count).timeUTC     = str2double(newline{1,1});
        message(count).pad1        = str2double(newline{1,2});
        message(count).satNum      = str2double(newline{1,3});
        message(count).losVec(1,1) = str2double(newline{1,4});
        message(count).losVec(1,2) = str2double(newline{1,5});
        message(count).losVec(1,3) = str2double(newline{1,6});
        message(count).satVec(1,1) = str2double(newline{1,7});
        message(count).satVec(1,2) = str2double(newline{1,8});
        message(count).satVec(1,3) = str2double(newline{1,9});
        message(count).radioInt  = str2double(newline{1,10});
        message(count).band      = str2double(newline{1,11});
        message(count).pad2      = str2double(newline{1,12});
    end

end
Cov9x9 = eye(9);

%%
% Convert messages into Our format
rows       = size(message,2);
M          = unique(track_data,'rows');
track_data = M;

saveName = sprintf('%s_los_data.txt',strrep(filename,'.txt',''));
save(fullfile(out_dir,saveName),'track_data','-ascii','-double');
save(fullfile(out_dir,'track_data.txt'),'track_data','-ascii','-double');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check Out Polynomial Fits   Sensor Orbits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%plot3(track_data(:,3), track_data(:,4),track_data(:,4))
Satellites = unique(track_data(:,1))
ChiSqSat      = [];
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
    iMid       = floor(numel(track_data(:,2))/2);
    tFit       = median(TrackOrbit(:,2))
    %tFit       = 0.;
    [state_vector_fit, fit, Iorder, ChiSq] = PolyFit(Object, TrackOrbit,tFit,Iorder,iType);
    plot3(state_vector_fit(1,:),state_vector_fit(2,:),state_vector_fit(3,:),'linewidth',2);
    hold on
    Diff = state_vector_fit(1:3,:) - TrackOrbit(:,3:5)';
    ChiSqSat = [ChiSqSat, sum(sum(Diff.*Diff))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check Out Polynomial Fits   Line of Sight Polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%plot3(track_data(:,3), track_data(:,4),track_data(:,4))
Satellites = unique(track_data(:,1))
ChiSqSat      = [];
for iSat=1:numel(Satellites)
    Index      = find(track_data(:,1) == Satellites(iSat));
    TrackOrbit = track_data(Index,:);
    plot3(TrackOrbit(:,6),TrackOrbit(:,7),TrackOrbit(:,8),'linewidth',2);
    hold on
    rows       = size(TrackOrbit,1);
    Iorder     = 4;
    %Iorder     = 1;
    Object     = 2; % for los
    iType      = 1;  % iType is hard-wired at the moment
    iMid       = floor(numel(track_data(:,2))/2);
    tFit       = median(TrackOrbit(:,2))
    %tFit       = 0.;
    [state_vector_fit, fitSat, Iorder, ChiSq] = PolyFit(Object, TrackOrbit,tFit,Iorder,iType);
    plot3(state_vector_fit(1,:),state_vector_fit(2,:),state_vector_fit(3,:),'linewidth',2);
    hold on
    Diff = state_vector_fit(1:3,:) - TrackOrbit(:,6:8)';
    ChiSqSat = [ChiSqSat, sum(sum(Diff.*Diff))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[KeplerFit, CovFit, tFit] = LeastSquaresAll(track_data, tFirst, 0, 0, zeros(7,1))
[KeplerFit402, CovFit402, tFit, SV7Fit] = LeastSquares(track_data, tFirst, 402, 0, zeros(7,1) )
%[KeplerFit401, CovFit401, tFit, SV7Fit] = LeastSquares(track_data, tFirst, 401, 0, zeros(7,1) )
%[KeplerFit401, CovFit401, tFit] = LeastSquares(track_data, tFirst, 401, 1, SV7Fit )
InvComb     = inv(CovFit401 + CovFit402)
KeplerF     = CovFit402*InvComb*KeplerFit401 + CovFit401*InvComb*KeplerFit402
CovFit1      = CovFit401*InvComb*CovFit402
CovFit2      = CovFit402*InvComb*CovFit401


    function [KeplerFit, CovFit, tFit, SV7Fit] = LeastSquares(track_data, tFirst, Isat, iFirst, SV7Initial)
    
        orbitparams = wgs84Constants;
        TU = orbitparams.TU;
        DU = orbitparams.DU;
        VU = orbitparams.VU;
        AU = orbitparams.AU;
        %mu = orbitparams.mu;
        mu = orbitparams.Grav;
        twopi  = 2.0*pi;         % 2*pi
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Index = find(track_data(:,1) == Isat);
        track_data = track_data(Index,:);
        rows       = size(track_data,1);
        iMid       = floor(numel(track_data(:,2))/2);
        %tFit       = track_data(iMid,2)
        tFit       = median(track_data(:,2))
        % Intialize the Least Squares algorithm
        % Method 1:
        % First do a simple position at closest approach of two line of
        % sight vectors
        % make initial guess using first two sightings -- assume same
        % point is described by both sightings
        SVLineCrossing = ClosestPointBetweenTwoLines(track_data,tFit);
        SVLineCrossing'
        time = tFit/TU;
        Rpos = SVLineCrossing(1:3)/DU;
        Rdot = SVLineCrossing(4:6)/VU;
        [e a Inclination omega Omega Mp Period Jacobian] = kepler(time, Rpos, Rdot);
        KeplerInitLC = [e a Inclination omega Omega Mp] % Canonical Units
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Method 2 -- construct two polynomial fits to the Satellite and the Target
        % Polynomial Interpolation for the Satellite Sensor
        % Polynomial Fit Tools
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Iorder = 10;
        Object = 1; % for Sensor
        iType  = 1;  % iType is hard-wired at the moment
        [state_vector_satArray, fitSat, Iorder, ChiSq] = PolyFit(Object, track_data,tFit,Iorder,iType);
        ChiSqRoot = sqrt(sum(ChiSq));
        figure
        iData = 1:numel(ChiSq);
        plot(iData,ChiSq)
        title('Sensor Polynomial');
        xlabel('Observation Number');
        ylabel('\chi^2 from polynomial');
        disp(['root(Summed Residials) for SAT in meters = ', num2str(sqrt(ChiSqRoot))])
        %state_vector_sat               = state_vector_satArray(:,iMid);
        state_vector_sat                = PolyExtrapToFitTime(fitSat, tFit);
        state_vector_sat'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Measurement Angles:
        %losM    = state_vector_losArray(1:3,ii);
        losM     = track_data(iMid,6:8)';
        losM     = losM/norm(losM);   % los
        rangeM   = norm(losM);
        thetaM   = acos(losM(3)/rangeM);
        phiM     = atan2(losM(2), losM(1));
        Observed =[thetaM; phiM];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Code for least squares fit to Kepler Orbit to Sensor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           
        elements                        = state_vector_sat(1:3);
        derivElements                   = state_vector_sat(4:6);
        % for ii = 1:rows
        %ii   = iMid;
        %time = track_data(ii,2)/TU;
        time  = tFit/TU;
        Rpos =  elements/DU;
        Rdot =  derivElements/VU;
        [e a Inclination omega Omega Mp Period Jacobian] = kepler(time, Rpos, Rdot);
        [e a Inclination omega Omega Mp Period]
        Kepler  = [e a Inclination omega Omega Mp]
        epsilon = 0.000001*ones(6,1);         
        KeplerSat  = LagrangePlanetary(Kepler);
        Sensor   = [Rpos+epsilon(1:3); zeros(6,1)];
        Observed = [zeros(3,1)];
        [Rinit,Vinit,Hessian,Gradient,Jacobian,JacobianLos,ParamList,...
         GravityCan,HessianLos,Phi,MLagrange,RVHessian] = KeplerSat.J6Gravity(time, Sensor, Observed);
        
        deltaPar     = zeros(6,1);
        ChiSqArray   = [];
        TangentArray = [];
        InvCovR = eye(3);
        deltaPar   = zeros(6,1);
        ChiSqArray = [];
        InvCovR = eye(3);
        for NR = 1:20
            ChiSq  = 0.0;
            betaNR  = zeros(6,1);
            hessNR  = zeros(6);
            derivChiSq   = zeros(3,1);
            hessianChiSq = zeros(3);
            KeplerSat  = LagrangePlanetary(Kepler);
            rSatKepler = [];
            for ii =1:rows               
                time = track_data(ii,2)/TU;
                Rdata  = [track_data(ii,3:5)']/DU;
                Sensor = [Rdata+epsilon(1:3); zeros(6,1)];
                % Measurement Angles:
                %losM    = state_vector_losArray(1:3,ii);
                losM     = track_data(ii,6:8)';
                losM     = losM/norm(losM);   % los
                rangeM   = norm(losM);
                thetaM   = acos(losM(3)/rangeM);
                phiM     = atan2(losM(2), losM(1));
                Observed =[thetaM; phiM];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [Rinit,Vinit,Hessian,Gradient,Jacobian,JacobianLos,ParamList,...
                    GravityCan,HessianLos,Phi,MLagrange,RVHessian] = KeplerSat.J6Gravity(time, Sensor, Observed);
                deltaR   = Rinit - Rdata;
                rSatKepler = [rSatKepler, Rinit];
                %deltaR   = Rinit - Rdata + Jacobian(1:3,2:7)'*delZ;
                ChiSq    = ChiSq  + 0.5*deltaR'*InvCovR*deltaR;
                betaNR   = betaNR + Jacobian(1:3,2:7)'*InvCovR*deltaR;
                hessPart = zeros(6);
                for k = 1:3
                    hessPart = hessPart + RVHessian{k}(2:7,2:7)*deltaR(k);
                end
                %hessNR  = hessNR + Jacobian(1:3,2:7)'*InvCovR*Jacobian(1:3,2:7);
                hessNR  = hessNR + hessPart + Jacobian(1:3,2:7)'*InvCovR*Jacobian(1:3,2:7);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            ChiSqArray  =  [ChiSqArray, log(ChiSq)];
            deltaPar    = -hessNR\betaNR
            TangentArray = [TangentArray, log(norm(betaNR))]
            Kepler      =  Kepler + deltaPar'  
        end
        
        %KeplerSat  = LagrangePlanetary(Kepler);
        
        figure
        plot([1:numel(ChiSqArray)], ChiSqArray)
        title('log(\chi^2) per iteration')
        ylabel('log(\chi^2)')
        xlabel('Newton-Raphson Iteration Number')
        
        figure
        plot([1:numel(TangentArray)], TangentArray)
        title('Log norm of the Derivatives of \chi^2 per iteration')
        ylabel('Log norm \chi^2 Derivatives')
        xlabel('Newton-Raphson Iteration Number')
        
        figure
        %plot3(track_data(:,3), track_data(:,4),track_data(:,4))
        Satellites = unique(track_data(:,1))
        ChiSqSat      = [];
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
            iMid       = floor(numel(track_data(:,2))/2);
            tFit       = median(TrackOrbit(:,2))
            %tFit       = 0.;
            [state_vector_fit, fit, Iorder, ChiSq] = PolyFit(Object, TrackOrbit,tFit,Iorder,iType);
            %plot3(state_vector_fit(1,:),state_vector_fit(2,:),state_vector_fit(3,:),'linewidth',2);
            %hold on
            plot3(DU*rSatKepler(1,:),DU*rSatKepler(2,:),DU*rSatKepler(3,:),'linewidth',2);
            hold on
            %Diff = state_vector_fit(1:3,:) - TrackOrbit(:,3:5)';
            Diff = DU*rSatKepler(1:3,:) - TrackOrbit(:,3:5)';
            ChiSqSat = [ChiSqSat, sum(Diff.*Diff)];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ChiSqRoot = sqrt(sum(ChiSqSat));
        figure
        plot([1:numel(ChiSqSat)],ChiSqSat)
        title('Kepler Fit to Satellite ');
        xlabel('Observation Number');
        ylabel('\chi^2 from residuals');
        disp(['root(Summed Residials) for Kepler Satellite = ', num2str(sqrt(ChiSqRoot))])
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  Polynomial Interpolation for the Line of Sight
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Iorder = 10;
        Object = 2; % for Target
        [state_vector_losArray, fitLos, Iorder, ChiSq] = PolyFit(Object, track_data, tFit, Iorder, iType);
        
        ChiSqRoot = sqrt(sum(ChiSq));
        figure
        iData = 1:numel(ChiSq);
        plot(iData,ChiSq)
        title('Line of Sight Polynomial');
        xlabel('Observation Number');
        ylabel('\chi^2 from polynomial');
        disp(['root(Summed Residials) for LOS = ', num2str(sqrt(ChiSqRoot))])
        
        %state_vector_los               = state_vector_losArray(:,iMid);
        state_vector_los                = PolyExtrapToFitTime(fitLos, tFit);
        state_vector_los'
        
        if iFirst == 0
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Method for Initial State -- Laplace's Angles-Only Algorithm
            StateVectorECI   = LaplaceInECI(DU, mu, state_vector_sat, state_vector_los);
            StateVectorECI'
            SVInitial        = [StateVectorECI(1:6);tFit];
            % SV7Initial       = [SVDiffCorrection(1:6); tFit + tFirst];
            %SVDiffCorrection =  DifferentialCorrectionToLaplace(tFit, track_data, StateVectorECI)
            %
            % toTime           = 71579.507
            % SV7Initial       = [SVDiffCorrection(1:6); tFit + tFirst];
            % SV7Prop          = propagateECI(SV7Initial, toTime);
            % Mprop            = ECR_to_ECI_PVA_Rotation_Matrix( -toTime );
            % StateVectorECEF = [Mprop(1:6,1:6)*SVDiffCorrection(1:6)/1000.0]'
            % M = ECR_to_ECI_PVA_Rotation_Matrix( -toTime );
            % %StateVector = [rVec*DU; vVec*VU; tFit]
            % %rVec = rho*Rlos + Rpos*DU
            % %[W T] = Earth_Rotation_Matrix(-tFit);
            % %T     = T';
            % %T*rVec/1000.0
            % % Move ECI 9-State Vector Back to ECEF
            % StateVectorECEF = [M(1:6,1:6)*SVDiffCorrection(1:6)/1000.0]'
            % %Mrotate    = ECR_to_ECI_PVA_Rotation_Matrix( -tFit );
            % %StateVectorECEF2 = ECR_to_ECI_PVA(StateVectorECI,-tFit)/1000.0
            % SV_For_Comparison = [-6418.325,  2835.917,  1254.261, 2.143745,  1.678883,  9.287525]
            % SV_Diff = StateVectorECEF - SV_For_Comparison
            % Fraction = norm(SV_Diff)/norm(SV_For_Comparison)
            % % SVDiffCorrection is the 7-state Initial Guess from the
            % % application of Differential Corrections (numerically) to
            % % Laplace's "Angles Only" Orbit Determination Method. This
            % % StateVector has been propagated to the last time of the
            % % current block of LOS history Data.
            %SV6init = propagateECI(SVInitial, tFit);
            SV6init  = SVInitial;
            %iFirst   = -1;
        end
        time    = tFit/TU
        Rpos    = SV6init(1:3)/DU;
        Rdot    = SV6init(4:6)/VU;
        % comment out 506 to see what happens. !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        [e a Inclination omega Omega Mp Period Jacobian] = kepler(time, Rpos, Rdot);
        KeplerInit = [e a Inclination omega Omega Mp] % Canonical Units
        %KeplerInit = [0.9525    0.5965    0.5637    3.3546    1.5967    3.1799]
        %KeplerInit = [0.6129    0.7078    1.5828    3.2409    1.7374    3.1735]
        %KeplerInit = KeplerInitLC % Canonical Units
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if iFirst == 1
        %     SV6init  = propagateECI(SVInitial, tFit + tFirst);
        %     %SV6init = SVDiffCorrection
        %     %SV6init = StateVectorECI(1:6)
        % end
        % time    = (tFit+tFirst)/TU;
        % Rpos    = SV6init(1:3)/DU;
        % Rdot    = SV6init(4:6)/VU;
        % %StandardGravity    = gravityJ4(SV7Mid);
        % [e a Inclination omega Omega Mp Period Jacobian] = kepler(time, Rpos, Rdot);
        %[J6 Jperturb Jhom] = J6Gravity(Rpos);
        %[EscobalGravity Phi] = Escobal(Rpos*DU);
        %KeplerInit = [e a Inclination omega Omega Mp] % Canonical Units
        %lna       = log(a);
        meanMotion = exp(-1.5*log(a));
        KLagrange  = LagrangePlanetary(KeplerInit);
        %J6Check   = KLagrange.J6Gravity(time, J6);
        %J6Check   = KLagrange.J6Gravity(time, EscobalGravity*TU^2/DU);
        %[J6Check] = KLagrange.J6Gravity(time);
        %[Rextrap, Vextrap, Hessian, Gradient, Jacobian] = KLagrange.J6Gravity(time)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %J6Check      = KLagrange.J6Gravity(tTest, EscobalGravity*TU^2/DU)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Zhat           = zeros(6,1);
        delta          = zeros(6,1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %      Determine the Kepler Chief's parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        %KeplerIteration = KeplerInit;   % Reset the Gravity Model Reference Point
        %KeplerIteration = KeplerInit + Zhat';   % Reset the Gravity Model Reference Point
        %KLagrange       = LagrangePlanetary(KeplerIteration);
        Sensor     = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU];
        %Sensor = [track_data(1,3:5)'/DU; zeros(6,1)];
        % Measurement Angles:
        %losM    = state_vector_losArray(1:3,ii);
        losM     = track_data(iMid,6:8)';
        losM     = losM/norm(losM);   % los
        rangeM   = norm(losM);
        thetaM   = acos(losM(3)/rangeM);
        phiM     = atan2(losM(2), losM(1));
        Observed =[thetaM; phiM];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %[Rinit,Vinit,Hessian,Gradient,Jacobian, ~, ~, ~,~,~,~] = KLagrange.J6Gravity(time, Sensor);
        [Rinit0,Vinit0,Hessian0,Gradient0,Jacobian0,JacobianLos0,ParamList0,...
            GravityCan0,HessianLos0,Phi0,MLagrange0, RVHessian0] = KLagrange.J6Gravity(time, Sensor, Observed);
        
        [Zdot0, Inhomogeneous0, Perturbation0, DMatrix0, LagrangeTerms0, LagrangeGrads0] = KLagrange.EquationsOfMotion(Zhat);
        SVstn0 = [Rinit0;Vinit0;time];
                [EscobalGravity Phi] = Escobal(Rpos*DU);
        EscobalGravity*TU^2/DU
        GravityCan0

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %      Begin Finite Differences
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 
        % GradientCheck = [];
        % JacobianCheck = [];
        % LagrangeCheck = [];
        % MLagrangeCheck = cell(1,6);
        % HessLosCheck   = cell(1,6);
        % 
        % for LineCheck = 1:130
        %     [Jacob, JacobFinite, Hess, HessFinite] = KLagrange.DigitalJustice(LineCheck, Sensor, Observed);
        %     Jacob
        %     JacobFinite
        %     Hess
        %     HessFinite
        %     Sum = sum(sum(Hess - HessFinite));
        %     disp([' Line = ', num2str(LineCheck),'  sum = ', num2str(Sum)])
        % end
        % 
        % epsilon     = 0.0000001
        % KeplerDelta = epsilon*eye(6);
        % Kepler1     = KeplerIteration + KeplerDelta(1,:)
        % Kdelta1     = LagrangePlanetary(Kepler1);
        % [Rinit1,Vinit1,Hessian1,Gradient1,Jacobian1,JacobianLos1,ParamList1,...
        %     GravityCan1,HessianLos1,Phi1,MLagrange1, RVHessian1] = Kdelta1.J6Gravity(time, Sensor);
        % [Zdot1, Inhomogeneous1, Perturbation1, DMatrix1, LagrangeTerms1, LagrangeGrads1] = Kdelta1.EquationsOfMotion(Zhat);
        % %MLagrangeCheck{1} = (MLagrange1 - MLagrange0)/epsilon;
        % MLagrangeCheck{1} = (MLagrange1*Gradient1 - MLagrange0*Gradient0)/epsilon
        % GradientCheck = [GradientCheck; (Phi1-Phi0)/epsilon];
        % LagrangeCheck = [LagrangeCheck, (LagrangeTerms1-LagrangeTerms0)/epsilon];
        % SVstn1 = [Rinit1;Vinit1;time];
        % JacobianCheck = [JacobianCheck; (SVstn1-SVstn0)'/epsilon];
        % 
        % Kepler2     = KeplerIteration + KeplerDelta(2,:)
        % Kdelta2     = LagrangePlanetary(Kepler2);
        % [Rinit2,Vinit2,Hessian2,Gradient2,Jacobian2,JacobianLos2,ParamList2,...
        %     GravityCan2,HessianLos2,Phi2,MLagrange2, RVHessian2] = Kdelta2.J6Gravity(time, Sensor);
        % [Zdot2, Inhomogeneous2, Perturbation2, DMatrix2, LagrangeTerms2, LagrangeGrads2] = Kdelta2.EquationsOfMotion(Zhat);
        % %MLagrangeCheck{2} = (MLagrange2 - MLagrange0)/epsilon
        % MLagrangeCheck{2} = (MLagrange2*Gradient2 - MLagrange0*Gradient0)/epsilon;
        % GradientCheck = [GradientCheck; (Phi2-Phi0)/epsilon];
        % LagrangeCheck = [LagrangeCheck, (LagrangeTerms2-LagrangeTerms0)/epsilon];
        % SVstn2 = [Rinit2;Vinit2;time];
        % JacobianCheck = [JacobianCheck; (SVstn2-SVstn0)'/epsilon];
        % 
        % Kepler3     = KeplerIteration + KeplerDelta(3,:)
        % Kdelta3     = LagrangePlanetary(Kepler3);
        % [Rinit3,Vinit3,Hessian3,Gradient3,Jacobian3,JacobianLos3,ParamList3,...
        %     GravityCan3,HessianLos3,Phi3,MLagrange3, RVHessian3] = Kdelta3.J6Gravity(time, Sensor);
        % [Zdot3, Inhomogeneous3, Perturbation3, DMatrix3, LagrangeTerms3, LagrangeGrads3] = Kdelta3.EquationsOfMotion(Zhat);
        % %MLagrangeCheck{3} = (MLagrange3 - MLagrange0)/epsilon;
        % MLagrangeCheck{3} = (MLagrange3*Gradient3 - MLagrange0*Gradient0)/epsilon;
        % GradientCheck = [GradientCheck; (Phi3-Phi0)/epsilon];
        % LagrangeCheck = [LagrangeCheck, (LagrangeTerms3-LagrangeTerms0)/epsilon];
        % SVstn3 = [Rinit3;Vinit3;time];
        % JacobianCheck = [JacobianCheck; (SVstn3-SVstn0)'/epsilon];
        % 
        % Kepler4     = KeplerIteration + KeplerDelta(4,:)
        % Kdelta4     = LagrangePlanetary(Kepler4);
        % [Rinit4,Vinit4,Hessian,Gradient4,Jacobian4,JacobianLos4,ParamList4,...
        %     GravityCan4,HessianLos4,Phi4,MLagrange4,RVHessian4] = Kdelta4.J6Gravity(time, Sensor);
        % [Zdot4, Inhomogeneous4, Perturbation4, DMatrix4, LagrangeTerms4, LagrangeGrads4] = Kdelta4.EquationsOfMotion(Zhat);
        % %MLagrangeCheck{4} = (MLagrange1 - MLagrange0)/epsilon;
        % MLagrangeCheck{4} = (MLagrange4*Gradient4 - MLagrange0*Gradient0)/epsilon;
        % GradientCheck = [GradientCheck; (Phi4-Phi0)/epsilon];
        % LagrangeCheck = [LagrangeCheck, (LagrangeTerms4-LagrangeTerms0)/epsilon];
        % SVstn4 = [Rinit4;Vinit4;time];
        % JacobianCheck = [JacobianCheck; (SVstn4-SVstn0)'/epsilon];
        % 
        % Kepler5     = KeplerIteration + KeplerDelta(5,:)
        % Kdelta5     = LagrangePlanetary(Kepler5);
        % [Rinit5,Vinit5,Hessian5,Gradient5,Jacobian5,JacobianLos5,ParamList5,...
        %     GravityCan5,HessianLos5,Phi5,MLagrange5, RVHessian5] = Kdelta5.J6Gravity(time, Sensor);
        % [Zdot5, Inhomogeneous5, Perturbation5, DMatrix5, LagrangeTerms5, LagrangeGrads5] = Kdelta5.EquationsOfMotion(Zhat);
        % %MLagrangeCheck{5} = (MLagrange1 - MLagrange0)/epsilon;
        % MLagrangeCheck{5} = (MLagrange5*Gradient5 - MLagrange0*Gradient0)/epsilon;
        % GradientCheck = [GradientCheck; (Phi5-Phi0)/epsilon];
        % SVstn5 = [Rinit5;Vinit5;time];
        % LagrangeCheck = [LagrangeCheck, (LagrangeTerms5-LagrangeTerms0)/epsilon];
        % JacobianCheck = [JacobianCheck; (SVstn5-SVstn0)'/epsilon];
        % 
        % Kepler6     = KeplerIteration + KeplerDelta(6,:)
        % Kdelta6     = LagrangePlanetary(Kepler6);
        % [Rinit6,Vinit6,Hessian6,Gradient6,Jacobian6,JacobianLos6,ParamList6,...
        %     GravityCan6,HessianLos6,Phi6,MLagrange6, RVHessian6] = Kdelta6.J6Gravity(time, Sensor);
        % [Zdot6, Inhomogeneous6, Perturbation6, DMatrix6, LagrangeTerms6, LagrangeGrads6] = Kdelta6.EquationsOfMotion(Zhat);
        % %MLagrangeCheck{6} = (MLagrange1 - MLagrange0)/epsilon;
        % MLagrangeCheck{6} = (MLagrange6*Gradient6 - MLagrange0*Gradient0)/epsilon;
        % GradientCheck = [GradientCheck; (Phi6-Phi0)/epsilon];
        % SVstn6 = [Rinit6;Vinit6;time];
        % LagrangeCheck = [LagrangeCheck, (LagrangeTerms6-LagrangeTerms0)/epsilon];
        % JacobianCheck = [JacobianCheck; (SVstn6-SVstn0)'/epsilon];
        % 
        % GradientCompare = [GradientCheck, Gradient0, GradientCheck - Gradient0]
        % JacobianCheck(:,1:6)'
        % Jacobian0(:,2:7)
        % JacobianCheck(:,1:6)' - Jacobian0(:,2:7)
        % 
        % LagrangeCheck
        % LagrangeGrads0'
        % 
        % MLagrangeCheck{1} - Inhomogeneous0 
        % 
        % MLagrangeCheck{2} - Inhomogeneous0
        % 
        % MLagrangeCheck{3} - Inhomogeneous0
        % 
        % MLagrangeCheck{4} - Inhomogeneous0
        % 
        % MLagrangeCheck{5} - Inhomogeneous0
        % 
        % MLagrangeCheck{6} - Inhomogeneous0
        % 
        % for iCell = 1:6
        % HessLosCheck{iCell} = [(JacobianLos1(iCell,2:7) - JacobianLos0(iCell,2:7));...
        %                       (JacobianLos2(iCell,2:7) - JacobianLos0(iCell,2:7));...
        %                       (JacobianLos3(iCell,2:7) - JacobianLos0(iCell,2:7));...
        %                       (JacobianLos4(iCell,2:7) - JacobianLos0(iCell,2:7));...
        %                       (JacobianLos5(iCell,2:7) - JacobianLos0(iCell,2:7));...
        %                       (JacobianLos6(iCell,2:7) - JacobianLos0(iCell,2:7))]/epsilon
        % end
        % 
        % 
        % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %      End Finite Difference Checking    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 0)     Begin Monte Carlo Search for the Minimum
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        ChiSqMin  = 10^100;;
        ChiSq     = 0.0;
        TangentArray   = [];
        ZhatArray      = [];
        DelZhatArray   = [];
        
        ChiSqArray     = [];
        DelChiSqArray  = [];
        
        DChiSqArray    = [];
        DiffChiSqArray = [];

%         for NR = 1:10000
%             ChiSq     = 0.0;
%             R  = rand(6,1);
%             Kepler = [R(1), 2*R(2), pi*R(3), twopi*R(4), twopi*R(5), twopi*R(6)]; 
%             KeplerObj    = LagrangePlanetary(Kepler);            
%             for ii =1:rows               
%                 tRecord = track_data(ii,2)/TU;                
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%                 %losM    = state_vector_losArray(1:3,ii);
%                 losM     = track_data(ii,6:8)';
%                 losM     = losM/norm(losM);   % los
%                 rangeM   = norm(losM);
%                 thetaM   = acos(losM(3)/rangeM);
%                 phiM     = atan2(losM(2), losM(1));
%                 Observed = [thetaM; phiM];
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%                 % %Sat Position at measured position
%                 Rsensor  = [track_data(ii,3:5)']/DU;
%                 Sensor   = [Rsensor; zeros(6,1)];                
%                 % %Sat Position from Polynomial at tFit
%                 % Sensor = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU]
%                 
%                 % %Sat Position from Orbit Fit to Kepler elements:
%                 %[Rsensor, Vsensor, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessSat] = KeplerSat.J6Gravity(tRecord, Sensor, Observed);
%                 %Sensor   = [Rsensor; Vsensor; zeros(3,1)];
% 
%                 [Rinit,Vinit,Hessian,Gradient,Jacobian,JacobianLos,ParamList,...
%                     GravityCan,HessianLos,Phi,MLagrange,RVHessian] = KeplerObj.J6Gravity(tRecord, Sensor, Observed);
%                 los       = [ParamList(7); ParamList(8); ParamList(9)];
%                 los       = los/norm(los);
% 
%                 delta     = los - losM;
%                 deltaSQ   = delta'*delta;
%                 ChiSq    = ChiSq + deltaSQ;                                
%             end
%             if ChiSq < ChiSqMin
%                 ChiSqMin  = ChiSq
%                 KeplerMin = Kepler
%             end
%         end
%         
%         figure
%         plot([1:numel(ChiSqArray)], ChiSqArray)
%         title('log(\chi^2) per iteration')
%         ylabel('log(\chi^2)')
%         xlabel('Newton-Raphson Iteration Number')
%         
%         figure
%         plot([1:numel(TangentArray)], TangentArray)
%         title('Log norm of the Derivatives of \chi^2 per iteration')
%         ylabel('Log norm \chi^2 Derivatives')
%         xlabel('Newton-Raphson Iteration Number')        
%                 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % 1)     Begin Application of Hessian and Second Order fit    
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%         deltaPar       = zeros(6,1);
%         
%         TangentArray   = [];
%         ZhatArray      = [];
%         DelZhatArray   = [];
% 
%         ChiSqArray     = [];
%         DelChiSqArray  = [];
%         
%         DChiSqArray    = [];
%         DiffChiSqArray = [];
% 
%         PhiFundamental = eye(6);
%         InvCovMeas   = eye(2);
%         InvCovM      = eye(2);
%         deltaPar     = zeros(6,1);
%         epsilon      = 0.000000001*ones(6,1);         
%         Kepler       = KeplerInit
%         %Kepler       = KeplerInitLC
%         %Kepler = [   0.834837539441724   1.855977574696516   1.563985062845812   4.396785028396810   1.476396527817562   0.288360418410242 ]    
%         angle         = Kepler(4)
%         while(angle <   0.0)
%             angle = angle + twopi;
%         end
%         while(angle > twopi)
%             angle = angle - twopi;
%         end
%         Kepler(4) = angle
%         
%         angle         = Kepler(5)
%         while(angle <   0.0)
%             angle = angle + twopi;
%         end
%         while(angle > twopi)
%             angle = angle - twopi;
%         end
%         Kepler(5) = angle
%         
%         for NR = 1:100
%             
%             ChiSq        = 0.0;
%             BetaInc      = zeros(6,1);
%             HessianSyst  = zeros(6);
%             KeplerObj    = LagrangePlanetary(Kepler);
%             
%             for ii =1:rows               
%                 tRecord = track_data(ii,2)/TU;
%                 
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%                 %losM    = state_vector_losArray(1:3,ii);
%                 losM     = track_data(ii,6:8)';
%                 losM     = losM/norm(losM);   % los
%                 rangeM   = norm(losM);
%                 thetaM   = acos(losM(3)/rangeM);
%                 phiM     = atan2(losM(2), losM(1));
%                 Observed = [thetaM; phiM];
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%                 
%                 % %Sat Position at measured position
%                 Rsensor  = [track_data(ii,3:5)']/DU;
%                 Sensor   = [Rsensor; zeros(6,1)];
%                 
%                 % %Sat Position from Polynomial at tFit
%                 % Sensor = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU]
%                 
%                 % %Sat Position from Orbit Fit to Kepler elements:
%                 %[Rsensor, Vsensor, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessSat] = KeplerSat.J6Gravity(tRecord, Sensor, Observed);
%                 %Sensor   = [Rsensor; Vsensor; zeros(3,1)];
% 
%                 [Rinit,Vinit,Hessian,Gradient,Jacobian,JacobianLos,ParamList,...
%                     GravityCan,HessianLos,Phi,MLagrange,RVHessian] = KeplerObj.J6Gravity(tRecord, Sensor, Observed);
%                 los       = [ParamList(7); ParamList(8); ParamList(9)];
%                 los       = los/norm(los);
% 
%                 delta     = los - losM
%                 deltaSQ   = delta'*delta
%                 JacobPred = JacobianLos(7:9,2:7);                
%                 BetaInc   = BetaInc + JacobPred'*delta;
%                 
%                 Jac2nd    = JacobPred'*JacobPred;
%                 Hess      = zeros(6);
%                 for ik = 1:3
%                     Hess      = Hess + delta(ik)*HessianLos{6+ik}(2:7,2:7);
%                 end
%                 HessianSyst = HessianSyst + Hess + Jac2nd;
%                 
%                 ChiSq    = ChiSq + deltaSQ;
%                                 
%             end
%             ChiSqArray   = [ChiSqArray, log(ChiSq)];
%             delta        = -(HessianSyst)\BetaInc;
%             %deltaPar     = -hessNR\betaNR
%             %TangentArray = [TangentArray, log(norm(betaNR))]
%             %Kepler       =  Kepler + deltaPar'
%             TangentArray = [TangentArray, log(norm(BetaInc))]
%             Kepler       =  Kepler + delta'
%             % angle         = Kepler(4);
%             % while(angle <   0.0)
%             %     angle = angle + twopi;
%             % end
%             % while(angle > twopi)
%             %     angle = angle - twopi;
%             % end
%             % Kepler(4) = angle;
%             %
%             % angle         = Kepler(5);
%             % while(angle <   0.0)
%             %     angle = angle + twopi;
%             % end
%             % while(angle > twopi)
%             %     angle = angle - twopi;
%             % end
%             % Kepler(5) = angle;
%             
%         end
%         
%         figure
%         plot([1:numel(ChiSqArray)], ChiSqArray)
%         title('log(\chi^2) per iteration')
%         ylabel('log(\chi^2)')
%         xlabel('Newton-Raphson Iteration Number')
%         
%         figure
%         plot([1:numel(TangentArray)], TangentArray)
%         title('Log norm of the Derivatives of \chi^2 per iteration')
%         ylabel('Log norm \chi^2 Derivatives')
%         xlabel('Newton-Raphson Iteration Number')        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2)     Begin Application of Hessian and Second Order fit    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        deltaPar       = zeros(6,1);
        
        TangentArray   = [];
        ZhatArray      = [];
        DelZhatArray   = [];

        ChiSqArray     = [];
        DelChiSqArray  = [];
        
        DChiSqArray    = [];
        DiffChiSqArray = [];

        PhiFundamental = eye(6);
        InvCovMeas   = eye(2);
        InvCovM      = eye(2);
        deltaPar     = zeros(6,1);
        epsilon      = 0.000000001*ones(6,1);         
        Kepler       = KeplerInit
        %Kepler       = KeplerInitLC
                
        angle         = Kepler(4)
        while(angle <   0.0)
            angle = angle + twopi;
        end
        while(angle > twopi)
            angle = angle - twopi;
        end
        Kepler(4) = angle
        
        angle         = Kepler(5)
        while(angle <   0.0)
            angle = angle + twopi;
        end
        while(angle > twopi)
            angle = angle - twopi;
        end
        Kepler(5) = angle
        
        for NR = 1:50
            
            ChiSq        = 0.0;
            BetaInc      = zeros(6,1);
            HessianSyst  = zeros(6);
            KeplerObj    = LagrangePlanetary(Kepler);
            
            for ii =1:rows               
                tRecord = track_data(ii,2)/TU;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
                %losM    = state_vector_losArray(1:3,ii);
                losM     = track_data(ii,6:8)';
                losM     = losM/norm(losM);   % los
                rangeM   = norm(losM);
                thetaM   = acos(losM(3)/rangeM);
                phiM     = atan2(losM(2), losM(1));
                Observed = [thetaM; phiM];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
                
                % %Sat Position at measured position
                % Rsensor  = [track_data(ii,3:5)']/DU;
                % Sensor   = [Rsensor; zeros(6,1)];
                
                % %Sat Position from Polynomial at tFit
                % Sensor = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU]
                
                % %Sat Position from Orbit Fit to Kepler elements:
                [Rsensor, Vsensor, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessSat] = KeplerSat.J6Gravity(tRecord, Sensor, Observed);
                Sensor   = [Rsensor; Vsensor; zeros(3,1)];

                [Rinit,Vinit,Hessian,Gradient,Jacobian,JacobianLos,ParamList,...
                    GravityCan,HessianLos,Phi,MLagrange,RVHessian] = KeplerObj.J6Gravity(tRecord, Sensor, Observed);
                los       = [ParamList(7); ParamList(8); ParamList(9)];
                los       = los/norm(los);

                delta     = dot(los,losM) - 1.0;
                JacobPred = JacobianLos(7:9,2:7);
                JacFac    = JacobPred'*losM;
                BetaInc   = BetaInc + JacFac*delta;
                
                Jac2nd    = JacFac*JacFac';
                Hess      = zeros(6);
                for ik = 1:3
                    Hess      = Hess + losM(ik)*HessianLos{6+ik}(2:7,2:7);
                end
                Hess = Hess*delta;
                HessianSyst = HessianSyst + Hess + Jac2nd;
                
                ChiSq    = ChiSq + delta^2;
                                
            end
            ChiSqArray   = [ChiSqArray, log(ChiSq)];
            delta        = -(HessianSyst)\BetaInc;
            %deltaPar     = -hessNR\betaNR
            %TangentArray = [TangentArray, log(norm(betaNR))]
            %Kepler       =  Kepler + deltaPar'
            TangentArray = [TangentArray, log(norm(BetaInc))]
            Kepler       =  Kepler + delta'
            angle         = Kepler(4);
            while(angle <   0.0)
                angle = angle + twopi;
            end
            while(angle > twopi)
                angle = angle - twopi;
            end
            Kepler(4) = angle;
            
            angle         = Kepler(5);
            while(angle <   0.0)
                angle = angle + twopi;
            end
            while(angle > twopi)
                angle = angle - twopi;
            end
            Kepler(5) = angle;
            
        end
        
        figure
        plot([1:numel(ChiSqArray)], ChiSqArray)
        title('log(\chi^2) per iteration')
        ylabel('log(\chi^2)')
        xlabel('Newton-Raphson Iteration Number')
        
        figure
        plot([1:numel(TangentArray)], TangentArray)
        title('Log norm of the Derivatives of \chi^2 per iteration')
        ylabel('Log norm \chi^2 Derivatives')
        xlabel('Newton-Raphson Iteration Number')        
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % 3)     Begin Application of Hessian and Second Order fit    
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%         deltaPar       = zeros(6,1);
%         
%         TangentArray   = [];
%         ZhatArray      = [];
%         DelZhatArray   = [];
% 
%         ChiSqArray     = [];
%         DelChiSqArray  = [];
%         
%         DChiSqArray    = [];
%         DiffChiSqArray = [];
% 
%         PhiFundamental = eye(6);
%         %InvCovM     = eye(3);
%         InvCovMeas   = eye(2);
%         InvCovM      = eye(2);
%         deltaPar     = zeros(6,1);
%         epsilon      = 0.000000001*ones(6,1);         
%         Kepler       = KeplerInit
%         %Kepler       = KeplerInitLC
%         for NR = 1:100
%             
%             ChiSq        = 0.0;
%             InvCov       = zeros(6);
%             HessianSyst  = zeros(6);
%             hessNR       = zeros(6);
%             hessianChiSq = zeros(6);
% 
%             BetaInc      = zeros(6,1);
%             betaNR       = zeros(6,1);
%             derivChiSq   = zeros(6,1);
%             
%             KeplerObj    = LagrangePlanetary(Kepler);
%             for ii =1:rows               
%                 tRecord = track_data(ii,2)/TU;
%                 
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%                 % Measurement Angles: 
%                 %losM    = state_vector_losArray(1:3,ii);
%                 losM     = track_data(ii,6:8)';
%                 losM     = losM/norm(losM);   % los
%                 rangeM   = norm(losM);
%                 thetaM   = acos(losM(3)/rangeM);
%                 phiM     = atan2(losM(2), losM(1));
%                 Observed =[thetaM; phiM];
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%                 
%                 % %Sat Position from Polynomial at tFit
%                 % Sensor = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU]
%                 
%                 % %Sat Position from Polynomial at t_ii
%                 % Sensor = [state_vector_satArray(1:3,ii)/DU; state_vector_satArray(4:6,ii)/VU; state_vector_satArray(7:9,ii)/AU]
%                 
%                 % %Sat Position at measured position
%                 Rdata  = [track_data(ii,3:5)']/DU;
%                 Sensor = [Rdata; zeros(6,1)];
%                 %Sensor = [Rdata+epsilon(1:3); zeros(6,1)];
%                 
%                 % %Sat Position from Orbit Fit to Kepler elements: 
%                 %[Rsensor, Vsensor, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessSat] = KeplerSat.J6Gravity(tRecord, Sensor);
%                 %Sensor   = [Rsensor; Vsensor; zeros(3,1)];
%                 
%                 [Rinit,Vinit,Hessian,Gradient,Jacobian,JacobianLos,ParamList,...
%                     GravityCan,HessianLos,Phi,MLagrange,RVHessian] = KeplerObj.J6Gravity(tRecord, Sensor, Observed);
%                 los       = [ParamList(7); ParamList(8); ParamList(9)];
%                 los       = los/norm(los);
%                 theta     = ParamList(2);
%                 phi       = ParamList(3);
%                 Predicted = [theta; phi];
%                 % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %[Jacob, JacobFinite, Hess, HessFinite] = KeplerObj.DigitalJustice(108,Sensor)
%                 %JacobianLos(2,:) = Jacob
%                 %HessianLos{2}    = HessFinite
%                 %[Jacob, JacobFinite, Hess, HessFinite] = KeplerObj.DigitalJustice(109,Sensor)
%                 %JacobianLos(3,:) = Jacob
%                 %HessianLos{3}    = HessFinite
%                 % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % JacobLos = [];
%                 % %
%                 % F       = zeros(8,1);
%                 % F(7)    = 1;     % Theta
%                 % f       = losDerivs(los,F);
%                 % JacobLos = [JacobLos; f(1), f(2), f(3)];
%                 % %
%                 % F       = zeros(8,1);
%                 % F(8)    = 1;   % Phi
%                 % f       = losDerivs(los,F);
%                 % JacobLos = [JacobLos; f(1), f(2), f(3)];
%                 % JacRVToK  = Jacobian(1:3,2:7);
%                 % JacobAlphaToK = JacobLos*JacRVToK
%                 % JacobianLos(2:3,2:7)
%                 %
%                 %F       = zeros(8,1);
%                 %F(5)    = 1;      % Range
%                 %f       = losDerivs(los,F);
%                 %JacobLos = [JacobLos; f(1), f(2), f(3)];
%                 %
%                 % JacobMeas = inv(JacobLos);
%                 % % Drop the range column
%                 % JacobMeas = JacobMeas(:,1:2);
%                 % JacToZ    = Jacobian(1:3,2:7);
%                 % CombJac   = JacToZ'*JacobMeas;
%                 % InvCov    = InvCov + CombJac'*InvCovMeas*CombJac;
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %deltaA  =  Predicted - Observed;
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %   Theta_diff and Phi_diff
%                 thetaDiff =  theta - thetaM;
%                 phiDiff   =  phi   - phiM;
%                 while(phiDiff <   -pi)
%                     phiDiff = PhiDiff + twopi;
%                 end
%                 while(phiDiff > pi)
%                     phiDiff = phiDiff - twopi;
%                 end
%                 
%                 crossP = cross([cos(phiM), sin(phiM), 0], [cos(phi),sin(phi), 0] );
%                 dotP   =   dot([cos(phiM), sin(phiM), 0], [cos(phi),sin(phi), 0] );
%                 phiDiff = atan2(crossP,dotP);
%                 deltaA  =  [thetaDiff; phiDiff(3)]
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % Turn off Gravity !!!!!!
%                 %CombJac   = JacobianLos(2:3,2:7);
%                 % Include Gravity Perturbation
%                 CombJac   = JacobianLos(2:3,2:7)*PhiFundamental;
%                 %CombJac   = JacobianLos0(2:3,2:7)*PhiFundamental;
%                 InvCov     = InvCov + CombJac'*InvCovMeas*CombJac;
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % line of sight components:
%                 %Mdelta   = [ cos(theta)*cos(phi) -cos(theta)*sin(phi)  -sin(theta);...
%                 %            -sin(theta)*sin(phi)  sin(theta)*cos(phi)      0      ];
%                 %deltaA   = Mdelta*[los-losM];
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 
%                 CovVec   = InvCovMeas*deltaA;
%                 dataVec  = CombJac'*CovVec;
%                 BetaInc  = BetaInc + dataVec;
%                 ChiSq    = ChiSq  + 0.5*deltaA'*InvCovMeas*deltaA;
%                 
%                 Matrix   = zeros(6);
%                 % for icol = 1:2
%                 %     %Matrix = Matrix + HessianLos{icol}(2:7,2:7)*PhiFundamental*CovVec(icol)*PhiFundamental;
%                 %     Matrix = Matrix + PhiFundamental'*HessianLos{icol}(2:7,2:7)*PhiFundamental*CovVec(icol);
%                 %     %Matrix = Matrix + PhiFundamental'*HessianLos0{icol}(2:7,2:7)*PhiFundamental*CovVec(icol);
%                 % end
%                 
%                 % JRS 
%                 % The following cell addresses are WRONG!!! 1->2 and 2->3
%                 % But the code doesn't converge yet. 
%                 
%                 Matrix = Matrix + PhiFundamental'*HessianLos{2}(2:7,2:7)*PhiFundamental*CovVec(1);
%                 Matrix = Matrix + PhiFundamental'*HessianLos{3}(2:7,2:7)*PhiFundamental*CovVec(2);
%                 %HessianSyst = HessianSyst + PhiFundamental*PhiFundamental*Matrix;
%                 %HessianSyst = HessianSyst + PhiFundamental*Matrix*PhiFundamental;
%                 HessianSyst = HessianSyst + Matrix;
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %   End Theta_diff and Phi_diff
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 
%                 % deltaR   = Predicted - Observed
%                 % %deltaR   = Rinit - Rdata + Jacobian(1:3,2:7)'*delZ;
%                 % ChiSq    = ChiSq  + 0.5*deltaR'*InvCovM*deltaR;
%                 %
%                 % [Jacob, JacobFinite, Hess, HessFinite] = KLagrange.DigitalJustice(LineCheck,Sensor);
%                 %
%                 % betaNR   = betaNR + JacobianLos(2:3,2:7)'*InvCovM*deltaR;
%                 % hessPart = zeros(6);
%                 % for k = 2:3
%                 %     hessPart = hessPart + HessianLos{k}(2:7,2:7)*deltaR(k-1);
%                 % end
%                 % hessNR  = hessPart + JacobianLos(2:3,2:7)'*InvCovM*JacobianLos(2:3,2:7);
%                 % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % %   U & V views:
%                 % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % deltaA  =  [ParamList(10); ParamList(11)];
%                 % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % % Turn off Gravity !!!!!!
%                 % %CombJac   = JacobianLos(2:3,2:7);
%                 % % Include Gravity Perturbation
%                 % CombJac   = JacobianLos(10:11,2:7)*PhiFundamental;
%                 % %CombJac   = JacobianLos0(2:3,2:7)*PhiFundamental;
%                 % %InvCov     = InvCov + CombJac'*InvCovMeas*CombJac;                
%                 % CovVec   = InvCovMeas*deltaA;
%                 % dataVec  = CombJac'*CovVec;
%                 % %Jacobian(7:8,2:7)'*InvCovMeas*deltaA
%                 % BetaInc  = BetaInc + dataVec;
%                 % ChiSq    = ChiSq  + 0.5*deltaA'*InvCovMeas*deltaA;
%                 % Matrix   = zeros(6);
%                 % %for icol = 1:2
%                 %     %Matrix = Matrix + HessianLos{icol}(2:7,2:7)*PhiFundamental*CovVec(icol)*PhiFundamental;
%                 %     %Matrix = Matrix + PhiFundamental'*HessianLos0{icol}(2:7,2:7)*PhiFundamental*CovVec(icol);
%                 % %end
%                 % Matrix = Matrix + PhiFundamental'*HessianLos{10}(2:7,2:7)*PhiFundamental*CovVec(1);
%                 % Matrix = Matrix + PhiFundamental'*HessianLos{11}(2:7,2:7)*PhiFundamental*CovVec(2);
%                 % %HessianSyst = HessianSyst + PhiFundamental*PhiFundamental*Matrix;
%                 % %HessianSyst = HessianSyst + PhiFundamental*Matrix*PhiFundamental;
%                 % %HessianSyst = HessianSyst + CombJac'*InvCovMeas*CombJac;
%                 % HessianSyst = HessianSyst + Matrix + CombJac'*InvCovMeas*CombJac;
%                 % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % %   End of U & V views
%                 % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %[Jacob1, JacobFinite1, Hess1, HessFinite1] = KeplerObj.DigitalJustice(127, Sensor, Observed);
%                 %[Jacob2, JacobFinite2, Hess2, HessFinite2] = KeplerObj.DigitalJustice(128, Sensor, Observed);
%                 % Jacob1
%                 % JacobianLos(7,:)
%                 % Jacob2
%                 % JacobianLos(8,:)
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 
%             end
%             ChiSqArray   = [ChiSqArray, log(ChiSq)];
%             delta        = -(HessianSyst+InvCov)\BetaInc;
%             %delta        = -(HessianSyst)\BetaInc;
%             %deltaPar     = -hessNR\betaNR
%             %TangentArray = [TangentArray, log(norm(betaNR))]
%             %Kepler       =  Kepler + deltaPar'
%             TangentArray = [TangentArray, log(norm(BetaInc))]
%             Kepler       =  Kepler + delta'
%         end
%         
%         figure
%         plot([1:numel(ChiSqArray)], ChiSqArray)
%         title('log(\chi^2) per iteration')
%         ylabel('log(\chi^2)')
%         xlabel('Newton-Raphson Iteration Number')
%         
%         figure
%         plot([1:numel(TangentArray)], TangentArray)
%         title('Log norm of the Derivatives of \chi^2 per iteration')
%         ylabel('Log norm \chi^2 Derivatives')
%         xlabel('Newton-Raphson Iteration Number')        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 4)     Begin Application of Hessian and Second Order fit    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TangentArray   = [];
        ZhatArray      = [];
        DelZhatArray   = [];

        ChiSqArray     = [];
        DelChiSqArray  = [];
        
        DChiSqArray    = [];
        DiffChiSqArray = [];
        
        Zhat           = zeros(6,1);
        delta          = zeros(6,1);
        PhiFundamental = eye(6);
        %KeplerIteration  = KeplerInit
        KeplerIteration  = Kepler
        for NR = 1:400
            KeplerIteration  = KeplerIteration - delta'
            % angle            = KeplerIteration(4);
            % while(angle <   0.0)
            %     angle = angle + twopi;
            % end
            % while(angle > twopi)
            %     angle = angle - twopi;
            % end
            % KeplerIteration(4) = angle;
            %
            % angle         = KeplerIteration(5);
            % while(angle <   0.0)
            %     angle = angle + twopi;
            % end
            % while(angle > twopi)
            %     angle = angle - twopi;
            % end
            % KeplerIteration(5) = angle;
            
            KLagrange       = LagrangePlanetary(KeplerIteration);

            if NR == 1
                KLagrange       = LagrangePlanetary(KeplerIteration);
                %Sensor     = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU];
                Sensor = [track_data(1,3:5)'/DU; zeros(6,1)];
                % Measurement Angles:
                %losM    = state_vector_losArray(1:3,ii);
                losM     = track_data(1,6:8)';
                losM     = losM/norm(losM);   % los
                rangeM   = norm(losM);
                thetaM   = acos(losM(3)/rangeM);
                phiM     = atan2(losM(2), losM(1));
                Observed =[thetaM; phiM];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %[Rinit,Vinit,Hessian,Gradient,Jacobian, ~, ~, ~,~,~,~] = KLagrange.J6Gravity(time, Sensor);
                [Rinit,Vinit,Hessian,Gradient,Jacobian,JacobianLos,ParamList,GravityCan,HessianLos,Phi,MLagrange]...
                    = KLagrange.J6Gravity(time, Sensor, Observed);
                [Rpos, Rdot]
                [(Rpos-Rinit), (Rdot-Vinit)]
                [Zdot, Inhomogeneous, Perturbation, DMatrix, LagrangeTerms, LagrangeGrads]= KLagrange.EquationsOfMotion(zeros(6,1));
                %     SV6init  = StateVectorECI(1:6)
                %     SV6check = [Rinit*DU; Vinit*VU]
            end
            
            ChiSq       = 0.0;
            InvCov      = zeros(6);
            BetaInc     = zeros(6,1);
            HessianSyst = zeros(6);

            for ii = 1:rows
                tRecord   =  track_data(ii,2)/TU;
                t = tRecord;
                tInterval    = tRecord - time;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Measurement
                %losM    = state_vector_losArray(1:3,ii);
                losM    = track_data(ii,6:8)';
                losM    = losM/norm(losM);   % los
                rangeM  = norm(losM);
                thetaM  = acos(losM(3)/rangeM)
                phiM    = atan2(losM(2), losM(1))
                Observed =[thetaM; phiM];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Sat Position from Polynomial at tFit
                %Sensor = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU]
                % Sat Position from Polynomial at t_ii
                Sensor = [state_vector_satArray(1:3,ii)/DU; state_vector_satArray(4:6,ii)/VU; state_vector_satArray(7:9,ii)/AU];
                % Sat Position from Data Measurements
                % Sensor = [track_data(ii,3:5)'/DU; zeros(6,1)];
                % Sat Position from Orbit Fit to Kepler elements: 
                %[Rsensor, Vsensor, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessSat] = KeplerSat.J6Gravity(tRecord, Sensor, Observed);
                %Sensor = [Rsensor; Vsensor; zeros(3,1)];
                %[Rextrap, Vextrap, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan] = KLagrangeProp.J6Gravity(tRecord, Sensor);
                [Rextrap, Vextrap, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessLos] = KLagrange.J6Gravity(tRecord, Sensor, Observed);
                %[Rextrap, Vextrap, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessLos] = KLagrangeProp.J6Gravity(tRecord, Sensor, Observed);
                theta = ParamList(2)
                phi   = ParamList(3)
                InvCovMeas = eye(2)*10^(-2);
                los   = [ParamList(7); ParamList(8); ParamList(9)];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lo
                deltaA  = [theta - thetaM; phi - phiM];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lo
                %CombJac   = JacobianLos(10:11,2:7)*PhiFundamental;
                %InvCov     = InvCov + CombJac'*InvCovMeas*CombJac;                
                %Mdelta   = [ cos(theta)*cos(phi) cos(theta)*sin(phi)  -sin(theta);...
                %            -sin(phi)  cos(phi)      0      ];
                %deltaA   = Mdelta*[los-losM]
                %CovVec     = InvCovMeas*deltaA;
                %Matrix     = zeros(6);
                %Matrix = Matrix + PhiFundamental'*HessianLos{10}(2:7,2:7)*PhiFundamental*CovVec(1);
                %Matrix = Matrix + PhiFundamental'*HessianLos{11}(2:7,2:7)*PhiFundamental*CovVec(2);
                %HessianSyst = HessianSyst + Matrix;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lo
                Mdelta    = [ cos(thetaM)*cos(phiM) cos(thetaM)*sin(phiM)  -sin(thetaM);...
                              -sin(phiM)  cos(phiM)      0      ];
                InvCovM  = Mdelta'*Mdelta;
                %InvCovM  = eye(3);
                CombJac  = JacobianLos(7:9,2:7);
                deltaA   = los-losM;
                dataVec  = CombJac'*InvCovM*deltaA;
                CovVec    = InvCovM*deltaA;
                %Jac2nd    = CombJac'*InvCovM*CombJac;
                Matrix    = zeros(6);
                Matrix = Matrix + HessianLos{7}(2:7,2:7)*CovVec(1);
                Matrix = Matrix + HessianLos{8}(2:7,2:7)*CovVec(2);
                Matrix = Matrix + HessianLos{9}(2:7,2:7)*CovVec(3);
                HessianSyst = HessianSyst + Matrix;
                InvCov      = InvCov      + CombJac'*InvCovM*CombJac;               
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lo
                ChiSqInc   = 0.5*deltaA'*InvCovM*deltaA;
                ChiSq      = ChiSq   + ChiSqInc;
                
                BetaInc    = BetaInc + dataVec;
            end
            delta      = (HessianSyst+InvCov)\BetaInc
            Zhat     = Zhat - delta;
            TangentArray = [TangentArray, log(norm(BetaInc))]

            ParamComp = [delta, Zhat, KeplerInit' + Zhat]
            
            if NR > 1
                ZhatPrev       = ZhatArray(:,NR-1) - Zhat;
                DelZhatArray   = [DelZhatArray, ZhatPrev];
                ChiSqPrev      = ChiSqArray(NR-1);
                deltaChiSq     = ChiSqPrev - log(ChiSq);
                DelChiSqArray  = [DelChiSqArray, deltaChiSq];
                diffChiSq      =  DChiSqArray(NR-1) - log(BetaInc);
                DiffChiSqArray = [DiffChiSqArray, diffChiSq];
            end
            DChiSqArray = [DChiSqArray, log(BetaInc)];
            ChiSqArray  = [ChiSqArray, log(ChiSq)];
            ZhatArray   = [ZhatArray, Zhat];
        end
                        
        KeplerFit = KeplerIteration;
        CovFit    = inv(InvCov);
        
        figure
        plot([1:numel(ChiSqArray)], ChiSqArray)
        title('Log \chi^2 per iteration')
        ylabel('\chi^2')
        xlabel('Newton-Raphson Iteration Number')

        figure
        plot([1:size(DelZhatArray,2)], DelZhatArray)
        title('Parameter Differences per iteration')
        ylabel('\Delta Parameter Set')
        xlabel('Newton-Raphson Iteration Number')

        figure
        plot([1:size(DiffChiSqArray,2)], DiffChiSqArray)
        title('Change in \chi^2 derivative per iteration')
        ylabel('\Delta \chi^2 derivatives')
        xlabel('Newton-Raphson Iteration Number')

        figure
        plot([1:size(DelChiSqArray,2)], DelChiSqArray)
        title('Log \chi^2 differences per iteration')
        ylabel('\Delta \chi^2')
        xlabel('Newton-Raphson Iteration Number')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 5)     Begin Application of Hessian and Second Order fit    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TangentArray   = [];
        ZhatArray      = [];
        DelZhatArray   = [];

        ChiSqArray     = [];
        DelChiSqArray  = [];
        
        DChiSqArray    = [];
        DiffChiSqArray = [];
        
        Zhat           = zeros(6,1);
        delta          = zeros(6,1);
        PhiFundamental = eye(6); 
        for NR = 1:100
            %ii = 1            
            % Update the Reference Orbit Chief -> New Chief:
            KeplerIteration  = KeplerInit
            %KeplerIteration  = Kepler
            %KLagrange        = LagrangePlanetary(KeplerIteration);
            %KeplerIteration = KeplerInit + Zhat';   % Reset the Gravity Model Reference Point

            %KeplerIteration  = Kepler + Zhat;
            %Zhat           = zeros(6,1);
            KLagrange       = LagrangePlanetary(KeplerIteration);

            if NR == 1
                KLagrange       = LagrangePlanetary(KeplerIteration);
                %Sensor     = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU];
                Sensor = [track_data(1,3:5)'/DU; zeros(6,1)];
                % Measurement Angles:
                %losM    = state_vector_losArray(1:3,ii);
                losM     = track_data(1,6:8)';
                losM     = losM/norm(losM);   % los
                rangeM   = norm(losM);
                thetaM   = acos(losM(3)/rangeM);
                phiM     = atan2(losM(2), losM(1));
                Observed =[thetaM; phiM];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %[Rinit,Vinit,Hessian,Gradient,Jacobian, ~, ~, ~,~,~,~] = KLagrange.J6Gravity(time, Sensor);
                [Rinit,Vinit,Hessian,Gradient,Jacobian,JacobianLos,ParamList,GravityCan,HessianLos,Phi,MLagrange]...
                    = KLagrange.J6Gravity(time, Sensor, Observed);
                [Rpos, Rdot]
                [(Rpos-Rinit), (Rdot-Vinit)]
                [Zdot, Inhomogeneous, Perturbation, DMatrix, LagrangeTerms, LagrangeGrads]= KLagrange.EquationsOfMotion(zeros(6,1));
                %     SV6init  = StateVectorECI(1:6)
                %     SV6check = [Rinit*DU; Vinit*VU]
            end
            
            ChiSq       = 0.0;
            InvCov      = zeros(6);
            BetaInc     = zeros(6,1);
            HessianSyst = zeros(6);

            for ii = 1:rows
                tRecord   =  track_data(ii,2)/TU;
                t = tRecord;
                tInterval    = tRecord - time;
                %[Zdot, Inhomogeneous, Perturbation] = KLagrange.EquationsOfMotion(KeplerInit');
                %[Zdot, Inhomogeneous, Perturbation] = KLagrange.EquationsOfMotion(KeplerIteration');
                %[Zdot, Inhomogeneous, Perturbation] = KLagrange.EquationsOfMotion(zeros(6,1));
                %[Zdot, Inhomogeneous, Perturbation, DMatrix, LagrangeTerms, LagrangeGrads]= KLagrange.EquationsOfMotion(zeros(6,1));
                %[Zdot, Inhomogeneous, Perturbation, DMatrix, LagrangeTerms, LagrangeGrads]= KLagrange.EquationsOfMotion(delta);
                %t = time;
                Nbins            = 600;
                delt             = tInterval/Nbins;
                tRiem            = -delt/2.0;
                PhiInvIntegrated = zeros(6,1);
                for i = 1:Nbins
                    tRiem = tRiem + delt;
                    PhiInvIntegrated = PhiInvIntegrated + expm(-Perturbation*tRiem)*Inhomogeneous;
                    %PhiInvIntegrated = PhiInvIntegrated + ComplementarySoln(-Perturbation*tRiem);
                end
                PhiInvIntegrated = PhiInvIntegrated*delt;
                PhiFundamental   = expm(Perturbation*tInterval);
                %PhiFundamental   = ComplementarySoln(Perturbation*tInterval)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % %Set One
                % delZ             = PhiFundamental*(PhiInvIntegrated*Inhomogeneous);
                % delZ'
                % KeplerCurrent      = KeplerInit      + delZ';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Set Two
                %delZ            = PhiFundamental*(Zhat + PhiInvIntegrated*Inhomogeneous);
                %delZ             = PhiFundamental*(PhiInvIntegrated*Inhomogeneous);
                %delZ             = PhiFundamental*PhiInvIntegrated;
                delZ            = PhiFundamental*(Zhat + PhiInvIntegrated);
                delZ'
                %KeplerCurrent     = KeplerInit      + delZ' + Zhat'
                KeplerCurrent     = KeplerIteration + delZ';
                %KeplerCurrent      = KeplerIteration;
                %KeplerCurrent      = KeplerInit      + delZ';
                %KeplerCurrent      = KeplerInit;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                KLagrangeProp      = LagrangePlanetary(KeplerCurrent);
                % Error in following call: check out order of calls:
                %[Zdot, Inhomogeneous, Perturbation] = KLagrangeProp.EquationsOfMotion(KeplerCurrent');
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % %KeplerProp = [e a Inclination omega Omega Mp]  % Canonical Units
                % M                    = meanMotion*t + Mp;
                % Kepler               = [e a Inclination omega Omega M] + delZ';
                % KLagrange            = LagrangePlanetary(Kepler);
                % %t = 0.0
                % [Rextrap, Vextrap]   = KLagrange.extrapolate(0)
                % [Rextrap, Vextrap, Hessian, Gradient, Jacobian] = KLagrange.J6Gravity(0)
                % Rextrap
                % Vextrap
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %t = 0.0
                %[Rextrap, Vextrap]   = KLagrange.extrapolate(0);
                %Rextrap
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Measurement Angles:
                %losM    = state_vector_losArray(1:3,ii);
                losM     = track_data(ii,6:8)';
                losM     = losM/norm(losM);   % los
                rangeM   = norm(losM);
                thetaM   = acos(losM(3)/rangeM);
                phiM     = atan2(losM(2), losM(1));
                Observed =[thetaM; phiM];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Sat Position from Polynomial at tFit
                %Sensor = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU]
                % Sat Position from Polynomial at t_ii
                %Sensor = [state_vector_satArray(1:3,ii)/DU; state_vector_satArray(4:6,ii)/VU; state_vector_satArray(7:9,ii)/AU]
                % Sat Position from Data Measurements
                % Sensor = [track_data(ii,3:5)'/DU; zeros(6,1)];
                % Sat Position from Orbit Fit to Kepler elements: 
                [Rsensor, Vsensor, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessSat] = KeplerSat.J6Gravity(tRecord, Sensor, Observed);
                Sensor = [Rsensor; Vsensor; zeros(3,1)];
                %[Rextrap, Vextrap, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan] = KLagrangeProp.J6Gravity(tRecord, Sensor);
                [Rextrap, Vextrap, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessLos] = KLagrangeProp.J6Gravity(tRecord, Sensor, Observed);
                theta = ParamList(2)
                phi   = ParamList(3)
                InvCovMeas = eye(2)*10^(-2);
                los   = [ParamList(7); ParamList(8); ParamList(9)];
                %InvCovMeas = eye(2)*10^(6);
                %InvCovMeas = eye(2);
                % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % % %Delta = Jacobian(:,2:7)*delZ;
                % % %Rextrap + Delta(1:3);
                % % %Vextrap + Delta(4:6);
                % % %los     = Rextrap - state_vector_satArray(1:3,ii)/DU;
                % los      = Rextrap -            track_data(ii,3:5)'/DU
                % % % On the fly Backwards Differentiation:
                % % % Prediction:
                % % % lx = los(1)    Line 1
                % % % ly = los(2)    Line 2
                % % % lz = los(3)    Line 3
                % % %range   = norm(los)
                % % %theta   = acos(los(3)/range)
                % % %phi     = atan2(los(2),los(1))
                % rangeSq  = los(1)^2 + los(2)^2 + los(3)^2;   % line 4
                % range    = sqrt(rangeSq)                    % line 5
                % u        = los(3)/range;                     % line 6
                % theta   = acos(u)                           % line 7
                % phi     = atan2(los(2),los(1))              % line 8
                % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % JacobLos = [];
                % %
                % F       = zeros(8,1);
                % F(7)    = 1;     % Theta
                % f       = losDerivs(los,F);
                % JacobLos = [JacobLos; f(1), f(2), f(3)];
                % %
                % F       = zeros(8,1);
                % F(8)    = 1;   % Phi
                % f       = losDerivs(los,F);
                % JacobLos = [JacobLos; f(1), f(2), f(3)];
                % JacRVToK  = Jacobian(1:3,2:7);
                % JacobAlphaToK = JacobLos*JacRVToK
                % JacobianLos(2:3,2:7)
                %
                %F       = zeros(8,1);
                %F(5)    = 1;      % Range
                %f       = losDerivs(los,F);
                %JacobLos = [JacobLos; f(1), f(2), f(3)];
                %
                % JacobMeas = inv(JacobLos);
                % % Drop the range column
                % JacobMeas = JacobMeas(:,1:2);
                % JacToZ    = Jacobian(1:3,2:7);
                % CombJac   = JacToZ'*JacobMeas;
                % InvCov    = InvCov + CombJac'*InvCovMeas*CombJac;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Turn off Gravity !!!!!!
                % CombJac   = JacobianLos(2:3,2:7);
                % %Include Gravity Perturbation 
                % CombJac   = JacobianLos(2:3,2:7)*PhiFundamental;
                % %CombJac   = JacobianLos0(2:3,2:7)*PhiFundamental;
                % InvCov     = InvCov + CombJac'*InvCovMeas*CombJac;                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % %Measurement
                % %losM    = state_vector_losArray(1:3,ii);
                % %losM    = track_data(ii,6:8)';
                % %losM    = losM/norm(losM);   % los
                % rangeM  = norm(losM);
                % thetaM  = acos(losM(3)/rangeM)
                % phiM    = atan2(losM(2), losM(1))
                % deltaA  = [theta - thetaM; phi - phiM];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                CombJac   = JacobianLos(10:11,2:7)*PhiFundamental;
                % %CombJac   = JacobianLos0(2:3,2:7)*PhiFundamental;
                InvCov     = InvCov + CombJac'*InvCovMeas*CombJac;                
                % line of sight components: JRS
                Mdelta   = [ cos(theta)*cos(phi) cos(theta)*sin(phi)  -sin(theta);...
                            -sin(phi)  cos(phi)      0      ];
                deltaA   = Mdelta*[los-losM]
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lo
                ChiSqInc   = 0.5*deltaA'*InvCovMeas*deltaA;
                ChiSq      = ChiSq   + ChiSqInc;
                
                CovVec     = InvCovMeas*deltaA;
                dataVec    = CombJac'*CovVec;
                %BetaInc    = BetaInc + dataVec;
                %SymCheck   = [deltaA'*InvCovMeas*CombJac]'
                BetaInc    = BetaInc + dataVec;
                Matrix     = zeros(6);
                Matrix = Matrix + PhiFundamental'*HessianLos{10}(2:7,2:7)*PhiFundamental*CovVec(1);
                Matrix = Matrix + PhiFundamental'*HessianLos{11}(2:7,2:7)*PhiFundamental*CovVec(2);
                % for icol   = 1:2
                %     %Matrix = Matrix + HessianLos{icol}(2:7,2:7)*PhiFundamental*CovVec(icol)*PhiFundamental;
                %     Matrix = Matrix + PhiFundamental'*HessianLos{icol}(2:7,2:7)*PhiFundamental*CovVec(icol);
                %     %Matrix = Matrix + PhiFundamental'*HessianLos0{icol}(2:7,2:7)*PhiFundamental*CovVec(icol);
                % end
                %HessianSyst = HessianSyst + PhiFundamental*PhiFundamental*Matrix;
                %HessianSyst = HessianSyst + PhiFundamental*Matrix*PhiFundamental;
                HessianSyst = HessianSyst + Matrix;
            end
            delta      = (HessianSyst+InvCov)\BetaInc
            %delParam   = InvCov\BetaInc
            %KeplerInit = KeplerInit + delParam';
            %KeplerInit'
            %delZhat  = -delParam - Zhat; 
            %Zhat     = -delParam
            Zhat     = Zhat - delta;
            
            ParamComp = [delta, Zhat, KeplerInit' + Zhat]
            
            if NR > 1
                ZhatPrev       = ZhatArray(:,NR-1) - Zhat;
                DelZhatArray   = [DelZhatArray, ZhatPrev];
                ChiSqPrev      = ChiSqArray(NR-1);
                deltaChiSq     = ChiSqPrev - log(ChiSq);
                DelChiSqArray  = [DelChiSqArray, deltaChiSq];
                diffChiSq      =  DChiSqArray(NR-1) - log(BetaInc);
                DiffChiSqArray = [DiffChiSqArray, diffChiSq];
            end
            DChiSqArray = [DChiSqArray, log(BetaInc)];
            ChiSqArray  = [ChiSqArray, log(ChiSq)];
            ZhatArray   = [ZhatArray, Zhat];
            %KeplerProp = [e a Inclination omega Omega Mp]  % Canonical Units
            %M               = meanMotion*t + Mp;
            %Kepler          = [e a Inclination omega Omega M];
            %KeplerInit       = [e a Inclination omega Omega Mp]+ Zhat'
            %KLagrange        = LagrangePlanetary(KeplerInit);
        end
                        
        KeplerFit = KeplerIteration;
        CovFit    = inv(InvCov);
        
        figure
        plot([1:numel(ChiSqArray)], ChiSqArray)
        title('Log \chi^2 per iteration')
        ylabel('\chi^2')
        xlabel('Newton-Raphson Iteration Number')

        figure
        plot([1:size(DelZhatArray,2)], DelZhatArray)
        title('Parameter Differences per iteration')
        ylabel('\Delta Parameter Set')
        xlabel('Newton-Raphson Iteration Number')

        figure
        plot([1:size(DiffChiSqArray,2)], DiffChiSqArray)
        title('Change in \chi^2 derivative per iteration')
        ylabel('\delta \chi^2 derivatives')
        xlabel('Newton-Raphson Iteration Number')

        figure
        plot([1:size(DelChiSqArray,2)], DelChiSqArray)
        title('\chi^2 differences per iteration')
        ylabel('\Delta \chi^2')
        xlabel('Newton-Raphson Iteration Number')
        %hold off
        % Update the Reference Orbit Chief -> New Chief:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        losM    = state_vector_losArray(1:3,ii);
        %losM     = track_data(ii,6:8)';
        losM     = losM/norm(losM);   % los
        rangeM   = norm(losM);
        thetaM   = acos(losM(3)/rangeM);
        phiM     = atan2(losM(2), losM(1));
        Observed = [thetaM; phiM];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        KLagrange  = LagrangePlanetary(KeplerIteration+Zhat');
        Sensor     = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU];
        [Rpos, Rdot, Hessian, Gradient, Jacobian, ~, ~] = KLagrange.J6Gravity(time, Sensor, Observed);
        [e a Inclination omega Omega Mp Period Jacobian] = kepler(time, Rpos, Rdot)

        %tFit        = tFit;
        %toTime     = track_data(end,2)
        %toTime    = track_data(end,2) + tFirst
        %toTime    = 71577.703  
        toTime     = 71579.507
        SV7Fit     = [Rpos*DU; Rdot*VU; tFit + tFirst];
        SV7Prop    = propagateECI(SV7Fit, toTime);

        SV7Fit'
        SV7Prop'
        Mrotate           = ECR_to_ECI_PVA_Rotation_Matrix( -toTime );
        StateVectorECEF   = [Mrotate(1:6,1:6)*SV7Prop(1:6)/1000.0]'
        SV_For_Comparison = [-6418.325,  2835.917,  1254.261, 2.143745,  1.678883,  9.287525]

        % StateVectorECI(1:6)
        % Mrotate    = ECR_to_ECI_PVA_Rotation_Matrix( -tFit );
        %
        % toTime           = 71579.507
        % SV7Initial       = [SVDiffCorrection(1:6); tFit + tFirst];
        % SV7Prop          = propagateECI(SV7Initial, toTime);
        % Mprop            = ECR_to_ECI_PVA_Rotation_Matrix( -toTime );
        % StateVectorECEF = [Mprop(1:6,1:6)*SVDiffCorrection(1:6)/1000.0]'
        % SV_For_Comparison = [-6418.325,  2835.917,  1254.261, 2.143745,  1.678883,  9.287525]'
        %
        % M = ECR_to_ECI_PVA_Rotation_Matrix( -tFit );
        % %StateVector = [rVec*DU; vVec*VU; tFit]
        % %rVec = rho*Rlos + Rpos*DU
        % %[W T] = Earth_Rotation_Matrix(-tFit);
        % %T     = T';
        % %T*rVec/1000.0
        % % Move ECI 9-State Vector Back to ECEF
        % StateVectorECEF = M(1:6,1:6)*SVDiffCorrection(1:6)/1000.0

        %figure
        %histogram(los_difference)

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Try to Integrate the Equations of Motion for the Classical Elements
        % Y            = zeros(6,1);
        % options = odeset('RelTol', 1e-10, 'AbsTol', 1e-13);
        % [Zdot, Inhomogeneous, Perturbation] = KLagrange.EquationsOfMotion(KeplerInit');
        % Y            = zeros(6,1);
        % % Propagate from fit Reference time back to dataTime(ii)
        % [Tout, Yout] = ode45(@(tint,y) KLagrange.EquationsOfMotion(y), [time, tRecord], Y, options);
        % Yout(end,:)
        % %KeplerProp = [e a Inclination omega Omega Mp]  % Canonical Units
        % M = meanMotion*t + Mp;
        % %Kepler   = [e a Inclination omega Omega M];  % Canonical Units
        % %Kepler   = [e a Inclination omega Omega M] + Yout(end,:);
        % Kepler    = [e a Inclination omega Omega M] + delZ';
        % KLagrange = LagrangePlanetary(Kepler);
        % %t = 0.0
        % [Rextrap, Vextrap]   = KLagrange.extrapolate(0);
        % SV7     = [StateVectorECI(1:6); tFit];
        % SV7Prop = propagateECI(SV7, t*DU);
        % Rpos    = SV7Prop(1:3)/DU;
        % Rdot    = SV7Prop(4:6)/VU;
        % [e a Inclination omega Omega Mp Period Jacobian] = kepler(t, Rpos, Rdot);
        % [e a Inclination omega Omega Mp]
        % Kepler
        % [Rpos, Rdot]
        % [Rextrap; Vextrap]'
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % % Attempt to Integrate the non-linear Lagrange Planetary Equations
        % % Y            = KeplerInit';
        % % KeplerDot    = FirstOrder(time,Y)
        % % options      = odeset('RelTol', 1e-8, 'AbsTol', 1e-13);
        % % [Tout, Yout] = ode45(@(tint,y) FirstOrder(tint,y), [time, tRecord], Y, options);
        % % Yout(end,:)'
        % % KeplerProp         = Yout(end,:);
        % % KLagrange          = LagrangePlanetary(KeplerProp);
        % % [Rextrap, Vextrap] = KLagrange.extrapolate(t);
        % % [Rextrap; Vextrap]'
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    end

    function [KeplerFit, CovFit, tFit, SV7Fit] = LeastSquaresAll(track_data, tFirst, Isat, iFirst, SV7Initial)
    
        orbitparams = wgs84Constants;
        TU = orbitparams.TU;
        DU = orbitparams.DU;
        VU = orbitparams.VU;
        AU = orbitparams.AU;
        %mu = orbitparams.mu;
        mu = orbitparams.Grav;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Index = find(track_data(:,1) == Isat);
        %track_data = track_data(Index,:);
        rows       = size(track_data,1);
        iMid       = floor(numel(track_data(:,2))/2);
        tFit       = median(track_data(:,2))
        
        %KeplerInit = [0.9541    0.5972    0.3180    3.5156    1.4248    3.1796]
        %KeplerInit = [0.0001    0.5    0.3    3.5    1.4    3.]
        KeplerInit = [0.6129    0.7078    1.5828    3.2409    1.7374    3.1735]
        KLagrange  = LagrangePlanetary(KeplerInit);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %      Begin Application of Hessian and Second Order fit    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        deltaPar       = zeros(6,1);
        
        TangentArray   = [];
        ZhatArray      = [];
        DelZhatArray   = [];

        ChiSqArray     = [];
        DelChiSqArray  = [];
        
        DChiSqArray    = [];
        DiffChiSqArray = [];

        PhiFundamental = eye(6);
        InvCovMeas   = eye(2);
        InvCovM      = eye(2);
        deltaPar     = zeros(6,1);
        epsilon      = 0.000000001*ones(6,1);         
        Kepler       = KeplerInit
        %Kepler       = KeplerInitLC
        for NR = 1:100
            
            ChiSq        = 0.0;
            BetaInc      = zeros(6,1);
            HessianSyst  = zeros(6);
            KeplerObj    = LagrangePlanetary(Kepler);
            
            for ii =1:rows               
                tRecord = track_data(ii,2)/TU;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
                %losM    = state_vector_losArray(1:3,ii);
                losM     = track_data(ii,6:8)';
                losM     = losM/norm(losM);   % los
                rangeM   = norm(losM);
                thetaM   = acos(losM(3)/rangeM);
                phiM     = atan2(losM(2), losM(1));
                Observed = [thetaM; phiM];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
                
                % %Sat Position at measured position
                Rsensor  = [track_data(ii,3:5)']/DU;
                Sensor   = [Rsensor; zeros(6,1)];
                
                % %Sat Position from Polynomial at tFit
                % Sensor = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU]
                
                % %Sat Position from Orbit Fit to Kepler elements:
                %[Rsensor, Vsensor, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessSat] = KeplerSat.J6Gravity(tRecord, Sensor, Observed);
                %Sensor   = [Rsensor; Vsensor; zeros(3,1)];

                [Rinit,Vinit,Hessian,Gradient,Jacobian,JacobianLos,ParamList,...
                    GravityCan,HessianLos,Phi,MLagrange,RVHessian] = KeplerObj.J6Gravity(tRecord, Sensor, Observed);
                los       = [ParamList(7); ParamList(8); ParamList(9)];
                los       = los/norm(los);

                delta     = dot(los,losM) - 1.0;
                JacobPred = JacobianLos(7:9,2:7);
                JacFac    = JacobPred'*losM;
                BetaInc   = BetaInc + JacFac*delta;
                
                Jac2nd    = JacFac*JacFac';
                Hess      = zeros(6);
                for ik = 1:3
                    Hess      = Hess + losM(ik)*HessianLos{6+ik}(2:7,2:7);
                end
                Hess = Hess*delta;
                HessianSyst = HessianSyst + Hess + Jac2nd;
                
                ChiSq    = ChiSq + delta^2;
                                
            end
            ChiSqArray   = [ChiSqArray, log(ChiSq)];
            delta        = -(HessianSyst)\BetaInc;
            %deltaPar     = -hessNR\betaNR
            %TangentArray = [TangentArray, log(norm(betaNR))]
            %Kepler       =  Kepler + deltaPar'
            TangentArray = [TangentArray, log(norm(BetaInc))]
            Kepler       =  Kepler + delta'
        end
        
        figure
        plot([1:numel(ChiSqArray)], ChiSqArray)
        title('log(\chi^2) per iteration')
        ylabel('log(\chi^2)')
        xlabel('Newton-Raphson Iteration Number')
        
        figure
        plot([1:numel(TangentArray)], TangentArray)
        title('Log norm of the Derivatives of \chi^2 per iteration')
        ylabel('Log norm \chi^2 Derivatives')
        xlabel('Newton-Raphson Iteration Number')        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        Zhat           = zeros(6,1);
        delta          = zeros(6,1);
        
        ZhatArray      = [];
        DelZhatArray   = [];

        ChiSqArray     = [];
        DelChiSqArray  = [];
        
        DChiSqArray    = [];
        DiffChiSqArray = [];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %      Begin Application of Hessian and Second Order fit    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        deltaPar     = zeros(6,1);
        ChiSqArray   = [];
        TangentArray = [];
        PhiFundamental = eye(6);
        %InvCovM     = eye(3);
        InvCovMeas   = eye(2);
        InvCovM      = eye(2);
        deltaPar     = zeros(6,1);
        epsilon      = 0.000000001*ones(6,1);         
        ChiSqArray   = [];
        Kepler       = KeplerInit
        %Kepler       = KeplerInitLC
        for NR = 1:600
            
            ChiSq        = 0.0;
            InvCov       = zeros(6);
            HessianSyst  = zeros(6);
            hessNR       = zeros(6);
            hessianChiSq = zeros(6);

            BetaInc      = zeros(6,1);
            betaNR       = zeros(6,1);
            derivChiSq   = zeros(6,1);
            
            KeplerObj    = LagrangePlanetary(Kepler);
            for ii =1:rows               
                tRecord = track_data(ii,2)/TU;
                
                % %Sat Position from Polynomial at tFit
                % Sensor = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU]
                
                % %Sat Position from Polynomial at t_ii
                % Sensor = [state_vector_satArray(1:3,ii)/DU; state_vector_satArray(4:6,ii)/VU; state_vector_satArray(7:9,ii)/AU]
                
                % %Sat Position at measured position
                Rdata  = [track_data(ii,3:5)']/DU;
                Sensor = [Rdata; zeros(6,1)];
                %Sensor = [Rdata+epsilon(1:3); zeros(6,1)];
                
                % %Sat Position from Orbit Fit to Kepler elements: 
                %[Rsensor, Vsensor, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessSat] = KeplerSat.J6Gravity(tRecord, Sensor);
                %Sensor   = [Rsensor; Vsensor; zeros(3,1)];
                
                [Rinit,Vinit,Hessian,Gradient,Jacobian,JacobianLos,ParamList,...
                    GravityCan,HessianLos,Phi,MLagrange,RVHessian] = KeplerObj.J6Gravity(tRecord, Sensor);
                theta     = ParamList(2);
                phi       = ParamList(3);
                Predicted = [theta; phi];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Turn off Gravity !!!!!!
                %CombJac   = JacobianLos(2:3,2:7);
                % Include Gravity Perturbation 
                CombJac   = JacobianLos(2:3,2:7)*PhiFundamental;
                %CombJac   = JacobianLos0(2:3,2:7)*PhiFundamental;
                InvCov     = InvCov + CombJac'*InvCovMeas*CombJac;                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
                
                % Measurement Angles: 
                %losM    = state_vector_losArray(1:3,ii);
                losM     = track_data(ii,6:8)';
                losM     = losM/norm(losM);   % los
                rangeM   = norm(losM);
                thetaM   = acos(losM(3)/rangeM);
                phiM     = atan2(losM(2), losM(1));
                Observed =[thetaM; phiM];
                
                deltaA   = Predicted - Observed;
                CovVec   = InvCovMeas*deltaA;
                dataVec  = CombJac'*CovVec;
                BetaInc  = BetaInc + dataVec;
                ChiSq    = ChiSq  + 0.5*deltaA'*InvCovMeas*deltaA;

                Matrix   = zeros(6);
                for icol = 1:2
                    %Matrix = Matrix + HessianLos{icol}(2:7,2:7)*PhiFundamental*CovVec(icol)*PhiFundamental;
                    Matrix = Matrix + PhiFundamental'*HessianLos{icol}(2:7,2:7)*PhiFundamental*CovVec(icol);
                    %Matrix = Matrix + PhiFundamental'*HessianLos0{icol}(2:7,2:7)*PhiFundamental*CovVec(icol);
                end
                %HessianSyst = HessianSyst + PhiFundamental*PhiFundamental*Matrix;
                %HessianSyst = HessianSyst + PhiFundamental*Matrix*PhiFundamental;
                HessianSyst = HessianSyst + Matrix;
                
                % deltaR   = Predicted - Observed
                % %deltaR   = Rinit - Rdata + Jacobian(1:3,2:7)'*delZ;
                % ChiSq    = ChiSq  + 0.5*deltaR'*InvCovM*deltaR;
                % 
                % [Jacob, JacobFinite, Hess, HessFinite] = KLagrange.DigitalJustice(LineCheck,Sensor);
                % 
                % betaNR   = betaNR + JacobianLos(2:3,2:7)'*InvCovM*deltaR;
                % hessPart = zeros(6);
                % for k = 2:3
                %     hessPart = hessPart + HessianLos{k}(2:7,2:7)*deltaR(k-1);
                % end
                % hessNR  = hessPart + JacobianLos(2:3,2:7)'*InvCovM*JacobianLos(2:3,2:7);                
                
            end
            ChiSqArray   = [ChiSqArray, log(ChiSq)];        
            delta        = -(HessianSyst+InvCov)\BetaInc;           
            %deltaPar     = -hessNR\betaNR
            %TangentArray = [TangentArray, log(norm(betaNR))]
            %Kepler       =  Kepler + deltaPar'  
            TangentArray = [TangentArray, log(norm(BetaInc))]
            Kepler       =  Kepler + delta'  
        end
        
        figure
        plot([1:numel(ChiSqArray)], ChiSqArray)
        title('log(\chi^2) per iteration')
        ylabel('log(\chi^2)')
        xlabel('Newton-Raphson Iteration Number')
        
        figure
        plot([1:numel(TangentArray)], TangentArray)
        title('Log norm of the Derivatives of \chi^2 per iteration')
        ylabel('Log norm \chi^2 Derivatives')
        xlabel('Newton-Raphson Iteration Number')        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    end

    function [expAt] = ComplementarySoln(A)
             dim = size(A,2);
             Inc = eye(dim);
             expAt = zeros(6);
             for i = 1:20
                 expAt = expAt + Inc; 
                 Inc   = Inc*A/i;
             end
    end
    
    function [F] = losDerivs(los,f)
        % On the fly Backwards Differentiation: 
        % Prediction: 
        % lx = los(1)    Line 1
        % ly = los(2)    Line 2
        % lz = los(3)    Line 3
        %range   = norm(los)
        %theta   = acos(los(3)/range)
        %phi     = atan2(los(2),los(1))
        rangeSq  = los(1)^2 + los(2)^2 + los(3)^2;   % line 4  
        range    = sqrt(rangeSq);                    % line 5
        u        = los(3)/range;                     % line 6
        theta   = acos(u);                           % line 7
        phi     = atan2(los(2),los(1));              % line 8
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        denom   = los(1)^2 + los(2)^2 ;
        % from line 8: phi     = atan2(los(2),los(1))
        f(1)   = f(1) - f(8)*los(2)/denom; 
        f(2)   = f(2) + f(8)*los(1)/denom ;
        % from line 7: theta   = acos(u)
        f(6)   = f(6) - f(7)/sqrt(1 - u^2); 
        %from line 6: u        = los(3)/range
        f(3)   = f(3) + f(6)/range; 
        f(5)   = f(5) - f(6)*u/range;
        % from line 5: range    = sqrt(rangeSq)
        f(4)   = f(4) + 0.5*f(5)/range;
        % from line 4: rangeSq  = los(1)^2 + los(2)^2 + los(3)^2
        f(1)  = f(1) + 2.0*f(4)*los(1); 
        f(2)  = f(2) + 2.0*f(4)*los(2); 
        f(3)  = f(3) + 2.0*f(4)*los(3);
        F = f; 
    end


    function [SVFit] = ClosestPointBetweenTwoLines(track_data,tFit)
        % Method 1: 
        % First do a simple position at closest approach of two line of
        % sight vectors
        % make initial guess using first two sightings -- assume same
        % point is described by both sightings
        r1 = track_data(1,3:5)';
        u1 = track_data(1,6:8)';
        t1 = track_data(1,2);

        r2 = track_data(2,3:5)';
        u2 = track_data(2,6:8)';
        t2 = track_data(2,2);

        r3 = track_data(3,3:5)';
        u3 = track_data(3,6:8)';
        t3 = track_data(3,2);

        CoefMatrix = [u1'*u1,  -u1'*u2;...
            -u2'*u1,  u2'*u2];
        bPoint     = [(r2 - r1)'*u1;...
            (r1 - r2)'*u2];
        tValues    = CoefMatrix\bPoint;
        l1         = r1 + tValues(1)*u1;
        l2         = r2 + tValues(2)*u2;
        v1         = (l2 - l1)/(t2 - t1);

        CoefMatrix = [u2'*u2,  -u2'*u3;...
            -u3'*u2,  u3'*u3];
        bPoint     = [(r3 - r2)'*u2;...
            (r2 - r3)'*u3];
        tValues    = CoefMatrix\bPoint;
        l2a        = r2 + tValues(1)*u2;
        l3         = r3 + tValues(2)*u3;
        v2         = (l3 - l2a)/(t3 - t2);

        SV1        = [ (l2 + l2)/2.0; v1; (t1+t2)/2.0];
        SV2        = [ (l3 + l2)/2.0; v2; (t3+t2)/2.0];
        toTime     = (t1 + t2)/2.0;
        SV2Prop    = propagateECI(SV2, toTime);            
        SVStart    = (SV1 + SV2Prop)/2.0;
        SVFit    = propagateECI(SVStart, tFit);
    end

    function [state_vector, fit] = PolyFitOld(Object, track_data, tFit, Iorder, iType)
        rows       = size(track_data,1);

        % Modfied Polynomial Fit Tools
        MatRow = Iorder*3;
        MatCol = MatRow + 1;
        IncMat = zeros(MatRow, MatRow);
        Minc   = zeros(MatRow, MatRow);
        Matrix = zeros(MatRow, MatCol);

        for ii = 1:rows
            trackID   = track_data(ii,1);
            tinc      = track_data(ii,2) - tFit;
            if Object == 1
                states_used  = track_data(ii,3:5)';   % Satellite
            else
                states_used  = track_data(ii,6:8)';   % Line of Site
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Polynomial Interpolation for the Sensor
            X           = [];   % Incidence/Design Matrix
            Y           = [];
            Column      = eye(3);
            Covariance  = eye(3);
            CL          = states_used;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Polynomials in time
            if iType == 1
                for j      = 1:Iorder
                    X      = [X, Column];
                    Y      = [Y; CL];
                    Column = Column*tinc;
                    CL     = CL*tinc;
                end
            end
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Chebyshev Polynomials of the First Kind
            % if iType == 2
            %     tscaled = tinc/tMag;
            %     T       = [1, tscaled];
            %     for j   = 1:Iorder
            %         if j >= 3
            %             T  = [T, (2*tscaled*T(j-1) - T(j-2))];
            %         end
            %         Column = Column*T(j);
            %         CL     = CL*T(j);
            %         X      = [X, Column];
            %         Y      = [Y; CL];
            %     end
            % end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            IncMat    = X'*Covariance*X;
            Minc      = Minc   + IncMat;
            Matrix    = Matrix + [IncMat, Y];
        end

        MatrixRank  = rank(Minc);
        if MatrixRank < size(Matrix,1)
            Iorder      = Iorder - 1;
            MatrixRank  = rank(Minc);
            disp(' Matrix rank is too low')
        else
            A      = Matrix(:,end);
            fit    = Minc\A;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Reshape into matrix form
            fit     = reshape(fit,3,Iorder);
        end
        [state_vector, ChiSq] = PolyExtrapArray(Object,track_data,fit, tFit);
    end

    function [state_vector, ChiSq] = PolyExtrapArray(Object, track_data, fit, tFit)
    
        rows        = size(track_data,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elements         = zeros(3,rows);
        derivElements    = zeros(3,rows);
        secDerivElements = zeros(3,rows);
        %elements      = zeros(3,1);
        %derivElements = zeros(3,1);
        %secDerivElements = zeros(3,1);
        state_vector_sat = [];
        %tIndex = find(track_data(:,2) == tFit);
        
        ChiSq = [];
        
        for ii = 1:rows
            if Object == 1
                states_used  = track_data(ii,3:5)';   % Satellite
            else
                states_used  = track_data(ii,6:8)';   % Line of Site
            end
            %ii = 1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Tinc     = 0.0;
            Tinc     = track_data(ii,2) - tFit;
            %timeFac  = ones(1,rows);
            %derivFac = ones(1,rows);
            timeFac     = 1;
            derivFac    = 1;
            secDerivFac = 1;
            
            for i = 1:size(fit,2)
                %elements = elements + fit(:,i)*Tinc.^(i-1);
                deriv1st    = i-1;
                deriv2nd    = (i-2)*deriv1st;
                elements(:,ii) = elements(:,ii) + fit(:,i)*timeFac;
                timeFac        = timeFac.*Tinc;
                if i > 1
                    derivElements(:,ii) = derivElements(:,ii) + fit(:,i)*deriv1st*derivFac;
                    derivFac      = derivFac.*Tinc;
                end
                if i > 2
                    secDerivElements(:,ii) = secDerivElements(:,ii) + fit(:,i)*deriv2nd*secDerivFac;
                    secDerivFac      = secDerivFac.*Tinc;
                end
            end
            
            Diff     = elements(:,ii) - states_used(:); 
            ChiSqInc = sum(Diff.*Diff);
            ChiSq    = [ChiSq; ChiSqInc];
        end
        %ChiSq
        %%%%  Satellite 9-State Vector at the median time
        %state_vector_sat = [elements/DU; derivElements/VU; secDerivElements/AU]
        state_vector = [elements; derivElements; secDerivElements];
    end

    function [state_vector, fit, Iorder, ChiSq] = PolyFit(Object, track_data, tFit, Iorder, iType)
        ChiSq  = [];
        [Matrix, Minc]   =  BuildMatrix(Object, track_data, tFit, Iorder, iType);
        MatrixRank  = rank(Minc);
        while MatrixRank < size(Matrix,1)
            Iorder      = Iorder - 1;
            [Matrix, Minc]   =  BuildMatrix(Object, track_data, tFit, Iorder, iType);
            MatrixRank  = rank(Minc);
        end
        if MatrixRank < size(Matrix,1)
            disp(' Polynomial Fit Fails: Matrix rank is too low')
        else
            A      = Matrix(:,end);
            fit    = Minc\A;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Reshape into matrix form
            fit     = reshape(fit,3,Iorder);
        end

        [state_vector, ChiSq] = PolyExtrapArray(Object,track_data,fit, tFit);
        
    end

    function [Matrix, Minc] = BuildMatrix(Object, track_data, tFit, Iorder, iType)
        rows       = size(track_data,1);

        % Modfied Polynomial Fit Tools
        MatRow = Iorder*3;
        MatCol = MatRow + 1;
        IncMat = zeros(MatRow, MatRow);
        Minc   = zeros(MatRow, MatRow);
        Matrix = zeros(MatRow, MatCol);

        for ii = 1:rows
            trackID   = track_data(ii,1);
            tinc      = track_data(ii,2) - tFit;
            if Object == 1
                states_used  = track_data(ii,3:5)';   % Satellite
            else
                states_used  = track_data(ii,6:8)';   % Line of Site
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Polynomial Interpolation for the Sensor
            X           = [];   % Incidence/Design Matrix
            Y           = [];
            Column      = eye(3);
            Covariance  = eye(3);
            CL          = states_used;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Polynomials in time
            if iType == 1
                for j      = 1:Iorder
                    X      = [X, Column];
                    Y      = [Y; CL];
                    Column = Column*tinc;
                    CL     = CL*tinc;
                end
            end
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Chebyshev Polynomials of the First Kind
            % if iType == 2
            %     tscaled = tinc/tMag;
            %     T       = [1, tscaled];
            %     for j   = 1:Iorder
            %         if j >= 3
            %             T  = [T, (2*tscaled*T(j-1) - T(j-2))];
            %         end
            %         Column = Column*T(j);
            %         CL     = CL*T(j);
            %         X      = [X, Column];
            %         Y      = [Y; CL];
            %     end
            % end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            IncMat    = X'*Covariance*X;
            Minc      = Minc   + IncMat;
            Matrix    = Matrix + [IncMat, Y];
        end

        MatrixRank  = rank(Minc);

    end

    function [state_vector] = PolyExtrapToFitTime(fit, tFit)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %elements         = zeros(3,rows);
        %derivElements    = zeros(3,rows);
        %secDerivElements = zeros(3,rows);
        elements      = zeros(3,1);
        derivElements = zeros(3,1);
        secDerivElements = zeros(3,1);
        state_vector_sat = [];
        %tIndex = find(track_data(:,2) == tFit);
        %for ii = 1:rows
            ii = 1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Tinc     = 0.0;
            %Tinc     = track_data(ii,2) - tFit;
            %timeFac  = ones(1,rows);
            %derivFac = ones(1,rows);
            timeFac     = 1;
            derivFac    = 1;
            secDerivFac = 1;
            for i = 1:size(fit,2)
                %elements = elements + fit(:,i)*Tinc.^(i-1);
                deriv1st    = i-1;
                deriv2nd    = (i-2)*deriv1st;
                elements(:,ii) = elements(:,ii) + fit(:,i)*timeFac;
                timeFac        = timeFac.*Tinc;
                if i > 1
                    derivElements(:,ii) = derivElements(:,ii) + fit(:,i)*deriv1st*derivFac;
                    derivFac      = derivFac.*Tinc;
                end
                if i > 2
                    secDerivElements(:,ii) = secDerivElements(:,ii) + fit(:,i)*deriv2nd*secDerivFac;
                    secDerivFac      = secDerivFac.*Tinc;
                end
            end
        %end
        %%%%  Satellite 9-State Vector at the median time
        %state_vector_sat = [elements/DU; derivElements/VU; secDerivElements/AU]
        state_vector = [elements; derivElements; secDerivElements];
    end

    function [StateVectorECI] = LaplaceInECI(DU,mu, state_vector_sat, state_vector_los)

        % Method 2 -- Laplace's Angles-Only Algorithm
        Rpos = state_vector_sat(1:3);
        Vpos = state_vector_sat(4:6);
        Apos = state_vector_sat(7:9);

        Rlos = state_vector_los(1:3);
        Vlos = state_vector_los(4:6);
        Alos = state_vector_los(7:9);

        D  = det([Rlos, 2*Vlos, Alos]);
        D1 = det([Rlos,   Vlos, Apos]);
        D2 = det([Rlos,   Vlos, Rpos]);
        D3 = det([Rlos,   Apos, Alos]);
        D4 = det([Rlos,   Rpos, Alos]);

        %delr = 0.00001
        %r    = 1.155;
        delr = 0.01*DU;
        r    = 1.0*DU;
        for i = 1:1000
            r = r + delr;
            %Col3 = Alos + Rlos/r^3;
            %D    = det([Rlos, 2*Vlos, Col3]);
            rho  = -2.0*(D1 + mu*D2/r^3)/D;
            %F    = rho^2 + (2.0*rho*state_vector_los(1:3)' + state_vector_sat(1:3)')*state_vector_sat(1:3) - r^2
            F    = rho^2 - r^2 + 2.0*rho*dot(Rlos,Rpos) + dot(Rpos,Rpos);
            if F < 0
                r = r - delr;
                delr = delr/2;
            end
        end
        rVec   = rho*Rlos + Rpos;
        rMag   = norm(rVec);
        rhodot = -(D3 + mu*D4/r^3)/D;
        vVec   = rhodot*Rlos + rho*Vlos + Vpos;
        %aVec   = -mu*rVec/rMag^3;
        [aVec Phi] = Escobal(rVec);
        % Estimated Initial 9-State Vector from Laplace Method in ECI
        StateVectorECI  = [rVec; vVec; aVec]
    end

    function [SVDiffCorrection] =  DifferentialCorrectionToLaplace(tFit, track_data, StateVectorECI)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Try to do a numerical differential correction:
        rows      = size(track_data,1);

        epsilon   = 0.01;
        col(:,1)  = epsilon*[1,0,0,0,0,0,0]';
        col(:,2)  = epsilon*[0,1,0,0,0,0,0]';
        col(:,3)  = epsilon*[0,0,1,0,0,0,0]';
        col(:,4)  = epsilon*[0,0,0,1,0,0,0]';
        col(:,5)  = epsilon*[0,0,0,0,1,0,0]';
        col(:,6)  = epsilon*[0,0,0,0,0,1,0]';
        SV7  = [StateVectorECI(1:3); StateVectorECI(4:6); tFit];
        % Want to apply Differential Corrections to
        %  LOS_MEASURED = LOS_PREDICTED|i + (partial_LOS_Predicted/partial_pk)|i * Delta
        % LOS Measured  = los_state
        % LOS PREDICTED = target_predicted - sensor_location
        for i = 1:100
            X = zeros(6,6);
            Y = zeros(6,1);
            for ii = 1:rows
                toTime          = track_data(ii,2);
                % Get the Measured Target LOS in ECI
                sat_state       = track_data(ii,3:5)';
                los_state       = track_data(ii,6:8)';   % Line of Site
                % Get the Measured Sensor Position in ECI
                % Propagate the target position to current time
                SV7Prop         = propagateECI(SV7, toTime);
                Rtarget         = SV7Prop(1:3);
                RLOS            = Rtarget - sat_state;
                Runit           = RLOS/norm(RLOS);   % LOS_PREDICTED
                Jacobian  = [];
                for irow = 1:6
                    SVTemp           = SV7 + col(:,irow);
                    SVTemp           = propagateECI(SVTemp, toTime);
                    RTemp            = SVTemp(1:3);
                    RTemp            = RTemp - sat_state;
                    UTemp            = RTemp/norm(RTemp);
                    delCol           = (UTemp - Runit)/epsilon;
                    Jacobian        = [Jacobian, [delCol(1:3)]];
                end
                Jacobian  = Jacobian(1:3,:);
                % Compute the Propagated LOS vector
                % Measured LOS Unit Vector
                ECI_LOS_difference = los_state - Runit;
                X = X + Jacobian'*Jacobian;
                Y = Y + Jacobian'*ECI_LOS_difference;
                %
            end
            diff = [X\Y;0];
            %norm(diff);
            SV7      = SV7 + diff;
            %SV7Prop  = propagateECI(SV7, tFit)
            %M = ECR_to_ECI_PVA_Rotation_Matrix( -tFit );
            %StateVectorECEF = M(1:3,1:3)*SV7Prop(1:3)/1000
            % Divide by 1000 to work in km, km/s, km/s^2 units
            %M*SV7(1:6)/1000.0
        end

        epsilon   = 0.001;
        col(:,1)  = epsilon*[1,0,0,0,0,0,0]';
        col(:,2)  = epsilon*[0,1,0,0,0,0,0]';
        col(:,3)  = epsilon*[0,0,1,0,0,0,0]';
        col(:,4)  = epsilon*[0,0,0,1,0,0,0]';
        col(:,5)  = epsilon*[0,0,0,0,1,0,0]';
        col(:,6)  = epsilon*[0,0,0,0,0,1,0]';
        SV7  = [StateVectorECI(1:3); StateVectorECI(4:6); tFit];
        % Want to apply Differential Corrections to
        %  LOS_MEASURED = LOS_PREDICTED|i + (partial_LOS_Predicted/partial_pk)|i * Delta
        % LOS Measured  = los_state
        % LOS PREDICTED = target_predicted - sensor_location
        for i = 1:100
            X = zeros(6,6);
            Y = zeros(6,1);
            for ii = 1:rows
                toTime          = track_data(ii,2);
                % Get the Measured Target LOS in ECI
                sat_state       = track_data(ii,3:5)';
                los_state       = track_data(ii,6:8)';   % Line of Site
                % Get the Measured Sensor Position in ECI
                % Propagate the target position to current time
                SV7Prop         = propagateECI(SV7, toTime);
                Rtarget         = SV7Prop(1:3);
                RLOS            = Rtarget - sat_state;
                Runit           = RLOS/norm(RLOS);   % LOS_PREDICTED
                Jacobian  = [];
                for irow = 1:6
                    SVTemp           = SV7 + col(:,irow);
                    SVTemp           = propagateECI(SVTemp, toTime);
                    RTemp            = SVTemp(1:3);
                    RTemp            = RTemp - sat_state;
                    UTemp            = RTemp/norm(RTemp);
                    delCol           = (UTemp - Runit)/epsilon;
                    Jacobian        = [Jacobian, [delCol(1:3)]];
                end
                Jacobian  = Jacobian(1:3,:);
                % Compute the Propagated LOS vector
                % Measured LOS Unit Vector
                ECI_LOS_difference = los_state - Runit;
                X = X + Jacobian'*Jacobian;
                Y = Y + Jacobian'*ECI_LOS_difference;
                %
            end
            diff = [X\Y;0];
            %norm(diff);
            SV7      = SV7 + diff;
            %SV7Prop  = propagateECI(SV7, tFit)
            %M = ECR_to_ECI_PVA_Rotation_Matrix( -tFit );
            %StateVectorECEF = M(1:3,1:3)*SV7Prop(1:3)/1000
            % Divide by 1000 to work in  km, km/s, km/s^2 units
            %M*SV7(1:6)/1000.0
        end

        % epsilon   = 0.0001;
        % col(:,1)  = epsilon*[1,0,0,0,0,0,0]';
        % col(:,2)  = epsilon*[0,1,0,0,0,0,0]';
        % col(:,3)  = epsilon*[0,0,1,0,0,0,0]';
        % col(:,4)  = epsilon*[0,0,0,1,0,0,0]';
        % col(:,5)  = epsilon*[0,0,0,0,1,0,0]';
        % col(:,6)  = epsilon*[0,0,0,0,0,1,0]';
        % SV7  = [StateVectorECI(1:3); StateVectorECI(4:6); tFit];
        % % Want to apply Differential Corrections to
        % %  LOS_MEASURED = LOS_PREDICTED|i + (partial_LOS_Predicted/partial_pk)|i * Delta
        % % LOS Measured  = los_state
        % % LOS PREDICTED = target_predicted - sensor_location
        % for i = 1:100
        %     X = zeros(6,6);
        %     Y = zeros(6,1);
        %     for ii = 1:rows
        %         toTime          = track_data(ii,2);
        %         % Get the Measured Target LOS in ECI
        %         sat_state       = track_data(ii,3:5)';
        %         los_state       = track_data(ii,6:8)';   % Line of Site
        %         % Get the Measured Sensor Position in ECI
        %         % Propagate the target position to current time
        %         SV7Prop         = propagateECI(SV7, toTime);
        %         Rtarget         = SV7Prop(1:3);
        %         RLOS            = Rtarget - sat_state;
        %         Runit           = RLOS/norm(RLOS);   % LOS_PREDICTED
        %         Jacobian  = [];
        %         for irow = 1:6
        %             SVTemp           = SV7 + col(:,irow);
        %             SVTemp           = propagateECI(SVTemp, toTime);
        %             RTemp            = SVTemp(1:3);
        %             RTemp            = RTemp - sat_state;
        %             UTemp            = RTemp/norm(RTemp);
        %             delCol           = (UTemp - Runit)/epsilon;
        %             Jacobian        = [Jacobian, [delCol(1:3)]];
        %         end
        %         Jacobian  = Jacobian(1:3,:);
        %         % Compute the Propagated LOS vector
        %         % Measured LOS Unit Vector
        %         ECI_LOS_difference = los_state - Runit;
        %         X = X + Jacobian'*Jacobian;
        %         Y = Y + Jacobian'*ECI_LOS_difference;
        %         %
        %     end
        %     diff = [X\Y;0];
        %     %norm(diff);
        %     SV7      = SV7 + diff;
        %     %SV7Prop  = propagateECI(SV7, tFit)
        %     %M = ECR_to_ECI_PVA_Rotation_Matrix( -tFit );
        %     %StateVectorECEF = M(1:3,1:3)*SV7Prop(1:3)/1000
        %     % Divide by 1000 to work in  km, km/s, km/s^2 units
        %     %M*SV7(1:6)/1000.0
        % end


        % Report Findings at the End of the Propagation
        %M = ECR_to_ECI_PVA_Rotation_Matrix( -tFit );
        %StateVectorECEF = M(1:6,1:6)*SV7Prop(1:6)/1000
        SVDiffCorrection = SV7Prop;
    end

    function [KeplerDot] = FirstOrder(t,K)             
        Kepler     = LagrangePlanetary(K);
        %J6Check    = Kepler.J6Gravity(time, EscobalGravity*TU^2/DU);
        J6Check    = Kepler.J6Gravity(t);
        KeplerDot  = Kepler.FirstOrder(t);
    end
