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
        Iorder = 6
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
        [Rinit,Vinit,ParamList,Jacobian,Hessian,GravityCan] = KeplerSat.OrbitDerivatives(time, Sensor, Observed);
        
        deltaPar     = zeros(6,1);
        ChiSqArray   = [];
        TangentArray = [];
        InvCovR = eye(3);
        deltaPar   = zeros(6,1);
        ChiSqArray = [];
        InvCovR = eye(3);
        
        Mind  = 6
        oparm = zeros(Mind,1);
        
        for NR = 1:20
            
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
                [Rinit,Vinit,ParamList,Jacobian,Hessian,GravityCan] = KeplerSat.OrbitDerivatives(time, Sensor, Observed);
                deltaR   = Rinit - Rdata;
                rSatKepler = [rSatKepler, Rinit];
                %deltaR   = Rinit - Rdata + Jacobian(1:3,2:7)'*delZ;
                ChiSq    = ChiSq  + 0.5*deltaR'*InvCovR*deltaR;
                betaNR   = betaNR + Jacobian(1:3,2:7)'*InvCovR*deltaR;
                hessPart = zeros(6);
                for k = 1:3
                    hessPart = hessPart + Hessian{k}(2:7,2:7)*deltaR(k);
                end
                %hessNR  = hessNR + Jacobian(1:3,2:7)'*InvCovR*Jacobian(1:3,2:7);
                hessNR  = hessNR + hessPart + Jacobian(1:3,2:7)'*InvCovR*Jacobian(1:3,2:7);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            ChiSqArray  =  [ChiSqArray, log(ChiSq)];
            %deltaPar    = -hessNR\betaNR
            
            AugmentedMatrix = [hessNR,-betaNR];
            AugmentedMatrix = AugmentedMatrix(IndexRows, IndexCols);
            hitzero   = 0;
            opzero = 0.0000000000001;
            disp("Entering NRMin ");
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
            
            TangentArray = [TangentArray, log(norm(betaNR))]
            %Kepler      =  Kepler + deltaPar'
            Kepler      =  Kepler + delp
        end
        KeplerSat  = LagrangePlanetary(Kepler);
        
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
            Iorder     = 3;
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
            StateVectorCell   = LaplaceInECI(state_vector_sat, state_vector_los);
            ChiSqMin = 10^100
            for icell = 1:size(StateVectorCell,2)
                StateVector = StateVectorCell{icell};
                time        = tFit/TU
                Rpos        = StateVector(1:3)/DU;
                Rdot        = StateVector(4:6)/VU;
                ChiSq = 0.0
                % comment out 506 to see what happens. !!!!!!!!!!!!!!!!!!!!!!!!!!!!
                [e a Inclination omega Omega Mp Period Jacobian] = kepler(time, Rpos, Rdot);
                KeplerTest = [e a Inclination omega Omega Mp] % Canonical Units
                
                Imaginary  = numel(find(imag(KeplerTest)))
                
                if Imaginary == 0
                    KeplerObj  = LagrangePlanetary(KeplerTest);
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
                        [Rsensor, Vsensor, ParamList, Jacobian, Hessian, GravityCan] = KeplerSat.OrbitDerivatives(tRecord, Sensor, Observed);
                        Sensor   = [Rsensor; Vsensor; zeros(3,1)];
                        [Rinit,Vinit,ParamList,Jacobian,Hessian,GravityCan] = KeplerObj.OrbitDerivatives(tRecord, Sensor, Observed);
                        los       = [ParamList(7); ParamList(8); ParamList(9)];
                        los       = los/norm(los);
                        AngleDiff = 1 - dot(los,losM)
                        ChiSq     = ChiSq + AngleDiff;
                    end
                    if ChiSq < ChiSqMin
                        ChiSqMin = ChiSq
                        StateVectorECI   = StateVector;
                    end
                end
            end
            
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
        %[EscobalGravity Phi] = Escobal(Rpos, mu);
        %KeplerInit = [e a Inclination omega Omega Mp] % Canonical Units
        %lna       = log(a);
        meanMotion = exp(-1.5*log(a));
        KLagrange  = LagrangePlanetary(KeplerInit);
        %J6Check   = KLagrange.OrbitDerivatives(time, J6);
        %J6Check   = KLagrange.OrbitDerivatives(time, EscobalGravity*TU^2/DU);
        %[J6Check] = KLagrange.OrbitDerivatives(time);
        %[Rextrap, Vextrap, Hessian, Gradient, Jacobian] = KLagrange.OrbitDerivatives(time)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %J6Check      = KLagrange.OrbitDerivatives(tTest, EscobalGravity*TU^2/DU)
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
        % Check on J6 Gravity Perturbations:
        %[Rinit,Vinit,Hessian,Gradient,Jacobian, ~, ~, ~,~,~,~] = KLagrange.OrbitDerivatives(time, Sensor);
        [Rinit0,Vinit0,ParamList0,Jacobian0Hessian0,Gradient0,GravityCan0] = KLagrange.OrbitDerivatives(time, Sensor, Observed);
        
        [Zdot0, Inhomogeneous0, Perturbation0, DMatrix0, LagrangeTerms0, LagrangeGrads0] = KLagrange.EquationsOfMotion(Zhat);
        SVstn0 = [Rinit0;Vinit0;time];
                [EscobalGravity Phi] = Escobal(Rpos, mu);
        EscobalGravity'*TU^2/DU
        GravityCan0'
        [J6 Jperturb Jhom] = J6Gravity(Rpos);
        Jperturb
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %      Begin Finite Differences Check
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for LineCheck = 1:130
        %     [Jacob, JacobFinite, Hess, HessFinite] = DigitalJustice(KLagrange, LineCheck, Sensor, Observed);
        %     Jacob
        %     JacobFinite
        %     Hess
        %     HessFinite
        %     Sum = sum(sum(Hess - HessFinite));
        %     disp([' Line = ', num2str(LineCheck),'  sum = ', num2str(Sum)])
        % end
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
%                 %[Rsensor, Vsensor, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessSat] = KeplerSat.OrbitDerivatives(tRecord, Sensor, Observed);
%                 %Sensor   = [Rsensor; Vsensor; zeros(3,1)];
% 
%                 [Rinit,Vinit,Hessian,Gradient,Jacobian,JacobianLos,ParamList,...
%                     GravityCan,HessianLos,Phi,MLagrange,RVHessian] = KeplerObj.OrbitDerivatives(tRecord, Sensor, Observed);
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
        deltaPar   = zeros(6,1);
        delp       = zeros(6,1);
        
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
        
        
        Mind  = 6
        oparm = zeros(Mind,1);
        
        for NR = 1:30
            
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
                [Rsensor, Vsensor, ParamList, Jacobian, Hessian, GravityCan] = KeplerSat.OrbitDerivatives(tRecord, Sensor, Observed);
                Sensor   = [Rsensor; Vsensor; zeros(3,1)];

                [Rinit,Vinit,ParamList,Jacobian,Hessian,GravityCan] = KeplerObj.OrbitDerivatives(tRecord, Sensor, Observed);
                los       = [ParamList(7); ParamList(8); ParamList(9)];
                los       = los/norm(los);

                deltaM    =  1 - dot(los,losM);
                JacobPred =  Jacobian(14:16,2:7);
                JacFac    = -JacobPred'*losM;
                BetaInc   =  BetaInc + JacFac;
                
                %Jac2nd    = JacFac*JacFac';
                Hess      = zeros(6);
                for ik = 1:3
                    Hess      = Hess -losM(ik)*Hessian{13+ik}(2:7,2:7);
                end
                HessianSyst = HessianSyst + Hess;
                %Hess = Hess*delta;
                %HessianSyst = HessianSyst + Hess + Jac2nd;
                
                ChiSq    = ChiSq + deltaM;
                                
            end
            ChiSqArray   = [ChiSqArray, log(ChiSq)];
            delta        = -(HessianSyst)\BetaInc
            %deltaPar     = -hessNR\betaNR
            %TangentArray = [TangentArray, log(norm(betaNR))]
            %Kepler       =  Kepler + deltaPar'
            %Kepler       =  Kepler + delta';
                                    
            AugmentedMatrix = [HessianSyst,-BetaInc];
            AugmentedMatrix = AugmentedMatrix(IndexRows, IndexCols);
            hitzero   = 0;
            opzero = 0.0000000000001;
            disp("Entering NRMin ");
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
            delp
            Kepler       =  Kepler + delp';

            
            TangentArray = [TangentArray, log(norm(BetaInc))]
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
%                 %[Rsensor, Vsensor, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessSat] = KeplerSat.OrbitDerivatives(tRecord, Sensor);
%                 %Sensor   = [Rsensor; Vsensor; zeros(3,1)];
%                 
%                 [Rinit,Vinit,Hessian,Gradient,Jacobian,JacobianLos,ParamList,...
%                     GravityCan,HessianLos,Phi,MLagrange,RVHessian] = KeplerObj.OrbitDerivatives(tRecord, Sensor, Observed);
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
        for NR = 1:800
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
                %[Rinit,Vinit,Hessian,Gradient,Jacobian, ~, ~, ~,~,~,~] = KLagrange.OrbitDerivatives(time, Sensor);
                [Rinit,Vinit,Hessian,Gradient,Jacobian,JacobianLos,ParamList,GravityCan,HessianLos,Phi,MLagrange]...
                    = KLagrange.OrbitDerivatives(time, Sensor, Observed);
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
                %[Rsensor, Vsensor, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessSat] = KeplerSat.OrbitDerivatives(tRecord, Sensor, Observed);
                %Sensor = [Rsensor; Vsensor; zeros(3,1)];
                %[Rextrap, Vextrap, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan] = KLagrangeProp.OrbitDerivatives(tRecord, Sensor);
                [Rextrap, Vextrap, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessLos] = KLagrange.OrbitDerivatives(tRecord, Sensor, Observed);
                %[Rextrap, Vextrap, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessLos] = KLagrangeProp.OrbitDerivatives(tRecord, Sensor, Observed);
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
        KeplerIteration  = KeplerFit
        %KLagrange        = LagrangePlanetary(KeplerIteration);
        %KeplerIteration = KeplerInit + Zhat';   % Reset the Gravity Model Reference Point
        
        %KeplerIteration  = Kepler + Zhat;
        %Zhat           = zeros(6,1);
        KLagrange       = LagrangePlanetary(KeplerIteration);
        [Rinit,Vinit,Hessian,Gradient,Jacobian,JacobianLos,ParamList,GravityCan,HessianLos,Phi,MLagrange]...
            = KLagrange.OrbitDerivatives(time, Sensor, Observed);
        [Rpos, Rdot]
        [(Rpos-Rinit), (Rdot-Vinit)]
        [Zdot, Inhomogeneous, Perturbation, DMatrix, LagrangeTerms, LagrangeGrads]= KLagrange.EquationsOfMotion(zeros(6,1));
        
        for NR = 1:100

            ChiSq       = 0.0;
            InvCov      = zeros(6);
            BetaInc     = zeros(6,1);
            HessianSyst = zeros(6);

            for ii = 1:rows
                [Zdot, Inhomogeneous, Perturbation, DMatrix, LagrangeTerms, LagrangeGrads]= KLagrange.EquationsOfMotion(Zhat);
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
                % [Rextrap, Vextrap, Hessian, Gradient, Jacobian] = KLagrange.OrbitDerivatives(0)
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
                Sensor = [track_data(ii,3:5)'/DU; zeros(6,1)];
                % Sat Position from Orbit Fit to Kepler elements: 
                %[Rsensor, Vsensor, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessSat] = KeplerSat.OrbitDerivatives(tRecord, Sensor, Observed);
                %Sensor = [Rsensor; Vsensor; zeros(3,1)];
                %[Rextrap, Vextrap, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan] = KLagrangeProp.OrbitDerivatives(tRecord, Sensor);
                [Rextrap, Vextrap, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessLos] = KLagrangeProp.OrbitDerivatives(tRecord, Sensor, Observed);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Alternative Method of correcting for gravitational perturbations:
                %Delta = Jacobian(:,2:7)*delZ;
                %Rextrap + Delta(1:3);
                %Vextrap + Delta(4:6);
                %los     = Rextrap - state_vector_satArray(1:3,ii)/DU;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % CombJac   = JacobianLos(10:11,2:7)*PhiFundamental;
                % % %CombJac   = JacobianLos0(2:3,2:7)*PhiFundamental;
                % InvCov     = InvCov + CombJac'*InvCovMeas*CombJac;
                % % line of sight components: JRS
                % Mdelta   = [ cos(theta)*cos(phi) cos(theta)*sin(phi)  -sin(theta);...
                %             -sin(phi)  cos(phi)      0      ];
                % deltaA   = Mdelta*[los-losM]
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lo
                % ChiSqInc   = 0.5*deltaA'*InvCovMeas*deltaA;
                % ChiSq      = ChiSq   + ChiSqInc;
                %
                % CovVec     = InvCovMeas*deltaA;
                % dataVec    = CombJac'*CovVec;
                % %BetaInc    = BetaInc + dataVec;
                % %SymCheck   = [deltaA'*InvCovMeas*CombJac]'
                % BetaInc    = BetaInc + dataVec;
                % Matrix     = zeros(6);
                % Matrix = Matrix + PhiFundamental'*HessianLos{10}(2:7,2:7)*PhiFundamental*CovVec(1);
                % Matrix = Matrix + PhiFundamental'*HessianLos{11}(2:7,2:7)*PhiFundamental*CovVec(2);
                % % for icol   = 1:2
                % %     %Matrix = Matrix + HessianLos{icol}(2:7,2:7)*PhiFundamental*CovVec(icol)*PhiFundamental;
                % %     Matrix = Matrix + PhiFundamental'*HessianLos{icol}(2:7,2:7)*PhiFundamental*CovVec(icol);
                % %     %Matrix = Matrix + PhiFundamental'*HessianLos0{icol}(2:7,2:7)*PhiFundamental*CovVec(icol);
                % % end
                % %HessianSyst = HessianSyst + PhiFundamental*PhiFundamental*Matrix;
                % %HessianSyst = HessianSyst + PhiFundamental*Matrix*PhiFundamental;
                % HessianSyst = HessianSyst + Matrix;
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
                CombJac  = JacobianLos(7:9,2:7)*PhiFundamental;
                deltaA   = los-losM;
                dataVec  = CombJac'*InvCovM*deltaA;
                CovVec    = InvCovM*deltaA;
                %Jac2nd    = CombJac'*InvCovM*CombJac;
                Matrix    = zeros(6);
                Matrix = Matrix + PhiFundamental'*HessianLos{7}(2:7,2:7)*PhiFundamental*CovVec(1);
                Matrix = Matrix + PhiFundamental'*HessianLos{8}(2:7,2:7)*PhiFundamental*CovVec(2);
                Matrix = Matrix + PhiFundamental'*HessianLos{9}(2:7,2:7)*PhiFundamental*CovVec(3);
                HessianSyst = HessianSyst + Matrix;
                InvCov      = InvCov      + CombJac'*InvCovM*CombJac;               
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lo
                ChiSqInc   = 0.5*deltaA'*InvCovM*deltaA;
                ChiSq      = ChiSq   + ChiSqInc;
                BetaInc    = BetaInc + dataVec;         
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
        [Rpos, Rdot, Hessian, Gradient, Jacobian, ~, ~] = KLagrange.OrbitDerivatives(time, Sensor, Observed);
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