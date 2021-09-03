    function [KeplerFit, CovFit, TFit, SV7Fit] = LeastSquaresAllSatellites(track_data, TFit, Isat, iFirst, InvCovUV, units, InitialState)
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        % Initial State         
        Iorder        = InitialState.Iorder;
        %TEnd          = InitialState.TEnd;
        %TInitialState = InitialState.TInitialState;
        SVInitial     = InitialState.SV10Initial;
        TPolyFit      = InitialState.TPolyFit;
        fitInitial    = InitialState.fit;
        CovInitial    = InitialState.CovFit;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        rows       = size(track_data,1);
        iFit       = floor(numel(track_data(:,2))/2);
        %TFit       = track_data(iFit,2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % For Choice of Units
        TU = units.TU;
        DU = units.DU;
        VU = units.VU;
        AU = units.AU;
        mu = units.mu;
        sqrtmu = units.sqrtmu;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        losM     = track_data(iFit,6:8)';
        losM     = losM/norm(losM);   % los
        rangeM   = norm(losM);
        thetaM   = acos(losM(3)/rangeM);
        phiM     = atan2(losM(2), losM(1));
        Observed = [thetaM; phiM];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % %Sat Position at measured position
        Rsensor  = [track_data(iFit,3:5)']/DU;
        Sensor   = [Rsensor; zeros(6,1)];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %tFit = TFit/TU;
        tFit = TFit/TU;
        timeArray = track_data(:,2)/TU;
        tBeg = min(timeArray);
        TEnd = track_data(end,2); 
        tEnd = TEnd/TU;
        
        tSVInit  = SVInitial(10)/TU;
        Rpos     = SVInitial(1:3)/DU;
        Rdot     = SVInitial(4:6)/VU;
        [e a Inclination omega Omega Mp Period Jacobian] = kepler(tSVInit, Rpos, Rdot, mu);
        KeplerInit    = [e a Inclination omega Omega Mp] % Canonical Units
        KeplerIngress = LagrangePlanetary(KeplerInit, units);
        [Rinit,Vinit,ParamList,Jacobian,Hessian,GravityCan] = KeplerIngress.OrbitDerivatives(tSVInit, Sensor, losM, InvCovUV);
        [Rinit, Rpos, Rinit-Rpos]
        [Vinit, Rdot, Vinit-Rdot]
        KeplerFit = KeplerInit;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % 2A)     Attempt to incorporate BFGS-type quasi-Newton method  
        %[KeplerFit, ChiSq, CovFit] = CollinearityBFGS(track_data, KeplerInit, InvCovUV, units)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2B)     Begin Application of Hessian and Second Order fit    
        %[KeplerFit, ChiSq, CovFit] = Collinearity(track_data, KeplerFit, InvCovUV, units)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 4)     Begin Application of Hessian and Second Order fit    
        %[KeplerFit, ChiSq, CovFit] = NewtonStepUV(track_data, KeplerFit, InvCovUV, units)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 5)     Begin Application of Hessian and Second Order fit 
        [KeplerFit, ChiSq, CovFit, SV7Final, CovFinal] = NewtonStepUVWithPerturbations(track_data, KeplerFit, tFit, InvCovUV, units)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % % Try to Integrate the Equations of Motion for the Classical Elements
        % Y             = zeros(6,1);
        % KLagrange = LagrangePlanetary(KeplerFit, mu);
        % options       = odeset('RelTol', 1e-10, 'AbsTol', 1e-13);
        % dJdtCartesian = perturbations_odefun(KLagrange,tFit,Y);
        % Y            = zeros(6,1);
        % % Propagate from fit Reference time tFit back to tBeg
        % [Tout, Yout] = ode45(@(t,y) perturbations_odefun(KLagrange, t, y), [tFit, tBeg], Y, options);
        % Yout(end,:)
        % Perturbed = [Tout, Yout];
        % Perturbed = sortrows(Perturbed,1)
        % Y = Perturbed(end,2:7)
        % % Propagate from fit Reference time tFit back to tEnd
        % [Tout, Yout] = ode45(@(t,y) perturbations_odefun(KLagrange, t, y), [tFit, tEnd], Y, options);
        % Yout(end,:)
        % Perturbed(end,:) = [];
        % Perturbed = [Perturbed; Tout, Yout];
        % Perturbed = sortrows(Perturbed,1)
        % KeplerInterpolated = interp1(Perturbed(:,1),Perturbed(:,2:7),timeArray);
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % 3)     Begin Application of Hessian and Second Order fit
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % deltaPar       = zeros(6,1);
        %
        % TangentArray   = [];
        % ZhatArray      = [];
        % DelZhatArray   = [];
        %
        % ChiSqArray     = [];
        % DelChiSqArray  = [];
        %
        % DChiSqArray    = [];
        % DiffChiSqArray = [];
        %
        % PhiFundamental = eye(6);
        % %InvCovM     = eye(3);
        % InvCovMeas   = eye(2);
        % InvCovM      = eye(2);
        % deltaPar     = zeros(6,1);
        % epsilon      = 0.000000001*ones(6,1);
        % Kepler       = KeplerInit
        % %Kepler       = KeplerInitLC
        % for NR = 1:100
        %
        %     ChiSq        = 0.0;
        %     InvCov       = zeros(6);
        %     HessianSyst  = zeros(6);
        %     hessNR       = zeros(6);
        %     hessianChiSq = zeros(6);
        %
        %     BetaInc      = zeros(6,1);
        %     betaNR       = zeros(6,1);
        %     derivChiSq   = zeros(6,1);
        %
        %     KeplerObj    = LagrangePlanetary(Kepler, units);
        %     for ii =1:rows
        %         tRecord = track_data(ii,2)/TU;
        %
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % Measurement Angles:
        %         %losM    = state_vector_losArray(1:3,ii);
        %         losM     = track_data(ii,6:8)';
        %         losM     = losM/norm(losM);   % los
        %         rangeM   = norm(losM);
        %         thetaM   = acos(losM(3)/rangeM);
        %         phiM     = atan2(losM(2), losM(1));
        %         Observed =[thetaM; phiM];
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %         % %Sat Position from Polynomial at tFit
        %         % Sensor = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU]
        %
        %         % %Sat Position from Polynomial at t_ii
        %         % Sensor = [state_vector_satArray(1:3,ii)/DU; state_vector_satArray(4:6,ii)/VU; state_vector_satArray(7:9,ii)/AU]
        %
        %         % %Sat Position at measured position
        %         Rdata  = [track_data(ii,3:5)']/DU;
        %         Sensor = [Rdata; zeros(6,1)];
        %         %Sensor = [Rdata+epsilon(1:3); zeros(6,1)];
        %
        %         % %Sat Position from Orbit Fit to Kepler elements:
        %         %[Rsensor, Vsensor, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessSat] = KeplerSat.OrbitDerivatives(tRecord, Sensor);
        %         %Sensor   = [Rsensor; Vsensor; zeros(3,1)];
        %
        %         [Rinit,Vinit,Hessian,Gradient,Jacobian,JacobianLos,ParamList,...
        %             GravityCan,HessianLos,Phi,MLagrange,RVHessian] = KeplerObj.OrbitDerivatives(tRecord, Sensor, Observed);
        %         los       = [ParamList(7); ParamList(8); ParamList(9)];
        %         los       = los/norm(los);
        %         theta     = ParamList(2);
        %         phi       = ParamList(3);
        %         Predicted = [theta; phi];
        %         % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         %[Jacob, JacobFinite, Hess, HessFinite] = KeplerObj.DigitalJustice(108,Sensor)
        %         %JacobianLos(2,:) = Jacob
        %         %HessianLos{2}    = HessFinite
        %         %[Jacob, JacobFinite, Hess, HessFinite] = KeplerObj.DigitalJustice(109,Sensor)
        %         %JacobianLos(3,:) = Jacob
        %         %HessianLos{3}    = HessFinite
        %         % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % JacobLos = [];
        %         % %
        %         % F       = zeros(8,1);
        %         % F(7)    = 1;     % Theta
        %         % f       = losDerivs(los,F);
        %         % JacobLos = [JacobLos; f(1), f(2), f(3)];
        %         % %
        %         % F       = zeros(8,1);
        %         % F(8)    = 1;   % Phi
        %         % f       = losDerivs(los,F);
        %         % JacobLos = [JacobLos; f(1), f(2), f(3)];
        %         % JacRVToK  = Jacobian(1:3,2:7);
        %         % JacobAlphaToK = JacobLos*JacRVToK
        %         % JacobianLos(2:3,2:7)
        %         %
        %         %F       = zeros(8,1);
        %         %F(5)    = 1;      % Range
        %         %f       = losDerivs(los,F);
        %         %JacobLos = [JacobLos; f(1), f(2), f(3)];
        %         %
        %         % JacobMeas = inv(JacobLos);
        %         % % Drop the range column
        %         % JacobMeas = JacobMeas(:,1:2);
        %         % JacToZ    = Jacobian(1:3,2:7);
        %         % CombJac   = JacToZ'*JacobMeas;
        %         % InvCov    = InvCov + CombJac'*InvCovMeas*CombJac;
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         %deltaA  =  Predicted - Observed;
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         %   Theta_diff and Phi_diff
        %         thetaDiff =  theta - thetaM;
        %         phiDiff   =  phi   - phiM;
        %         while(phiDiff <   -pi)
        %             phiDiff = PhiDiff + twopi;
        %         end
        %         while(phiDiff > pi)
        %             phiDiff = phiDiff - twopi;
        %         end
        %
        %         crossP = cross([cos(phiM), sin(phiM), 0], [cos(phi),sin(phi), 0] );
        %         dotP   =   dot([cos(phiM), sin(phiM), 0], [cos(phi),sin(phi), 0] );
        %         phiDiff = atan2(crossP,dotP);
        %         deltaA  =  [thetaDiff; phiDiff(3)]
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % Turn off Gravity !!!!!!
        %         %CombJac   = JacobianLos(2:3,2:7);
        %         % Include Gravity Perturbation
        %         CombJac   = JacobianLos(2:3,2:7)*PhiFundamental;
        %         %CombJac   = JacobianLos0(2:3,2:7)*PhiFundamental;
        %         InvCov     = InvCov + CombJac'*InvCovMeas*CombJac;
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % line of sight components:
        %         %Mdelta   = [ cos(theta)*cos(phi) -cos(theta)*sin(phi)  -sin(theta);...
        %         %            -sin(theta)*sin(phi)  sin(theta)*cos(phi)      0      ];
        %         %deltaA   = Mdelta*[los-losM];
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %         CovVec   = InvCovMeas*deltaA;
        %         dataVec  = CombJac'*CovVec;
        %         BetaInc  = BetaInc + dataVec;
        %         ChiSq    = ChiSq  + 0.5*deltaA'*InvCovMeas*deltaA;
        %
        %         Matrix   = zeros(6);
        %         % for icol = 1:2
        %         %     %Matrix = Matrix + HessianLos{icol}(2:7,2:7)*PhiFundamental*CovVec(icol)*PhiFundamental;
        %         %     Matrix = Matrix + PhiFundamental'*HessianLos{icol}(2:7,2:7)*PhiFundamental*CovVec(icol);
        %         %     %Matrix = Matrix + PhiFundamental'*HessianLos0{icol}(2:7,2:7)*PhiFundamental*CovVec(icol);
        %         % end
        %
        %         % JRS
        %         % The following cell addresses are WRONG!!! 1->2 and 2->3
        %         % But the code doesn't converge yet.
        %
        %         Matrix = Matrix + PhiFundamental'*HessianLos{2}(2:7,2:7)*PhiFundamental*CovVec(1);
        %         Matrix = Matrix + PhiFundamental'*HessianLos{3}(2:7,2:7)*PhiFundamental*CovVec(2);
        %         %HessianSyst = HessianSyst + PhiFundamental*PhiFundamental*Matrix;
        %         %HessianSyst = HessianSyst + PhiFundamental*Matrix*PhiFundamental;
        %         HessianSyst = HessianSyst + Matrix;
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         %   End Theta_diff and Phi_diff
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %         % deltaR   = Predicted - Observed
        %         % %deltaR   = Rinit - Rdata + Jacobian(1:3,2:7)'*delZ;
        %         % ChiSq    = ChiSq  + 0.5*deltaR'*InvCovM*deltaR;
        %         %
        %         % [Jacob, JacobFinite, Hess, HessFinite] = KLagrange.DigitalJustice(LineCheck,Sensor);
        %         %
        %         % betaNR   = betaNR + JacobianLos(2:3,2:7)'*InvCovM*deltaR;
        %         % hessPart = zeros(6);
        %         % for k = 2:3
        %         %     hessPart = hessPart + HessianLos{k}(2:7,2:7)*deltaR(k-1);
        %         % end
        %         % hessNR  = hessPart + JacobianLos(2:3,2:7)'*InvCovM*JacobianLos(2:3,2:7);
        %         % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % %   U & V views:
        %         % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % deltaA  =  [ParamList(10); ParamList(11)];
        %         % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % % Turn off Gravity !!!!!!
        %         % %CombJac   = JacobianLos(2:3,2:7);
        %         % % Include Gravity Perturbation
        %         % CombJac   = JacobianLos(10:11,2:7)*PhiFundamental;
        %         % %CombJac   = JacobianLos0(2:3,2:7)*PhiFundamental;
        %         % %InvCov     = InvCov + CombJac'*InvCovMeas*CombJac;
        %         % CovVec   = InvCovMeas*deltaA;
        %         % dataVec  = CombJac'*CovVec;
        %         % %Jacobian(7:8,2:7)'*InvCovMeas*deltaA
        %         % BetaInc  = BetaInc + dataVec;
        %         % ChiSq    = ChiSq  + 0.5*deltaA'*InvCovMeas*deltaA;
        %         % Matrix   = zeros(6);
        %         % %for icol = 1:2
        %         %     %Matrix = Matrix + HessianLos{icol}(2:7,2:7)*PhiFundamental*CovVec(icol)*PhiFundamental;
        %         %     %Matrix = Matrix + PhiFundamental'*HessianLos0{icol}(2:7,2:7)*PhiFundamental*CovVec(icol);
        %         % %end
        %         % Matrix = Matrix + PhiFundamental'*HessianLos{10}(2:7,2:7)*PhiFundamental*CovVec(1);
        %         % Matrix = Matrix + PhiFundamental'*HessianLos{11}(2:7,2:7)*PhiFundamental*CovVec(2);
        %         % %HessianSyst = HessianSyst + PhiFundamental*PhiFundamental*Matrix;
        %         % %HessianSyst = HessianSyst + PhiFundamental*Matrix*PhiFundamental;
        %         % %HessianSyst = HessianSyst + CombJac'*InvCovMeas*CombJac;
        %         % HessianSyst = HessianSyst + Matrix + CombJac'*InvCovMeas*CombJac;
        %         % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % %   End of U & V views
        %         % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         %[Jacob1, JacobFinite1, Hess1, HessFinite1] = KeplerObj.DigitalJustice(127, Sensor, Observed);
        %         %[Jacob2, JacobFinite2, Hess2, HessFinite2] = KeplerObj.DigitalJustice(128, Sensor, Observed);
        %         % Jacob1
        %         % JacobianLos(7,:)
        %         % Jacob2
        %         % JacobianLos(8,:)
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %     end
        %     ChiSqArray   = [ChiSqArray, log(ChiSq)];
        %     delta        = -(HessianSyst+InvCov)\BetaInc;
        %     %delta        = -(HessianSyst)\BetaInc;
        %     %deltaPar     = -hessNR\betaNR
        %     %TangentArray = [TangentArray, log(norm(betaNR))]
        %     %Kepler       =  Kepler + deltaPar'
        %     TangentArray = [TangentArray, log(norm(BetaInc))]
        %     Kepler       =  Kepler + delta'
        % end
        %
        % figure
        % plot([1:numel(ChiSqArray)], ChiSqArray)
        % title('log(\chi^2) per iteration')
        % ylabel('log(\chi^2)')
        % xlabel('Newton-Raphson Iteration Number')
        %
        % figure
        % plot([1:numel(TangentArray)], TangentArray)
        % title('Log norm of the Derivatives of \chi^2 per iteration')
        % ylabel('Log norm \chi^2 Derivatives')
        % xlabel('Newton-Raphson Iteration Number')

        %
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % 5)     Begin Application of Hessian and Second Order fit
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TangentArray   = [];
        % ZhatArray      = [];
        % DelZhatArray   = [];
        %
        % ChiSqArray     = [];
        % DelChiSqArray  = [];
        %
        % DChiSqArray    = [];
        % DiffChiSqArray = [];
        %
        % Zhat           = zeros(6,1);
        % delta          = zeros(6,1);
        % PhiFundamental = eye(6);
        % KeplerIteration  = KeplerFit
        % %KLagrange        = LagrangePlanetary(KeplerIteration, units);
        % %KeplerIteration = KeplerInit + Zhat';   % Reset the Gravity Model Reference Point
        %
        % %KeplerIteration  = Kepler + Zhat;
        % %Zhat           = zeros(6,1);
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % losM     = track_data(iFit,6:8)';
        % losM     = losM/norm(losM);   % los
        % rangeM   = norm(losM);
        % thetaM   = acos(losM(3)/rangeM);
        % phiM     = atan2(losM(2), losM(1));
        % Observed = [thetaM; phiM];
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % %Sat Position at measured position
        % Rsensor  = [track_data(iFit,3:5)']/DU;
        % Sensor   = [Rsensor; zeros(6,1)];
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % KLagrange       = LagrangePlanetary(KeplerIteration, units);
        % [Rinit,Vinit,ParamList,Jacobian,Hessian,GravityCan] = KLagrange.OrbitDerivatives(tFit, Sensor, losM);
        % [Rpos, Rdot]
        % [(Rpos-Rinit), (Rdot-Vinit)]
        % [Zdot, Inhomogeneous, Perturbation, DMatrix, LagrangeTerms, LagrangeGrads]= KLagrange.EquationsOfMotion(zeros(6,1));
        %
        % for NR = 1:30
        %
        %     ChiSq       = 0.0;
        %     InvCov      = zeros(6);
        %     BetaInc     = zeros(6,1);
        %     HessianSyst = zeros(6);
        %
        %     for ii = 1:rows
        %         [Zdot, Inhomogeneous, Perturbation, DMatrix, LagrangeTerms, LagrangeGrads]= KLagrange.EquationsOfMotion(Zhat);
        %         tRecord   =  track_data(ii,2)/TU;
        %         t = tRecord;
        %         tInterval    = tRecord - tFit;
        %         %[Zdot, Inhomogeneous, Perturbation] = KLagrange.EquationsOfMotion(KeplerInit');
        %         %[Zdot, Inhomogeneous, Perturbation] = KLagrange.EquationsOfMotion(KeplerIteration');
        %         %[Zdot, Inhomogeneous, Perturbation] = KLagrange.EquationsOfMotion(zeros(6,1));
        %         %[Zdot, Inhomogeneous, Perturbation, DMatrix, LagrangeTerms, LagrangeGrads]= KLagrange.EquationsOfMotion(zeros(6,1));
        %         %[Zdot, Inhomogeneous, Perturbation, DMatrix, LagrangeTerms, LagrangeGrads]= KLagrange.EquationsOfMotion(delta);
        %         %t = tFit;
        %         Nbins            = 600;
        %         delt             = tInterval/Nbins;
        %         tRiem            = -delt/2.0;
        %         PhiInvIntegrated = zeros(6,1);
        %         for i = 1:Nbins
        %             tRiem = tRiem + delt;
        %             PhiInvIntegrated = PhiInvIntegrated + expm(-Perturbation*tRiem)*Inhomogeneous;
        %             %PhiInvIntegrated = PhiInvIntegrated + ComplementarySoln(-Perturbation*tRiem);
        %         end
        %         PhiInvIntegrated = PhiInvIntegrated*delt;
        %         PhiFundamental   = expm(Perturbation*tInterval);
        %         %PhiFundamental   = ComplementarySoln(Perturbation*tInterval)
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % %Set One
        %         % delZ             = PhiFundamental*(PhiInvIntegrated*Inhomogeneous);
        %         % delZ'
        %         % KeplerCurrent      = KeplerInit      + delZ';
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % Set Two
        %         %delZ            = PhiFundamental*(Zhat + PhiInvIntegrated*Inhomogeneous);
        %         %delZ             = PhiFundamental*(PhiInvIntegrated*Inhomogeneous);
        %         %delZ             = PhiFundamental*PhiInvIntegrated;
        %         delZ            = PhiFundamental*(Zhat + PhiInvIntegrated);
        %         delZ'
        %         %KeplerCurrent     = KeplerInit      + delZ' + Zhat'
        %         KeplerCurrent     = KeplerIteration + delZ';
        %         %KeplerCurrent      = KeplerIteration;
        %         %KeplerCurrent      = KeplerInit      + delZ';
        %         %KeplerCurrent      = KeplerInit;
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         KLagrangeProp      = LagrangePlanetary(KeplerCurrent, units);
        %         % Error in following call: check out order of calls:
        %         %[Zdot, Inhomogeneous, Perturbation] = KLagrangeProp.EquationsOfMotion(KeplerCurrent');
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % %KeplerProp = [e a Inclination omega Omega Mp]  % Canonical Units
        %         % M                    = meanMotion*t + Mp;
        %         % Kepler               = [e a Inclination omega Omega M] + delZ';
        %         % KLagrange            = LagrangePlanetary(Kepler, units);
        %         % %t = 0.0
        %         % [Rextrap, Vextrap]   = KLagrange.extrapolate(0)
        %         % [Rextrap, Vextrap, Hessian, Gradient, Jacobian] = KLagrange.OrbitDerivatives(0)
        %         % Rextrap
        %         % Vextrap
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         %t = 0.0
        %         %[Rextrap, Vextrap]   = KLagrange.extrapolate(0);
        %         %Rextrap
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % Measurement Angles:
        %         %losM    = state_vector_losArray(1:3,ii);
        %         losM     = track_data(ii,6:8)';
        %         losM     = losM/norm(losM);   % los
        %         rangeM   = norm(losM);
        %         thetaM   = acos(losM(3)/rangeM);
        %         phiM     = atan2(losM(2), losM(1));
        %         Observed =[thetaM; phiM];
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        %         % Sat Position from Polynomial at tFit
        %         %Sensor = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU]
        %         % Sat Position from Polynomial at t_ii
        %         %Sensor = [state_vector_satArray(1:3,ii)/DU; state_vector_satArray(4:6,ii)/VU; state_vector_satArray(7:9,ii)/AU]
        %         % Sat Position from Data Measurements
        %         Sensor = [track_data(ii,3:5)'/DU; zeros(6,1)];
        %         % Sat Position from Orbit Fit to Kepler elements:
        %         %[Rsensor, Vsensor, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessSat] = KeplerSat.OrbitDerivatives(tRecord, Sensor, Observed);
        %         %Sensor = [Rsensor; Vsensor; zeros(3,1)];
        %         %[Rextrap, Vextrap, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan] = KLagrangeProp.OrbitDerivatives(tRecord, Sensor);
        %         [Rextrap, Vextrap, ParamList, Jacobian, Hessian, GravityCan] = KLagrangeProp.OrbitDerivatives(tRecord, Sensor, losM);
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % Alternative Method of correcting for gravitational perturbations:
        %         %Delta = Jacobian(:,2:7)*delZ;
        %         %Rextrap + Delta(1:3);
        %         %Vextrap + Delta(4:6);
        %         %los     = Rextrap - state_vector_satArray(1:3,ii)/DU;
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         theta = ParamList(2)
        %         phi   = ParamList(3)
        %         InvCovMeas = eye(2)*10^(-2);
        %         los   = [ParamList(7); ParamList(8); ParamList(9)];
        %         %InvCovMeas = eye(2)*10^(6);
        %         %InvCovMeas = eye(2);
        %         % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % % %Delta = Jacobian(:,2:7)*delZ;
        %         % % %Rextrap + Delta(1:3);
        %         % % %Vextrap + Delta(4:6);
        %         % % %los     = Rextrap - state_vector_satArray(1:3,ii)/DU;
        %         % los      = Rextrap -            track_data(ii,3:5)'/DU
        %         % % % On the fly Backwards Differentiation:
        %         % % % Prediction:
        %         % % % lx = los(1)    Line 1
        %         % % % ly = los(2)    Line 2
        %         % % % lz = los(3)    Line 3
        %         % % %range   = norm(los)
        %         % % %theta   = acos(los(3)/range)
        %         % % %phi     = atan2(los(2),los(1))
        %         % rangeSq  = los(1)^2 + los(2)^2 + los(3)^2;   % line 4
        %         % range    = sqrt(rangeSq)                    % line 5
        %         % u        = los(3)/range;                     % line 6
        %         % theta   = acos(u)                           % line 7
        %         % phi     = atan2(los(2),los(1))              % line 8
        %         % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % JacobLos = [];
        %         % %
        %         % F       = zeros(8,1);
        %         % F(7)    = 1;     % Theta
        %         % f       = losDerivs(los,F);
        %         % JacobLos = [JacobLos; f(1), f(2), f(3)];
        %         % %
        %         % F       = zeros(8,1);
        %         % F(8)    = 1;   % Phi
        %         % f       = losDerivs(los,F);
        %         % JacobLos = [JacobLos; f(1), f(2), f(3)];
        %         % JacRVToK  = Jacobian(1:3,2:7);
        %         % JacobAlphaToK = JacobLos*JacRVToK
        %         % JacobianLos(2:3,2:7)
        %         %
        %         %F       = zeros(8,1);
        %         %F(5)    = 1;      % Range
        %         %f       = losDerivs(los,F);
        %         %JacobLos = [JacobLos; f(1), f(2), f(3)];
        %         %
        %         % JacobMeas = inv(JacobLos);
        %         % % Drop the range column
        %         % JacobMeas = JacobMeas(:,1:2);
        %         % JacToZ    = Jacobian(1:3,2:7);
        %         % CombJac   = JacToZ'*JacobMeas;
        %         % InvCov    = InvCov + CombJac'*InvCovMeas*CombJac;
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % Turn off Gravity !!!!!!
        %         % CombJac   = JacobianLos(2:3,2:7);
        %         % %Include Gravity Perturbation
        %         % CombJac   = JacobianLos(2:3,2:7)*PhiFundamental;
        %         % %CombJac   = JacobianLos0(2:3,2:7)*PhiFundamental;
        %         % InvCov     = InvCov + CombJac'*InvCovMeas*CombJac;
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % %Measurement
        %         % %losM    = state_vector_losArray(1:3,ii);
        %         % %losM    = track_data(ii,6:8)';
        %         % %losM    = losM/norm(losM);   % los
        %         % rangeM  = norm(losM);
        %         % thetaM  = acos(losM(3)/rangeM)
        %         % phiM    = atan2(losM(2), losM(1))
        %         % deltaA  = [theta - thetaM; phi - phiM];
        %         % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         % CombJac   = JacobianLos(10:11,2:7)*PhiFundamental;
        %         % % %CombJac   = JacobianLos0(2:3,2:7)*PhiFundamental;
        %         % InvCov     = InvCov + CombJac'*InvCovMeas*CombJac;
        %         % % line of sight components: JRS
        %         % Mdelta   = [ cos(theta)*cos(phi) cos(theta)*sin(phi)  -sin(theta);...
        %         %             -sin(phi)  cos(phi)      0      ];
        %         % deltaA   = Mdelta*[los-losM]
        %         % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lo
        %         % ChiSqInc   = 0.5*deltaA'*InvCovMeas*deltaA;
        %         % ChiSq      = ChiSq   + ChiSqInc;
        %         %
        %         % CovVec     = InvCovMeas*deltaA;
        %         % dataVec    = CombJac'*CovVec;
        %         % %BetaInc    = BetaInc + dataVec;
        %         % %SymCheck   = [deltaA'*InvCovMeas*CombJac]'
        %         % BetaInc    = BetaInc + dataVec;
        %         % Matrix     = zeros(6);
        %         % Matrix = Matrix + PhiFundamental'*HessianLos{10}(2:7,2:7)*PhiFundamental*CovVec(1);
        %         % Matrix = Matrix + PhiFundamental'*HessianLos{11}(2:7,2:7)*PhiFundamental*CovVec(2);
        %         % % for icol   = 1:2
        %         % %     %Matrix = Matrix + HessianLos{icol}(2:7,2:7)*PhiFundamental*CovVec(icol)*PhiFundamental;
        %         % %     Matrix = Matrix + PhiFundamental'*HessianLos{icol}(2:7,2:7)*PhiFundamental*CovVec(icol);
        %         % %     %Matrix = Matrix + PhiFundamental'*HessianLos0{icol}(2:7,2:7)*PhiFundamental*CovVec(icol);
        %         % % end
        %         % %HessianSyst = HessianSyst + PhiFundamental*PhiFundamental*Matrix;
        %         % %HessianSyst = HessianSyst + PhiFundamental*Matrix*PhiFundamental;
        %         % HessianSyst = HessianSyst + Matrix;
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lo
        %         %CombJac   = JacobianLos(10:11,2:7)*PhiFundamental;
        %         %InvCov     = InvCov + CombJac'*InvCovMeas*CombJac;
        %         %Mdelta   = [ cos(theta)*cos(phi) cos(theta)*sin(phi)  -sin(theta);...
        %         %            -sin(phi)  cos(phi)      0      ];
        %         %deltaA   = Mdelta*[los-losM]
        %         %CovVec     = InvCovMeas*deltaA;
        %         %Matrix     = zeros(6);
        %         %Matrix = Matrix + PhiFundamental'*HessianLos{10}(2:7,2:7)*PhiFundamental*CovVec(1);
        %         %Matrix = Matrix + PhiFundamental'*HessianLos{11}(2:7,2:7)*PhiFundamental*CovVec(2);
        %         %HessianSyst = HessianSyst + Matrix;
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lo
        %         Mdelta    = [ cos(thetaM)*cos(phiM) cos(thetaM)*sin(phiM)  -sin(thetaM);...
        %                       -sin(phiM)  cos(phiM)      0      ];
        %         InvCovM  = Mdelta'*Mdelta;
        %         %InvCovM  = eye(3);
        %         CombJac  = Jacobian(14:16,2:7)*PhiFundamental;
        %         deltaA   = los-losM;
        %         dataVec  = CombJac'*InvCovM*deltaA;
        %         CovVec    = InvCovM*deltaA;
        %         %Jac2nd    = CombJac'*InvCovM*CombJac;
        %         Matrix    = zeros(6);
        %         Matrix = Matrix + PhiFundamental'*Hessian{14}(2:7,2:7)*PhiFundamental*CovVec(1);
        %         Matrix = Matrix + PhiFundamental'*Hessian{15}(2:7,2:7)*PhiFundamental*CovVec(2);
        %         Matrix = Matrix + PhiFundamental'*Hessian{16}(2:7,2:7)*PhiFundamental*CovVec(3);
        %         HessianSyst = HessianSyst + Matrix;
        %         InvCov      = InvCov      + CombJac'*InvCovM*CombJac;
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lo
        %         ChiSqInc   = 0.5*deltaA'*InvCovM*deltaA;
        %         ChiSq      = ChiSq   + ChiSqInc;
        %         BetaInc    = BetaInc + dataVec;
        %     end
        %     delta      = (HessianSyst+InvCov)\BetaInc
        %     %delParam   = InvCov\BetaInc
        %     %KeplerInit = KeplerInit + delParam';
        %     %KeplerInit'
        %     %delZhat  = -delParam - Zhat;
        %     %Zhat     = -delParam
        %     Zhat     = Zhat - delta;
        %
        %     ParamComp = [delta, Zhat, KeplerInit' + Zhat]
        %
        %     if NR > 1
        %         ZhatPrev       = ZhatArray(:,NR-1) - Zhat;
        %         DelZhatArray   = [DelZhatArray, ZhatPrev];
        %         ChiSqPrev      = ChiSqArray(NR-1);
        %         deltaChiSq     = ChiSqPrev - log(ChiSq);
        %         DelChiSqArray  = [DelChiSqArray, deltaChiSq];
        %         diffChiSq      =  DChiSqArray(NR-1) - log(BetaInc);
        %         DiffChiSqArray = [DiffChiSqArray, diffChiSq];
        %     end
        %     DChiSqArray = [DChiSqArray, log(norm(BetaInc))];
        %     ChiSqArray  = [ChiSqArray, log(abs(ChiSq))];
        %     ZhatArray   = [ZhatArray, Zhat];
        %     %KeplerProp = [e a Inclination omega Omega Mp]  % Canonical Units
        %     %M               = meanMotion*t + Mp;
        %     %Kepler          = [e a Inclination omega Omega M];
        %     %KeplerInit       = [e a Inclination omega Omega Mp]+ Zhat'
        %     %KLagrange        = LagrangePlanetary(KeplerInit, units);
        % end
        %
        % KeplerFit = KeplerIteration;
        % CovFit    = inv(InvCov);
        %
        % figure
        % plot([1:numel(ChiSqArray)], ChiSqArray)
        % title('Log \chi^2 per iteration')
        % ylabel('\chi^2')
        % xlabel('Newton-Raphson Iteration Number')
        %
        % figure
        % plot([1:size(DelZhatArray,2)], DelZhatArray)
        % title('Parameter Differences per iteration')
        % ylabel('\Delta Parameter Set')
        % xlabel('Newton-Raphson Iteration Number')
        %
        % figure
        % plot([1:size(DiffChiSqArray,2)], DiffChiSqArray)
        % title('Change in \chi^2 derivative per iteration')
        % ylabel('\delta \chi^2 derivatives')
        % xlabel('Newton-Raphson Iteration Number')
        %
        % figure
        % plot([1:size(DelChiSqArray,2)], DelChiSqArray)
        % title('\chi^2 differences per iteration')
        % ylabel('\Delta \chi^2')
        % xlabel('Newton-Raphson Iteration Number')
        % %hold off
        % % Update the Reference Orbit Chief -> New Chief:
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %losM    = state_vector_losArray(1:3,ii);
        %losM     = track_data(ii,6:8)';
        %losM     = losM/norm(losM);   % los
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %KLagrange  = LagrangePlanetary(KeplerIteration+Zhat', units);
        KLagrange  = LagrangePlanetary(KeplerFit, units);
        %Sensor     = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU];
        losM     = track_data(iFit,6:8)';
        losM     = losM/norm(losM);
        Rsensor  = [track_data(iFit,3:5)']/DU;
        Sensor   = [Rsensor; zeros(6,1)];
        [Rpos, Rdot, ParamList, Jacobian, Hessian, GradientCan] = KLagrange.OrbitDerivatives(tFit, Sensor, losM, InvCovUV);
        SV7Fit = [Rpos*DU; Rdot*VU; tFit*TU]
        %[e a Inclination omega Omega Mp Period JacECI_2_Kepler] = kepler(tFit, Rpos, Rdot, mu);
        Jac        = eye(6);
        JacKepler_2_ECI = Jacobian(1:6,2:7);
        Jac        = JacKepler_2_ECI*Jac;        
        CovMKS = (Jac*CovFit*Jac')./units.COV6Ref
        %JacKepler_2_ECI*JacECI_2_Kepler
        KeplerObject = KeplerFromECI(tFit, Rpos, Rdot, units);
        InMotion = KeplerObject.Extrapolator(tEnd);
        %KeplerObject.ECIPropagator
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TStartPlot = -500.0;
        tStartPlot =  TStartPlot/TU;
        TEndPlot   =   1500.0;
        tEndPlot   =  TEndPlot/TU;
        timeArray = [(TStartPlot-0.5):(TEndPlot-0.5)]/TU;        
        % Try to Integrate the Equations of Motion for the Classical Elements
        Y             = zeros(6,1);
        options   = odeset('RelTol', 1e-10, 'AbsTol', 1e-13);
        dJdt      = perturbations_odefun(KLagrange,tFit,Y,units);
        Y         = zeros(6,1);
        % Propagate from fit Reference time tFit back to tBeg
        [Tout, Yout] = ode45(@(t,y) perturbations_odefun(KLagrange, t, y, units), [tFit, tStartPlot], Y, options);
        %Yout(end,:)
        Perturbed = [Tout, Yout];
        Perturbed = sortrows(Perturbed,1);
        Y = Perturbed(end,2:7);
        % Propagate from fit Reference time tFit back to tEnd
        if tFit < tEndPlot
            [Tout, Yout] = ode45(@(t,y) perturbations_odefun(KLagrange, t, y, units), [tFit, tEndPlot], Y, options);
            %Yout(end,:)
            Perturbed(end,:) = [];
            Perturbed = [Perturbed; Tout, Yout];
            Perturbed = sortrows(Perturbed,1);
        end
        KeplerInterpolated = interp1(Perturbed(:,1),Perturbed(:,2:7),timeArray);
        KeplerInterpolated(1,:) = [];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        Nbins          = 2000;
        % tOffset        = 200/TU;
        % tInterval      = tEnd - tBeg;
        % tRange         = tEnd + 3.0*tOffset - tBeg;
        % tStart         = tBeg - 1.75*tOffset;
        % tBins          = tRange/Nbins;
        % time           = tStart - tBins/2.0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        % Fit Results
        TimeArray      = [];
        CovArrayPos       = [];
        CovArrayVel       = [];
        %Time  = 200.0 - 0.5;
        Time  = TStartPlot - 0.5;
        for I = 1:Nbins
            Time       = Time + 1.0;
            time       = Time/TU;
            TimeArray  = [TimeArray, Time];
            
            Zhat      = KeplerInterpolated(I,:);
            KLagrange = LagrangePlanetary(KeplerFit+Zhat, units);
            Jac        = eye(6);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [Rpos, Rdot, ParamList, Jacobian, ~, ~] = OrbitDerivatives(KLagrange, time, Sensor, losM, InvCovUV);            
            %[Rpos, Rdot, ParamList, Jacobian, ~, ~] = KLagrange.OrbitDerivatives(time, Sensor, losM, InvCovUV);
            JacKepler_2_ECI = Jacobian(1:6,2:7);
            Jac        = JacKepler_2_ECI*Jac;
            %InMotion   = KeplerObject.Extrapolator(time);
            %Propagator = KeplerObject.ECIPropagator(1:6,1:6);
            %Propagator
            %Jac        = Propagator*Jac;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Final UNits Conversion back to MKS;
            %Jac        = units.JacRef6*Jac;
            % Extrapolated Covariance in MKS
            CovMKS      = (Jac*CovFit*Jac')./units.COV6Ref;
            %SizeCov    = sqrt(sum(eig(CovMKS(1:3,1:3))));
            SizeCov     = sqrt(trace(CovMKS(1:3,1:3)));
            CovArrayPos = [CovArrayPos, SizeCov];
            %SizeCov    = sqrt((eig(CovMKS(4:6,4:6))));
            SizeCov     = sqrt(trace(CovMKS(4:6,4:6)));
            CovArrayVel = [CovArrayVel, SizeCov];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        SizeInitialStatePos = [];
        SizeInitialStateVel = [];
        %TimeArray = [];
        %Time  = 200.0 - 0.5;
        Time  =  -500 - 0.5;
        for I = 1:Nbins
            Time       = Time + 1.0;
            %TimeArray = [TimeArray, Time];
            [state_vector, CovSV] = PolyCorrExtrapToFitTime(Iorder, Time, TPolyFit, fitInitial, CovInitial, units);
            sizeCov = sqrt(trace(CovSV(1:3,1:3)));
            SizeInitialStatePos = [SizeInitialStatePos, sizeCov];
            sizeCov = sqrt(trace(CovSV(4:6,4:6)));
            SizeInitialStateVel = [SizeInitialStateVel, sizeCov];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        % Simulated Data
        %Residuals     = [];
        EventTime      = [];
        TimeData      = [];
        DiffArrayPos   = [];
        DiffArrayVel   = [];
        RExtrapPos     = [];
        RExtrapVel     = [];
        for I = 1:rows
            time       = track_data(I,2)/TU;
            TimeData  = [TimeData, time*TU];
            
            Rdata      = track_data(I,10:12);
            %EventTime   = [EventTime, time];
            InMotion   = KeplerObject.Extrapolator(time);
            
            Rextrap    = KeplerObject.Rextrap*DU;
            RExtrapPos = [RExtrapPos, Rextrap];
            Diff       = [Rdata(1) - Rextrap(1);...
                          Rdata(2) - Rextrap(2);...
                          Rdata(3) - Rextrap(3)];
            Resids     = sqrt(sum(Diff.*Diff));
            DiffArrayPos  = [DiffArrayPos, Resids];
            
            Vdata      = track_data(I,13:15);
            Vextrap    = KeplerObject.Vextrap*VU;
            RExtrapVel = [RExtrapVel, Vextrap];
            Diff       = [Vdata(1) - Vextrap(1);...
                          Vdata(2) - Vextrap(2);...
                          Vdata(3) - Vextrap(3)];
            Resids     = sqrt(sum(Diff.*Diff));
            DiffArrayVel  = [DiffArrayVel, Resids];
                        
            %StateVector7 = KeplerObject.StateVector7
            %StateVector7(1:3)*DU
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        % Position Resolution
        figure 
        plot(TimeArray,CovArrayPos,'r','LineWidth',2)        
        hold on
        plot(TimeArray, SizeInitialStatePos,'b','LineWidth',2);
        hold on        
        plot(TimeData,DiffArrayPos,'g','LineWidth',2)
        title('Square-root of the trace of position-space Covariance Matrix' )
        ylim([0 1000]);
        ylabel('Meters')
        xlim([0 1400])
        xlabel('Time (seconds)');        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        % Velocity Resolution
        figure 
        plot(TimeArray,CovArrayVel,'r','LineWidth',2)        
        hold on
        plot(TimeArray, SizeInitialStateVel,'b','LineWidth',2);
        hold on        
        plot(TimeData,DiffArrayVel,'g','LineWidth',2)
        title('Square-root of the trace of Velocity-space Covariance Matrix' )
        ylim([0 50]);
        ylabel('Meters per Second')
        xlim([0 1400])
        xlabel('Time (seconds)');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        figure
        plot3(track_data(:,10), track_data(:,11), track_data(:,12));
        hold on
        plot3(RExtrapPos(1,:), RExtrapPos(2,:), RExtrapPos(3,:));        
        title('Data and Kepler Extrapolation')
        xlabel(' X meters')
        ylabel(' Y meters')
        zlabel(' Z meters')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        figure
        plot3(track_data(:,13), track_data(:,14), track_data(:,15));
        hold on
        plot3(RExtrapVel(1,:), RExtrapVel(2,:), RExtrapVel(3,:));        
        title('Data and Kepler Extrapolation')
        xlabel(' Vx meters per second')
        ylabel(' Vy meters per second')
        zlabel(' Vz meters per second')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        %tFit        = tFit;
        %toTime     = track_data(end,2)
        %toTime    = track_data(end,2) + tFirst
        %toTime    = 71577.703  
        %toTime     = 71579.507
        toTime     = tEnd*TU
        %SV7Fit     = [Rpos*DU; Rdot*VU; TFit + tFirst];
        [Rpos, Rdot, ParamList, Jacobian, Hessian, GradientCan] = KLagrange.OrbitDerivatives(tFit, Sensor, losM, InvCovUV);
        SV7Fit     = [Rpos*DU; Rdot*VU; TFit]
        SV7Prop    = propagateECI(SV7Fit, toTime)
        %SV7Final   
        
        KLagrange     = LagrangePlanetary(KeplerFit, units);
        [Rpos, Rdot, ParamList, Jacobian, Hessian, GradientCan] = KLagrange.OrbitDerivatives(tFit, Sensor, losM, InvCovUV);
        [Rpos*DU; Rdot*VU; tFit*TU]
        options       = odeset('RelTol', 1e-10, 'AbsTol', 1e-13);
        Y             = zeros(6,1);
        Perturbed     = [tFit, Y'];
        dJdt          = perturbations_odefun(KLagrange,tFit,Y,units);        
        % Propagate from fit Reference time tFit back to tBeg
        % [Tout, Yout] = ode45(@(t,y) perturbations_odefun(KLagrange, t, y, units), [tFit, tBeg], Y, options);
        % %Yout(end,:)
        % Perturbed = [Tout, Yout];
        % Perturbed = sortrows(Perturbed,1);
        % Y = Perturbed(end,2:7);
        % % Propagate from fit Reference time tFit back to tEnd
        if tFit < tEnd
            [Tout, Yout] = ode45(@(t,y) perturbations_odefun(KLagrange, t, y, units), [tFit, tEnd], Y, options);
            %Yout(end,:)
            Perturbed(end,:) = [];
            Perturbed = [Perturbed; Tout, Yout];
            Perturbed = sortrows(Perturbed,1);
        end        
        Zhat = Perturbed(end,2:7)
        KLagrange     = LagrangePlanetary(KeplerFit+Zhat, units);
        [Rpos, Rdot, ParamList, Jacobian, Hessian, GradientCan] = KLagrange.OrbitDerivatives(tEnd, Sensor, losM, InvCovUV);
        JacKepler_2_ECI = Jacobian(1:6,2:7);
        Jac        = JacKepler_2_ECI*Jac;
        Jac        = Jac;
        CovMKS     = (Jac*CovFit*Jac')./units.COV6Ref;        
        SizeCov = sqrt(trace(CovMKS(1:3,1:3)))
        %SizeCov    = sqrt(sum(eig(CovMKS(1:3,1:3))));
        %SizeCov    = sqrt(sum(eig(CovMKS(4:6,4:6))));
        SizeCov    = sqrt(trace(CovMKS(4:6,4:6)));
        [Rpos*DU; Rdot*VU; tEnd*TU]
        

        SV7Fit'
        SV7Prop'
        Mrotate           = ECR_to_ECI_PVA_Rotation_Matrix( -toTime );
        StateVectorECEF   = [Mrotate(1:6,1:6)*SV7Prop(1:6)/1000.0]'
        %SV_For_Comparison = [-6418.325,  2835.917,  1254.261, 2.143745,  1.678883,  9.287525]

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
        % KLagrange = LagrangePlanetary(Kepler, units);
        % %t = 0.0
        % [Rextrap, Vextrap]   = KLagrange.extrapolate(0);
        % SV7     = [StateVectorECI(1:6); tFit];
        % SV7Prop = propagateECI(SV7, t*DU);
        % Rpos    = SV7Prop(1:3)/DU;
        % Rdot    = SV7Prop(4:6)/VU;
        % [e a Inclination omega Omega Mp Period Jacobian] = kepler(t, Rpos, Rdot, units);
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
        % % KLagrange          = LagrangePlanetary(KeplerProp, units);
        % % [Rextrap, Vextrap] = KLagrange.extrapolate(t);
        % % [Rextrap; Vextrap]'
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    end