function [KeplerFit, ChiSq, CovFit, SV7MKSFinal, CovMKSFinal] = NewtonStepUVWithPerturbations(track_data, KeplerInit, tFit, InvCovUV, units)

    rows       = size(track_data,1);
    % For Choice of Units
    twopi = units.twopi;
    TU = units.TU;
    DU = units.DU;
    VU = units.VU;
    %AU = units.AU;
    mu = units.mu;
    %sqrtmu = units.sqrtmu;

    timeArray = track_data(:,2)/TU;
    tBeg = min(timeArray);
    tEnd = max(timeArray);


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
    KeplerIteration  = KeplerInit
    
    Mind  = 6
    oparm = zeros(Mind,1);
    oparm(1) = true;
    oparm(2) = true;
    oparm(3) = true;
    oparm(4) = true;
    oparm(5) = true;
    oparm(6) = true;

    for NR = 1:30
        
        nfree    = 0;
        
        IndexRows = find(oparm == 1);
        IndexCols = [IndexRows; Mind+1];

        for ipar=1:Mind
            if oparm(ipar) == true
                nfree = nfree+1;
            end
        end
        
        delp      = zeros(6,1);
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Try to Integrate the Equations of Motion for the Classical Elements
        Y             = zeros(6,1);
        KLagrange = LagrangePlanetary(KeplerIteration, units);
        losM = eye(3,1);
        losM = losM/norm(losM);
        Sensor = [losM; zeros(6,1)];
        [Rextrap, Vextrap, ParamList] = OrbitAtTime(KLagrange, tFit, Sensor, losM, InvCovUV);
        options   = odeset('RelTol', 1e-10, 'AbsTol', 1e-13);
        dJdt      = perturbations_odefun(KLagrange,tFit,Y,units);
        % Propagate from fit Reference time tFit back to tBeg
        %[Tout, Yout] = ode45(@(t,y) perturbations_odefun(KLagrange, t, y, units), [tFit, tBeg], Y, options);
        [Tout, Yout] = ode45(@(t,y) LagrangePlanetary_odefun(KeplerIteration, t, y, units), [tFit, tBeg], Y, options);
        %Yout(end,:)
        Perturbed = [Tout, Yout];
        Perturbed = sortrows(Perturbed,1);
        Y = Perturbed(end,2:7);
        % Propagate from fit Reference time tFit back to tEnd
        if tFit < tEnd
            %[Tout, Yout] = ode45(@(t,y) perturbations_odefun(KLagrange, t, y, units), [tFit, tEnd], Y, options);
            [Tout, Yout] = ode45(@(t,y) LagrangePlanetary_odefun(KeplerIteration, t, y, units), [tFit, tEnd], Y, options);
            %Yout(end,:)
            Perturbed(end,:) = [];
            Perturbed = [Perturbed; Tout, Yout];
            Perturbed = sortrows(Perturbed,1);
        end
        KeplerInterpolated = interp1(Perturbed(:,1),Perturbed(:,2:7),timeArray);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ChiSq            = 0.0;
        BetaInc          = zeros(6,1);
        HessianSyst      = zeros(6);
        InvCovarianceFit = zeros(6);
        InvCovFit        = zeros(6);

        for ii = 1:rows
            tRecord   =  track_data(ii,2)/TU;
            t = tRecord;
            Zhat = KeplerInterpolated(ii,:);
            KLagrange     = LagrangePlanetary(KeplerIteration+Zhat, units);
            %tInterval    = tRecord - time;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Measurement
            %losM    = state_vector_losArray(1:3,ii);
            losM    = track_data(ii,6:8)';
            rangeM  = norm(losM);
            losM    = losM/rangeM;   % los
            % thetaM  = acos(losM(3));
            % phiM    = atan2(losM(2), losM(1));
            % Observed =[thetaM; phiM];
            % Mdelta    = [ cos(thetaM)*cos(phiM) cos(thetaM)*sin(phiM)  -sin(thetaM);...
            %              -sin(phiM)  cos(phiM)      0      ];
            % InvCovM  = Mdelta'*Mdelta;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Sat Position from Polynomial at tFit
            %Sensor = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU]
            % Sat Position from Polynomial at t_ii
            % Sensor = [state_vector_satArray(1:3,ii)/DU; state_vector_satArray(4:6,ii)/VU; state_vector_satArray(7:9,ii)/AU];
            % Sat Position from Data Measurements
            Sensor = [track_data(ii,3:5)'/DU; zeros(6,1)];
            % Sat Position from Orbit Fit to Kepler elements:
            %[Rsensor, Vsensor, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessSat] = KeplerSat.OrbitDerivatives(tRecord, Sensor, Observed);
            %Sensor = [Rsensor; Vsensor; zeros(3,1)];
            %[Rextrap, Vextrap, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan] = KLagrangeProp.OrbitDerivatives(tRecord, Sensor);
            %[Rextrap, Vextrap, ParamList, Jacobian, Hessian, GravityCan] = KLagrange.OrbitDerivatives(tRecord, Sensor, losM);
            %[Rextrap, Vextrap, Hessian, Gradient, Jacobian,
            %JacobianLos, ParamList, GravityCan, HessianLos, RVHessLos] = KLagrangeProp.OrbitDerivatives(tRecord, Sensor, Observed);
            [Rextrap, Vextrap, ParamList] = OrbitAtTime(KLagrange, tRecord, Sensor, losM, InvCovUV);
            % if norm(Rextrap) < 1.0
            %     continue
            % end
            [H_I, Jacob, Hess]            = WorkOrder(KLagrange, 133);
            ChiSq       = ChiSq   + H_I;
            BetaInc     = BetaInc     + Jacob(2:7)';
            HessianSyst = HessianSyst + Hess(2:7,2:7);

            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Mdelta    = [ cos(thetaM)*cos(phiM) cos(thetaM)*sin(phiM)  -sin(thetaM);...
            %               -sin(phiM)  cos(phiM)      0      ];
            % InvCovM  = Mdelta'*Mdelta;
            % %InvCovM  = eye(3);
            % CombJac  = Jacobian(1:3,2:7);
            % deltaA   = los-losM;
            % dataVec  = CombJac'*InvCovM*deltaA;
            % CovVec    = InvCovM*deltaA;
            % %Jac2nd    = CombJac'*InvCovM*CombJac;
            % Matrix    = zeros(6);
            % Matrix = Matrix + Hessian{1}(2:7,2:7)*CovVec(1);
            % Matrix = Matrix + Hessian{2}(2:7,2:7)*CovVec(2);
            % Matrix = Matrix + Hessian{3}(2:7,2:7)*CovVec(3);
            % HessianSyst = HessianSyst + Matrix;
            % InvCov      = InvCov      + CombJac'*InvCovM*CombJac;
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lo
            % ChiSqInc   = 0.5*deltaA'*InvCovM*deltaA;
            % ChiSq      = ChiSq   + ChiSqInc;
            % BetaInc    = BetaInc + dataVec;
            
        end
        %Hessian   = HessianSyst + InvCov;
        Hessian   = HessianSyst;
        %delp       = -Hessian\BetaInc
        Zhat     = Zhat + delp;
        TangentArray = [TangentArray, log(norm(BetaInc))]
        %ParamComp = [delta, Zhat, KeplerInit' + Zhat]
        
        AugmentedMatrix = [Hessian,-BetaInc];
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
        
        % if abs(delp(4) - delp(6)) < 0.0001
        %     oparm(4) = false;
        %     oparm(6) = false;
        % end
        
        KeplerIteration  = KeplerIteration + delp'
        angle            = KeplerIteration(4);
        while(angle <   0.0)
            angle = angle + twopi;
        end
        while(angle > twopi)
            angle = angle - twopi;
        end
        KeplerIteration(4) = angle;
        
        angle         = KeplerIteration(5);
        while(angle <   0.0)
            angle = angle + twopi;
        end
        while(angle > twopi)
            angle = angle - twopi;
        end
        KeplerIteration(5) = angle;

        angle         = KeplerIteration(6);
        while(angle <   0.0)
            angle = angle + twopi;
        end
        while(angle > twopi)
            angle = angle - twopi;
        end
        KeplerIteration(6) = angle;

        if NR > 1
            ZhatPrev       = ZhatArray(:,NR-1) - Zhat;
            DelZhatArray   = [DelZhatArray, ZhatPrev];
            ChiSqPrev      = ChiSqArray(NR-1);
            deltaChiSq     = ChiSqPrev - log(abs(ChiSq));
            DelChiSqArray  = [DelChiSqArray, deltaChiSq];
            diffChiSq      =  DChiSqArray(NR-1) - log(norm(BetaInc));
            DiffChiSqArray = [DiffChiSqArray, diffChiSq];
        end
        DChiSqArray = [DChiSqArray, log(norm(BetaInc))];
        ChiSqArray  = [ChiSqArray, log(abs(ChiSq))];
        ZhatArray   = [ZhatArray, Zhat];
    end

    KeplerFit = KeplerIteration;
    %CovFit    = inv(InvCov);
    CovFit    = inv(Hessian);
    
    ZhatFinal = KeplerInterpolated(rows,:);
    tFinal    =  track_data(rows,2)/TU;
    KLagrange = LagrangePlanetary(KeplerIteration+Zhat, units);
    [Rextrap, Vextrap, ParamList, Jacobian, ~, ~] = OrbitDerivatives(KLagrange, tFinal, Sensor, losM, InvCovUV);
    Jac       = Jacobian(1:6,2:7);
    Jac       = units.JacRef6*Jac;
    SV7MKSFinal  = [Rextrap*DU; Vextrap*VU; tFinal*TU] 
    CovMKSFinal  = Jac*CovFit*Jac';
    
    CovPosArray  = [];
    CovVelArray  = [];
    TimeArray    = [];
    for ii = 1:rows
        tRecord   = track_data(ii,2);
        time      = tRecord/TU;
        TimeArray = [TimeArray, tRecord];
        Zhat      = KeplerInterpolated(ii,:);
        KLagrange = LagrangePlanetary(KeplerIteration+Zhat, units);
        [Rextrap, Vextrap, ParamList, Jacobian, ~, ~] = OrbitDerivatives(KLagrange, time, Sensor, losM, InvCovUV);
        Jac       = Jacobian(1:6,2:7);
        Jac       = units.JacRef6*Jac;
        SV7MKSFinal  = [Rextrap*DU; Vextrap*VU; tRecord];
        CovMKSFinal  = Jac*CovFit*Jac';
        sizeCovPos = sqrt(trace(CovMKSFinal(1:3,1:3)));
        sizeCovVel = sqrt(trace(CovMKSFinal(4:6,4:6)));
        CovPosArray = [CovPosArray, sizeCovPos];
        CovVelArray = [CovVelArray, sizeCovVel];
    end
    
    % Time at the End: 
    disp([' Covariance traces at TEnd: ',num2str(tRecord),' seconds,  ',num2str(sizeCovPos),' m,   ',num2str(sizeCovVel),' m/s']);
    
    figure
    plot(TimeArray, CovPosArray,'r','LineWidth',2)
    title('Square-root of the trace of position-space Covariance Matrix' )
    ylim([0 1000]);
    ylabel('Meters')
    xlim([0 1000])
    xlabel('Time (seconds)');
    
    figure
    plot(TimeArray, CovVelArray,'r','LineWidth',2)
    title('Square-root of the trace of Velocity-space Covariance Matrix' )
    ylim([0 5]);
    ylabel('Meters per Second')
    xlim([0 1000])
    xlabel('Time (seconds)');
    
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
end