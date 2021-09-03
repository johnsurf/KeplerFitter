function [KeplerFit, ChiSq, CovFit] = NewtonStepUV(track_data, KeplerInit, InvCovUV, units)
    rows       = size(track_data,1);

    % For Choice of Units
    twopi = units.twopi;
    TU = units.TU;
    DU = units.DU;
    %VU = units.VU;
    %AU = units.AU;
    mu = units.mu;
    %sqrtmu = units.sqrtmu;
    
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
    for NR = 1:30
        KeplerIteration  = KeplerIteration - delta'
        
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

        KLagrange       = LagrangePlanetary(KeplerIteration, units);

        % if NR == 1
        %     KLagrange       = LagrangePlanetary(KeplerIteration);
        %     %Sensor     = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU];
        %     Sensor = [track_data(1,3:5)'/DU; zeros(6,1)];
        %     % Measurement Angles:
        %     %losM    = state_vector_losArray(1:3,ii);
        %     losM     = track_data(1,6:8)';
        %     losM     = losM/norm(losM);   % los
        %     rangeM   = norm(losM);
        %     thetaM   = acos(losM(3)/rangeM);
        %     phiM     = atan2(losM(2), losM(1));
        %     Observed =[thetaM; phiM];
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     %[Rinit,Vinit,Hessian,Gradient,Jacobian, ~, ~, ~,~,~,~] = KLagrange.OrbitDerivatives(time, Sensor);
        %     [Rinit,Vinit,ParamList,Jacobian,Hessian,GravityCan] = KLagrange.OrbitDerivatives(time, Sensor, Observed);
        %     [Rpos, Rdot]
        %     [(Rpos-Rinit), (Rdot-Vinit)]
        %     [Zdot, Inhomogeneous, Perturbation, DMatrix, LagrangeTerms, LagrangeGrads]= KLagrange.EquationsOfMotion(zeros(6,1));
        %     %     SV6init  = StateVectorECI(1:6)
        %     %     SV6check = [Rinit*DU; Vinit*VU]
        % end

        ChiSq       = 0.0;
        BetaInc     = zeros(6,1);
        HessianSyst = zeros(6);
        InvCov      = zeros(6);

        for ii = 1:rows
            tRecord   =  track_data(ii,2)/TU;
            t = tRecord;
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
            [H_I, Jacob, Hess]            = WorkOrder(KLagrange, 133);
            ChiSq       = ChiSq   + H_I;
            BetaInc     = BetaInc     + Jacob(2:7)';
            HessianSyst = HessianSyst + Hess(2:7,2:7);

            % Jacobian = [];
            % %LineCheck    = 129;                % xLos
            % [H_I, Jacob, Hess] = WorkOrder(KLagrange, 129);
            % Jacobian      = [Jacobian; Jacob];   % row 1
            % %LineCheck    = 130;                % yLos
            % [H_I, Jacob, Hess] = WorkOrder(KLagrange, 130);
            % Jacobian      = [Jacobian; Jacob];   % row 2
            % %LineCheck    = 107;                % zLos
            % [H_I, Jacob, Hess] = WorkOrder(KLagrange, 107);
            % Jacobian      = [Jacobian; Jacob];  % row 3
            % InvCov      = InvCov      + Jacobian(:,2:7)'*InvCovM*Jacobian(:,2:7);

            % theta = ParamList(2);
            % phi   = ParamList(3);
            % %InvCovMeas = eye(2)*10^(-2);
            % los   = [ParamList(7); ParamList(8); ParamList(9)];
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lo
            % %deltaA  = [theta - thetaM; phi - phiM];
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lo
            % %CombJac   = JacobianLos(10:11,2:7)*PhiFundamental;
            % %InvCov     = InvCov + CombJac'*InvCovMeas*CombJac;
            % %Mdelta   = [ cos(theta)*cos(phi) cos(theta)*sin(phi)  -sin(theta);...
            % %            -sin(phi)  cos(phi)      0      ];
            % %deltaA   = Mdelta*[los-losM]
            % %CovVec     = InvCovMeas*deltaA;
            % %Matrix     = zeros(6);
            % %Matrix = Matrix + PhiFundamental'*HessianLos{10}(2:7,2:7)*PhiFundamental*CovVec(1);
            % %Matrix = Matrix + PhiFundamental'*HessianLos{11}(2:7,2:7)*PhiFundamental*CovVec(2);
            % %HessianSyst = HessianSyst + Matrix;
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Hessian = cell(7, 7);
            % Jacobian = [];
            % %LineCheck    = 129;                % xLos
            % [H_I, Jacob, Hess] = WorkOrder(KLagrange, 129);
            % Hessian{1}  = Hess;
            % Jacobian      = [Jacobian; Jacob];   % row 1
            % %LineCheck    = 130;                % yLos
            % [H_I, Jacob, Hess] = WorkOrder(KLagrange, 130);
            % Hessian{2}  = Hess;
            % Jacobian      = [Jacobian; Jacob];   % row 2
            % %LineCheck    = 107;                % zLos
            % [H_I, Jacob, Hess] = WorkOrder(KLagrange, 107);
            % Hessian{3}  = Hess;
            % Jacobian      = [Jacobian; Jacob];  % row 3
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
        delta      = Hessian\BetaInc
        Zhat     = Zhat - delta;
        TangentArray = [TangentArray, log(norm(BetaInc))]
        %ParamComp = [delta, Zhat, KeplerInit' + Zhat]

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