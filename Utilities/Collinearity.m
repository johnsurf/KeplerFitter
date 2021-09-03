function [KeplerFit, ChiSq, CovFit] = Collinearity(track_data, KeplerInit, InvCovUV, units)

    rows       = size(track_data,1);

% For Choice of Units
twopi = units.twopi;
TU = units.TU;
DU = units.DU;
%VU = units.VU;
%AU = units.AU;
mu = units.mu;
%sqrtmu = units.sqrtmu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    delta      = zeros(6,1);
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

    %           deltaPar = zeros(6,1);
    %           delta     = zeros(6,1);
        delp      = zeros(6,1);

        IndexRows = find(oparm == 1);
        IndexCols = [IndexRows; Mind+1];

        for ipar=1:Mind
            if oparm(ipar) == true
                nfree = nfree+1;
            end
        end

        ChiSq        = 0.0;
        BetaInc      = zeros(6,1);
        %HessianSyst  = eye(6);
        HessianSyst  = zeros(6);
        %HessianSyst  = 0.0000001*eye(6);
        KeplerObj    = LagrangePlanetary(Kepler, units);

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
            % [Rsensor, Vsensor, ParamList, Jacobian, Hessian, GravityCan] = KeplerSat.OrbitDerivatives(tRecord, Sensor, Observed);
            % Sensor   = [Rsensor; Vsensor; zeros(3,1)];
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % [Rinit,Vinit,ParamList,Jacobian,Hessian,GravityCan] = KeplerObj.OrbitDerivatives(tRecord, Sensor, losM);
            % los       = [ParamList(7); ParamList(8); ParamList(9)];
            % los       = los/norm(los);
            % %
            % deltaM     =  1 - dot(los,losM);
            % JacobPred =  Jacobian(14:16,2:7);
            % JacFac    = -JacobPred'*losM;
            % BetaInc   =  BetaInc + JacFac;
            % %
            % % %Jac2nd    = JacFac*JacFac';
            % Hess      = zeros(6);
            % for ik = 1:3
            %     Hess      = Hess -losM(ik)*Hessian{13+ik}(2:7,2:7);
            % end
            % %Hess
            % HessianSyst = HessianSyst + Hess;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [Rinit,Vinit,ParamList] = OrbitAtTime(KeplerObj, tRecord, Sensor, losM, InvCovUV);
            % LineCheck    = 131 Collinearity
            [deltaM, JacFac, Hess] = WorkOrder(KeplerObj, 131);
            %[deltaM, deltaM1]
            %JacFac
            %JacFac1(2:7)'
            %Hess
            %Hess1(2:7,2:7)
            BetaInc     =  BetaInc    + JacFac(2:7)';
            HessianSyst = HessianSyst + Hess(2:7,2:7);
             
            %Hess = Hess*delta;
            %HessianSyst = HessianSyst + Hess + Jac2nd;

            % deltaRange    = ParamList(1) - track_data(ii,9)/DU;
            % if NR > 5
            %     deltaRSq      = 0.5*deltaRange^2;
            %     JacobianRange = Jacobian(8,2:7)*deltaRange;
            %     HessianRange  = Hessian{8}(2:7,2:7)*deltaRange +...
            %         Jacobian(8,2:7)'*Jacobian(8,2:7);
            %     BetaInc       = BetaInc + JacobianRange;
            %     HessianSyst   = HessianSyst + HessianRange;
            % end

            ChiSq    = ChiSq + abs(deltaM);

        end
        ChiSqArray   = [ChiSqArray, log(ChiSq)];
        %delta        = -(HessianSyst)\BetaInc
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
        % Put in a couple of rounds of Gradient Descent before the full
        % Newton-Raphson 2nd Order Method kicks in.
        if NR < 3
            delp = -0.15*BetaInc/rows;
        end
        delp
        % if norm(delp)<10^(-12)
        %     break
        % end
        Kepler =  Kepler + delp'

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
    
    KeplerFit = Kepler;
    CovFit    = inv(HessianSyst);
    % Adjust for Spec Values for cov_UV = 10^(-8)*eye*2(
    CovFit    = 10^(-8)*CovFit;

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
end
        