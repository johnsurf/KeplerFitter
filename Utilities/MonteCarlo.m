function KeplerMin = MonteCarlo(track_data) 
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

    for NR = 1:10000
        ChiSq     = 0.0;
        R  = rand(6,1);
        Kepler = [R(1), 2*R(2), pi*R(3), twopi*R(4), twopi*R(5), twopi*R(6)]; 
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
            %[Rsensor, Vsensor, Hessian, Gradient, Jacobian, JacobianLos, ParamList, GravityCan, HessianLos, RVHessSat] = KeplerSat.OrbitDerivatives(tRecord, Sensor, Observed);
            %Sensor   = [Rsensor; Vsensor; zeros(3,1)];

            [Rinit,Vinit,Hessian,Gradient,Jacobian,JacobianLos,ParamList,...
                GravityCan,HessianLos,Phi,MLagrange,RVHessian] = KeplerObj.OrbitDerivatives(tRecord, Sensor, Observed);
            los       = [ParamList(7); ParamList(8); ParamList(9)];
            los       = los/norm(los);

            delta     = los - losM;
            deltaSQ   = delta'*delta;
            ChiSq    = ChiSq + deltaSQ;                                
        end
        if ChiSq < ChiSqMin
            ChiSqMin  = ChiSq
            KeplerMin = Kepler
        end
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

end