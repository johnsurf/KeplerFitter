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