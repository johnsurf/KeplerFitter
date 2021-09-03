        function [Jacob, JacobFinite, Hess, HessFinite] = DigitalJustice(this, I, Sensor, Observed)
        
            time          = this.time;
            e             = this.e;
            a             = this.a;
            Inclination   = this.Inclination;
            omega         = this.omega;
            Omega         = this.Omega;
            Mp            = this.Mp;
            Kepler        = [e, a, Inclination, omega, Omega, Mp];
            Kdelta        = LagrangePlanetary(Kepler);
            
            [H0, Jacob, Hess] = WorkOrder(this,I);
            %Hess
            
            % Do time Increment Separately
            
            epsilon     = 0.000000001;
            KeplerDelta = epsilon*eye(6);
            [~,~,~,~,~,~] = Kdelta.OrbitDerivatives(time+epsilon, Sensor, Observed);
            [HT, JacobT, HessT] = Kdelta.WorkOrder(I);
            
            Kepler1     = Kepler + KeplerDelta(1,:);
            Kdelta1     = LagrangePlanetary(Kepler1);
            [~,~,~,~,~,~] = Kdelta1.OrbitDerivatives(time, Sensor, Observed);
            [H1, Jacob1, Hess1] = Kdelta1.WorkOrder(I);
            
            Kepler2     = Kepler + KeplerDelta(2,:);
            Kdelta2     = LagrangePlanetary(Kepler2);
            [~,~,~,~,~,~] = Kdelta2.OrbitDerivatives(time, Sensor, Observed);
            [H2, Jacob2, Hess2] = Kdelta2.WorkOrder(I);
            
            Kepler3     = Kepler + KeplerDelta(3,:);
            Kdelta3     = LagrangePlanetary(Kepler3);
            [~,~,~,~,~,~] = Kdelta3.OrbitDerivatives(time, Sensor, Observed);
            [H3, Jacob3, Hess3] = Kdelta3.WorkOrder(I);
            
            Kepler4     = Kepler + KeplerDelta(4,:);
            Kdelta4     = LagrangePlanetary(Kepler4);
            [~,~,~,~,~,~] = Kdelta4.OrbitDerivatives(time, Sensor, Observed);
            [H4, Jacob4, Hess4] = Kdelta4.WorkOrder(I);
            
            Kepler5     = Kepler + KeplerDelta(5,:);
            Kdelta5     = LagrangePlanetary(Kepler5);
            [~,~,~,~,~,~] = Kdelta5.OrbitDerivatives(time, Sensor, Observed);
            [H5, Jacob5, Hess5] = Kdelta5.WorkOrder(I);
            
            Kepler6     = Kepler + KeplerDelta(6,:);
            Kdelta6     = LagrangePlanetary(Kepler6);
            [~,~,~,~,~,~] = Kdelta6.OrbitDerivatives(time, Sensor, Observed);
            [H6, Jacob6, Hess6] = Kdelta6.WorkOrder(I);
            
            JacobFinite = [HT-H0,H1-H0, H2-H0, H3-H0, H4-H0, H5-H0, H6-H0]/epsilon; 

            HessFinite = [(JacobT-Jacob);...
                          (Jacob1-Jacob);...
                          (Jacob2-Jacob);...
                          (Jacob3-Jacob);...
                          (Jacob4-Jacob);...
                          (Jacob5-Jacob);...
                          (Jacob6-Jacob)]/epsilon;
                      
        %GradientCheck = [];
        %JacobianCheck = [];
        %LagrangeCheck = [];
        %MLagrangeCheck = cell(1,6);
        %HessLosCheck   = cell(1,6);        
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
                      
        end