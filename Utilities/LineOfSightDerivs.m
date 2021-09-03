        function [LineOfSight, GradientLos] = LineOfSightDerivs(this, tBegin, tEnd, SatArray)
            %t            = time;
            e             = this.e;
            a             = this.a;
            Inclination   = this.Inclination;
            omega         = this.omega; 
            Mp            = this.Mp; 
            meanMotion    = this.meanMotion;
            Perturbation  = this.Pertubation;
            Inhomogeneous = this.Inhomogeneous;
            
            [Rextrap, Vextrap, Hessian, Gradient, Jacobian] = J6Gravity(0)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            Nbins            = 500;
            tInterval        = tEnd  - tBegin
            delt             = tInterval/Nbins;
            tRiem            = tBegin - delt/2.0;
            PhiInvIntegrated = zeros(6); 
            for i = 1:Nbins
                tRiem = tRiem + delt;
                PhiInvIntegrated = PhiInvIntegrated + expm(-Perturbation*tRiem);
            end
            PhiInvIntegrated  = PhiInvIntegrated*delt;
            PhiFundamental    = expm(Perturbation*tInterval);
            KPropagated       = (PhiFundamental*(PhiInvIntegrated*Inhomogeneous));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %KeplerProp = [e a Inclination omega Omega Mp + n*t]  % Canonical Units
            %M                    = meanMotion*tEnd + Mp;
            %Kepler               = [e a Inclination omega Omega M] + KPropagated';
            %KLagrangeNew         = LagrangePlanetary(Kepler);
            %[Rextrap, Vextrap, Hessian, Gradient, Jacobian] = KLagrangeNew.OrbitDerivatives(0)
            %Rextrap
            %Vextrap
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %KeplerProp = [e a Inclination omega Omega Mp + n*t]  % Canonical Units
            M                    = meanMotion*t + Mp;
            Kepler               = [e a Inclination omega Omega M];
            KLagrangeNew         = LagrangePlanetary(Kepler);
            %[Rextrap, Vextrap]   = KLagrange.extrapolate(0);
            %Rextrap
            [Rextrap, Vextrap, Hessian, Gradient, Jacobian] = KLagrangeNew.OrbitDerivatives(0);
            Delta = Jacobian(:,2:7)*KPropagated;           
            Rextrap + Delta(1:3)
            Vextrap + Delta(4:6)
        end