    function [StateVectorECI, RootsLaplace] = LaplaceInECI(state_vector_sat, state_vector_los)
        
        orbitparams = wgs84Constants;
        TU = orbitparams.TU;
        DU = orbitparams.DU;
        VU = orbitparams.VU;
        AU = orbitparams.AU;
        %mu = orbitparams.mu;
        mu = orbitparams.Grav;

        % Method 2 -- Laplace's Angles-Only Algorithm
        Rpos = state_vector_sat(1:3);
        Vpos = state_vector_sat(4:6);
        Apos = state_vector_sat(7:9);
        Tpos = state_vector_sat(10); 

        Rlos = state_vector_los(1:3);
        Vlos = state_vector_los(4:6);
        Alos = state_vector_los(7:9);

        D  = det([Rlos, 2*Vlos, Alos]);
        D1 = det([Rlos,   Vlos, Apos]);
        D2 = det([Rlos,   Vlos, Rpos]);
        D3 = det([Rlos,   Apos, Alos]);
        D4 = det([Rlos,   Rpos, Alos]);
        
        LdotR = dot(Rlos, Rpos);
        RSq   = dot(Rpos, Rpos); 
        
        c6    = (4.0*D1*(D*LdotR - D1) - D^2*RSq)/D^2;
        c3    =  4.0*mu*D2*(D*LdotR - 2.0*D1)/D^2;
        c1    = -(2*mu*D2/D)^2;
        C     = [1,0,c6,0,0,c3,0,0,c1];
        RootsLaplace = roots(C);
        %Index = 1:8;
        Index = find(imag(RootsLaplace) == 0);
        %Index = find(real(RootsLaplace(Index) > 0));
        StateVectorECI = cell(1);
        for i = 1:numel(Index)
            r    = RootsLaplace(Index(i));
            %r    = real(RootsLaplace(Index(i)))
            rho  = -2.0*(D1 + mu*D2/r^3)/D;
            rVec   = rho*Rlos + Rpos;
            rMag   = norm(rVec);
            rhodot = -(D3 + mu*D4/r^3)/D;
            vVec   = rhodot*Rlos + rho*Vlos + Vpos;
            aVec   = -mu*rVec/rMag^3;
            %[aVec ~, ~, ~] = Escobal(rVec);
            % Estimated Initial 9-State Vector from Laplace Method in ECI
            SV6  = [rVec; vVec; aVec; Tpos];
            StateVectorECI{i}  = SV6;
        end
        
        % %delr = 0.00001
        % %r    = 1.155;
        % delr = 0.01*DU;
        % r    = 1.0*DU;
        % for i = 1:1000
        %     r = r + delr;
        %     %Col3 = Alos + Rlos/r^3;
        %     %D    = det([Rlos, 2*Vlos, Col3]);
        %     rho  = -2.0*(D1 + mu*D2/r^3)/D;
        %     %F    = rho^2 + (2.0*rho*state_vector_los(1:3)' + state_vector_sat(1:3)')*state_vector_sat(1:3) - r^2
        %     F    = rho^2 - r^2 + 2.0*rho*LdotR + RSq;
        %     if F < 0
        %         r = r - delr;
        %         delr = delr/2;
        %     end
        % end
        % %r
        % rho  = -2.0*(D1 + mu*D2/r^3)/D;
        % rVec   = rho*Rlos + Rpos;
        % rMag   = norm(rVec);
        % rhodot = -(D3 + mu*D4/r^3)/D;
        % vVec   = rhodot*Rlos + rho*Vlos + Vpos;
        % %aVec   = -mu*rVec/rMag^3;
        % [aVec, ~, ~, ~] = Escobal(rVec);
        % % Estimated Initial 9-State Vector from Laplace Method in ECI
        % %LoopSearch = [rVec; vVec; aVec];
        
    end