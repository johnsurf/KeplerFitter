classdef KeplerOrbit < handle
    
    properties (Constant = true)
        % Constants
        %extract conversion factors and ititialize to 1 in case they are not used
        %twopi  = 2.0*pi;         % 2*pi
        %XMNPDA = 1440.0;         % 1440 minutes per day
        %XKE    = 0.0743669161;   % Hoot's ke Gravitational constant for Earth in RE^(1.5)/minutes
    end
    
    properties (SetAccess = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        time;
        Rpos;
        Rdot;
        Derivs;
        ClassicalElements;
        EquinoctiaElements;
        JECI2Kepler;
        JECI2Equinoctial;
        Hessian;
        DataList;
        JacobianAll;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    methods
        
        % Constructor:
        function this = KeplerOrbit(time, Rpos, Rdot)
            
            wgs84     = wgs84Constants();
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Canonical Coordinates assumed throughout -- length in DU, time in TU
            % Modifications can be made to the choice of units in wgs84Constants!
            % Orbital elements e, p, I, omega, Omega, PTime
            % or lines [17, 24, 64, 50, 47, 62]
            % and the parameters rx, ry, rz, vx, vy, vz
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Constants
            orbitparams = wgs84Constants;
            twopi       = orbitparams.twopi;
            TU          = orbitparams.TU;  % Canonical Time UnitTime Unit in Seconds
            DU          = orbitparams.DU;  % Canonical Distance Unit
            VU          = orbitparams.VU;  % Canonical Velocity Unit
            AU          = wgs84.AU;        % Canonical Acceleration Unit
            mu          = orbitparams.mu;
            sqrtmu      = sqrt(mu);
            
            Ivec  = [1.0 0.0 0.0];
            Jvec  = [0.0 1.0 0.0];
            Kvec  = [0.0 0.0 1.0];
            
            this.time = time;
            this.Rpos = Rpos;
            this.Rdot = Rdot;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            t    = time;                        %  1) t
            rx   = Rpos(1);                     %  2) Rx
            ry   = Rpos(2);                     %  3) Ry
            rz   = Rpos(3);                     %  4) Rz
            vx   = Rdot(1);                     %  5) Vx
            vy   = Rdot(2);                     %  6) Vy
            vz   = Rdot(3);                     %  7) Vz
            
            rmag = sqrt(rx*rx + ry*ry + rz*rz); % 10) rmag = norm(Rpos);
            vsq  = vx*vx + vy*vy + vz*vz;       % 11) vsq  = dot(Rdot,Rdot);
            er   = vsq - mu/rmag;               % 12) er   = v^2 - mu/rmag;
            ev   = rx*vx + ry*vy + rz*vz;       % 13) ev   = dot(Rpos,Rdot);
            % Evec = er*Rpos - ev*Vpos;
            ex   = (er*rx - ev*vx)/mu;               % 14) Evec(1)
            ey   = (er*ry - ev*vy)/mu;               % 15) Evec(2)
            ez   = (er*rz - ev*vz)/mu;               % 16) Evec(3)
            Evec = [ex ey ez];
            % Kepler 1:
            e    = sqrt(ex*ex + ey*ey + ez*ez); % 17) e    = norm(Evec);
            % Pvec = Evec/e;
            Px   = ex/e;                        % 18) Pvec(1)
            Py   = ey/e;                        % 19) Pvec(2)
            Pz   = ez/e;                        % 20) Pvec(3)
            Pvec = [Px  Py  Pz];
            % hVec = cross(Rpos, Rdot);
            hx   = ry*vz - rz*vy;               % 21) hVec(1)
            hy   = rz*vx - rx*vz;               % 22) hVec(2)
            hz   = rx*vy - ry*vx;               % 23) hVec(3)
            hVec = [hx hy hz];
            % p    = dot(hVec,hVec)/mu;
            p    = (hx*hx + hy*hy + hz*hz)/mu;  % 24) p = hVec^2/mu
            hmag = sqrt(mu*p);
            % Wvec = hVec/norm(hVec);
            Wx   = hx/hmag;                  % 25) Wvec(1)
            Wy   = hy/hmag;                  % 26) Wvec(2)
            Wz   = hz/hmag;                  % 27) Wvec(3)
            Wvec = [Wx Wy Wz];
            % Qvec = cross(Wvec, Pvec);
            Qx   = Wy*Pz - Wz*Py;               % 28) Qvec(1)
            Qy   = Wz*Px - Wx*Pz;               % 29) Qvec(2)
            Qz   = Wx*Py - Wy*Px;               % 30) Qvec(3)
            Qvec = [Qx Qy Qz];
            % Dot1 = Pvec'*Qvec
            % Dot2 = Pvec'*Wvec
            % Dot2 = Qvec'*Wvec
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Other Vectors
            rUnitx = rx/rmag;                   % 31) rUnit(1)
            rUnity = ry/rmag;                   % 32) rUnit(2)
            rUnitz = rz/rmag;                   % 33) rUnit(3)
            rUnit  = [rUnitx rUnity rUnitz];
            Ax     = -Wy;                       % 34)
            Ay     =  Wx;                       % 35)
            normN  = sqrt(Ax*Ax + Ay*Ay);       % 36)
            Nx     = Ax/normN;                  % 37) Nvec(1)
            Ny     = Ay/normN;                  % 38) Nvec(2)
            Nvec   = [Nx  Ny 0.0];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Semi-Major Axis:
            onePe = 1.0 + e;                    % 39)
            oneMe = 1.0 - e;                    % 40)
            fac   = onePe*oneMe;                % 41)
            if (fac > 0.0)
                rootfac = sqrt(fac);            % 42)
            else
                rootfac = 0;                    %
            end
            a = p/fac;                          % 43)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Period
            PeriodTU   = twopi*a^(1.5)/sqrtmu;  % 44) Period in TU's
            % Period     = PeriodTU*time_unit;    % 45) Period in Seconds
            Period     = PeriodTU;                % 45) Work in TU's for now
            meanMotion = twopi/Period;            % 46) Radians per TU.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Kepler 3) Omega
            % Longitude of the Ascending Node
            Omega = atan2(Ny, Nx);              % 47) Omega = atan(Ny/Nx);
            sinOmega = sin(Omega);
            cosOmega = cos(Omega);
            while Omega<0.0
                Omega = Omega + twopi;
            end
            while Omega>twopi
                Omega = Omega - twopi;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Kepler 4) omega
            % Argument of Periapsis
            cosPeriapsis = dot(Nvec, Pvec);
            cosP =Nx*Px + Ny*Py;                           % 48) cos(omega)
            sinPeriapsis = dot(Wvec, cross(Nvec, Pvec));
            sinP =(Wx*Ny - Wy*Nx)*Pz + (Nx*Py - Ny*Px)*Wz; % 49) sin(omega)
            
            omega=atan2(sinP, cosP);                       % 50) omega = atan(sinP/cosP)
            while omega<0.0
                omega = omega + twopi;
            end
            while omega>twopi
                omega = omega - twopi;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % True Anomaly at Epoch
            cosnu1   = dot(rUnit,Pvec);
            cosnu    = rUnitx*Px + rUnity*Py + rUnitz*Pz;  % 51) cos(nu)
            sinnu1   = dot(Wvec,(cross(Pvec,rUnit)));
            sinnu    = Wx*(Py*rUnitz - Pz*rUnity) + ...
                Wy*(Pz*rUnitx - Px*rUnitz) + ...
                Wz*(Px*rUnity - Py*rUnitx);         % 52) sin(nu)
            % rUnit   = rUnitx^2 + rUnity^2 + rUnitz^2
            % Pmag    = Px^2 + Py^2 + Pz^2
            % trig    = cosnu*cosnu + sinnu*sinnu
            nuEpoch = atan2(sinnu, cosnu);                 % 53) nu = atan(sinnu/cosnu)
            while(nuEpoch <   0.0)
                nuEpoch = nuEpoch + twopi;
            end
            while(nuEpoch > twopi)
                nuEpoch = nuEpoch - twopi;
            end
            % nuEpoch
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Eccentric Anomaly at Epoch
            %     cosE    = ( e + cosnu )/(1.0 + e*cosnu);      % 54) cos(E)
            %     sinE    = sinnu*rootfac/(1.0 + e*cosnu);      % 55) sin(E)
            cosT    = e + cosnu;                          % 54) ~cos(E)
            sinT    = sinnu*rootfac;                      % 55) ~sin(E)
            EccentricAnomalyEpoch = atan2(sinT, cosT);    % 56) E = atan(sinE/cosE)
            cosE    = cos(EccentricAnomalyEpoch);         % 57  (old 54)
            sinE    = sin(EccentricAnomalyEpoch);         % 58  (old 55)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Kepler 5: Periapsis Time
            % Find the Time at Periapsis
            % Solve Kepler's Equation for the MeanAnomaly at Epoch
            % -- this side is easy!
            MeanAnomalyEpoch                = EccentricAnomalyEpoch - e*sinE;% 59)
            TimeSincePeriapsis              = MeanAnomalyEpoch/meanMotion;   % 60)
            % PTime                         = TimeSincePeriapsis;
            DeltaTime                       = t - TimeSincePeriapsis;     % 61)
            %FractionOfPeriodSincePeriapsis = TimeSincePeriapsis/Period;
            % DeltaTime is in the intervale [-T/2, T/2]
            % If DeltaTime > 0 then DeltaTime is in [0, T/2], then do nothing
            % If DeltaTime in [-T/2, 0] then add a Period DeltaTime -> [T/2,T]
            % If you see PTime > Period/2, then above if-statement was executed.
            PTime = DeltaTime;  % DEFAULT                     % 62) in TUs
            if DeltaTime < 0
                % The OVERWRITE PTime HERE:
                PTime = DeltaTime + Period;
            end
            %Tnext     = PTime*time_unit;                    % 63 in Seconds
            Mp = MeanAnomalyEpoch - meanMotion*t;            % 63 in radians
            while Mp < 0
                Mp = Mp + twopi;
            end
            while Mp > twopi
                Mp = Mp - twopi;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Kepler 6: Inclination
            Inclination = acos(Wz);                          % 64)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Test Section
            % Okay given the TimeSincePeriapsis -- can we go back to the true Anomaly?
            % // 1)   Find Eccentric Anomaly from Mean Anomaly
            % //      Solve Kepler's Equation:  M = E - e*sin(E);
            % //      write as E = M + e*sin(E) and iterate 10 times.
            %
            % //  Brute Force iteration (good way to seed Newton's method which follows)
            % Here we will insert the time and proceed to calculate
            % time > MeanAnomaly > Eccentric Anomaly > True Anomaly > Position & Velocity
            % and the calculate the Terms in Kepler, J2-J6 Perturbations using the
            % linearization obtained from the Taylor expansion of the Gravitational
            % Potential from Escobal pages 48-52
            % The time variable is set up so that time=0 corresponds to the Epoch Time
            % We can always use Time to propagate before (time < 0) or after (time > 0)
            % the Epoch Time -- but this requires solving Kepler's Equation.
            % JRS TEST Use for Extrapolation plus time t:
            %M = meanMotion*t + MeanAnomalyEpoch;                             % 65)
            % JRS TEST: Uses the Current MeanAnomaly at this event EPOCH.
            M = MeanAnomalyEpoch;                                             % 65)
            %  meanAnomaly M to eccentric Anomaly EM
            EM = M + e;
            if (M > pi)
                EM = M - e;
            elseif (-pi < M && M < 0)
                EM = M - e;
            end
            %         //cout << " Input meanAnomaly = " << meanAnomaly << endl;
            for i=1:1:10
                EM  = M + e*sin(EM);
                %//cout << " Mean Anomaly Solution " << E << endl;
            end
            % //      10 rounds of Newton's root finding method based on the above "seed".
            for i=1:1:10
                Eprime      = 1.0 - e*cos(EM);
                EM          = EM + (M - EM + e*sin(EM))/Eprime;
            end
            % Solve Kepler's Equation for E:  M = EM - e*sin(EM);   % 66) EM = E(e,M)
            % Need to differentiate Implicitly here!
            eDenom   = Eprime*Eprime*Eprime;
            
            % KeplerInv = EM;
            % //cout  << " Newton-Raphson Kepler's Equation Eccentric Anomaly =" << E << endl;
            % eDenom   = Eprime*Eprime*Eprime;
            
            %rmag = p/(1 + e*cosnu);
            %rRecon = rmag*(cosnu*Pvec + sinnu*Qvec);
            %QReco = rmag*cosnu*Rdot + (1.0/sqrt(p))*sinnu*Rpos;
            %QReco  = (1.0/norm(QReco))*QReco;
            
            cosK  = cos(EM);                                        % 67)
            sinK  = sin(EM);                                        % 68)
            %E     = atan2(sinK,cosK)
            %%%%%%%%%%%%%%%%%%% Implicit Function Theorem Derivatives
            dE_dM    = 1.0/Eprime;
            dE_de    = sinK/Eprime;
            d2E_dMdM = -e*sinK/eDenom;
            %// Mike Cain Corrections!
            d2E_dMde = (cosK - e)/eDenom;
            d2E_dedM = (cosK - e)/eDenom;
            d2E_dede = ((2.0 - e*cosK)*cosK*sinK - e*sinK)/eDenom;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % kfac  = 1.0 - e*cosK;                                 % 69)
            % coss  = (cosK - e)/kfac;                              % 70)
            % sins  = rootfac*sinK/kfac;                            % 71)
            % 	s     = atan2(sins, coss);                          % 72
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   % for convenience -- you can eliminate kfac;
            tanX  = cosK - e;                                     % line 70
            tanY  = rootfac*sinK;                                 % line 71
            nu     = atan2(tanY, tanX);                           % line 72
            coss  = cos(nu);                                      % line 73
            sins  = sin(nu);                                      % line 74
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Position and Unit Vector along the Position
            %Kvec  = a*kfac*(cosnu*Pvec + sinnu*Qvec);
            %Runit = Kvec/norm(Kvec);
            
            %   Unit Vector along the Velocity
            %Vunit =(-sinnu*Pvec + (e + cosnu)*Qvec);
            %Vunit = Vunit/norm(Vunit);
            
            %   Unit Vector out of the R-V plane
            %Wlocal = cross(Runit, Vunit);
            %Wlocal = Wlocal/norm(Wlocal);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Let us propagate this Orbit to time = t
            rorbit   = p/(1.0 + e*coss);                          % line 75
            rorbitx  = rorbit*(coss*Px + sins*Qx);                % line 76
            rorbity  = rorbit*(coss*Py + sins*Qy);                % line 77
            rorbitz  = rorbit*(coss*Pz + sins*Qz);                % line 78
            rtpinv   = sqrt(mu/p);
            vorbitx  = rtpinv*(-sins*Px + (e + coss)*Qx);         % line 79
            vorbity  = rtpinv*(-sins*Py + (e + coss)*Qy);         % line 80
            vorbitz  = rtpinv*(-sins*Pz + (e + coss)*Qz);         % line 81
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Rextrap  = [rorbitx rorbity rorbitz];
            Vextrap  = [vorbitx vorbity vorbitz];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Aeq      = omega + Omega;                             % line 82
            cosAeq   = cos(Aeq);                                  % line 83
            sinAeq   = sin(Aeq);                                  % line 84
            tanHalf  = tan(Inclination/2.0);                      % line 85
            Feq      = e*cosAeq;                                  % line 86
            Geq      = e*sinAeq;                                  % line 87
            Heq      = tanHalf*cosOmega;                          % line 88
            Keq      = tanHalf*sinOmega;                          % line 89
            Leq      = Aeq + nuEpoch;                             % line 90
            CosL     = cos(Leq);
            SinL     = sin(Leq);
            
            alphaSq = Heq^2 - Keq^2;
            Seq     = 1.0 + Heq^2 + Keq^2;
            Weq     = 1.0 + Feq*CosL + Geq*SinL;
            Req     = p/Weq;
            RovS    = Req/Seq;
            srtpinv = rtpinv/Seq;
            
            Xeq     = RovS*((1.+alphaSq)*CosL + 2.0*Heq*Keq*SinL);
            Yeq     = RovS*((1.-alphaSq)*SinL + 2.0*Heq*Keq*CosL);
            Zeq     = 2.0*RovS*(Heq*SinL - Keq*CosL);
            Req     = [Xeq, Yeq, Zeq];
            Rpos';
            
            VXeq     = -srtpinv*((1.+alphaSq)*SinL - 2.0*Heq*Keq*CosL + Geq - 2.0*Feq*Heq*Keq+alphaSq*Geq);
            VYeq     = -srtpinv*((alphaSq-1.0)*CosL + 2.0*Heq*Keq*SinL - Feq + 2.0*Geq*Heq*Keq+alphaSq*Feq);
            VZeq     = 2.0*srtpinv*(Heq*CosL + Keq*SinL + Feq*Heq + Geq*Keq);
            Veq      = [VXeq, VYeq, VZeq];
            Rdot';
            
            ClassicalElements       = [e, a, Inclination, omega, Omega, Mp];
            EquinoctiaElements      = [p, Feq, Geq, Heq, Keq, Leq];
            this.ClassicalElements  = ClassicalElements;
            this.EquinoctiaElements = EquinoctiaElements;
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Checking First Derivatives
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % format short
            % iStart = 1;
            % CombinedDifferences = 0.0;
            % for LineCheck = iStart:90
            %     ForDerivs  = [];
            %     %BackDerivs = [];
            %     %Derivs     = [];
            %     %Combined   = [];
            %
            %     F            = zeros(90,1);
            %     F(LineCheck) = 1;
            %     F            = backdiff(F, LineCheck);
            %     % BackDerivs   = [F(1) F(2) F(3) F(4) F(5) F(6) F(7)];
            %     % Combined     = [Index F];
            %     % Kate's Combined'
            %
            %     forward      = [];
            %     for jk = 1:90 % check range of derivatives wrt forwards differentiation
            %         %jk
            %         Q        = zeros(90,1);
            %         Q(jk)     = 1;
            %         Q        = fordiff(Q, LineCheck);
            %         % Q(LineCheck) = dK_LineCheck/dLj
            %         %forward  = [forward; j Q(LineCheck)];
            %         forward  = [forward; jk+1-iStart Q(LineCheck)];
            %         %Derivs   = [Derivs; Q(LineCheck)];
            %         %Combined = [Combined Q];
            %     end
            %     ForDerivs    = [ForDerivs; forward];
            %     Derivs       = [ForDerivs F (F-ForDerivs(:,2))];
            %     Derivs(iStart:LineCheck,:)
            %     Difference   = sum(F-ForDerivs(:,2));
            %     CombinedDifferences = CombinedDifferences + Difference;
            %     disp([' Line = ', num2str(LineCheck),'  Sum of Differences = ',num2str(Difference)]);
            %     % Kate's ForDerivs'
            %     % Derivs      = [F Derivs (F-Derivs)];
            %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %     %     ForDerivs    = [ForDerivs; forward];
            %     %     Derivs = [ForDerivs F (F-ForDerivs(:,2))];
            %     %     Derivs = Derivs(1:LineCheck,:)
            %     %     %Derivs = [F Derivs (F-Derivs)];
            %     % 'Debug Point'
            % end
            % disp([' Combined Sum of Differences = ',num2str(CombinedDifferences)]);
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            JacobianAll = [];
            for I = 1:90
                LineCheck    = I;                % e
                F            = zeros(90,1);
                F(LineCheck) = 1;
                F            = backdiff(F, LineCheck);
                JacobianAll  = [JacobianAll; F(2) F(3) F(4) F(5) F(6) F(7)];
            end
            this.JacobianAll = JacobianAll;
            % We will want to construct the 6 by 6 Jacobian of the Kepler
            % Orbital elements e, p, I, omega, Omega, PTime
            % or lines [17, 24, 64, 50, 47, 62]
            % and the parameters rx, ry, rz, vx, vy, vz
            
            JECI2Kepler = [];
            LineCheck    = 17;                % e
            F            = zeros(90,1);
            F(LineCheck) = 1;
            F            = backdiff(F, LineCheck);
            JECI2Kepler     = [JECI2Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 1
            
            LineCheck    = 43;                % a
            F            = zeros(90,1);
            F(LineCheck) = 1;
            F            = backdiff(F, LineCheck);
            JECI2Kepler     = [JECI2Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 2
            
            LineCheck    = 64;                % I
            F            = zeros(90,1);
            F(LineCheck) = 1;
            F            = backdiff(F, LineCheck);
            JECI2Kepler     = [JECI2Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 3
            
            LineCheck    = 50;                % omega
            F            = zeros(90,1);
            F(LineCheck) = 1;
            F            = backdiff(F, LineCheck);
            JECI2Kepler     = [JECI2Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 4
            
            LineCheck    = 47;                % Omega
            F            = zeros(90,1);
            F(LineCheck) = 1;
            F            = backdiff(F, LineCheck);
            JECI2Kepler     = [JECI2Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 5
            
            LineCheck    = 63;                % Mp in Radians
            F            = zeros(90,1);
            F(LineCheck) = 1;
            F            = backdiff(F, LineCheck);
            JECI2Kepler     = [JECI2Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 6
            
            this.JECI2Kepler = JECI2Kepler;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            JECI2Equinoctial  = [];
            LineCheck    = 24;                % p
            F            = zeros(90,1);
            F(LineCheck) = 1;
            F            = backdiff(F, LineCheck);
            JECI2Equinoctial     = [JECI2Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 2
            
            LineCheck    = 86;                % f
            F            = zeros(90,1);
            F(LineCheck) = 1;
            F            = backdiff(F, LineCheck);
            JECI2Equinoctial     = [JECI2Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 2
            
            LineCheck    = 87;                % g
            F            = zeros(90,1);
            F(LineCheck) = 1;
            F            = backdiff(F, LineCheck);
            JECI2Equinoctial     = [JECI2Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 3
            
            LineCheck    = 88;                % h
            F            = zeros(90,1);
            F(LineCheck) = 1;
            F            = backdiff(F, LineCheck);
            JECI2Equinoctial     = [JECI2Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 4
            
            LineCheck    = 89;                % k
            F            = zeros(90,1);
            F(LineCheck) = 1;
            F            = backdiff(F, LineCheck);
            JECI2Equinoctial     = [JECI2Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 5
            
            LineCheck    = 90;                % MLin Radians
            F            = zeros(90,1);
            F(LineCheck) = 1;
            F            = backdiff(F, LineCheck);
            JECI2Equinoctial     = [JECI2Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 6
            
            this.JECI2Equinoctial = JECI2Equinoctial;
            % OrbitalCov   = Jacobian*Cov*Jacobian';  % transformed Covariances
            % 'Debug Point'
            
            this.DataList    = zeros(90,1);
            
            this.DataList( 1) = t;                      %  1) t
            this.DataList( 2) = rx;                     %  2) Rx
            this.DataList( 3) = ry;                     %  3) Ry
            this.DataList( 4) = rz;                     %  4) Rz
            this.DataList( 5) = vx;                     %  5) Vx
            this.DataList( 6) = vy;                     %  6) Vy
            this.DataList( 7) = vz;                     %  7) Vz
            
            this.DataList(10) = rmag;                   % 10) rmag
            this.DataList(11) = vsq;                    % 11) vsq
            this.DataList(12) = er;                     % 12) er
            this.DataList(13) = ev;                     % 13) ev
            this.DataList(14) = ex;                     % 14) Evec(1)
            this.DataList(15) = ey;                     % 15) Evec(2)
            this.DataList(16) = ez;                     % 16) Evec(3)
            this.DataList(17) = e;                      % 17) e;
            this.DataList(18) = Px;                     % 18) Pvec(1)
            this.DataList(19) = Py;                     % 19) Pvec(2)
            this.DataList(20) = Pz;                     % 20) Pvec(3)
            this.DataList(21) = hx;                     % 21) hVec(1)
            this.DataList(22) = hy;                     % 22) hVec(2)
            this.DataList(23) = hz;                     % 23) hVec(3)
            this.DataList(24) = p;;                     % 24) p = hVec^2/mu
            this.DataList(25) = Wx;                     % 25) Wvec(1)
            this.DataList(26) = Wy;                     % 26) Wvec(2)
            this.DataList(27) = Wz;                     % 27) Wvec(3)
            this.DataList(28) = Qx;                     % 28) Qvec(1)
            this.DataList(29) = Qy;                     % 29) Qvec(2)
            this.DataList(30) = Qz;                     % 30) Qvec(3)
            this.DataList(31) = rUnitx;                 % 31) rUnit(1)
            this.DataList(32) = rUnity;                 % 32) rUnit(2)
            this.DataList(33) = rUnitz;                 % 33) rUnit(3)
            this.DataList(34) = Ax;                     % 34) Ax
            this.DataList(35) = Ay;                     % 35) Ay
            this.DataList(36) = normN;                  % 36) normN
            this.DataList(37) = Nx;                     % 37) Nvec(1)
            this.DataList(38) = Ny;                     % 38) Nvec(2)
            this.DataList(39) = onePe;                  % 39)
            this.DataList(40) = oneMe;                  % 40)
            this.DataList(41) = fac;                    % 41)
            this.DataList(42) = rootfac;                % 42)
            this.DataList(43) = a;                      % 43)
            this.DataList(44) = PeriodTU;               % 44) Period in TU's
            this.DataList(45) = Period;                 % 45) Work in TU's for now
            this.DataList(46) = meanMotion;             % 46) Radians per TU.
            this.DataList(47) = Omega;                  % 47) Omega = atan(Ny/Nx);
            this.DataList(48) = cosP;                   % 48) cos(omega)
            this.DataList(49) = sinP;                   % 49) sin(omega)
            this.DataList(50) = omega;                  % 50) omega = atan(sinP/cosP)
            this.DataList(51) = cosnu;                  % 51) cos(nu)
            this.DataList(52) = sinnu;                  % 52) sin(nu)
            this.DataList(53) = nuEpoch;                % 53) nu = atan(sinnu/cosnu)
            this.DataList(54) = cosT;                   % 54) ~cos(E)
            this.DataList(55) = sinT;                   % 55) ~sin(E)
            this.DataList(56) = EccentricAnomalyEpoch;  % 56) E = atan(sinE/cosE)
            this.DataList(57) = cosE;                   % 57  (old 54)
            this.DataList(58) = sinE;                   % 58  (old 55)
            this.DataList(59) = MeanAnomalyEpoch;       % 59)
            this.DataList(60) = TimeSincePeriapsis;     % 60)
            this.DataList(61) = DeltaTime;              % 61)
            this.DataList(62) = PTime;  % DEFAULT       % 62) in TUs
            this.DataList(63) = Mp;                     % 63 in radians
            this.DataList(64) = Inclination;            % 64)
            this.DataList(65) = M;                      % 65)
            this.DataList(66) = EM;                     % 66) EM = E(e,M)
            this.DataList(67) = cosK;                   % 67)
            this.DataList(68) = sinK;                   % 68)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            this.DataList(70) = tanX;                   % 70)
            this.DataList(71) = tanY;                   % 71)
            this.DataList(72) = nu;                     % 72)
            this.DataList(73) = coss;                   % 73)
            this.DataList(74) = sins;                   % 74)
            this.DataList(75) = rorbit;                 % 75)
            this.DataList(76) = rorbitx;                % 76)
            this.DataList(77) = rorbity;                % 77)
            this.DataList(78) = rorbitz;                % 78)
            this.DataList(79) = vorbitx;                % 79)
            this.DataList(80) = vorbity;                % 80)
            this.DataList(81) = vorbitz;                % 81)
            this.DataList(81) = Aeq;                    % line 82
            this.DataList(81) = cosAeq;                 % line 83
            this.DataList(81) = sinAeq;                 % line 84
            this.DataList(81) = tanHalf;                % line 85
            this.DataList(81) = Feq;                    % line 86
            this.DataList(81) = Geq;                    % line 87
            this.DataList(81) = Heq;                    % line 88
            this.DataList(81) = Keq;                    % line 89
            this.DataList(81) = Leq;                    % line 90
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            function F = backdiff( f, LineCheck)
                % eliminate some repetitions:
                % rtpinv  = sqrt(1.0/p);
                %//////////////////////////////////////////////////////////////////////
                if LineCheck > 64
                    % from line 90: Leq      = Aeq + nuEpoch;         % line 90
                    f(53) = f(53) + f(90);
                    f(82) = f(82) + f(90);
                    % from line 89: Keq      = tanHalf*sin(Omega);    % line 89
                    f(85) = f(85) + f(89)*sinOmega;
                    f(47) = f(47) + f(89)*tanHalf*cosOmega;
                    % from line 88: Heq      = tanHalf*cos(Omega);    % line 88
                    f(85) = f(85) + f(88)*cosOmega;
                    f(47) = f(47) - f(88)*tanHalf*sinOmega;
                    % from line 87: Geq      = e*sinAeq;              % line 87
                    f(17) = f(17) + f(87)*sinAeq;
                    f(84) = f(84) + f(87)*e;
                    % from line 86: Feq      = e*cosAeq;              % line 86
                    f(17) = f(17) + f(86)*cosAeq;
                    f(83) = f(83) + f(86)*e;
                    % from line 85: tanHalf  = tan(Inclination/2.0);  % line 85
                    f(64) = f(64) + 0.5*f(85)*(1.0 + tanHalf^2);
                    % from line 84: sinAeq   = sin(Aeq);              % line 84
                    f(82) = f(82) + f(84)*cosAeq;
                    % from line 83: cosAeq   = cos(Aeq);              % line 83
                    f(82) = f(82) - f(83)*sinAeq;
                    % from line 82: Aeq      = omega + Omega;         % line 82
                    f(47) = f(47) + f(82);
                    f(50) = f(50) + f(82);
                    % Computation of backwards derivative starting at the bottom of the algorithm
                    %////////////  /////Block 81-75////////////////////////////////////////
                    % from line 81: double vorbitz  = sqrt(mu/p)*(-sins*Pz + (e + coss)*Qz);
                    f(74) = f(74) - f(81)*rtpinv*Pz;                 % sins    #74
                    f(73) = f(73) + f(81)*rtpinv*Qz;                 % coss    #73
                    f(30) = f(30) + f(81)*rtpinv*(e + coss);         % Qz      #30
                    f(24) = f(24) - f(81)*rtpinv*(0.5/p)*(-sins*Pz + (e + coss)*Qz); % p #24
                    f(20) = f(20) - f(81)*rtpinv*sins;               % Pz      #20
                    f(17) = f(17) + f(81)*rtpinv*Qz;                 % e       #17
                    % from line 80: double vorbity  = sqrt(mu/p)*(-sins*Py + (e +coss)*Qy);
                    f(74) = f(74) - f(80)*rtpinv*Py;                 % sins    #74
                    f(73) = f(73) + f(80)*rtpinv*Qy;                 % coss    #73
                    f(29) = f(29) + f(80)*rtpinv*(e + coss);         % Qy      #29
                    f(24) = f(24) - f(80)*rtpinv*(0.5/p)*(-sins*Py + (e + coss)*Qy); % p #24
                    f(19) = f(19) - f(80)*rtpinv*sins;               % Py      #19
                    f(17) = f(17) + f(80)*rtpinv*Qy;                 % e       #17
                    % from line 79: double vorbitx  = sqrt(mu/p)*(-sins*Px + (e + coss)*Qx);
                    f(74) = f(74) - f(79)*rtpinv*Px;                 % sins    #74
                    f(73) = f(73) + f(79)*rtpinv*Qx;                 % coss    #73
                    f(28) = f(28) + f(79)*rtpinv*(e + coss);         % Qx      #28
                    f(24) = f(24) - f(79)*rtpinv*(0.5/p)*(-sins*Px + (e + coss)*Qx); % p #24
                    f(18) = f(18) - f(79)*rtpinv*sins;               % Px      #18
                    f(17) = f(17) + f(79)*rtpinv*Qx;                 % e       #17
                    % from line 78: double rorbitz  = rorbit*(coss*Pz + sins*Qz);
                    f(75) = f(75) + f(78)*(coss*Pz + sins*Qz);            % rorbit  #75
                    f(74) = f(74) + f(78)*rorbit*Qz;                      % sins    #74
                    f(73) = f(73) + f(78)*rorbit*Pz;                      % coss    #73
                    f(30) = f(30) + f(78)*rorbit*sins;                    % Qz      #30
                    f(20) = f(20) + f(78)*rorbit*coss;                    % Pz      #20
                    % from line 77: double rorbity  = rorbit*(coss*Py + sins*Qy);
                    f(75) = f(75) + f(77)*(coss*Py + sins*Qy);            % rorbit  #75
                    f(74) = f(74) + f(77)*rorbit*Qy;                      % sins    #74
                    f(73) = f(73) + f(77)*rorbit*Py;                      % coss    #73
                    f(29) = f(29) + f(77)*rorbit*sins;                    % Qy      #29
                    f(19) = f(19) + f(77)*rorbit*coss;                    % Py      #19
                    % from line 76: double rorbitx  = rorbit*(coss*Px + sins*Qx);
                    f(75) = f(75) + f(76)*(coss*Px + sins*Qx);            % rorbit  #75
                    f(74) = f(74) + f(76)*rorbit*Qx;                      % sins    #74
                    f(73) = f(73) + f(76)*rorbit*Px;                      % coss    #73
                    f(28) = f(28) + f(76)*rorbit*sins;                    % Qx      #28
                    f(18) = f(18) + f(76)*rorbit*coss;                    % Px      #18
                    % from line 75: double rorbit   = p/(1.0 + e*coss);
                    denom = 1.0/(1.0 + e*coss);
                    f(73) = f(73) - f(75)*p*e*denom*denom;     % coss  #73
                    f(24) = f(24) + f(75)*denom;               % p     #24
                    f(17) = f(17) - f(75)*p*coss*denom*denom;  % e     #17
                    %//////////////////////////////////////////////////////////////////////
                    %     Block 74-72
                    % from line 74: sins  = sin(nu);
                    f(72) = f(72) + f(74)*coss;                           % nu      #72
                    % from line: 73: coss  = cos(nu);
                    f(72) = f(72) - f(73)*sins;                           % nu      #72
                    %//////////////////////////////////////////////////////////////////////
                    % from line: 72: nu     = atan2(tanY, tanX);
                    % JRS: Comment out following two lines to take
                    % (nu, rx, ry, rx, vx, vy, vz) as independent variables
                    f(71) = f(71) + f(72)*tanX/(tanX*tanX + tanY*tanY);   % tanY    #71
                    f(70) = f(70) - f(72)*tanY/(tanX*tanX + tanY*tanY);   % tanX    #70
                    %/////////////////////Block 71-65//////////////////////////////////////
                    % from line 71: sins  = rootfac*sinK/kfac;
                    %  f(69) = f(69) - f(71)*rootfac*sinK/(kfac*kfac);       % kfac    #69
                    %  f(68) = f(68) + f(71)*rootfac/kfac;                   % sinK    #68
                    %  f(42) = f(42) + f(71)*sinK/kfac;                      % rootfac #42
                    % from line 70: coss = (cosK - e)/kfac;
                    %  f(69) = f(69) - f(70)*(cosK - e)/(kfac*kfac);         % kfac    #69
                    %  f(67) = f(67) + f(70)/kfac;                           % cosK    #67
                    %  f(17) = f(17) - f(70)/kfac;                           % e       #17
                    % from line 71: tanY  = rootfac*sinK;
                    f(68) = f(68) + f(71)*rootfac;                         % sinK    #68
                    f(42) = f(42) + f(71)*sinK;                            % rootfac #42
                    % from line 70: tanX  = cosK - e;
                    f(67) = f(67) + f(70);                                 % cosK    #67
                    f(17) = f(17) - f(70);                                 % e       #17
                    % from line 69: kfac  = 1.0 - e*cosK;
                    %  f(67) = f(67) - f(69)*e;                              % cosK    #67
                    %  f(17) = f(17) - f(69)*cosK;                           % e       #17
                    % from line 68: sinK  = sin(EM);
                    f(66) = f(66) + f(68)*cosK;                           % EM      #66
                    % from line 67: cosK  = cos(EM);
                    f(66) = f(66) - f(67)*sinK;                           % EM      #66
                    % from line 66: Implicitly E = E(e, meanAnomaly)
                    % eccentricity line 17, meanAnomaly line 65
                    f(65) = f(65) + f(66)*dE_dM;                          % M       #65
                    f(17) = f(17) + f(66)*dE_de;                          % e       #17
                    % from line 65: M = meanMotion*t + MeanAnomalyEpoch;
                    f(59) = f(59) + f(65);                                % MeanAnomalyEpoch #59
                    f(46) = f(46) + f(65)*t;                              % meanMotion       #46
                    f( 1) = f( 1) + f(65)*meanMotion;                     % t       # 1
                end
                %/////////////////////Block 64-54//////////////////////////////////////
                % from line 64: Inclination = acos(Wz);
                f(27) = f(27) - f(64)/sqrt(1.0 - Wz*Wz);              % Wz      #27
                % from line 63: Tnext     = PTime*time_unit;
                % f(62) = f(62) + f(63)*time_unit;                      % PTime   # 62
                % from line 63: Mp = MeanAnomalyEpoch - meanMotion*t;
                f(59) = f(59) + f(63);                                % MeanAnomalyEpoch #59
                f(46) = f(46) - f(63)*t;                              % meanMotion  #46
                f( 1) = f( 1) - f(63)*meanMotion;                     % t     # 1
                %//////////////////////////////////////////////////////////////////////
                % from line 62 : if DeltaTime > 0 PTime = DeltaTime;
                %                if DeltaTime < 0 PTime = DeltaTime + Period
                %  PTime = DeltaTime  % DEFAULT                       % 62)
                %  if DeltaTime < 0
                %    Then OVERWRITE PTime HERE:
                %    PTime = DeltaTime + Period;
                %  end
                f(61) = f(61) + f(62);                                % DeltaTime  #61
                if DeltaTime < 0
                    f(45) = f(45) + f(62);                            % Period  #45
                end
                % from line 61: DeltaTime = t - TimeSincePeriapsis;
                f(60) = f(60) - f(61);                                % TimeSincePeriapsis #60
                f( 1) = f( 1) + f(61);                                % t        # 1
                % from line 60: TimeSincePeriapsis = MeanAnomalyEpoch/meanMotion;
                f(59) = f(59) + f(60)/meanMotion;                     % MeanAnomalyEpoch #59
                f(46) = f(46) - f(60)*MeanAnomalyEpoch/(meanMotion*meanMotion);  % meanMotion #46
                % from line 59: MeanAnomalyEpoch   = EccentricAnomalyEpoch - e*sinE;
                f(56) = f(56) + f(59);                                % EccentricAnomalyEpoch #56
                f(58) = f(58) - f(59)*e;                              % sinE    #58
                f(17) = f(17) - f(59)*sinE;                           % e       #17
                % from line 58: sinE    = sin(EccentricAnomalyEpoch);
                f(56) = f(56) + f(58)*cosE;                           % EccentricAnomalyEpoch #56
                % from line 57: cosE    = cos(EccentricAnomalyEpoch);
                f(56) = f(56) - f(57)*sinE;                           % EccentricAnomaly      #56
                % Eccentric Anomaly at Epoch
                % from line 56:  EccentricAnomalyEpoch = atan2(sinT, cosT);
                f(55) = f(55) + f(56)*cosT/(sinT*sinT + cosT*cosT);   % sinE    #55
                f(54) = f(54) - f(56)*sinT/(sinT*sinT + cosT*cosT);   % cosE    #54
                % old line 55: sinE    = sinnu*rootfac/(1.0 + e*cosnu);
                %  f(52) = f(52) + f(55)*rootfac/(1.0 + e*cosnu);       % sinnu   #52
                %  f(51) = f(51) - f(55)*sinnu*rootfac*e/((1.0 + e*cosnu)*(1.0 + e*cosnu));     % cosnu # 51
                %  f(42) = f(42) + f(55)*sinnu/(1.0 + e*cosnu);         % rootfac #42
                %  f(17) = f(17) - f(55)*sinnu*rootfac*cosnu/((1.0 + e*cosnu)*(1.0 + e*cosnu)); % e # 17
                % new line 55: sinT    = sinnu*rootfac;
                f(52) = f(52) + f(55)*rootfac;                        % sinnu   #52
                f(42) = f(42) + f(55)*sinnu;                          % rootfac #42
                % old line 54: cosE    = (e + cosnu)/(1.0 + e*cosnu);
                % f(51) = f(51) + f(54)*fac/((1.0 + e*cosnu)*(1.0 + e*cosnu));                  % cosnu #51
                % f(17) = f(17) + f(54)*sinnu*sinnu/((1.0 + e*cosnu)*(1.0 + e*cosnu));          % e     #17
                % new line 54: cosT    = e + cosnu;
                f(51) = f(51) + f(54);                                % cosnu   #51
                f(17) = f(17) + f(54);                                % e       #17
                %/////////////////////Block 53-51//////////////////////////////////////
                % True Anomaly at Epoch
                % from line 53: nuEpoch = atan2(sinnu, cosnu);
                f(52) = f(52) + f(53)*cosnu/(sinnu*sinnu + cosnu*cosnu);                      % sinnu #52
                f(51) = f(51) - f(53)*sinnu/(sinnu*sinnu + cosnu*cosnu);                      % cosnu #51
                % from line 52: sinnu1 = Wx*(Py*rz - Pz*ry) + Wy*(Pz*rx - Px*rz) + Wz*(Px*ry - Py*rx);
                %         sinnu    = Wx*(Py*rUnitz - Pz*rUnity) + ...
                %                    Wy*(Pz*rUnitx - Px*rUnitz) + ...
                %                    Wz*(Px*rUnity - Py*rUnitx);
                %      wrt Wvec
                f(27) = f(27) + f(52)*(Px*rUnity - Py*rUnitx);        % Wz      #27
                f(26) = f(26) + f(52)*(Pz*rUnitx - Px*rUnitz);        % Wy      #26
                f(25) = f(25) + f(52)*(Py*rUnitz - Pz*rUnity);        % Wx      #25
                %      wrt Pvec
                f(20) = f(20) + f(52)*(Wy*rUnitx - Wx*rUnity);        % Pz      #20
                f(19) = f(19) + f(52)*(Wx*rUnitz - Wz*rUnitx);        % Py      #19
                f(18) = f(18) + f(52)*(Wz*rUnity - Wy*rUnitz);        % Px      #18
                %      wrt Runit
                f(33) = f(33) + f(52)*(Wx*Py - Wy*Px);                % rUnitz  #33
                f(32) = f(32) + f(52)*(Wz*Px - Wx*Pz);                % rUnity  #32
                f(31) = f(31) + f(52)*(Wy*Pz - Wz*Py);                % rUnitx  #31
                
                % from line 51: cosnu    = rUnitx*Px + rUnity*Py + rUnitz*Pz;
                f(33) = f(33) + f(51)*Pz;                             % rUnitz  #33
                f(32) = f(32) + f(51)*Py;                             % rUnity  #32
                f(31) = f(31) + f(51)*Px;                             % rUnitx  #31
                
                f(20) = f(20) + f(51)*rUnitz;                         % Pz      #20
                f(19) = f(19) + f(51)*rUnity;                         % Py      #19
                f(18) = f(18) + f(51)*rUnitx;                         % Px      #18
                %/////////////////////Block 50-47//////////////////////////////////////
                % Kepler 4) omega
                % Argument of Periapsis
                % from line 50: omega=atan2(sinP, cosP);
                f(49) = f(49) + f(50)*cosP/(sinP*sinP + cosP*cosP);   % sinP    #49
                f(48) = f(48) - f(50)*sinP/(sinP*sinP + cosP*cosP);   % cosP    #48
                
                % from line 49: sinP =(Wx*Ny - Wy*Nx)*Pz + (Nx*Py - Ny*Px)*Wz;
                f(38) = f(38) + f(49)*(Wx*Pz - Wz*Px);                % Ny      #38
                f(37) = f(37) + f(49)*(Wz*Py - Wy*Pz);                % Nx      #37
                
                f(27) = f(27) + f(49)*(Nx*Py - Ny*Px);                % Wz      #27
                f(26) = f(26) + f(49)*(      - Nx*Pz);                % Wy      #26
                f(25) = f(25) + f(49)*(Ny*Pz        );                % Wx      #25
                
                f(20) = f(20) + f(49)*(Wx*Ny - Wy*Nx);                % Pz      #20
                f(19) = f(19) + f(49)*(Wz*Nx        );                % Py      #19
                f(18) = f(18) + f(49)*(      - Wz*Ny);                % Px      #18
                
                % from line 48: cosP =Nx*Px + Ny*Py;
                f(38) = f(38) + f(48)*Py;                             % Ny      #38
                f(37) = f(37) + f(48)*Px;                             % Nx      #37
                
                f(19) = f(19) + f(48)*Ny;                             % Py      #19
                f(18) = f(18) + f(48)*Nx;                             % Px      #18
                % Kepler 3) Omega
                % Longitude of the Ascending Node
                % from line 47:  Omega = atan2(Ny, Nx);
                f(38) = f(38) + f(47)*Nx/(Nx*Nx + Ny*Ny);             % Ny      #38
                f(37) = f(37) - f(47)*Ny/(Nx*Nx + Ny*Ny);             % Nx      #37
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %/////////////////////Block 46-39//////////////////////////////////////
                % from line 46: meanMotion = twopi/Period;              % Radians per Second
                f(45) = f(45) - f(46)*meanMotion/Period;              % Period  #45
                % Period
                % from line: 45 Period   = PeriodTU*time_unit;          % Period in Seconds
                %  f(44) = f(44) + f(45)*time_unit;                     % Period in Seconds
                % from line: 45 Period   = PeriodTU;                    % Period in TU's
                f(44) = f(44) + f(45);                                % Period in TU's #44
                % from line: 44 PeriodTU = twopi*a^(1.5)/sqrtmu;        % Period in TU's
                f(43) = f(43) + f(44)*PeriodTU*1.5/a;                 % a       #43
                % Semi-Major Axis:
                % from line 43: a = p/fac;
                f(41) = f(41) - f(43)*a/fac;                          % fac     #41
                f(24) = f(24) + f(43)/fac;                            % p       #24
                if (fac > 0.0)
                    % from line 42: rootfac = sqrt(fac);
                    f(41) = f(41) + f(42)*(0.5)/rootfac;               % fac     #41
                else
                    rootfac = 0;
                end
                % from line 41: fac   = onePe*oneMe;
                f(40) = f(40) + f(41)*onePe;                          % oneMe   #40
                f(39) = f(39) + f(41)*oneMe;                          % onePe   #39
                % from line 40: oneMe = 1.0 - e;
                f(17) = f(17) - f(40);                                % e       #17
                % from line 39: onePe = 1.0 + e;
                f(17) = f(17) + f(39);                                % e       #17
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %/////////////////////Block 38-1//////////////////////////////////////
                % from line 38: Ny     = Ay/normN;
                f(36) = f(36) - f(38)*Ny/normN;                       % normN   #36
                f(35) = f(35) + f(38)/normN;                          % Ay      #35
                % from line 37: Nx     = Ax/normN;
                f(36) = f(36) - f(37)*Nx/normN;                       % normN   #36
                f(34) = f(34) + f(37)/normN;                          % Ax      #34
                % from line 36: normN  = sqrt(Ax*Ax + Ay*Ay);
                f(35) = f(35) + f(36)*Ay/normN;                       % Ay      #35
                f(34) = f(34) + f(36)*Ax/normN;                       % Ax      #34
                % from line 35: Ay     =  Wx;
                f(25) = f(25) + f(35);                                % Wx      #25
                % from line 34: Ax     = -Wy;
                f(26) = f(26) - f(34);                                % Wy      #26
                % from line 33: rUnitz = rz/rmag;
                f(10) = f(10) - f(33)*rUnitz/rmag;                    % rmag    #10
                f( 4) = f( 4) + f(33)/rmag;                           % rz      # 4
                % from line 32: rUnity = ry/rmag;
                f(10) = f(10) - f(32)*rUnity/rmag;                    % rmag    #10
                f( 3) = f( 3) + f(32)/rmag;                           % ry      # 3
                % from line 31: rUnitx = rx/rmag;
                f(10) = f(10) - f(31)*rUnitx/rmag;                    % rmag    #10
                f( 2) = f( 2) + f(31)/rmag;                           % rx      # 2
                % from line 30: Qz   = Wx*Py - Wy*Px;
                f(26) = f(26) - f(30)*Px;                             % Wy      #26
                f(25) = f(25) + f(30)*Py;                             % Wx      #25
                f(19) = f(19) + f(30)*Wx;                             % Py      #19
                f(18) = f(18) - f(30)*Wy;                             % Px      #18
                % from line 29: Qy   = Wz*Px - Wx*Pz;
                f(27) = f(27) + f(29)*Px;                             % Wz      #27
                f(25) = f(25) - f(29)*Pz;                             % Wx      #25
                f(20) = f(20) - f(29)*Wx;                             % Pz      #20
                f(18) = f(18) + f(29)*Wz;                             % Px      #18
                % from line 28: Qx   = Wy*Pz - Wz*Py;
                f(27) = f(27) - f(28)*Py;                             % Wz      #27
                f(26) = f(26) + f(28)*Pz;                             % Wy      #26
                f(20) = f(20) + f(28)*Wy;                             % Pz      #20
                f(19) = f(19) - f(28)*Wz;                             % Py      #19
                % from line 27: Wz   = hz/sqrt(mu*p);  sqrt(mu*p) = hmag
                f(24) = f(24) - f(27)*Wz*0.5/p;                       % p       #24
                f(23) = f(23) + f(27)/hmag;                           % hz      #23
                % from line 26: Wy   = hy/sqrt(mu*p);  sqrt(mu*p) = hmag
                f(24) = f(24) - f(26)*Wy*0.5/p;                       % p       #24
                f(22) = f(22) + f(26)/hmag;                           % hy      #22
                % from line 25: Wx   = hx/sqrt(mu*p);  sqrt(mu*p) = hmag
                f(24) = f(24) - f(25)*Wx*0.5/p;                       % p       #24
                f(21) = f(21) + f(25)/hmag;                           % hx      #21
                % from line 24: p = (hx*hx + hy*hy + hz*hz)/mu;
                f(23) = f(23) + f(24)*2.0*hz/mu;                      % hz      #23
                f(22) = f(22) + f(24)*2.0*hy/mu;                      % hy      #22
                f(21) = f(21) + f(24)*2.0*hx/mu;                      % hx      #21
                % from line 23: hz   = rx*vy - ry*vx;
                f( 6) = f( 6) + f(23)*rx;                             % vy      # 6
                f( 5) = f( 5) - f(23)*ry;                             % vx      # 5
                f( 3) = f( 3) - f(23)*vx;                             % ry      # 3
                f( 2) = f( 2) + f(23)*vy;                             % rx      # 2
                % from line 22: hy   = rz*vx - rx*vz;
                f( 7) = f( 7) - f(22)*rx;                             % vz      # 7
                f( 5) = f( 5) + f(22)*rz;                             % vx      # 5
                f( 4) = f( 4) + f(22)*vx;                             % rz      # 4
                f( 2) = f( 2) - f(22)*vz;                             % rx      # 2
                % from line 21: hx   = ry*vz - rz*vy;
                f( 7) = f( 7) + f(21)*ry;                             % vz      # 7
                f( 6) = f( 6) - f(21)*rz;                             % vy      # 6
                f( 4) = f( 4) - f(21)*vy;                             % rz      # 4
                f( 3) = f( 3) + f(21)*vz;                             % ry      # 3
                % from line 20: Pz   = ez/e;
                f(17) = f(17) - f(20)*Pz/e;                           % e       #17
                f(16) = f(16) + f(20)/e;                              % ez      #16
                % from line 19: Py   = ey/e;
                f(17) = f(17) - f(19)*Py/e;                           % e       #17
                f(15) = f(15) + f(19)/e;                              % ey      #15
                % from line 18: Px   = ex/e
                f(17) = f(17) - f(18)*Px/e;                           % e       #17
                f(14) = f(14) + f(18)/e;                              % ex      #14
                % from line 17: e    = sqrt(ex*ex + ey*ey + ez*ez);
                f(16) = f(16) + f(17)*ez/e;                           % ez      #16
                f(15) = f(15) + f(17)*ey/e;                           % ey      #15
                f(14) = f(14) + f(17)*ex/e;                           % ex      #14
                % from line 16: ez   = (er*rz - ev*vz)/mu;
                f(13) = f(13) - f(16)*vz/mu;                             % ev      #13
                f(12) = f(12) + f(16)*rz/mu;                             % er      #12
                f( 7) = f( 7) - f(16)*ev/mu;                             % vz      # 7
                f( 4) = f( 4) + f(16)*er/mu;                             % rz      # 4
                % from line 15: ey   = (er*ry - ev*vy)/mu;
                f(13) = f(13) - f(15)*vy/mu;                             % ev      #13
                f(12) = f(12) + f(15)*ry/mu;                             % er      #12
                f( 6) = f( 6) - f(15)*ev/mu;                             % vy      # 6
                f( 3) = f( 3) + f(15)*er/mu;                             % ry      # 3
                % from line 14: ex   = (er*rx - ev*vx)/mu;
                f(13) = f(13) - f(14)*vx/mu;                             % ev      #13
                f(12) = f(12) + f(14)*rx/mu;                             % er      #12
                f( 5) = f( 5) - f(14)*ev/mu;                             % vx      # 5
                f( 2) = f( 2) + f(14)*er/mu;                             % rx      # 2
                % from line 13: ev   = rx*vx + ry*vy + rz*vz;
                f( 7) = f( 7) + f(13)*rz;                             % vz      # 7
                f( 6) = f( 6) + f(13)*ry;                             % vy      # 6
                f( 5) = f( 5) + f(13)*rx;                             % vx      # 5
                f( 4) = f( 4) + f(13)*vz;                             % rz      # 4
                f( 3) = f( 3) + f(13)*vy;                             % ry      # 3
                f( 2) = f( 2) + f(13)*vx;                             % rx      # 2
                % from line 12: er   = vsq - mu/rmag;
                f(11) = f(11) + f(12);                                % vsq     #11
                f(10) = f(10) + f(12)*mu/rmag^2;                      % rmag    #10
                % from line 11: vsq  = vx*vx + vy*vy + vz*vz;
                f( 7) = f( 7) + f(11)*2.0*vz;                         % vz      # 7
                f( 6) = f( 6) + f(11)*2.0*vy;                         % vy      # 6
                f( 5) = f( 5) + f(11)*2.0*vx;                         % vx      # 5
                % from line 10: rmag = sqrt(rx*rx + ry*ry + rz*rz);
                f( 4) = f( 4) + f(10)*rz/rmag;                        % rz      # 4
                f( 3) = f( 3) + f(10)*ry/rmag;                        % ry      # 3
                f( 2) = f( 2) + f(10)*rx/rmag;                        % rx      # 2
                
                %         vz   = Rdot(3);                                   % Vz      # 7
                %         vy   = Rdot(2);                                   % Vy      # 6
                %         vx   = Rdot(1);                                   % Vx      # 5
                %         rz   = Rpos(3);                                   % Rz      # 4
                %         ry   = Rpos(2);                                   % Ry      # 3
                %         rx   = Rpos(1);                                   % Rx      # 2
                %         t    = time;                                      % t       # 1
                F = f;
            end
            
            function Q = fordiff(q, LineCheck)
                %     eliminate some repetitions:
                %     rtpinv   = sqrt(1.0/p);
                %     Parameters 1-7:
                %     t    = time;                        %  1) t
                %     rx   = Rpos(1);                     %  2) Rx
                %     ry   = Rpos(2);                     %  3) Ry
                %     rz   = Rpos(3);                     %  4) Rz
                %     vx   = Rdot(1);                     %  5) Vx
                %     vy   = Rdot(2);                     %  6) Vy
                %     vz   = Rdot(3);                     %  7) Vz
                
                %///////////////////////////////////////Block 07-27/////////////////////////////////////////
                % Computation of forwards derivative starting at the top of the algorithm
                % from line 10: rmag = sqrt(rx*rx + ry*ry + rz*rz);
                q(10) = q(10) + ( q(2)*rx + q(3)*ry + q(4)*rz )/rmag;
                % from line 11: vsq  = vx*vx + vy*vy + vz*vz;
                q(11) = q(11) + q(5)*2.0*vx + q(6)*2.0*vy + q(7)*2.0*vz;
                % from line 12: er   = vsq - mu/rmag;
                q(12) = q(12) + q(11) + q(10)*mu/rmag^2;
                % from line 13: ev   = rx*vx + ry*vy + rz*vz;
                q(13) = q(13) + q(2)*vx + q(3)*vy + q(4)*vz;
                q(13) = q(13) + q(5)*rx + q(6)*ry + q(7)*rz;
                %//////////////////////////////////////////////////////////////////////////
                % 	// from line 14: ex   = (er*rx - ev*vx)/mu;
                q(14) = q(14) + (q( 2)*er - q( 5)*ev + q(12)*rx - q(13)*vx)/mu;
                % 	// from line 15: ey   = (er*ry - ev*vy)/mu;
                q(15) = q(15) + (q( 3)*er - q( 6)*ev + q(12)*ry - q(13)*vy)/mu;
                % 	// from line 16: ez   = (er*rz - ev*vz)/mu;
                q(16) = q(16) + (q( 4)*er - q( 7)*ev + q(12)*rz - q(13)*vz)/mu;
                % 	// from line 17: e    = sqrt(ex*ex + ey*ey + ez*ez);
                q(17) = q(17)  + (q(14)*ex + q(15)*ey + q(16)*ez)/e;
                % 	// from line 18: Px   = ex/e;
                q(18) = q(18) + (q(14) - q(17)*Px)/e;
                % 	// from line 19: Py   = ey/e;
                q(19) = q(19) + (q(15) - q(17)*Py)/e;
                % 	// from line 20: Pz   = ez/e;
                q(20) = q(20) + (q(16) - q(17)*Pz)/e;
                % 	// from line 21: hx   = ry*vz - rz*vy;
                q(21) = q(21) + q( 3)*vz - q( 4)*vy - q( 6)*rz + q( 7)*ry;
                % 	// from line 22: hy   = rz*vx - rx*vz;
                q(22) = q(22) - q( 2)*vz + q( 4)*vx + q( 5)*rz  - q( 7)*rx;
                % 	// from line 23: hz   = rx*vy - ry*vx;
                q(23) = q(23) + q( 2)*vy - q( 3)*vx - q( 5)*ry + q( 6)*rx;
                % 	// from line 24: p    = hx*hx + hy*hy + hz*hz;
                q(24) = q(24) + 2.0*(q(21)*hx + q(22)*hy + q(23)*hz)/mu;
                % 	// from line 25: Wx   = hx/sqrt(mu*p); sqrt(mu*p) = hmag
                q(25) = q(25) + q(21)/hmag - q(24)*0.5*Wx/p;
                % 	// from line 26: Wy   = hy/sqrt(mu*p); sqrt(mu*p) = hmag
                q(26) = q(26) + q(22)/hmag - q(24)*0.5*Wy/p;
                % 	// from line 27: Wz   = hz/sqrt(mu*p); sqrt(mu*p) = hmag
                q(27) = q(27) + q(23)/hmag - q(24)*0.5*Wz/p;
                % 	// from line 28: Qx   = Wy*Pz - Wz*Py;
                q(28) = q(28) - q(19)*Wz + q(20)*Wy + q(26)*Pz - q(27)*Py;
                % 	// from line 29: Qy   = Wz*Px - Wx*Pz;               % 29) Qvec(2)
                q(29) = q(29) + q(18)*Wz - q(20)*Wx - q(25)*Pz + q(27)*Px;
                % 	// from line 30: Qz   = Wx*Py - Wy*Px;
                q(30) = q(30) - q(18)*Wy + q(19)*Wx + q(25)*Py - q(26)*Px;
                % Other Vectors
                % 	// from line 31: rUnitx = rx/rmag;
                q(31) = q(31) + (q( 2) - q(10)*rUnitx)/rmag;
                % 	// from line 32: rUnity = ry/rmag;
                q(32) = q(32) + (q( 3) - q(10)*rUnity)/rmag;
                % 	// from line 33: rUnitz = rz/rmag;
                q(33) = q(33) + (q( 4) - q(10)*rUnitz)/rmag;
                % 	// from line 34: Ax     = -Wy;
                q(34) = q(34) - q(26);
                % 	// from line 35: Ay     =  Wx;
                q(35) = q(35) + q(25);
                % 	// from line 36: normN  = sqrt(Ax*Ax + Ay*Ay);
                q(36) = q(36) + (q(34)*Ax + q(35)*Ay)/normN;
                % 	// from line 37: Nx     = Ax/normN;
                q(37) = q(37) + (q(34) - q(36)*Nx)/normN;
                % 	// from line 38: Ny     = Ay/normN;
                q(38) = q(38) + (q(35) - q(36)*Ny)/normN;
                % Semi-Major Axis:
                % 	// from line 39: onePe = 1.0 + e;
                q(39) = q(39) + q(17);
                % 	// from line 40: oneMe = 1.0 - e;
                q(40) = q(40) - q(17);
                % 	// from line 41: fac   = onePe*oneMe;
                q(41) = q(41) + q(39)*oneMe + q(40)*onePe;
                % 	// from line 42: if(fac > 0) rootfac = sqrt(fac);
                if (fac > 0.0)
                    q(42) = q(42) + q(41)*0.5/rootfac;
                end
                % 	// from line 43: a = p/fac;
                q(43) = q(43) + q(24)/fac - q(41)*a/fac;
                %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %     % Period
                % 	// from line 44: PeriodTU = twopi*a^(1.5)/sqrt(mu);
                q(44) = q(44) + q(43)*1.5*PeriodTU/a;
                % 	// from line 45: Period   = PeriodTU*time_unit;
                %      q(45) = q(45) + q(44)*time_unit;
                % 	// from line 45: Period   = PeriodTU;
                q(45) = q(45) + q(44);
                % 	// from line 46: meanMotion = twopi/Period;
                q(46) = q(46) - q(45)*meanMotion/Period;
                %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %     % Kepler 3) Omega
                %     % Longitude of the Ascending Node
                % 	// from line 47: Omega = atan2(Ny, Nx);
                q(47) = q(47) - q(37)*Ny/(Nx*Nx + Ny*Ny) + q(38)*Nx/(Nx*Nx + Ny*Ny);
                %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %     % Kepler 4) omega
                %     % Argument of Periapsis
                % 	// from line 48: cosP =Nx*Px + Ny*Py;
                q(48) = q(48) + q(18)*Nx + q(19)*Ny + q(37)*Px + q(38)*Py;
                % 	// from line 49: sinP =(Wx*Ny - Wy*Nx)*Pz + (Nx*Py - Ny*Px)*Wz;
                q(49) = q(49) - q(18)*Ny*Wz + q(19)*Nx*Wz + q(20)*(Wx*Ny - Wy*Nx);
                q(49) = q(49) + q(25)*Ny*Pz - q(26)*Nx*Pz + q(27)*(Nx*Py - Ny*Px);
                q(49) = q(49) + q(37)*(Py*Wz - Wy*Pz) + q(38)*(Wx*Pz - Px*Wz);
                % 	// from line 50: omega=atan2(sinP, cosP);
                q(50) = q(50) + (q(49)*cosP - q(48)*sinP)/(cosP*cosP + sinP*sinP);
                % True Anomaly at Epoch
                % 	// from line 51: cosnu    = rUnitx*Px + rUnity*Py + rUnitz*Pz;
                q(51) = q(51) + q(18)*rUnitx + q(19)*rUnity + q(20)*rUnitz ...
                    + q(31)*Px     + q(32)*Py     + q(33)*Pz;
                % 	// from line 52: sinnu    = Wx*(Py*rUnitz - Pz*rUnity) + ...
                %                               Wy*(Pz*rUnitx - Px*rUnitz) + ...
                %                               Wz*(Px*rUnity - Py*rUnitx);
                q(52)  = q(52) + ...
                    q(18)*(Wz*rUnity - Wy*rUnitz) + q(19)*(Wx*rUnitz - Wz*rUnitx) + q(20)*(Wy*rUnitx - Wx*rUnity) + ...
                    q(25)*(Py*rUnitz - Pz*rUnity) + q(26)*(Pz*rUnitx - Px*rUnitz) + q(27)*(Px*rUnity - Py*rUnitx) + ...
                    q(31)*(Wy*Pz - Wz*Py)         + q(32)*(Wz*Px - Wx*Pz)         + q(33)*(Wx*Py - Wy*Px);
                %     trig    = cosnu*cosnu + sinnu*sinnu
                % 	// from line 53: nuEpoch = atan2(sinnu, cosnu);
                q(53) = q(53) + (q(52)*cosnu - q(51)*sinnu)/(cosnu*cosnu + sinnu*sinnu);
                %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %     % Eccentric Anomaly at Epoch
                % 	// from line 54: cosT    = e + cosnu;
                q(54) = q(54) + q(17) + q(51);
                % 	// from line 55: sinT    = sinnu*rootfac;
                q(55) = q(55) + q(42)*sinnu + q(52)*rootfac;
                % 	// from line 56: EccentricAnomalyEpoch = atan2(sinT, cosT);
                q(56) = q(56) + (q(55)*cosT - q(54)*sinT)/(cosT*cosT + sinT*sinT);
                %   // from line 57: cosE    = cos(EccentricAnomalyEpoch);
                q(57) = q(57) - q(56)*sinE;
                %   // from line 58: sinE    = sin(EccentricAnomalyEpoch);
                q(58) = q(58) + q(56)*cosE;
                % Kepler 5: Periapsis Time
                % Find the Time at Periapsis
                % Solve Kepler's Equation for the MeanAnomaly at Epoch
                % -- this side is easy!
                % // from line 59: MeanAnomalyEpoch = EccentricAnomalyEpoch - e*sinE;
                q(59) = q(59) - q(17)*sinE - q(58)*e + q(56);
                % // from line 60: TimeSincePeriapsis = MeanAnomalyEpoch/meanMotion;
                q(60) = q(60) + q(59)/meanMotion - q(46)*TimeSincePeriapsis/meanMotion;
                % PTime                  = TimeSincePeriapsis;
                % // from line 61: DeltaTime = t - TimeSincePeriapsis;
                q(61) = q(61) + q( 1) - q(60);
                % // from line 62: PTime:
                if DeltaTime > 0
                    q(62) = q(62) + q(61);
                else
                    q(62) = q(62) + q(61) + q(45);
                end
                % from line 63: Tnext     = PTime*time_unit;
                %  q(63) = q(63) + q(62)*time_unit;                     % PTime   # 62
                % from line 63: Mp = MeanAnomalyEpoch - meanMotion*t;
                q(63) = q(63) + q(59);                                % MeanAnomalyEpoch #59
                q(63) = q(63) - q(46)*t;                              % meanMotion  #46
                q(63) = q(63) - q( 1)*meanMotion;                     % t     # 1
                % Kepler 6: Inclination
                % // from line 64: Inclination = acos(Wz);
                q(64) = q(64) - q(27)/sqrt(1.0 - Wz*Wz);
                if LineCheck > 64
                    % // from line 65: M = meanMotion*t + MeanAnomalyEpoch;
                    q(65) = q(65) + q( 1)*meanMotion + q(46)*t + q(59);
                    % Solve Kepler's Equation for EM:  M = E - e*sin(E);   % 66) EM = E(e,M)
                    q(66) = q(66) + q(17)*dE_de + q(65)*dE_dM;
                    % from line 67: cosK  = cos(EM);
                    q(67) = q(67) - q(66)*sinK;
                    % from line 68: sinK  = sin(EM);
                    q(68) = q(68) + q(66)*cosK;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %   % for convenience -- you can eliminate kfac;
                    % kfac  = 1.0 - e*cosK;                                 % 69)
                    % // from line 70: tanX  = cosK - e;
                    q(70) = q(70) - q(17) + q(67);
                    % // from line 71: tanY  = rootfac*sinK;
                    q(71) = q(71) + q(42)*sinK + q(68)*rootfac;
                    % JRS: Comment out following two lines to take
                    %  Comment out to use (nu, rx, ry, rz, vx, vy, vz) as independent variables
                    % // from line 72: nu     = atan2(tanY, tanX);
                    q(72) = q(72) + (q(71)*tanX - q(70)*tanY)/(tanX*tanX + tanY*tanY);
                    % // from line 73: coss  = cos(nu);
                    q(73) = q(73) - q(72)*sins;
                    % // from line 74: sins  = sin(nu);
                    q(74) = q(74) + q(72)*coss;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Let us propagate this Orbit to time = t
                    % // from line 75: rorbit   = p/(1.0 + e*coss);
                    denom = 1.0/(1.0 + e*coss);
                    %      q(75) = q(75) + q(72)*p*e*sins*denom*denom;        % line 72 # nu
                    q(75) = q(75) - q(73)*p*e*denom*denom;              % line 73 # coss
                    q(75) = q(75) + q(24)*denom;
                    q(75) = q(75) - q(17)*p*coss*denom*denom;
                    % // from line 76: rorbitx  = rorbit*(coss*Px + sins*Qx);
                    %      q(76) = q(76)  + q(72)*rorbit*(-sins*Px + coss*Qx);% line 72 # nu
                    q(76) = q(76)  + q(75)*(coss*Px + sins*Qx);
                    q(76) = q(76)  + q(74)*rorbit*Qx;                   % line 74 # sins
                    q(76) = q(76)  + q(73)*rorbit*Px;                   % line 73 # coss
                    q(76) = q(76)  + q(28)*rorbit*sins;
                    q(76) = q(76)  + q(18)*rorbit*coss;
                    % // from line 77: rorbity  = rorbit*(coss*Py + sins*Qy);
                    %      q(77) = q(77) + q(72)*rorbit*(-sins*Py + coss*Qy); % line 72 # nu
                    q(77) = q(77) + q(75)*(coss*Py + sins*Qy);
                    q(77) = q(77) + q(74)*rorbit*Qy;                    % line 74 # sins
                    q(77) = q(77) + q(73)*rorbit*Py;                    % line 73 # coss
                    q(77) = q(77) + q(29)*rorbit*sins;
                    q(77) = q(77) + q(19)*rorbit*coss;
                    % // from line 78: rorbitz  = rorbit*(coss*Pz + sins*Qz);
                    %      q(78) = q(78) + q(72)*rorbit*(-sins*Pz + coss*Qz);  %line 72 # nu
                    q(78) = q(78) + q(75)*(coss*Pz + sins*Qz);
                    q(78) = q(78) + q(74)*rorbit*Qz;                     %line 74 # sins
                    q(78) = q(78) + q(73)*rorbit*Pz;                     %line 73 # coss
                    q(78) = q(78) + q(30)*rorbit*sins;
                    q(78) = q(78) + q(20)*rorbit*coss;
                    % // from line 79: vorbitx  = sqrt(1.0/p)*(-sins*Px + (e + coss)*Qx);
                    %      q(79) = q(79) - q(72)*rtpinv*(coss*Px + sins*Qx);%line 72 # nu
                    q(79) = q(79) - q(74)*rtpinv*Px;                  %line 74 # sins
                    q(79) = q(79) + q(73)*rtpinv*Qx;                  %line 73 # coss
                    q(79) = q(79) + q(28)*rtpinv*(e + coss);
                    q(79) = q(79) - q(24)*(0.5/(p*sqrt(p)))*(-sins*Px + (e + coss)*Qx);
                    q(79) = q(79) - q(18)*rtpinv*sins;
                    q(79) = q(79) + q(17)*rtpinv*Qx;
                    % // from line 80: vorbity  = sqrt(1.0/p)*(-sins*Py + (e + coss)*Qy);
                    %      q(80) = q(80) - q(72)*rtpinv*(coss*Py + sins*Qy);%line 72 # nu
                    q(80) = q(80) - q(74)*rtpinv*Py;                  %line 74 # sins
                    q(80) = q(80) + q(73)*rtpinv*Qy;                  %line 73 # coss
                    q(80) = q(80) + q(29)*rtpinv*(e + coss);
                    q(80) = q(80) - q(24)*(0.5/(p*sqrt(p)))*(-sins*Py + (e + coss)*Qy);
                    q(80) = q(80) - q(19)*rtpinv*sins;
                    q(80) = q(80) + q(17)*rtpinv*Qy;
                    % // from line 81: vorbitz  = sqrt(1.0/p)*(-sins*Pz + (e + coss)*Qz);
                    %      q(81) = q(81) - q(72)*rtpinv*(coss*Pz + sins*Qz);%line 72 # nu
                    q(81) = q(81) - q(74)*rtpinv*Pz;                  %line 74 # sins
                    q(81) = q(81) + q(73)*rtpinv*Qz;                  %line 73 # coss
                    q(81) = q(81) + q(30)*rtpinv*(e + coss);
                    q(81) = q(81) - q(24)*(0.5/(p*sqrt(p)))*(-sins*Pz + (e + coss)*Qz);
                    q(81) = q(81) - q(20)*rtpinv*sins;
                    q(81) = q(81) + q(17)*rtpinv*Qz;
                    
                    % from line 82: Aeq      = omega + Omega;         % line 82
                    q(82) = q(82) + q(47) + q(50);
                    % from line 83: cosAeq   = cos(Aeq);              % line 83
                    q(83) = q(83) - q(82)*sinAeq;
                    % from line 84: sinAeq   = sin(Aeq);              % line 84
                    q(84) = q(84) + q(82)*cosAeq;
                    % from line 85: tanHalf  = tan(Inclination/2.0);  % line 85
                    q(85) = q(85) + 0.5*q(64)*(1.0 + tanHalf^2);
                    % from line 86: Feq      = e*cosAeq;              % line 86
                    q(86) = q(86) + q(17)*cosAeq + q(83)*e;
                    % from line 87: Geq      = e*sinAeq;              % line 87
                    q(87) = q(87) + q(17)*sinAeq + q(84)*e;
                    % from line 88: Heq      = tanHalf*cos(Omega);    % line 88
                    q(88) = q(88) + q(85)*cosOmega - q(47)*tanHalf*sinOmega;
                    % from line 89: Keq      = tanHalf*sin(Omega);    % line 89
                    q(89) = q(89) + q(85)*sinOmega + q(47)*tanHalf*cosOmega;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % from line 90: Leq      = Aeq + nuEpoch;         % line 90
                    q(90) = q(90) + q(53) + q(82);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
                
                Q = q;
            end
            
            % function [H_I, Jacob, Hess] = WorkOrder(I)
            %
            %     H_I          = DataList(I);
            %     Nvar         = 7;
            %     LineCheck    = I;
            %     F            = zeros(90,1);
            %     F(LineCheck) = 1;
            %     F = backdiff(F, LineCheck);
            %     Jacob        = [F(1) F(2) F(3) F(4) F(5) F(6) F(7)];
            %     % e,a,...         t    e    a    I    w    W    Mp
            %     Hess   = [];
            %     % for j = 1: Nvar
            %     %     Q = zeros(133,1);
            %     %     S = ze/ros(133,1);
            %     %     %Q(ipoint(j)) = 1.0;
            %     %     Q(j) = 1.0;
            %     %     Q = this.fordiff(Q, LineCheck);
            %     %     S = this.secdiff(F, Q, S, LineCheck);
            %     %     S = this.backdiff(S, LineCheck);
            %     %     row = [];
            %     %     for k = 1 : Nvar
            %     %         row = [row, S(k)];
            %     %     end
            %     %     Hess  = [Hess; row];
            %     % end
            %
            %     this.Derivs = Jacob;
            %     this.H_I    = H_I
            % end
            
        end
        
        function [Jacob, JacobFinite, Hess, HessFinite] = DigitalJustice(this,I)
            
            time          = this.time;
            Rpos          = this.Rpos;
            Rdot          = this.Rdot;
            Kepler        = [Rpos(1), Rpos(2), Rpos(3), Rdot(1), Rdot(2), Rdot(3)];
            
            Kdelta        = KeplerOrbit(time, Rpos, Rdot);
            [H0, Jacob, Hess] = Kdelta.WorkOrder(I);
            %Hess
            
            % Do time Increment Separately
            
            epsilon     = 0.000000001;
            KeplerDelta = epsilon*eye(6);
            Kdelta = KeplerOrbit(time+epsilon, Kepler);
            [HT, JacobT, HessT] = Kdelta.WorkOrder(I);
            
            Kepler1     = Kepler + KeplerDelta(1,:);
            Kdelta1 = KeplerOrbit(time, Kepler1);
            [H1, Jacob1, Hess1] = Kdelta1.WorkOrder(I);
            
            Kepler2     = Kepler + KeplerDelta(2,:);
            Kdelta2 = KeplerOrbit(time, Kepler2);
            [H2, Jacob2, Hess2] = Kdelta2.WorkOrder(I);
            
            Kepler3     = Kepler + KeplerDelta(3,:);
            Kdelta3 = KeplerOrbit(time, Kepler3);
            [H3, Jacob3, Hess3] = Kdelta3.WorkOrder(I);
            
            Kepler4     = Kepler + KeplerDelta(4,:);
            Kdelta4 = KeplerOrbit(time, Kepler4);
            [H4, Jacob4, Hess4] = Kdelta4.WorkOrder(I);
            
            [~,~,~,~,~,~] = Kdelta5.OrbitDerivatives(time, Kepler);
            Kdelta5 = KeplerOrbit(time, Kepler5);
            [H5, Jacob5, Hess5] = Kdelta5.WorkOrder(I);
            
            Kepler6     = Kepler + KeplerDelta(6,:);
            Kdelta6 = KeplerOrbit(time, Kepler6);
            [H6, Jacob6, Hess6] = Kdelta6.WorkOrder(I);
            
            JacobFinite = [HT-H0,H1-H0, H2-H0, H3-H0, H4-H0, H5-H0, H6-H0]/epsilon;
            
            Hess = [];
            HessFinite = [];
        end
        
        % function [H_I, Jacob, Hess] = WorkOrder(this,I)
        %
        %     H_I          = this.DataList(I);
        %
        %     Nvar         = 7;
        %     LineCheck    = I;
        %     F            = zeros(90,1);
        %     F(LineCheck) = 1;
        %     F = backdiff(F, LineCheck);
        %     Jacob        = [F(1) F(2) F(3) F(4) F(5) F(6) F(7)];
        %     % e,a,...         t    e    a    I    w    W    Mp
        %     Hess   = [];
        %     % for j = 1: Nvar
        %     %     Q = zeros(133,1);
        %     %     S = ze/ros(133,1);
        %     %     %Q(ipoint(j)) = 1.0;
        %     %     Q(j) = 1.0;
        %     %     Q = this.fordiff(Q, LineCheck);
        %     %     S = this.secdiff(F, Q, S, LineCheck);
        %     %     S = this.backdiff(S, LineCheck);
        %     %     row = [];
        %     %     for k = 1 : Nvar
        %     %         row = [row, S(k)];
        %     %     end
        %     %     Hess  = [Hess; row];
        %     % end
        %
        % end
        
    end% end of member methods
end