classdef KeplerFromECI < handle
    
    properties (Constant = true)
        % Constants
        %extract conversion factors and ititialize to 1 in case they are not used
        %twopi  = 2.0*pi;         % 2*pi
        %XMNPDA = 1440.0;         % 1440 minutes per day
        %XKE    = 0.0743669161;   % Hoot's ke Gravitational constant for Earth in RE^(1.5)/minutes
    end
    
    properties (SetAccess = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        units;
        SignOrbit
        time;
        Rpos;
        Rdot;
        Derivs;
        RVeq;
        VVeq;
        StateVector7;
        ClassicalElements;
        EquinoctialElements;
        ExtrapolatedClass;
        KeplerToECI;
        
        JacECI_2_Kepler;
        JacKepler_2_ECI;
        
        JacECI_2_Equinoctial;
        JacEquinoctial_2_ECI;
        
        JacEquinoctial_2_Kepler
        JacKepler_2_Equinoctial;
        
        KeplerPropagator
        ECIPropagator
        EquinPropagator
        
        Hessian;
        DataList;
        Jacobian;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    properties (SetAccess = private)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        twopi;                  % twopi
        % TU;                     % TU  Canonical Time UnitTime Unit in Seconds
        % DU;                     % DU  Canonical Distance Unit
        % VU;                     % VU  Canonical Velocity Unit
        % AU;                     % AU  Canonical Acceleration Unit
        mu                     % mu  Canonical Gravitational Constant
        sqrtmu                 % sqrt(mu);
        hmag                   % sqrt(mu/p);
        
        %t;                                                    % line  1 t
        tEpoch;
        rx;                                                   % line  2 Rx
        ry;                                                   % line  3 Ry
        rz;                                                   % line  4 Rz
        vx;                                                   % line  5 Vx
        vy;                                                   % line  6 Vy
        vz;                                                   % line  7 Vz
        tExtrap;
        
        rmag;                                                 % line 10 rmag
        vsq;                                                  % line 11 vsq
        er;                                                   % line 12 er
        ev;                                                   % line 13 ev
        ex;                                                   % line 14 ex
        ey;                                                   % line 15 ey
        ez;                                                   % line 16 ez
        Evec;
        % Kepler 1:
        e;                                                    % line 17 e
        Px;                                                   % line 18 Px
        Py;                                                   % line 19 Py
        Pz;                                                   % line 20 Pz
        Pvec;
        % hVec = cross(Rpos, Rdot);
        hx;                                                   % line 21 hx
        hy;                                                   % line 22 hy
        hz;                                                   % line 23 hz
        hVec;
        % p    = dot(hVec,hVec)/mu;
        p;                                                    % line 24 p
        % Wvec = hVec/norm(hVec);
        Wx;                                                   % line 25 Wx
        Wy;                                                   % line 26 Wy
        Wz;                                                   % line 27 Wz
        Wvec;
        % Qvec = cross(Wvec, Pvec);
        Qx;                                                   % line 28 Qx
        Qy;                                                   % line 29 Qy
        Qz;                                                   % line 30 Qz
        Qvec;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rUnitx;                                               % line 31 rUnitx
        rUnity;                                               % line 32 rUnity
        rUnitz;                                               % line 33 rUnitz
        rUnit;
        Ax;                                                   % line 34 Ax
        Ay;                                                   % line 35 Ay
        normN;                                                % line 36 Az
        Nx;                                                   % line 37 Nx
        Ny;                                                   % line 38 Ny
        Nvec;
        onePe;                                                % line 39 onePe
        oneMe;                                                % line 40 oneMe
        fac;                                                  % line 41 fac                % 41)
        rootfac;                                              % line 42 rootfac
        % Kepler 2) a
        a;                                                    % line 43 a
        PeriodTU;                                             % line 44 PeriodTU
        Period;                                               % line 45
        meanMotion;                                           % line 46
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Kepler 5) Omega
        Omega;                                                % line 47 Omega
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Kepler 4) omega
        %cosPeriapsis; = dot(Nvec, Pvec);
        cosP;                                                 % line 48 cosP
        %sinPeriapsis; = dot(Wvec, cross(Nvec, Pvec));
        sinP;                                                 % line 49 sinP
        omega;                                                % line 50 omega
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % True Anomaly at Epoch
        %cosnu1;   = dot(rUnit,Pvec);
        %sinnu1;   = dot(Wvec,(cross(Pvec,rUnit)));
        cosnu;                                                % line 51 cos(nu)
        sinnu;                                                % line 52 sin(nu)
        nuEpoch;                                              % line 53 nu
        cosT;                                                 % line 54 cosT
        sinT;                                                 % line 55 sinT
        EccentricAnomalyEpoch;                                % line 56 E
        cosE;                                                 % line 57 cosE
        sinE;                                                 % line 58 sinE
        MeanAnomalyEpoch;                                     % line 59 MeanAnomalyEpoch
        TimeSincePeriapsis;                                   % line 60 TimeSincePeriapsis
        DeltaTime;                                            % line 61 DeltaTime
        PTime;                                                % line 62 PTime
        Mp;                                                   % line 63 Mp
        % Kepler 6: Inclination
        Inclination;                                          % line 64 Inclination
        cosI;
        sinI;
        cosom;
        sinom;
        cosO;
        sinO;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        M;                                                    % line 65 M
        EM;                                                   % line 66 EM = E(e,M)
        cosK;                                                 % line 67 cosK
        sinK;                                                 % line 68 sinK
        Eprime;
        eDenom;
        dE_dM;
        dE_de;
        d2E_dMdM;
        %// Mike Cain Corrections!
        d2E_dMde;
        d2E_dedM;
        d2E_dede;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tanX;                                                 % line 70 tanX
        tanY;                                                 % line 71 tanY
        % Kepler 4) omega
        nu;                                                   % line 72 nu
        coss;                                                 % line 73 coss
        sins;                                                 % line 74 sins
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Let us propagate this Orbit to time = t
        rorbit;                                               % line 75 rorbit
        rorbitx;                                              % line 76 rorbitx
        rorbity;                                              % line 77 rorbity
        rorbitz;                                              % line 78 rorbitz
        rtpinv;
        vorbitx;                                              % line 79 vorbitx
        vorbity;                                              % line 80 vorbity
        vorbitz;                                              % line 81 vorbitz
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Rextrap;
        Vextrap;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Aeq;                                                  % line 82 Aeq
        cosAeq;                                               % line 83 cosAeq
        sinAeq;                                               % line 84 sinAeq
        tanHalf;                                              % line 85 tanHalf
        Feq;                                                  % line 86 Feq
        Geq;                                                  % line 87 Geq
        Heq;                                                  % line 88 Heq
        Keq;                                                  % line 89 Keq
        Leq;                                                  % line 90 Leq
        CosL;                                                 % line 91 CosL
        SinL;                                                 % line 92 SinL
        
        alphaSq;                                              % line 93 alphaSq
        Seq;                                                     % line 94
        Weq;                                                     % line 95
        Req;                                                     % line 96
        RovS;                                                    % line 97
        srtpinv;                                                 % line 98
        HK;                                                      % line 99
        OnePalphaSqCosL;                                         % line 100
        OneMalphaSqCosL;                                         % line 101
        OnePalphaSqSinL;                                         % line 102
        OneMalphaSqSinL;                                         % line 103
        Xfac;                                                    % line 104
        Yfac;                                                    % line 105
        Zfac;                                                    % line 106
        VXfac;                                                   % line 107
        VYfac;                                                   % line 108
        VZfac;                                                   % line 109
        Xeq;                                                     % line 110
        Yeq;                                                     % line 111
        Zeq;                                                     % line 112
        VXeq;                                                    % line 113
        VYeq;                                                    % line 114
        VZeq;                                                    % line 115
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    methods
        
        % Constructor:
        function this = KeplerFromECI(tEpoch, Rpos, Rdot, units)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Canonical Coordinates assumed throughout -- length in DU, time in TU
            % Modifications can be made to the choice of units in wgs84Constants!
            % Orbital elements e, p, I, omega, Omega, PTime
            % or lines [17, 24, 64, 50, 47, 62]
            % and the parameters rx, ry, rz, vx, vy, vz
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Constants
            % For Choice of Units
            this.units = units; 
            twopi = units.twopi;
            this.twopi = twopi; 
            TU = units.TU;
            DU = units.DU;
            VU = units.VU;
            AU = units.AU;
            mu = units.mu;
            sqrtmu      = units.sqrtmu;
            this.mu     = mu;
            this.sqrtmu = sqrtmu;
            %wgs84       = wgs84Constants();
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Ivec  = [1.0 0.0 0.0];
            %Jvec  = [0.0 1.0 0.0];
            %Kvec  = [0.0 0.0 1.0];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            SignOrbit = 1;
            %t    = tEpoch;                      %  1) t
            rx   = Rpos(1);                     %  2) Rx
            ry   = Rpos(2);                     %  3) Ry
            rz   = Rpos(3);                     %  4) Rz
            vx   = Rdot(1);                     %  5) Vx
            vy   = Rdot(2);                     %  6) Vy
            vz   = Rdot(3);                     %  7) Vz
            tExtrap = tEpoch;
            
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
            this.hmag = hmag;
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
            Nvec   = [Nx,  Ny, 0.0];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Semi-Major Axis:
            onePe = 1.0 + e;                    % 39)
            oneMe = 1.0 - e;                    % 40)
            fac   = onePe*oneMe;                % 41)
            if (e < 1.0)
                SignOrbit =  1;
            else
                if e == 1
                    SignOrbit = 0;    % Reserve for Parabolic Orbits
                else
                    SignOrbit = -1;
                end
            end
            rootfac = sqrt(SignOrbit*fac);      % 42)
            a = p/fac;                          % 43)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Period
            PeriodTU   = twopi*(SignOrbit*a)^(1.5)/sqrtmu;  % 44) Period in TU's
            % Period     = PeriodTU*time_unit;    % 45) Period in Seconds
            Period     = PeriodTU;                % 45) Work in TU's for now
            meanMotion = twopi/Period;            % 46) Radians per TU.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Kepler 3) Omega
            % Longitude of the Ascending Node
            Omega = atan2(Ny, Nx);              % 47) Omega = atan(Ny/Nx);
            while Omega<0.0
                Omega = Omega + twopi;
            end
            while Omega>twopi
                Omega = Omega - twopi;
            end
            cosO = cos(Omega);
            sinO = sin(Omega);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Kepler 4) omega
            % Argument of Periapsis
            %cosPeriapsis = dot(Nvec, Pvec);
            % Note that Nz = 0!!!
            cosP =Nx*Px + Ny*Py;                           % 48) cos(omega)
            %sinPeriapsis = dot(Wvec, cross(Nvec, Pvec));
            sinP =(Wx*Ny - Wy*Nx)*Pz + (Nx*Py - Ny*Px)*Wz; % 49) sin(omega)
            
            omega=atan2(sinP, cosP);                       % 50) omega = atan(sinP/cosP)
            while omega<0.0
                omega = omega + twopi;
            end
            while omega>twopi
                omega = omega - twopi;
            end
            %omega
            cosom       = cos(omega);                  % line  19
            sinom       = sin(omega);                  % line  20
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % True Anomaly at Epoch
            %cosnu1   = dot(rUnit,Pvec);
            cosnu    = rUnitx*Px + rUnity*Py + rUnitz*Pz;  % 51) cos(nu)
            %sinnu1   = dot(Wvec,(cross(Pvec,rUnit)));
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
            %///////////////////////////////////////////////////////////////
            % TimeSincePeriapsis              = MeanAnomalyEpoch/meanMotion;   % 60)
            % PTime                         = TimeSincePeriapsis;
            % DeltaTime                       = tEpoch - TimeSincePeriapsis;     % 61)
            %FractionOfPeriodSincePeriapsis = TimeSincePeriapsis/Period;
            % DeltaTime is in the intervale [-T/2, T/2]
            % If DeltaTime > 0 then DeltaTime is in [0, T/2], then do nothing
            % If DeltaTime in [-T/2, 0] then add a Period DeltaTime -> [T/2,T]
            % If you see PTime > Period/2, then above if-statement was executed.
            % PTime = DeltaTime;  % DEFAULT                     % 62) in TUs
            % if DeltaTime < 0
            %     % The OVERWRITE PTime HERE:
            %     PTime = DeltaTime + Period;
            % end
            %Tnext     = PTime*time_unit;                    % 63 in Seconds
            %///////////////////////////////////////////////////////////////
            Mp = MeanAnomalyEpoch - meanMotion*tEpoch;            % 63 in radians
            while Mp < 0
                Mp = Mp + twopi;
            end
            while Mp > twopi
                Mp = Mp - twopi;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Kepler 6: Inclination
            Inclination = acos(Wz);                          % 64)
            cosI        = cos(Inclination);
            sinI        = sin(Inclination);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %M = meanMotion*t + MeanAnomalyEpoch;             % Orig     65)
            % JRS TEST: Uses the Current MeanAnomaly at this event EPOCH.
            %M = MeanAnomalyEpoch;                          % Test line  65)
            M = Mp + meanMotion*tExtrap;                      % New line 65
            while M < 0
                M = M + twopi;
            end
            while M > twopi
                M = M - twopi;
            end
            if e < 1
                % Elliptical Case:
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
                % Solve Kepler's Equation for E:  M = EM - e*sin(EM) % line  66 EM = E(e,M)
                % Need to differentiate Implicitly here!
                eDenom   = Eprime*Eprime*Eprime;
                this.Eprime = Eprime;
                this.eDenom = eDenom;
                
                % KeplerInv = EM;
                % //cout  << " Newton-Raphson Kepler's Equation Eccentric Anomaly =" << E << endl;
                % eDenom   = Eprime*Eprime*Eprime;
                
                %rmag = p/(1 + e*cosnu);
                %rRecon = rmag*(cosnu*Pvec + sinnu*Qvec);
                %QReco = rmag*cosnu*Rdot + (1.0/sqrt(p))*sinnu*Rpos;
                %QReco  = (1.0/norm(QReco))*QReco;
                
                cosK  = cos(EM);                                     % line  67
                sinK  = sin(EM);                                     % line  68
                %E     = atan2(sinK,cosK)
                %%%%%%%%%%%%%%%%%%% Implicit Function Theorem Derivatives
                dE_dM    = 1.0/Eprime;
                dE_de    = sinK/Eprime;
                d2E_dMdM = -e*sinK/eDenom;
                %// Mike Cain Corrections!
                d2E_dMde = (cosK - e)/eDenom;
                d2E_dedM = (cosK - e)/eDenom;
                d2E_dede = ((2.0 - e*cosK)*cosK*sinK - e*sinK)/eDenom;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                this.eDenom   = eDenom;
                this.Eprime   = Eprime;
                this.dE_dM    = dE_dM;
                this.dE_de    = dE_de;
                this.d2E_dMdM = d2E_dMdM;
                this.d2E_dMde = d2E_dMde;
                this.d2E_dedM = d2E_dedM;
                this.d2E_dede = d2E_dede;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % kfac  = 1.0 - e*cosK;                              % line  69
                % coss  = (cosK - e)/kfac;                           % line  70
                % sins  = rootfac*sinK/kfac;                         % line  71
                % 	s     = atan2(sins, coss);                       % line  72
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   % for convenience -- you can eliminate kfac;
                tanX  = cosK - e;                                    % line  70
                tanY  = rootfac*sinK;                                % line  71
            else
                if e == 1
                    % Reserved for Parabolic Case
                else
                    if e < 1.6
                        if ((-pi < M) && (M<0)) || M>pi
                            EM = M -e;
                        else
                            EM = M + e;
                        end
                    else
                        if (e < 3.6) &&  abs(M) >pi
                            EM = M - sign(M)*e
                        else
                            EM = M/(e - 1.);
                        end
                    end
                    
                    %  meanAnomaly M to eccentric Anomaly EM
                    %         //cout << " Input meanAnomaly = " << meanAnomaly << endl;
                    for i=1:1:10
                        EM  = e*sinh(EM) - M;
                        %//cout << " Mean Anomaly Solution " << E << endl;
                    end
                    % //      10 rounds of Newton's root finding method based on the above "seed".
                    for i=1:1:10
                        Eprime      = e*cosh(EM) - 1.0;
                        EM          = EM + (M + EM -e*sinh(EM))/Eprime;
                    end
                    % Solve Kepler's Equation for E:  M = EM - e*sin(EM) % line  66 EM = E(e,M)
                    % Need to differentiate Implicitly here!
                    eDenom   = Eprime*Eprime*Eprime;
                    this.Eprime = Eprime;
                    this.eDenom = eDenom;
                    
                    % KeplerInv = EM;
                    % //cout  << " Newton-Raphson Kepler's Equation Eccentric Anomaly =" << E << endl;
                    % eDenom   = Eprime*Eprime*Eprime;
                    
                    %rmag = p/(1 + e*cosnu);
                    %rRecon = rmag*(cosnu*Pvec + sinnu*Qvec);
                    %QReco = rmag*cosnu*Rdot + (1.0/sqrt(p))*sinnu*Rpos;
                    %QReco  = (1.0/norm(QReco))*QReco;
                    
                    %N.B. We're using the same variable names for Hyperbolic Functions!!!
                    cosK  = cosh(EM);    %Hyperbolic cos         % line  67
                    sinK  = sinh(EM);    %Hyperbolic sin         % line  68
                    %E     = atan2(sinK,cosK)
                    %%%%%%%%%%%%%%%%%%% Implicit Function Theorem Derivatives
                    dE_dM    =  1.0/Eprime;
                    dE_de    = -sinK/Eprime;
                    d2E_dMdM = -e*sinK/eDenom;
                    %// Mike Cain Corrections!
                    d2E_dMde = (cosK - e)/eDenom;
                    d2E_dedM = (cosK - e)/eDenom;
                    d2E_dede = ((e*cosK -2.0)*cosK*sinK + e*sinK)/eDenom;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    this.eDenom   = eDenom;
                    this.Eprime   = Eprime;
                    this.dE_dM    = dE_dM;
                    this.dE_de    = dE_de;
                    this.d2E_dMdM = d2E_dMdM;
                    this.d2E_dMde = d2E_dMde;
                    this.d2E_dedM = d2E_dedM;
                    this.d2E_dede = d2E_dede;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % kfac  = 1.0 - e*cosK;                              % line  69
                    % coss  = (cosK - e)/kfac;                           % line  70
                    % sins  = rootfac*sinK/kfac;                         % line  71
                    % 	s     = atan2(sins, coss);                       % line  72
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %   % for convenience -- you can eliminate kfac;
                    tanX  = cosK - e;                                    % line  70
                    tanY  = -rootfac*sinK;                                % line  71
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nu     = atan2(tanY, tanX);                          % line  72
            coss  = cos(nu);                                     % line  73
            sins  = sin(nu);                                     % line  74
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Position and Unit Vector along the Position
            %Kvec  = a*kfac*(cosnu*Pvec + sinnu*Qvec);
            %Runit = Kvec/norm(Kvec);
            
            %   Unit Vector along the Velocity
            %Vunit =(-sinnu*Pvec + (e + cosnu)*Qvec);
            %Vunit = Vunit/norm(Vunit);
            
            %   Unit Vector out of the R-V plane
            %Wlocal = cross(Runit, Vunit);
            %Wlocal = Wlocal/norm(Wlocal);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Let us propagate this Orbit to time = tExtrap
            %rorbit   = p/(1.0 + e*coss);                        % line  75
            rorbit   = a*(1.0 - e*cosK);
            
            rorbitx  = rorbit*(coss*Px + sins*Qx);               % line  76
            rorbity  = rorbit*(coss*Py + sins*Qy);               % line  77
            rorbitz  = rorbit*(coss*Pz + sins*Qz);               % line  78
            rtpinv   = sqrt(mu/p);                               % line  69
            vorbitx  = rtpinv*(-sins*Px + (e + coss)*Qx);        % line  79
            vorbity  = rtpinv*(-sins*Py + (e + coss)*Qy);        % line  80
            vorbitz  = rtpinv*(-sins*Pz + (e + coss)*Qz);        % line  81
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Rextrap  = [rorbitx; rorbity; rorbitz];
            Vextrap  = [vorbitx; vorbity; vorbitz];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Treat Equinoctial Variables:
            Aeq      = omega + Omega;                            % line  82
            while Aeq < 0
                Aeq = Aeq + twopi;
            end
            while Aeq > twopi
                Aeq = Aeq - twopi;
            end
            
            cosAeq   = cos(Aeq);                                 % line  83
            sinAeq   = sin(Aeq);                                 % line  84
            tanHalf  = tan(Inclination/2.0);                     % line  85
            Feq      = e*cosAeq;                                 % line  86
            Geq      = e*sinAeq;                                 % line  87
            Heq      = tanHalf*cosO;                             % line  88
            Keq      = tanHalf*sinO;                             % line  89
            Leq      = Aeq + nu;                                 % line  90
            while Leq < 0
                Leq = Leq + twopi;
            end
            while Leq > twopi
                Leq = Leq - twopi;
            end
                        
            CosL     = cos(Leq);                                 % line  91
            SinL     = sin(Leq);                                 % line  92
            alphaSq = Heq^2 - Keq^2;                             % line  93
            Seq     = 1.0 + Heq^2 + Keq^2;                       % line  94
            Weq     = 1.0 + Feq*CosL + Geq*SinL;                 % line  95
            Req     = p/Weq;                                     % line  96
            RovS    = Req/Seq;                                   % line  97
            srtpinv = rtpinv/Seq;                                % line  98
            HK      = Heq*Keq;                                   % line  99
            OnePalphaSqCosL = (1.+alphaSq)*CosL;                 % line 100
            OneMalphaSqCosL = (1.-alphaSq)*CosL;                 % line 101
            OnePalphaSqSinL = (1.+alphaSq)*SinL;                 % line 102
            OneMalphaSqSinL = (1.-alphaSq)*SinL;                 % line 103
            Xfac    = OnePalphaSqCosL + 2.0*HK*SinL;             % line 104
            Yfac    = OneMalphaSqSinL + 2.0*HK*CosL;             % line 105
            Zfac    = Heq*SinL - Keq*CosL;                       % line 106
            VXfac   =  OnePalphaSqSinL - 2.0*HK*(Feq + CosL) + Geq*(1. + alphaSq);% line 107
            VYfac   = -OneMalphaSqCosL + 2.0*HK*(Geq + SinL) + Feq*(alphaSq - 1.);% line 108
            VZfac   =  Heq*(Feq+CosL) + Keq*(Geq+SinL);          % line 109
            Xeq     =     RovS*Xfac;                             % line 110
            Yeq     =     RovS*Yfac;                             % line 111
            Zeq     = 2.0*RovS*Zfac;                             % line 112
            VXeq    =    -srtpinv*VXfac;                         % line 113
            VYeq    =    -srtpinv*VYfac;                         % line 114
            VZeq    = 2.0*srtpinv*VZfac;                         % line 115
            
            RVeq    = [Xeq; Yeq; Zeq];
            %Rpos'
            VVeq    = [VXeq; VYeq; VZeq];
            %Rdot'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Kepler      = [e, a, Inclination, omega, Omega, Mp];
            KeplerToECI = LagrangePlanetary(Kepler, this.units);
            this.KeplerToECI = KeplerToECI;
            epsilon     = 0.0001*ones(6,1);
            % Check DummySensor and DummyLos
            DummySensor = [Rpos + 10000*epsilon(1:3); zeros(6,1)];
            DummyLos    = [Rpos  - ones(3,1)];
            %ummyLos    = ones(1,3);
            DummyLos    = DummyLos/norm(DummyLos);
            InvCovUV    = eye(2);
            [Rinit,Vinit,ParamList,JacKepler_2_ECI,Hessian,GravityCan] = KeplerToECI.OrbitDerivatives(tEpoch, DummySensor, DummyLos, InvCovUV);
            [Rinit,Vinit,ParamList] = OrbitAtTime(KeplerToECI, tEpoch, DummySensor, DummyLos, InvCovUV);
            this.JacKepler_2_ECI = KeplerToECI.Jacobian;
            this.Hessian         = KeplerToECI.Hessian;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             % Test Section
            %             % Okay given the TimeSincePeriapsis -- can we go back to the true Anomaly?
            %             % // 1)   Find Eccentric Anomaly from Mean Anomaly
            %             % //      Solve Kepler's Equation:  M = E - e*sin(E);
            %             % //      write as E = M + e*sin(E) and iterate 10 times.
            %             %
            %             % //  Brute Force iteration (good way to seed Newton's method which follows)
            %             % Here we will insert the time and proceed to calculate
            %             % time > MeanAnomaly > Eccentric Anomaly > True Anomaly > Position & Velocity
            %             % and the calculate the Terms in Kepler, J2-J6 Perturbations using the
            %             % linearization obtained from the Taylor expansion of the Gravitational
            %             % Potential from Escobal pages 48-52
            %             % The time variable is set up so that time=0 corresponds to the Epoch Time
            %             % We can always use Time to propagate before (time < 0) or after (time > 0)
            %             % the Epoch Time -- but this requires solving Kepler's Equation.
            %             % JRS TEST Use for Extrapolation plus time t:
            %             %M = meanMotion*t + MeanAnomalyEpoch;                             % 65)
            %             % JRS TEST: Uses the Current MeanAnomaly at this event EPOCH.
            %             M = MeanAnomalyEpoch;                                             % 65)
            %             %  meanAnomaly M to eccentric Anomaly EM
            %             EM = M + e;
            %             if (M > pi)
            %                 EM = M - e;
            %             elseif (-pi < M && M < 0)
            %                 EM = M - e;
            %             end
            %             %         //cout << " Input meanAnomaly = " << meanAnomaly << endl;
            %             for i=1:1:10
            %                 EM  = M + e*sin(EM);
            %                 %//cout << " Mean Anomaly Solution " << E << endl;
            %             end
            %             % //      10 rounds of Newton's root finding method based on the above "seed".
            %             for i=1:1:10
            %                 Eprime      = 1.0 - e*cos(EM);
            %                 EM          = EM + (M - EM + e*sin(EM))/Eprime;
            %             end
            %             % Solve Kepler's Equation for E:  M = EM - e*sin(EM);   % 66) EM = E(e,M)
            %             % Need to differentiate Implicitly here!
            %             eDenom   = Eprime*Eprime*Eprime;
            %             this.Eprime = Eprime;
            %             this.eDenom = eDenom;
            %
            %             % KeplerInv = EM;
            %             % //cout  << " Newton-Raphson Kepler's Equation Eccentric Anomaly =" << E << endl;
            %             % eDenom   = Eprime*Eprime*Eprime;
            %
            %             %rmag = p/(1 + e*cosnu);
            %             %rRecon = rmag*(cosnu*Pvec + sinnu*Qvec);
            %             %QReco = rmag*cosnu*Rdot + (1.0/sqrt(p))*sinnu*Rpos;
            %             %QReco  = (1.0/norm(QReco))*QReco;
            %
            %             cosK  = cos(EM);                                        % 67)
            %             sinK  = sin(EM);                                        % 68)
            %             %E     = atan2(sinK,cosK)
            %             %%%%%%%%%%%%%%%%%%% Implicit Function Theorem Derivatives
            %             dE_dM    = 1.0/Eprime;
            %             dE_de    = sinK/Eprime;
            %             d2E_dMdM = -e*sinK/eDenom;
            %             %// Mike Cain Corrections!
            %             d2E_dMde = (cosK - e)/eDenom;
            %             d2E_dedM = (cosK - e)/eDenom;
            %             d2E_dede = ((2.0 - e*cosK)*cosK*sinK - e*sinK)/eDenom;
            %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             this.eDenom   = eDenom;
            %             this.Eprime   = Eprime;
            %             this.dE_dM    = dE_dM;
            %             this.dE_de    = dE_de;
            %             this.d2E_dMdM = d2E_dMdM;
            %             this.d2E_dMde = d2E_dMde;
            %             this.d2E_dedM = d2E_dedM;
            %             this.d2E_dede = d2E_dede;
            %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             % kfac  = 1.0 - e*cosK;                                 % 69)
            %             % coss  = (cosK - e)/kfac;                              % 70)
            %             % sins  = rootfac*sinK/kfac;                            % 71)
            %             % 	s     = atan2(sins, coss);                          % 72
            %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             %   % for convenience -- you can eliminate kfac;
            %             tanX  = cosK - e;                                     % line 70
            %             tanY  = rootfac*sinK;                                 % line 71
            %             nu     = atan2(tanY, tanX);                           % line 72
            %             coss  = cos(nu);                                      % line 73
            %             sins  = sin(nu);                                      % line 74
            %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             %   Position and Unit Vector along the Position
            %             %Kvec  = a*kfac*(cosnu*Pvec + sinnu*Qvec);
            %             %Runit = Kvec/norm(Kvec);
            %
            %             %   Unit Vector along the Velocity
            %             %Vunit =(-sinnu*Pvec + (e + cosnu)*Qvec);
            %             %Vunit = Vunit/norm(Vunit);
            %
            %             %   Unit Vector out of the R-V plane
            %             %Wlocal = cross(Runit, Vunit);
            %             %Wlocal = Wlocal/norm(Wlocal);
            %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             % Let us propagate this Orbit to time = t
            %             rorbit   = p/(1.0 + e*coss);                          % line 75
            %             rorbitx  = rorbit*(coss*Px + sins*Qx);                % line 76
            %             rorbity  = rorbit*(coss*Py + sins*Qy);                % line 77
            %             rorbitz  = rorbit*(coss*Pz + sins*Qz);                % line 78
            %             rtpinv   = sqrt(mu/p);                                % line 69
            %             vorbitx  = rtpinv*(-sins*Px + (e + coss)*Qx);         % line 79
            %             vorbity  = rtpinv*(-sins*Py + (e + coss)*Qy);         % line 80
            %             vorbitz  = rtpinv*(-sins*Pz + (e + coss)*Qz);         % line 81
            %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             Rextrap  = [rorbitx; rorbity; rorbitz];
            %             Vextrap  = [vorbitx; vorbity; vorbitz];
            %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             % Treat Equinoctial Variables:
            %             Aeq      = omega + Omega;                             % line 82
            %             cosAeq   = cos(Aeq);                                  % line 83
            %             sinAeq   = sin(Aeq);                                  % line 84
            %             tanHalf  = tan(Inclination/2.0);                      % line 85
            %             Feq      = e*cosAeq;                                  % line 86
            %             Geq      = e*sinAeq;                                  % line 87
            %             Heq      = tanHalf*cosO;                          % line 88
            %             Keq      = tanHalf*sinO;                          % line 89
            %             Leq      = Aeq + nu;                                  % line 90
            %             CosL     = cos(Leq);                                  % line 91
            %             SinL     = sin(Leq);                                  % line 92
            %             alphaSq = Heq^2 - Keq^2;                              % line 93
            %             Seq     = 1.0 + Heq^2 + Keq^2;                        % line 94
            %             Weq     = 1.0 + Feq*CosL + Geq*SinL;                  % line 95
            %             Req     = p/Weq;                                      % line 96
            %             RovS    = Req/Seq;                                    % line 97
            %             srtpinv = rtpinv/Seq;                                 % line 98
            %             HK      = Heq*Keq;                                    % line 99
            %             OnePalphaSqCosL = (1.+alphaSq)*CosL;                 % line 100
            %             OneMalphaSqCosL = (1.-alphaSq)*CosL;                 % line 101
            %             OnePalphaSqSinL = (1.+alphaSq)*SinL;                 % line 102
            %             OneMalphaSqSinL = (1.-alphaSq)*SinL;                 % line 103
            %             Xfac    = OnePalphaSqCosL + 2.0*HK*SinL;             % line 104
            %             Yfac    = OneMalphaSqSinL + 2.0*HK*CosL;             % line 105
            %             Zfac    = Heq*SinL - Keq*CosL;                       % line 106
            %             VXfac   =  OnePalphaSqSinL - 2.0*HK*(Feq + CosL) + Geq*(1. + alphaSq);% line 107
            %             VYfac   = -OneMalphaSqCosL + 2.0*HK*(Geq + SinL) + Feq*(alphaSq - 1.);% line 108
            %             VZfac   =  Heq*(Feq+CosL) + Keq*(Geq+SinL);  % line 109
            %             Xeq     =     RovS*Xfac;                             % line 110
            %             Yeq     =     RovS*Yfac;                             % line 111
            %             Zeq     = 2.0*RovS*Zfac;                             % line 112
            %             VXeq    =    -srtpinv*VXfac;                         % line 113
            %             VYeq    =    -srtpinv*VYfac;                         % line 114
            %             VZeq    = 2.0*srtpinv*VZfac;                         % line 115
            %
            %             RVeq    = [Xeq; Yeq; Zeq];
            %             %Rpos'
            %             VVeq     = [VXeq; VYeq; VZeq];
            %             %Rdot'
            %
            %             StateVector            = [Rextrap, Vextrap; tExtrap];
            %             ClassicalElements       = [e; a; Inclination; omega; Omega; Mp; tExtrap];
            %             EquinoctialElements      = [p; Feq; Geq; Heq; Keq; Leq; tExtrap];
            %             this.ClassicalElements  = ClassicalElements;
            %             this.EquinoctialElements = EquinoctialElements;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            this.SignOrbit = SignOrbit;
            this.time   = tEpoch;
            this.Rpos   = Rpos;
            this.Rdot   = Rdot;
            
            %this.t       = tEpoch;                                         %  1) t
            this.tEpoch  = tEpoch;                                         %  1) t
            this.rx      = rx;                                        %  2) Rx
            this.ry      = ry;                                        %  3) Ry
            this.rz      = rz;                                        %  4) Rz
            this.vx      = vx;                                        %  5) Vx
            this.vy      = vy;                                        %  6) Vy
            this.vz      = vz;                                        %  7) Vz
            this.tExtrap = tEpoch;
            
            this.rmag = rmag;                                         % 10) rmag = norm(Rpos);
            this.vsq  = vsq;                                          % 11) vsq  = dot(Rdot,Rdot);
            this.er   = vsq - mu/rmag;                                % 12) er   = v^2 - mu/rmag;
            this.ev   = ev;                    % 13) ev   = dot(Rpos,Rdot);
            this.ex   = ex;                    % 14) Evec(1)
            this.ey   = ey;                    % 15) Evec(2)
            this.ez   = ez;                    % 16) Evec(3)
            this.Evec = [ex ey ez];
            this.e    = e;                     % 17) e    = norm(Evec);
            this.Px   = Px;                    % 18) Pvec(1)
            this.Py   = Py;                    % 19) Pvec(2)
            this.Pz   = Pz;                    % 20) Pvec(3)
            this.Pvec = [Px  Py  Pz];
            this.hx   = hx;                    % 21) hVec(1)
            this.hy   = hy;                    % 22) hVec(2)
            this.hz   = hz;                    % 23) hVec(3)
            this.hVec = [hx hy hz];
            this.p    = p;                     % 24) p = hVec^2/mu
            this.hmag = hmag;
            this.Wx   = hx/hmag;               % 25) Wvec(1)
            this.Wy   = hy/hmag;               % 26) Wvec(2)
            this.Wz   = hz/hmag;               % 27) Wvec(3)
            this.Wvec = [Wx Wy Wz];
            this.Qx   = Wy*Pz - Wz*Py;         % 28) Qvec(1)
            this.Qy   = Wz*Px - Wx*Pz;         % 29) Qvec(2)
            this.Qz   = Wx*Py - Wy*Px;         % 30) Qvec(3)
            this.Qvec = [Qx Qy Qz];
            this.rUnitx = rx/rmag;             % 31) rUnit(1)
            this.rUnity = ry/rmag;             % 32) rUnit(2)
            this.rUnitz = rz/rmag;             % 33) rUnit(3)
            this.rUnit      = rUnit;
            this.Ax         = Ax;              % 34)
            this.Ay         = Ay;              % 35)
            this.normN      = normN;           % 36)
            this.Nx         = Nx;              % 37) Nvec(1)
            this.Ny         = Ny;              % 38) Nvec(2)
            this.Nvec       = Nvec;
            this.onePe      = onePe;           % 39)
            this.oneMe      = oneMe;           % 40)
            this.fac        = fac;             % 41)
            this.rootfac    = rootfac;         % 42)
            this.a          = a;               % 43)
            this.PeriodTU   = PeriodTU;        % 44) Period in TU's
            this.Period     = Period;          % 45) Work in TU's for now
            this.meanMotion = meanMotion;      % 46) Radians per TU.
            this.Omega      = Omega;           % 47) Omega = atan(Ny/Nx)
            this.cosO       = cosO;
            this.sinO       = sinO;
            this.cosP       = cosP;            % 48) cos(omega)
            this.sinP       = sinP;            % 49) sin(omega)
            this.omega      = omega;           % 50) omega = atan(sinP/cosP)
            this.cosom      = cosom;           % line  19
            this.sinom      = sinom;           % line  20
            this.cosnu      = cosnu;           % 51) cos(nu)
            this.sinnu      = sinnu;           % 52) sin(nu)
            this.nuEpoch    = nuEpoch;         % 53) nu = atan(sinnu/cosnu)
            this.cosT            = cosT;            % 54) ~cos(E)
            this.sinT            = sinT;            % 55) ~sin(E)
            this.EccentricAnomalyEpoch = EccentricAnomalyEpoch;    % 56) E = atan(sinE/cosE)
            this.cosE    = cosE;               % 57  (old 54)
            this.sinE    = sinE;               % 58  (old 55)
            this.MeanAnomalyEpoch   = MeanAnomalyEpoch;  % 59)
            %this.TimeSincePeriapsis = TimeSincePeriapsis;% 60)
            %this.DeltaTime   = DeltaTime;      % 61)
            %this.PTime       = PTime;          % 62) in TUs
            this.Mp          = Mp;             % 63 in radians
            this.Inclination = Inclination;    % 64
            this.cosI        = cosI;           %
            this.sinI        = sinI;           %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            this.M        = M;                                   % line  65
            this.EM       = EM;                                  % line  66
            this.Eprime   = Eprime;
            this.eDenom   = eDenom;
            this.eDenom   = eDenom;
            this.Eprime   = Eprime;
            this.cosK     = cosK;                                % line  67)
            this.sinK     = sinK;                                % line  68
            this.dE_dM    = dE_dM;
            this.dE_de    = dE_de;
            this.d2E_dMdM = d2E_dMdM;
            this.d2E_dMde = d2E_dMde;
            this.d2E_dedM = d2E_dedM;
            this.d2E_dede = d2E_dede;
            this.tanX     = tanX;                                % line  70
            this.tanY     = tanY;                                % line  71
            this.nu       = nu;                                  % line  72
            this.coss     = coss;                                % line  73
            this.sins     = sins;                                % line  74
            this.rorbit   = rorbit;                              % line  75
            this.rorbitx  = rorbitx;                             % line  76
            this.rorbity  = rorbity;                             % line  77
            this.rorbitz  = rorbitz;                             % line  78
            this.rtpinv   = rtpinv;                              % line  69
            this.vorbitx  = vorbitx;                             % line  79
            this.vorbity  = vorbity;                             % line  80
            this.vorbitz  = vorbitz;                             % line  81
            this.Rextrap  = Rextrap;
            this.Vextrap  = Vextrap;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            this.Aeq      = Aeq;                                 % line  82
            this.cosAeq   = cosAeq;                              % line  83
            this.sinAeq   = sinAeq;                              % line  84
            this.tanHalf  = tanHalf;                             % line  85
            this.Feq      = Feq;                                 % line  86
            this.Geq      = Geq;                                 % line  87
            this.Heq      = Heq;                                 % line  88
            this.Keq      = Keq;                                 % line  89
            this.Leq      = Leq;                                 % line  90
            this.CosL     = CosL;                                % line  91
            this.SinL     = SinL;                                % line  92
            this.alphaSq = alphaSq;                              % line  93
            this.Seq     = Seq;                                  % line  94
            this.Weq     = Weq;                                  % line  95
            this.Req     = Req;                                  % line  96
            this.RovS    = RovS;                                 % line  97
            this.srtpinv = srtpinv;                              % line  98
            this.HK      = HK;                                   % line  99
            this.OnePalphaSqCosL = OnePalphaSqCosL;              % line 100
            this.OneMalphaSqCosL = OneMalphaSqCosL;              % line 101
            this.OnePalphaSqSinL = OnePalphaSqSinL;              % line 102
            this.OneMalphaSqSinL = OneMalphaSqSinL;              % line 103
            this.Xfac    = Xfac;                                 % line 104
            this.Yfac    = Yfac;                                 % line 105
            this.Zfac    = Zfac;                                 % line 106
            this.VXfac   = VXfac;                                % line 107
            this.VYfac   = VYfac;                                % line 108
            this.VZfac   = VZfac;                                % line 109
            this.Xeq     = Xeq;                                  % line 110
            this.Yeq     = Yeq;                                  % line 111
            this.Zeq     = Zeq;                                  % line 112
            this.VXeq    = VXeq;                                 % line 113
            this.VYeq    = VYeq;                                 % line 114
            this.VZeq    = VZeq;                                 % line 115
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            this.DataList    = zeros(115,1);
            
            this.DataList( 1) = tEpoch;                 %  1) t
            this.DataList( 2) = rx;                     %  2) Rx
            this.DataList( 3) = ry;                     %  3) Ry
            this.DataList( 4) = rz;                     %  4) Rz
            this.DataList( 5) = vx;                     %  5) Vx
            this.DataList( 6) = vy;                     %  6) Vy
            this.DataList( 7) = vz;                     %  7) Vz
            this.DataList( 8) = tExtrap;                %  8) tExtrap
            
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
            this.DataList(24) = p;                      % 24) p = hVec^2/mu
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
            %this.DataList(60) = TimeSincePeriapsis;     % 60)
            %this.DataList(61) = DeltaTime;              % 61)
            %this.DataList(62) = PTime;  % DEFAULT       % 62) in TUs
            this.DataList(63) = Mp;                     % 63 in radians
            this.DataList(64) = Inclination;            % 64)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            this.DataList(65) = M;                      % 65)
            this.DataList(66) = EM;                     % 66) EM = E(e,M)
            this.DataList(67) = cosK;                   % 67)
            this.DataList(68) = sinK;                   % 68)
            this.DataList(69) = rtpinv;                 % 69
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
            this.DataList(82) = Aeq;                    % line 82
            this.DataList(83) = cosAeq;                 % line 83
            this.DataList(84) = sinAeq;                 % line 84
            this.DataList(85) = tanHalf;                % line 85
            this.DataList(86) = Feq;                    % line 86
            this.DataList(87) = Geq;                    % line 87
            this.DataList(88) = Heq;                    % line 88
            this.DataList(89) = Keq;                    % line 89
            this.DataList(90) = Leq;                    % line 90
            this.DataList(91) = CosL;                   % line 91 CosL
            this.DataList(92) = SinL;                   % line 92 SinL
            this.DataList(93) = alphaSq;                % line 93 alphaSq
            this.DataList(94) = Seq;                    % line 94
            this.DataList(95) = Weq;                    % line 95
            this.DataList(96) = Req;                    % line 96
            this.DataList(97) = RovS;                   % line 97
            this.DataList(98) = srtpinv;                % line 98
            this.DataList(99) = HK;                     % line 99
            this.DataList(100) = OnePalphaSqCosL;       % line 100
            this.DataList(101) = OneMalphaSqCosL;       % line 101
            this.DataList(102) = OnePalphaSqSinL;       % line 102
            this.DataList(103) = OneMalphaSqSinL;       % line 103
            this.DataList(104) = Xfac;                  % line 104
            this.DataList(105) = Yfac;                  % line 105
            this.DataList(106) = Zfac;                  % line 106
            this.DataList(107) = VXfac;                 % line 107
            this.DataList(108) = VYfac;                 % line 108
            this.DataList(109) = VZfac;                 % line 109
            this.DataList(110) = Xeq;                   % line 110
            this.DataList(111) = Yeq;                   % line 111
            this.DataList(112) = Zeq;                   % line 112
            this.DataList(113) = VXeq;                  % line 113
            this.DataList(114) = VYeq;                  % line 114
            this.DataList(115) = VZeq;                  % line 115
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %             % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             % % Checking First Derivatives
            %             % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             % format short
            %             % iStart = 1;
            %             % Nlines = 115;
            %             % CombinedDifferences = 0.0;
            %             % for LineCheck = iStart:Nlines
            %             %     ForDerivs  = [];
            %             %     %BackDerivs = [];
            %             %     %Derivs     = [];
            %             %     %Combined   = [];
            %             %
            %             %     F            = zeros(Nlines,1);
            %             %     F(LineCheck) = 1;
            %             %     F            = backdiff(this, F, LineCheck);
            %             %     % BackDerivs   = [F(1) F(2) F(3) F(4) F(5) F(6) F(7)];
            %             %     % Combined     = [Index F];
            %             %     % Kate's Combined'
            %             %
            %             %     forward      = [];
            %             %     for jk = 1:Nlines % check range of derivatives wrt forwards differentiation
            %             %         %jk
            %             %         Q        = zeros(Nlines,1);
            %             %         Q(jk)     = 1;
            %             %         Q        = fordiff(this, Q, LineCheck);
            %             %         % Q(LineCheck) = dK_LineCheck/dLj
            %             %         %forward  = [forward; j Q(LineCheck)];
            %             %         forward  = [forward; jk+1-iStart Q(LineCheck)];
            %             %         %Derivs   = [Derivs; Q(LineCheck)];
            %             %         %Combined = [Combined Q];
            %             %     end
            %             %     ForDerivs    = [ForDerivs; forward];
            %             %     Derivs       = [ForDerivs F (F-ForDerivs(:,2))];
            %             %     Derivs(iStart:LineCheck,:)
            %             %     Difference   = sum(F-ForDerivs(:,2));
            %             %     CombinedDifferences = CombinedDifferences + Difference;
            %             %     disp([' Line = ', num2str(LineCheck),'  Sum of Differences = ',num2str(Difference)]);
            %             %     % Kate's ForDerivs'
            %             %     % Derivs      = [F Derivs (F-Derivs)];
            %             %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             %     %     ForDerivs    = [ForDerivs; forward];
            %             %     %     Derivs = [ForDerivs F (F-ForDerivs(:,2))];
            %             %     %     Derivs = Derivs(1:LineCheck,:)
            %             %     %     %Derivs = [F Derivs (F-Derivs)];
            %             %     % 'Debug Point'
            %             % end
            %             % disp([' Combined Sum of Differences = ',num2str(CombinedDifferences)]);
            %             % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             % %%% TEST SECOND DERIVATIVES:
            %             % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             % %Debug Second Derivatives by checking symmetry:
            %             % %Begin Checking Out Second Derivatives:
            %             % Nvar     = 7;
            %             % % LineCheck    = 56;  % Kepler Part
            %             % % LineCheck    = 90;  % Escobal Gravitation to J6
            %             % % F            = zeros(115,1);
            %             % % F(LineCheck) = 1;
            %             % % F = backdiff(this, F, LineCheck);
            %             % HessianDifferences = 0.0;
            %             % Nlines = 115;
            %             % format short
            %             % for LineCheck = 1:Nlines
            %             %     ListMax   = LineCheck + 1;
            %             %     %Nvar = LineCheck;
            %             %     Hessian      = [];
            %             %     F            = zeros(Nlines,1);
            %             %     F(LineCheck) = 1;
            %             %     F = backdiff(this, F, LineCheck);
            %             %     %F(ListMax:Nlines) = 0;
            %             %     % for j = 1:Nvar
            %             %     for j = 1:LineCheck
            %             %         Q = zeros(Nlines,1);
            %             %         S = zeros(Nlines,1);
            %             %         Q(j) = 1.0;
            %             %         Q = fordiff(this, Q, LineCheck);
            %             %         %Q(ListMax:Nlines) = 0;
            %             %         S = secdiff(this, F,Q,S, LineCheck);
            %             %         %S(ListMax:Nlines) = 0;
            %             %         S = backdiff(this, S, LineCheck);
            %             %         %S(ListMax:Nlines) = 0;
            %             %         row = [];
            %             %         %for k = 1:Nvar
            %             %         for k = 1:LineCheck
            %             %             row = [row S(k)];
            %             %         end
            %             %         Hessian  = [Hessian; row];
            %             %     end
            %             %     %disp([' Line = ', num2str(Nvar)]);
            %             %     maxElement = max(max(abs(Hessian)));
            %             %     HessianTest = Hessian;
            %             %     if maxElement > 0
            %             %         HessianTest = Hessian/maxElement;
            %             %     end
            %             %
            %             %     Symmetry = HessianTest - HessianTest';
            %             %     Differences = max(max(abs(Symmetry)));
            %             %     HessianDifferences = HessianDifferences + Differences;
            %             %     disp([' Line Check = ', num2str(LineCheck),' max Hessian = ', num2str(maxElement),'  Summed Differences = ', num2str(Differences)]);
            %             %     %disp('Debug Point')
            %             % end
            %             % disp(['  Total Difference with Hessians = ',num2str(HessianDifferences)]);
            %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %             % JacobianAll = [];
            %             % for I = 1:100
            %             %     LineCheck    = I;                % e
            %             %     F            = zeros(115,1);
            %             %     F(LineCheck) = 1;
            %             %     F            = backdiff(this, F, LineCheck);
            %             %     JacobianAll  = [JacobianAll; F(2) F(3) F(4) F(5) F(6) F(7)];
            %             % end
            %             % this.JacobianAll = JacobianAll;
            %             % We will want to construct the 6 by 6 Jacobian of the Kepler
            %             % Orbital elements e, p, I, omega, Omega, PTime
            %             % or lines [17, 24, 64, 50, 47, 62]
            %             % and the parameters rx, ry, rz, vx, vy, vz
            %
            %             JacECI_2_Kepler = [];
            %             LineCheck    = 17;                % e
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             JacECI_2_Kepler     = [JacECI_2_Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 1
            %
            %             LineCheck    = 43;                % a
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             JacECI_2_Kepler     = [JacECI_2_Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 2
            %
            %             LineCheck    = 64;                % I
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             JacECI_2_Kepler     = [JacECI_2_Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 3
            %
            %             LineCheck    = 50;                % omega
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             JacECI_2_Kepler     = [JacECI_2_Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 4
            %
            %             LineCheck    = 47;                % Omega
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             JacECI_2_Kepler     = [JacECI_2_Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 5
            %
            %             LineCheck    = 63;                % Mp in Radians
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             JacECI_2_Kepler     = [JacECI_2_Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 6
            %
            %             this.JacECI_2_Kepler = JacECI_2_Kepler;
            %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %             JacECI_2_Equinoctial  = [];
            %             LineCheck    = 24;                % p
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 2
            %
            %             LineCheck    = 86;                % f
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 2
            %
            %             LineCheck    = 87;                % g
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 3
            %
            %             LineCheck    = 88;                % h
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 4
            %
            %             LineCheck    = 89;                % k
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 5
            %
            %             LineCheck    = 90;                % MLin Radians
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 6
            %
            %             this.JacECI_2_Equinoctial = JacECI_2_Equinoctial;
            %             this.JacKepler_2_Equinoctial = JacECI_2_Equinoctial*inv(JacECI_2_Kepler);
            %             % OrbitalCov   = Jacobian*Cov*Jacobian';  % transformed Covariances
            %             % 'Debug Point'
            %
            %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Try to Integrate the Equations of Motion for the Classical Elements
        % Y             = zeros(6,1);
        % KLagrange = LagrangePlanetary(KeplerIteration, this.units);
        % options       = odeset('RelTol', 1e-10, 'AbsTol', 1e-13);
        % dJdtCartesian = perturbations_odefun(KLagrange,tFit,Y);
        % Y            = zeros(6,1);
        % % Propagate from fit Reference time tFit back to tBeg
        % [Tout, Yout] = ode45(@(t,y) perturbations_odefun(KLagrange, t, y), [tFit, tBeg], Y, options);
        % %Yout(end,:)
        % Perturbed = [Tout, Yout];
        % Perturbed = sortrows(Perturbed,1);
        % Y = Perturbed(end,2:7);
        % % Propagate from fit Reference time tFit back to tEnd
        % [Tout, Yout] = ode45(@(t,y) perturbations_odefun(KLagrange, t, y), [tFit, tEnd], Y, options);
        % %Yout(end,:)
        % Perturbed(end,:) = [];
        % Perturbed = [Perturbed; Tout, Yout];
        % Perturbed = sortrows(Perturbed,1);
        % KeplerInterpolated = interp1(Perturbed(:,1),Perturbed(:,2:7),timeArray);
        % %Implementation
        % Zhat = KeplerInterpolated(ii,:);
        % KLagrange     = LagrangePlanetary(KeplerIteration+Zhat, this.units);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function Okay = Extrapolator(this, tExtrap)
            Okay = 0;
            wgs84     = wgs84Constants();
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Canonical Coordinates assumed throughout -- length in DU, time in TU
            % Modifications can be made to the choice of units in wgs84Constants!
            % Orbital elements e, p, I, omega, Omega, PTime
            % or lines [17, 24, 64, 50, 47, 62]
            % and the parameters rx, ry, rz, vx, vy, vz
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Constants
            %units = wgs84Constants;
            twopi       = this.twopi;
            % TU          = units.TU;  % Canonical Time UnitTime Unit in Seconds
            % DU          = units.DU;  % Canonical Distance Unit
            % VU          = units.VU;  % Canonical Velocity Unit
            % AU          = units.AU;        % Canonical Acceleration Unit
            % mu          = units.mu;  % Canonical Gravitational Constant
            mu            = this.mu;
            sqrtmu        = this.sqrtmu;
            
            Ivec  = [1.0 0.0 0.0];
            Jvec  = [0.0 1.0 0.0];
            Kvec  = [0.0 0.0 1.0];
            
            sinO = this.sinO;
            cosO = this.cosO;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %this.Rpos = Rpos;
            %this.Rdot = Rdot;
            %t    = tExtrap;                                     % line   1) t
            SignOrbit   = this.SignOrbit; 
            rx          = this.DataList( 2); % = rx;             % line   2) Rx
            ry          = this.DataList( 3); % = ry;             % line   3) Ry
            rz          = this.DataList( 4); % = rz;             % line   4) Rz
            vx          = this.DataList( 5); % = vx;             % line   5) Vx
            vy          = this.DataList( 6); % = vy;             %  line  6) Vy
            vz          = this.DataList( 7); % = vz;             % line   7) Vz
            this.time         = tExtrap;
            %this.tExtrap      = tExtrap;
            this.DataList( 8) = tExtrap;
            
            rmag        = this.DataList(10); % = rmag;           % line  10) rmag
            vsq         = this.DataList(11); % = vsq;            % line  11) vsq
            er          = this.DataList(12); % = er;             % line  12) er
            ev          = this.DataList(13); % = ev;             % line  13) ev
            ex          = this.DataList(14); % = ex;             % line  14) Evec(1)
            ey          = this.DataList(15); % = ey;             % line  15) Evec(2)
            ez          = this.DataList(16); % = ez;             % line  16) Evec(3)
            e           = this.DataList(17); % = e;              % line  17) e;
            Px          = this.DataList(18); % = Px;             % line  18) Pvec(1)
            Py          = this.DataList(19); % = Py;             % line  19) Pvec(2)
            Pz          = this.DataList(20); % = Pz;             % line  20) Pvec(3)
            hx          = this.DataList(21); % = hx;             % line  21) hVec(1)
            hy          = this.DataList(22); % = hy;             % line  22) hVec(2)
            hz          = this.DataList(23); % = hz;             % line  23) hVec(3)
            p           = this.DataList(24); % = p;;             % line  24) p = hVec^2/mu
            Wx          = this.DataList(25); % = Wx;             % line  25) Wvec(1)
            Wy          = this.DataList(26); % = Wy;             % line  26) Wvec(2)
            Wz          = this.DataList(27); % = Wz;             % line  27) Wvec(3)
            Qx          = this.DataList(28); % = Qx;             % line  28) Qvec(1)
            Qy          = this.DataList(29); % = Qy;             % line  29) Qvec(2)
            Qz          = this.DataList(30); % = Qz;             % line  30) Qvec(3)
            rUnitx      = this.DataList(31); % = rUnitx;         % line  31) rUnit(1)
            rUnity      = this.DataList(32); % = rUnity;         % line  32) rUnit(2)
            rUnitz      = this.DataList(33); % = rUnitz;         % line  33) rUnit(3)
            Ax          = this.DataList(34); % = Ax;             % line  34) Ax
            Ay          = this.DataList(35); % = Ay;             % line  35) Ay
            normN       = this.DataList(36); % = normN;          % line  36) normN
            Nx          = this.DataList(37); % = Nx;             % line  37) Nvec(1)
            Ny          = this.DataList(38); % = Ny;             % line  38) Nvec(2)
            onePe       = this.DataList(39); % = onePe;          % line  39)
            oneMe       = this.DataList(40); % = oneMe;          % line  40)
            fac         = this.DataList(41); % = fac;            % line  41)
            rootfac     = this.DataList(42); % = rootfac;        % line  42)
            a           = this.DataList(43); % = a;              % line  43)
            PeriodTU    = this.DataList(44); % = PeriodTU;       % line  44) Period in TU's
            Period      = this.DataList(45); % = Period;         % line  45) Work in TU's for now
            meanMotion  = this.DataList(46); % = meanMotion;     % line  46) Radians per TU.
            Omega       = this.DataList(47); % = Omega;          % line  47) Omega = atan(Ny/Nx);
            cosP        = this.DataList(48); % = cosP;           % line  48 cos(omega)
            sinP        = this.DataList(49); % = sinP;           % line  49 sin(omega)
            omega       = this.DataList(50); % = omega;          % line  50 omega = atan(sinP/cosP)
            cosnu       = this.DataList(51); % = cosnu;          % line  51 cos(nu)
            sinnu       = this.DataList(52); % = sinnu;          % line  52 sin(nu)
            nuEpoch     = this.DataList(53); % = nuEpoch;        % line  53 nu = atan(sinnu/cosnu)
            cosT        = this.DataList(54); % = cosT;           % line  54 ~cos(E)
            sinT        = this.DataList(55); % = sinT;           % line  55 ~sin(E)
            EccentricAnomalyEpoch = this.DataList(56); %= EccentricAnomalyEpoch;  % 56) E = atan(sinE/cosE)
            cosE        = this.DataList(57); % = cosE;           %  line 57  (old 54)
            sinE        = this.DataList(58); % = sinE;           %  line 58  (old 55)
            MeanAnomalyEpoch   = this.DataList(59); % = MeanAnomalyEpoch;    % line  59)
            %TimeSincePeriapsis = this.DataList(60); % = TimeSincePeriapsis; % line  60)
            %DeltaTime          = this.DataList(61); % = DeltaTime;          % line  61)
            %PTime       = this.DataList(62);% = PTime;  % DEFAULT %line 62 in TUs
            Mp          = this.DataList(63); % = Mp;             %  line 63 in radians
            Inclination = this.DataList(64); % = Inclination;    %  line 64
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            %M = meanMotion*t + MeanAnomalyEpoch;             % Orig     65)
            % JRS TEST: Uses the Current MeanAnomaly at this event EPOCH.
            %M = MeanAnomalyEpoch;                          % Test line  65)
            M = Mp + meanMotion*tExtrap;                      % New line 65
            while M < 0
                M = M + twopi;
            end
            while M > twopi
                M = M - twopi;
            end
            %  meanAnomaly M to eccentric Anomaly EM
            if e < 1
                EM = M + e;
                if (M > pi)
                    EM = M - e;
                elseif ( (-pi < M) && (M < 0) )
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
                % Solve Kepler's Equation for E:  M = EM - e*sin(EM) % line  66 EM = E(e,M)
                % Need to differentiate Implicitly here!
                eDenom   = Eprime*Eprime*Eprime;
                this.Eprime = Eprime;
                this.eDenom = eDenom;
                
                % KeplerInv = EM;
                % //cout  << " Newton-Raphson Kepler's Equation Eccentric Anomaly =" << E << endl;
                % eDenom   = Eprime*Eprime*Eprime;
                
                %rmag = p/(1 + e*cosnu);
                %rRecon = rmag*(cosnu*Pvec + sinnu*Qvec);
                %QReco = rmag*cosnu*Rdot + (1.0/sqrt(p))*sinnu*Rpos;
                %QReco  = (1.0/norm(QReco))*QReco;
                
                cosK  = cos(EM);                                     % line  67
                sinK  = sin(EM);                                     % line  68
                %E     = atan2(sinK,cosK)
                %%%%%%%%%%%%%%%%%%% Implicit Function Theorem Derivatives
                dE_dM    = 1.0/Eprime;
                dE_de    = sinK/Eprime;
                d2E_dMdM = -e*sinK/eDenom;
                %// Mike Cain Corrections!
                d2E_dMde = (cosK - e)/eDenom;
                d2E_dedM = (cosK - e)/eDenom;
                d2E_dede = ((2.0 - e*cosK)*cosK*sinK - e*sinK)/eDenom;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                this.eDenom   = eDenom;
                this.Eprime   = Eprime;
                this.dE_dM    = dE_dM;
                this.dE_de    = dE_de;
                this.d2E_dMdM = d2E_dMdM;
                this.d2E_dMde = d2E_dMde;
                this.d2E_dedM = d2E_dedM;
                this.d2E_dede = d2E_dede;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % kfac  = 1.0 - e*cosK;                              % line  69
                % coss  = (cosK - e)/kfac;                           % line  70
                % sins  = rootfac*sinK/kfac;                         % line  71
                % 	s     = atan2(sins, coss);                       % line  72
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   % for convenience -- you can eliminate kfac;
                tanX  = cosK - e;                                    % line  70
                tanY  = rootfac*sinK;                                % line  71
            else
                if e == 1
                    % Reserved for Parabolic Case
                else
                    if e < 1.6
                        if ((-pi < M) && (M<0)) || M>pi
                            EM = M -e;
                        else
                            EM = M + e;
                        end
                    else
                        if (e < 3.6) &&  abs(M) >pi
                            EM = M - sign(M)*e
                        else
                            EM = M/(e - 1.);
                        end
                    end
                    
                    %  meanAnomaly M to eccentric Anomaly EM
                    %         //cout << " Input meanAnomaly = " << meanAnomaly << endl;
                    for i=1:1:10
                        EM  = e*sinh(EM) - M;
                        %//cout << " Mean Anomaly Solution " << E << endl;
                    end
                    % //      10 rounds of Newton's root finding method based on the above "seed".
                    for i=1:1:10
                        Eprime      = e*cosh(EM) - 1.0;
                        EM          = EM + (M + EM -e*sinh(EM))/Eprime;
                    end
                    % Solve Kepler's Equation for E:  M = EM - e*sin(EM) % line  66 EM = E(e,M)
                    % Need to differentiate Implicitly here!
                    eDenom   = Eprime*Eprime*Eprime;
                    this.Eprime = Eprime;
                    this.eDenom = eDenom;
                    
                    % KeplerInv = EM;
                    % //cout  << " Newton-Raphson Kepler's Equation Eccentric Anomaly =" << E << endl;
                    % eDenom   = Eprime*Eprime*Eprime;
                    
                    %rmag = p/(1 + e*cosnu);
                    %rRecon = rmag*(cosnu*Pvec + sinnu*Qvec);
                    %QReco = rmag*cosnu*Rdot + (1.0/sqrt(p))*sinnu*Rpos;
                    %QReco  = (1.0/norm(QReco))*QReco;
                    
                    cosK  = cosh(EM);                                     % line  67
                    sinK  = sinh(EM);                                     % line  68
                    %E     = atan2(sinK,cosK)
                    %%%%%%%%%%%%%%%%%%% Implicit Function Theorem Derivatives
                    dE_dM    =  1.0/Eprime;
                    dE_de    = -sinK/Eprime;
                    d2E_dMdM = -e*sinK/eDenom;
                    %// Mike Cain Corrections!
                    d2E_dMde = (cosK - e)/eDenom;
                    d2E_dedM = (cosK - e)/eDenom;
                    d2E_dede = ((e*cosK -2.0)*cosK*sinK + e*sinK)/eDenom;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    this.eDenom   = eDenom;
                    this.Eprime   = Eprime;
                    this.dE_dM    = dE_dM;
                    this.dE_de    = dE_de;
                    this.d2E_dMdM = d2E_dMdM;
                    this.d2E_dMde = d2E_dMde;
                    this.d2E_dedM = d2E_dedM;
                    this.d2E_dede = d2E_dede;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % kfac  = 1.0 - e*cosK;                              % line  69
                    % coss  = (cosK - e)/kfac;                           % line  70
                    % sins  = rootfac*sinK/kfac;                         % line  71
                    % 	s     = atan2(sins, coss);                       % line  72
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %   % for convenience -- you can eliminate kfac;
                    tanX  = cosK - e;                                    % line  70
                    tanY  = -rootfac*sinK;                                % line  71
                end
            end
            nu     = atan2(tanY, tanX);                          % line  72
            coss  = cos(nu);                                     % line  73
            sins  = sin(nu);                                     % line  74
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   Position and Unit Vector along the Position
            %Kvec  = a*kfac*(cosnu*Pvec + sinnu*Qvec);
            %Runit = Kvec/norm(Kvec);
            
            %   Unit Vector along the Velocity
            %Vunit =(-sinnu*Pvec + (e + cosnu)*Qvec);
            %Vunit = Vunit/norm(Vunit);
            
            %   Unit Vector out of the R-V plane
            %Wlocal = cross(Runit, Vunit);
            %Wlocal = Wlocal/norm(Wlocal);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Let us propagate this Orbit to time = tExtrap
            %rorbit   = p/(1.0 + e*coss);                        % line  75
            rorbit   = a*(1.0 - e*cosK);
            
            rorbitx  = rorbit*(coss*Px + sins*Qx);               % line  76
            rorbity  = rorbit*(coss*Py + sins*Qy);               % line  77
            rorbitz  = rorbit*(coss*Pz + sins*Qz);               % line  78
            rtpinv   = sqrt(mu/p);                               % line  69
            vorbitx  = rtpinv*(-sins*Px + (e + coss)*Qx);        % line  79
            vorbity  = rtpinv*(-sins*Py + (e + coss)*Qy);        % line  80
            vorbitz  = rtpinv*(-sins*Pz + (e + coss)*Qz);        % line  81
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Rextrap  = [rorbitx; rorbity; rorbitz];
            Vextrap  = [vorbitx; vorbity; vorbitz];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Treat Equinoctial Variables:
            Aeq      = omega + Omega;                            % line  82
            while Aeq < 0
               Aeq = Aeq + twopi;
            end
            while Aeq > twopi
               Aeq = Aeq - twopi;
            end
            cosAeq   = cos(Aeq);                                 % line  83
            sinAeq   = sin(Aeq);                                 % line  84
            tanHalf  = tan(Inclination/2.0);                     % line  85
            Feq      = e*cosAeq;                                 % line  86
            Geq      = e*sinAeq;                                 % line  87
            Heq      = tanHalf*cosO;                             % line  88
            Keq      = tanHalf*sinO;                             % line  89
            Leq      = Aeq + nu;                                 % line  90
            CosL     = cos(Leq);                                 % line  91
            SinL     = sin(Leq);                                 % line  92
            alphaSq = Heq^2 - Keq^2;                             % line  93
            Seq     = 1.0 + Heq^2 + Keq^2;                       % line  94
            Weq     = 1.0 + Feq*CosL + Geq*SinL;                 % line  95
            Req     = p/Weq;                                     % line  96
            RovS    = Req/Seq;                                   % line  97
            srtpinv = rtpinv/Seq;                                % line  98
            HK      = Heq*Keq;                                   % line  99
            OnePalphaSqCosL = (1.+alphaSq)*CosL;                 % line 100
            OneMalphaSqCosL = (1.-alphaSq)*CosL;                 % line 101
            OnePalphaSqSinL = (1.+alphaSq)*SinL;                 % line 102
            OneMalphaSqSinL = (1.-alphaSq)*SinL;                 % line 103
            Xfac    = OnePalphaSqCosL + 2.0*HK*SinL;             % line 104
            Yfac    = OneMalphaSqSinL + 2.0*HK*CosL;             % line 105
            Zfac    = Heq*SinL - Keq*CosL;                       % line 106
            VXfac   =  OnePalphaSqSinL - 2.0*HK*(Feq + CosL) + Geq*(1. + alphaSq);% line 107
            VYfac   = -OneMalphaSqCosL + 2.0*HK*(Geq + SinL) + Feq*(alphaSq - 1.);% line 108
            VZfac   =  Heq*(Feq+CosL) + Keq*(Geq+SinL);          % line 109
            Xeq     =     RovS*Xfac;                             % line 110
            Yeq     =     RovS*Yfac;                             % line 111
            Zeq     = 2.0*RovS*Zfac;                             % line 112
            VXeq    =    -srtpinv*VXfac;                         % line 113
            VYeq    =    -srtpinv*VYfac;                         % line 114
            VZeq    = 2.0*srtpinv*VZfac;                         % line 115
            
            RVeq    = [Xeq; Yeq; Zeq];
            %Rpos'
            VVeq    = [VXeq; VYeq; VZeq];
            %Rdot'
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            this.M        = M;                                   % line  65
            this.EM       = EM;                                  % line  66
            this.Eprime   = Eprime;
            this.eDenom   = eDenom;
            this.eDenom   = eDenom;
            this.Eprime   = Eprime;
            this.cosK     = cosK;                                % line  67)
            this.sinK     = sinK;                                % line  68
            this.dE_dM    = dE_dM;
            this.dE_de    = dE_de;
            this.d2E_dMdM = d2E_dMdM;
            this.d2E_dMde = d2E_dMde;
            this.d2E_dedM = d2E_dedM;
            this.d2E_dede = d2E_dede;
            this.tanX     = tanX;                                % line  70
            this.tanY     = tanY;                                % line  71
            this.nu       = nu;                                  % line  72
            this.coss     = coss;                                % line  73
            this.sins     = sins;                                % line  74
            this.rorbit   = rorbit;                              % line  75
            this.rorbitx  = rorbitx;                             % line  76
            this.rorbity  = rorbity;                             % line  77
            this.rorbitz  = rorbitz;                             % line  78
            this.rtpinv   = rtpinv;                              % line  69
            this.vorbitx  = vorbitx;                             % line  79
            this.vorbity  = vorbity;                             % line  80
            this.vorbitz  = vorbitz;                             % line  81
            this.Rextrap  = Rextrap;
            this.Vextrap  = Vextrap;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            this.Aeq      = Aeq;                                 % line  82
            this.cosAeq   = cosAeq;                              % line  83
            this.sinAeq   = sinAeq;                              % line  84
            this.tanHalf  = tanHalf;                             % line  85
            this.Feq      = Feq;                                 % line  86
            this.Geq      = Geq;                                 % line  87
            this.Heq      = Heq;                                 % line  88
            this.Keq      = Keq;                                 % line  89
            this.Leq      = Leq;                                 % line  90
            this.CosL     = CosL;                                % line  91
            this.SinL     = SinL;                                % line  92
            this.alphaSq = alphaSq;                              % line  93
            this.Seq     = Seq;                                  % line  94
            this.Weq     = Weq;                                  % line  95
            this.Req     = Req;                                  % line  96
            this.RovS    = RovS;                                 % line  97
            this.srtpinv = srtpinv;                              % line  98
            this.HK      = HK;                                   % line  99
            this.OnePalphaSqCosL = OnePalphaSqCosL;              % line 100
            this.OneMalphaSqCosL = OneMalphaSqCosL;              % line 101
            this.OnePalphaSqSinL = OnePalphaSqSinL;              % line 102
            this.OneMalphaSqSinL = OneMalphaSqSinL;              % line 103
            this.Xfac    = Xfac;                                 % line 104
            this.Yfac    = Yfac;                                 % line 105
            this.Zfac    = Zfac;                                 % line 106
            this.VXfac   = VXfac;                                % line 107
            this.VYfac   = VYfac;                                % line 108
            this.VZfac   = VZfac;                                % line 109
            this.Xeq     = Xeq;                                  % line 110
            this.Yeq     = Yeq;                                  % line 111
            this.Zeq     = Zeq;                                  % line 112
            this.VXeq    = VXeq;                                 % line 113
            this.VYeq    = VYeq;                                 % line 114
            this.VZeq    = VZeq;                                 % line 115
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            this.StateVector7        = [Rextrap; Vextrap; tExtrap];
            this.ClassicalElements   = [e; a; Inclination; omega; Omega; Mp; tExtrap];
            this.EquinoctialElements = [p; Feq; Geq; Heq; Keq; Leq; tExtrap];
            
            %[XEquin, XKepler, XECI] =  ECI_2_All(this.StateVector7, mu)
            %DataList1 = this.DataList;
            
            this.DataList(65) = M;                      % 65)
            this.DataList(66) = EM;                     % 66) EM = E(e,M)
            this.DataList(67) = cosK;                   % 67)
            this.DataList(68) = sinK;                   % 68)
            this.DataList(69) = rtpinv;                 % 69
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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            this.DataList(82) = Aeq;                    % line 82
            this.DataList(83) = cosAeq;                 % line 83
            this.DataList(84) = sinAeq;                 % line 84
            this.DataList(85) = tanHalf;                % line 85
            this.DataList(86) = Feq;                    % line 86
            this.DataList(87) = Geq;                    % line 87
            this.DataList(88) = Heq;                    % line 88
            this.DataList(89) = Keq;                    % line 89
            this.DataList(90) = Leq;                    % line 90
            this.DataList(91) = CosL;                   % line 91 CosL
            this.DataList(92) = SinL;                   % line 92 SinL
            this.DataList(93) = alphaSq;                % line 93 alphaSq
            this.DataList(94) = Seq;                    % line 94
            this.DataList(95) = Weq;                    % line 95
            this.DataList(96) = Req;                    % line 96
            this.DataList(97) = RovS;                   % line 97
            this.DataList(98) = srtpinv;                % line 98
            this.DataList(99) = HK;                     % line 99
            this.DataList(100) = OnePalphaSqCosL;       % line 100
            this.DataList(101) = OneMalphaSqCosL;       % line 101
            this.DataList(102) = OnePalphaSqSinL;       % line 102
            this.DataList(103) = OneMalphaSqSinL;       % line 103
            this.DataList(104) = Xfac;                  % line 104
            this.DataList(105) = Yfac;                  % line 105
            this.DataList(106) = Zfac;                  % line 106
            this.DataList(107) = VXfac;                 % line 107
            this.DataList(108) = VYfac;                 % line 108
            this.DataList(109) = VZfac;                 % line 109
            this.DataList(110) = Xeq;                   % line 110
            this.DataList(111) = Yeq;                   % line 111
            this.DataList(112) = Zeq;                   % line 112
            this.DataList(113) = VXeq;                  % line 113
            this.DataList(114) = VYeq;                  % line 114
            this.DataList(115) = VZeq;                  % line 115
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %DataList2 = this.DataList;
            
            
            ECIPropagator = [];
            LineCheck    = 76;                % Rx
            F            = zeros(115,1);
            F(LineCheck) = 1;
            F            = backdiff(this, F, LineCheck);
            ECIPropagator     = [ECIPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(8)];   % row 1
            
            LineCheck    = 77;                % Ry
            F            = zeros(115,1);
            F(LineCheck) = 1;
            F            = backdiff(this, F, LineCheck);
            ECIPropagator     = [ECIPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(8)];   % row 2
            
            LineCheck    = 78;                % Ry
            F            = zeros(115,1);
            F(LineCheck) = 1;
            F            = backdiff(this, F, LineCheck);
            ECIPropagator     = [ECIPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(8)];   % row 3
            
            LineCheck    = 79;                % Vx
            F            = zeros(115,1);
            F(LineCheck) = 1;
            F            = backdiff(this, F, LineCheck);
            ECIPropagator     = [ECIPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(8)];   % row 4
            
            LineCheck    = 80;                % Vy
            F            = zeros(115,1);
            F(LineCheck) = 1;
            F            = backdiff(this, F, LineCheck);
            ECIPropagator     = [ECIPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(8)];   % row 5
            
            LineCheck    = 81;                % Vz
            F            = zeros(115,1);
            F(LineCheck) = 1;
            F            = backdiff(this, F, LineCheck);
            ECIPropagator     = [ECIPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(1)];   % row 6
            
            LineCheck    =  8;                % Vz
            F            = zeros(115,1);
            F(LineCheck) = 1;
            F            = backdiff(this, F, LineCheck);
            ECIPropagator     = [ECIPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(1)];   % row 6
            %ECIPropagator
            this.ECIPropagator = ECIPropagator;
            
            % %
            % JacECI_2_Kepler = [];
            % LineCheck    = 17;                % e
            % F            = zeros(115,1);
            % F(LineCheck) = 1;
            % F            = backdiff(this, F, LineCheck);
            % JacECI_2_Kepler     = [JacECI_2_Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 1
            %
            % LineCheck    = 43;                % a
            % F            = zeros(115,1);
            % F(LineCheck) = 1;
            % F            = backdiff(this, F, LineCheck);
            % JacECI_2_Kepler     = [JacECI_2_Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 2
            %
            % LineCheck    = 64;                % I
            % F            = zeros(115,1);
            % F(LineCheck) = 1;
            % F            = backdiff(this, F, LineCheck);
            % JacECI_2_Kepler     = [JacECI_2_Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 3
            %
            % LineCheck    = 50;                % omega
            % F            = zeros(115,1);
            % F(LineCheck) = 1;
            % F            = backdiff(this, F, LineCheck);
            % JacECI_2_Kepler     = [JacECI_2_Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 4
            %
            % LineCheck    = 47;                % Omega
            % F            = zeros(115,1);
            % F(LineCheck) = 1;
            % F            = backdiff(this, F, LineCheck);
            % JacECI_2_Kepler     = [JacECI_2_Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 5
            %
            % LineCheck    = 63;                % Mp in Radians
            % F            = zeros(115,1);
            % F(LineCheck) = 1;
            % F            = backdiff(this, F, LineCheck);
            % JacECI_2_Kepler     = [JacECI_2_Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 6
            % JacECI_2_Kepler
            %
            % %this.JacECI_2_Kepler = JacECI_2_Kepler;
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % JacECI_2_Equinoctial  = [];
            % LineCheck    = 24;                % p
            % F            = zeros(115,1);
            % F(LineCheck) = 1;
            % F            = backdiff(this, F, LineCheck);
            % JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 2
            %
            % LineCheck    = 86;                % f
            % F            = zeros(115,1);
            % F(LineCheck) = 1;
            % F            = backdiff(this, F, LineCheck);
            % JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 2
            %
            % LineCheck    = 87;                % g
            % F            = zeros(115,1);
            % F(LineCheck) = 1;
            % F            = backdiff(this, F, LineCheck);
            % JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 3
            %
            % LineCheck    = 88;                % h
            % F            = zeros(115,1);
            % F(LineCheck) = 1;
            % F            = backdiff(this, F, LineCheck);
            % JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 4
            %
            % LineCheck    = 89;                % k
            % F            = zeros(115,1);
            % F(LineCheck) = 1;
            % F            = backdiff(this, F, LineCheck);
            % JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 5
            %
            % LineCheck    = 90;                % MLin Radians
            % F            = zeros(115,1);
            % F(LineCheck) = 1;
            % F            = backdiff(this, F, LineCheck);
            % JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 6
            % JacECI_2_Equinoctial
            %
            % %this.JacECI_2_Equinoctial    = JacECI_2_Equinoctial;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %             KeplerPropagator = [];
            %             LineCheck    = 2;                % e
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             KeplerPropagator = [KeplerPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(1)];   % row 1
            %
            %             LineCheck    = 3;                % a
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             KeplerPropagator     = [KeplerPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(1)];   % row 2
            %
            %             LineCheck    = 4;                % I
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             KeplerPropagator     = [KeplerPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(1)];   % row 3
            %
            %             LineCheck    = 5;                % omega
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             KeplerPropagator     = [KeplerPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(1)];   % row 4
            %
            %             LineCheck    = 6;                % Omega
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             KeplerPropagator     = [KeplerPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(1)];   % row 5
            %
            %             LineCheck    = 65;                % M
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             KeplerPropagator     = [KeplerPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(1)];   % row 6
            %
            %             LineCheck    = 8;                % time
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             KeplerPropagator     = [KeplerPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(1)];   % row 6
            %             this.KeplerPropagator = KeplerPropagator;
            %
            %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %             EquinPropagator = [];
            %             LineCheck    = 76;                % p
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             EquinPropagator     = [EquinPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(1)];   % row 1
            %
            %             LineCheck    = 77;                % y
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             EquinPropagator     = [EquinPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(1)];   % row 2
            %
            %             LineCheck    = 78;                % Ry
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             EquinPropagator     = [EquinPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(1)];   % row 3
            %
            %             LineCheck    = 79;                % Vx
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             EquinPropagator     = [EquinPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(1)];   % row 4
            %
            %             LineCheck    = 80;                % Vy
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             EquinPropagator     = [EquinPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(1)];   % row 5
            %
            %             LineCheck    = 81;                % Vz
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             EquinPropagator     = [EquinPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(1)];   % row 6
            %
            %             LineCheck    =  8;                % Vz
            %             F            = zeros(115,1);
            %             F(LineCheck) = 1;
            %             F            = backdiff(this, F, LineCheck);
            %             EquinPropagator     = [EquinPropagator; F(2) F(3) F(4) F(5) F(6) F(7), F(1)];   % row 6
            %             this.EquinPropagator = EquinPropagator;
            
            %DataList = this.DataList(1:81);
            %Propagation   = [DataList(76),DataList(77),DataList(78),DataList(79),DataList(80),DataList(81), DataList( 8);...
            %                 DataList( 2),DataList( 3),DataList( 4),DataList( 5),DataList( 6),DataList( 7), DataList( 1)]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            this.ExtrapolatedClass  = KeplerFromECI(tExtrap, Rextrap, Vextrap, this.units);
            this.ExtrapolatedClass.DataList(65:115)   = this.DataList(65:115);
            this.ExtrapolatedClass.StateVector7       = [Rextrap; Vextrap; tExtrap];
            this.ExtrapolatedClass.ClassicalElements  = [e; a; Inclination; omega; Omega; Mp; tExtrap];
            this.ExtrapolatedClass.EquinoctialElements = [p; Feq; Geq; Heq; Keq; Leq; tExtrap];
            this.ExtrapolatedClass.M        = M;                            % 65
            this.ExtrapolatedClass.EM       = EM;                           % 66
            this.ExtrapolatedClass.Eprime   = Eprime;
            this.ExtrapolatedClass.eDenom   = eDenom;
            this.ExtrapolatedClass.eDenom   = eDenom;
            this.ExtrapolatedClass.Eprime   = Eprime;
            this.ExtrapolatedClass.cosK     = cosK;                                        % 67)
            this.ExtrapolatedClass.sinK     = sinK;                                        % 68)
            this.ExtrapolatedClass.dE_dM    = dE_dM;
            this.ExtrapolatedClass.dE_de    = dE_de;
            this.ExtrapolatedClass.d2E_dMdM = d2E_dMdM;
            this.ExtrapolatedClass.d2E_dMde = d2E_dMde;
            this.ExtrapolatedClass.d2E_dedM = d2E_dedM;
            this.ExtrapolatedClass.d2E_dede = d2E_dede;
            this.ExtrapolatedClass.tanX     = tanX;                                     % line 70
            this.ExtrapolatedClass.tanY     = tanY;                                 % line 71
            this.ExtrapolatedClass.nu       = nu;                           % line 72
            this.ExtrapolatedClass.coss     = coss;                                      % line 73
            this.ExtrapolatedClass.sins     = sins;                                      % line 74
            this.ExtrapolatedClass.rorbit   = rorbit;                          % line 75
            this.ExtrapolatedClass.rorbitx  = rorbitx;                % line 76
            this.ExtrapolatedClass.rorbity  = rorbity;                % line 77
            this.ExtrapolatedClass.rorbitz  = rorbitz;                % line 78
            this.ExtrapolatedClass.rtpinv   = rtpinv;                                % line 69
            this.ExtrapolatedClass.vorbitx  = vorbitx;         % line 79
            this.ExtrapolatedClass.vorbity  = vorbity;         % line 80
            this.ExtrapolatedClass.vorbitz  = vorbitz;         % line 81
            this.ExtrapolatedClass.Rextrap  = Rextrap;
            this.ExtrapolatedClass.Vextrap  = Vextrap;
            this.ExtrapolatedClass.Aeq      = Aeq;                             % line 82
            this.ExtrapolatedClass.cosAeq   = cosAeq;                                  % line 83
            this.ExtrapolatedClass.sinAeq   = sinAeq;                                  % line 84
            this.ExtrapolatedClass.tanHalf  = tanHalf;                      % line 85
            this.ExtrapolatedClass.Feq      = Feq;                                  % line 86
            this.ExtrapolatedClass.Geq      = Geq;                                  % line 87
            this.ExtrapolatedClass.Heq      = Heq;                          % line 88
            this.ExtrapolatedClass.Keq      = Keq;                          % line 89
            this.ExtrapolatedClass.Leq      = Leq;                                  % line 90
            this.ExtrapolatedClass.CosL     = CosL;                                  % line 91
            this.ExtrapolatedClass.SinL     = SinL;                                  % line 92
            this.ExtrapolatedClass.alphaSq = alphaSq;                              % line 93
            this.ExtrapolatedClass.Seq     = Seq;                        % line 94
            this.ExtrapolatedClass.Weq     = Weq;                  % line 95
            this.ExtrapolatedClass.Req     = Req;                                      % line 96
            this.ExtrapolatedClass.RovS    = RovS;                                    % line 97
            this.ExtrapolatedClass.srtpinv = srtpinv;                                 % line 98
            this.ExtrapolatedClass.HK      = HK;                                    % line 99
            this.ExtrapolatedClass.OnePalphaSqCosL = OnePalphaSqCosL;                 % line 100
            this.ExtrapolatedClass.OneMalphaSqCosL = OneMalphaSqCosL;                 % line 101
            this.ExtrapolatedClass.OnePalphaSqSinL = OnePalphaSqSinL;                 % line 102
            this.ExtrapolatedClass.OneMalphaSqSinL = OneMalphaSqSinL;                 % line 103
            this.ExtrapolatedClass.Xfac    = Xfac;             % line 104
            this.ExtrapolatedClass.Yfac    = Yfac;             % line 105
            this.ExtrapolatedClass.Zfac    = Zfac;                       % line 106
            this.ExtrapolatedClass.VXfac   = VXfac;% line 107
            this.ExtrapolatedClass.VYfac   = VYfac;% line 108
            this.ExtrapolatedClass.VZfac   = VZfac;  % line 109
            this.ExtrapolatedClass.Xeq     = Xeq;                             % line 110
            this.ExtrapolatedClass.Yeq     = Yeq;                             % line 111
            this.ExtrapolatedClass.Zeq     = Zeq;                             % line 112
            this.ExtrapolatedClass.VXeq    = VXeq;                         % line 113
            this.ExtrapolatedClass.VYeq    = VYeq;                         % line 114
            this.ExtrapolatedClass.VZeq    = VZeq;                         % line 115
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            epsilon     = 0.000001*ones(6,1);
            DummySensor = [this.StateVector7(1:3) + epsilon(1:3); zeros(6,1)];
            DummyLos    = [this.StateVector7(1:3) - 10000.*epsilon(1:3)];
            DummyLos    = DummyLos/norm(DummyLos);
            InvCovUV    = eye(2);
            [Rinit,Vinit,ParamList,JacKepler_2_ECI,Hessian,GravityCan] = this.KeplerToECI.OrbitDerivatives(tExtrap, DummySensor, DummyLos, InvCovUV);
            %[Rinit,Vinit,ParamList] = OrbitAtTime(this.KeplerToECI, tExtrap, DummySensor, DummyLos, InvCovUV);
            % OrbitalCov   = Jacobian*Cov*Jacobian';  % transformed Covariances
            %JacKepler_2_ECI      = JacKepler_2_ECI(1:6,2:7)
            JacKepler_2_ECI = this.KeplerToECI.JacKepler_2_ECI;
            %JacKepler_2_ECI
            this.JacKepler_2_ECI  = JacKepler_2_ECI;
            
            %JacKepler_2_ECI*JacECI_2_Kepler
            
            
            %JacKepler_2_Equinoctial
            JacKepler_2_Equinoctial      = this.KeplerToECI.JacKepler_2_Equinoctial;
            %JacKepler_2_Equinoctial
            this.JacKepler_2_Equinoctial = JacKepler_2_Equinoctial;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Check that the Inverse Transformation should produce the inverse Jacobian matrix.
            %The following matrix produce should be the 6x6 Identity matrix:
            [Okay] = Derivatives(this.ExtrapolatedClass, tExtrap);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Direct Jacobians from explicit differentiation of the
            % transformation equations.
            % JacECI_2_Kepler = J(Kepler/ECI)
            this.JacECI_2_Kepler         = this.ExtrapolatedClass.JacECI_2_Kepler;
            this.JacECI_2_Equinoctial    = this.ExtrapolatedClass.JacECI_2_Equinoctial;
            %this.EquinoctialElements
            %this.ClassicalElements
            %this.StateVector7
            %[XEquin, XKepler, XECI, JacECI_2_Equinoctial] = ECI_2_Equinoctial(this.StateVector7, this.ClassicalElements, this.EquinoctialElements,  mu)
            %JacECI_2_Equinoctial
            %this.JacECI_2_Equinoctial
            
            this.JacKepler_2_ECI         = this.ExtrapolatedClass.JacKepler_2_ECI;
            this.JacKepler_2_Equinoctial = this.ExtrapolatedClass.JacKepler_2_Equinoctial;
            %this.JacKepler_2_Equinoctial = this.JacECI_2_Equinoctial*inv(this.JacECI_2_Kepler);
            JacECI_2_Equinoctial = this.JacKepler_2_Equinoctial*inv(this.JacKepler_2_ECI);
            
            % Indirect Jacobians -- obtained by matrix multiplication and inversion
            this.JacEquinoctial_2_ECI    = inv(this.ExtrapolatedClass.JacECI_2_Equinoctial);
            this.JacEquinoctial_2_Kepler = inv(this.JacKepler_2_Equinoctial);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Debug Section:
            % [XEquin, JacKepler_2_Equinoctial] =  Kepler_2_EquinoctialCanon(this.ClassicalElements);
            % [this.EquinoctialElements, XEquin]
            % JacKepler_2_Equinoctial - this.JacKepler_2_Equinoctial
            %
            % XEquin = this.EquinoctialElements;
            % [XECI, JacEquinoctial_2_ECI] =  Equinoctial_2_ECICanon(XEquin);
            % [this.StateVector7, XECI]
            % JacEquinoctial_2_ECI - this.JacEquinoctial_2_ECI
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % JacECI_2_Kepler         = J(Kepler/ECI)           %direct method
            % % JacECI_2_Equinoctial    = J(Equinoctial/ECI)      %direct method
            % % JacKepler_2_ECI         = J(ECI/Kepler)           %direct method
            % % JacKepler_2_Equinoctial    = J(Equinoctial/Kepler)
            %
            % % JacEquinoctial_2_ECI    = J(ECI/Equinoctial);
            % % JacEquinoctial_2_Kepler = J(Kepler/Equinoctial);
            %
            % Jac1 = this.JacECI_2_Kepler
            % Jac2 = this.JacKepler_2_ECI
            % Jac3 = this.JacECI_2_Equinoctial
            %
            % % Indirect Jacobians -- obtained by matrix multiplication and inversion
            % Jac4 = this.JacEquinoctial_2_ECI
            % Jac5 = this.JacKepler_2_Equinoctial
            % Jac6 = this.JacEquinoctial_2_Kepler
            %
            % rank1 = rank(Jac1);
            % rank2 = rank(Jac2);
            % rank3 = rank(Jac3);
            % rank4 = rank(Jac4);
            % rank5 = rank(Jac5);
            % rank6 = rank(Jac6);
            %
            % disp([' Rank of Jacobian Matrices = ',num2str(rank1),'  ',...
            % num2str(rank2),'  ',num2str(rank3),'  ',num2str(rank4),'  ',...
            % num2str(rank5),'  ',num2str(rank6)])
            %
            % Jac1*Jac2
            % Jac3*Jac4
            % Jac5*Jac6
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Okay = 1;
        end
        
        function [Okay] = Derivatives(this, tExtrap)
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Checking First Derivatives
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % format short
            % iStart = 1;
            % Nlines = 115;
            % CombinedDifferences = 0.0;
            % for LineCheck = iStart:Nlines
            %     ForDerivs  = [];
            %     %BackDerivs = [];
            %     %Derivs     = [];
            %     %Combined   = [];
            %
            %     F            = zeros(Nlines,1);
            %     F(LineCheck) = 1;
            %     F            = backdiff(this, F, LineCheck);
            %     % BackDerivs   = [F(1) F(2) F(3) F(4) F(5) F(6) F(7)];
            %     % Combined     = [Index F];
            %     % Kate's Combined'
            %
            %     forward      = [];
            %     for jk = 1:Nlines % check range of derivatives wrt forwards differentiation
            %         %jk
            %         Q        = zeros(Nlines,1);
            %         Q(jk)     = 1;
            %         Q        = fordiff(this, Q, LineCheck);
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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% TEST SECOND DERIVATIVES:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %Debug Second Derivatives by checking symmetry:
            % %Begin Checking Out Second Derivatives:
            % Nvar     = 7;
            % % LineCheck    = 56;  % Kepler Part
            % % LineCheck    = 90;  % Escobal Gravitation to J6
            % % F            = zeros(115,1);
            % % F(LineCheck) = 1;
            % % F = backdiff(this, F, LineCheck);
            % HessianDifferences = 0.0;
            % Nlines = 115;
            % format short
            % for LineCheck = 1:Nlines
            %     ListMax   = LineCheck + 1;
            %     %Nvar = LineCheck;
            %     Hessian      = [];
            %     F            = zeros(Nlines,1);
            %     F(LineCheck) = 1;
            %     F = backdiff(this, F, LineCheck);
            %     %F(ListMax:Nlines) = 0;
            %     % for j = 1:Nvar
            %     for j = 1:LineCheck
            %         Q = zeros(Nlines,1);
            %         S = zeros(Nlines,1);
            %         Q(j) = 1.0;
            %         Q = fordiff(this, Q, LineCheck);
            %         %Q(ListMax:Nlines) = 0;
            %         S = secdiff(this, F,Q,S, LineCheck);
            %         %S(ListMax:Nlines) = 0;
            %         S = backdiff(this, S, LineCheck);
            %         %S(ListMax:Nlines) = 0;
            %         row = [];
            %         %for k = 1:Nvar
            %         for k = 1:LineCheck
            %             row = [row S(k)];
            %         end
            %         Hessian  = [Hessian; row];
            %     end
            %     %disp([' Line = ', num2str(Nvar)]);
            %     maxElement = max(max(abs(Hessian)));
            %     HessianTest = Hessian;
            %     if maxElement > 0
            %         HessianTest = Hessian/maxElement;
            %     end
            %
            %     Symmetry = HessianTest - HessianTest';
            %     Differences = max(max(abs(Symmetry)));
            %     Index = numel(find(abs(Symmetry) > 0.001));
            %     if Index > 0
            %         find (abs(Symmetry) > 0.001)
            %         for i = 1:LineCheck
            %             for j = i+1:LineCheck
            %                 if (abs(Symmetry(i,j)) > 0.001)
            %                     itrouble = i
            %                     jtrouble = j
            %                 end
            %             end
            %         end
            %         itrouble
            %         jtrouble
            %     end
            %
            %     HessianDifferences = HessianDifferences + Differences;
            %     disp([' Line Check = ', num2str(LineCheck),' max Hessian = ', num2str(maxElement),'  Summed Differences = ', num2str(Differences)]);
            %     %disp('Debug Point')
            % end
            % disp(['  Total Difference with Hessians = ',num2str(HessianDifferences)]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % JacobianAll = [];
            % for I = 1:100
            %     LineCheck    = I;                % e
            %     F            = zeros(115,1);
            %     F(LineCheck) = 1;
            %     F            = backdiff(this, F, LineCheck);
            %     JacobianAll  = [JacobianAll; F(2) F(3) F(4) F(5) F(6) F(7)];
            % end
            % this.JacobianAll = JacobianAll;
            % We will want to construct the 6 by 6 Jacobian of the Kepler
            % Orbital elements e, p, I, omega, Omega, PTime
            % or lines [17, 24, 64, 50, 47, 62]
            % and the parameters rx, ry, rz, vx, vy, vz
            %
            JacECI_2_Kepler = [];
            LineCheck    = 17;                % e
            F            = zeros(115,1);
            F(LineCheck) = 1;
            F            = backdiff(this, F, LineCheck);
            JacECI_2_Kepler     = [JacECI_2_Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 1
            
            LineCheck    = 43;                % a
            F            = zeros(115,1);
            F(LineCheck) = 1;
            F            = backdiff(this, F, LineCheck);
            JacECI_2_Kepler     = [JacECI_2_Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 2
            
            LineCheck    = 64;                % I
            F            = zeros(115,1);
            F(LineCheck) = 1;
            F            = backdiff(this, F, LineCheck);
            JacECI_2_Kepler     = [JacECI_2_Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 3
            
            LineCheck    = 50;                % omega
            F            = zeros(115,1);
            F(LineCheck) = 1;
            F            = backdiff(this, F, LineCheck);
            JacECI_2_Kepler     = [JacECI_2_Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 4
            
            LineCheck    = 47;                % Omega
            F            = zeros(115,1);
            F(LineCheck) = 1;
            F            = backdiff(this, F, LineCheck);
            JacECI_2_Kepler     = [JacECI_2_Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 5
            
            LineCheck    = 63;                % Mp in Radians
            F            = zeros(115,1);
            F(LineCheck) = 1;
            F            = backdiff(this, F, LineCheck);
            JacECI_2_Kepler     = [JacECI_2_Kepler; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 6
            
            this.JacECI_2_Kepler = JacECI_2_Kepler;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            JacECI_2_Equinoctial  = [];
            LineCheck    = 24;                % p
            F            = zeros(115,1);
            F(LineCheck) = 1;
            F            = backdiff(this, F, LineCheck);
            JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 2
            
            LineCheck    = 86;                % f
            F            = zeros(115,1);
            F(LineCheck) = 1;
            F            = backdiff(this, F, LineCheck);
            JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 2
            
            LineCheck    = 87;                % g
            F            = zeros(115,1);
            F(LineCheck) = 1;
            F            = backdiff(this, F, LineCheck);
            JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 3
            
            LineCheck    = 88;                % h
            F            = zeros(115,1);
            F(LineCheck) = 1;
            F            = backdiff(this, F, LineCheck);
            JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 4
            
            LineCheck    = 89;                % k
            F            = zeros(115,1);
            F(LineCheck) = 1;
            F            = backdiff(this, F, LineCheck);
            JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 5
            
            LineCheck    = 90;                % MLin Radians
            F            = zeros(115,1);
            F(LineCheck) = 1;
            F            = backdiff(this, F, LineCheck);
            JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(7)];   % row 6
            
            this.JacECI_2_Equinoctial    = JacECI_2_Equinoctial;
            
            % JacECI_2_Kepler = [];
            % LineCheck    = 17;                % e
            % [H_I, Jacob, Hess] = WorkOrder(this,LineCheck);
            % JacECI_2_Kepler     = [JacECI_2_Kepler; Jacob(2:7)];   % row 1
            %
            % LineCheck    = 43;                % a
            % [H_I, Jacob, Hess] = WorkOrder(this,LineCheck);
            % JacECI_2_Kepler     = [JacECI_2_Kepler; Jacob(2:7)];   % row 1
            %
            % LineCheck    = 64;                % I
            % [H_I, Jacob, Hess] = WorkOrder(this,LineCheck);
            % JacECI_2_Kepler     = [JacECI_2_Kepler; Jacob(2:7)];   % row 1
            %
            % LineCheck    = 50;                % omega
            % [H_I, Jacob, Hess] = WorkOrder(this,LineCheck);
            % JacECI_2_Kepler     = [JacECI_2_Kepler; Jacob(2:7)];   % row 1
            %
            % LineCheck    = 47;                % Omega
            % [H_I, Jacob, Hess] = WorkOrder(this,LineCheck);
            % JacECI_2_Kepler     = [JacECI_2_Kepler; Jacob(2:7)];   % row 1
            %
            % LineCheck    = 63;                % Mp in Radians
            % [H_I, Jacob, Hess] = WorkOrder(this,LineCheck);
            % JacECI_2_Kepler     = [JacECI_2_Kepler; Jacob(2:7)];   % row 1
            %
            % this.JacECI_2_Kepler = JacECI_2_Kepler;
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % JacECI_2_Equinoctial  = [];
            %
            % LineCheck    = 24;                % p
            % [H_I, Jacob, Hess] = WorkOrder(this,LineCheck);
            % JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; Jacob(2:7)];   % row 2
            %
            % LineCheck    = 86;                % f
            % [H_I, Jacob, Hess] = WorkOrder(this,LineCheck);
            % JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; Jacob(2:7)];   % row 2
            %
            % LineCheck    = 87;                % g
            % [H_I, Jacob, Hess] = WorkOrder(this,LineCheck);
            % JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; Jacob(2:7)];   % row 2
            %
            % LineCheck    = 88;                % h
            % [H_I, Jacob, Hess] = WorkOrder(this,LineCheck);
            % JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; Jacob(2:7)];   % row 2
            %
            % LineCheck    = 89;                % k
            % [H_I, Jacob, Hess] = WorkOrder(this,LineCheck);
            % JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; Jacob(2:7)];   % row 2
            %
            % LineCheck    = 90;                % Lin Radians
            % [H_I, Jacob, Hess] = WorkOrder(this,LineCheck);
            % JacECI_2_Equinoctial     = [JacECI_2_Equinoctial; Jacob(2:7)];   % row 2
            
            %JacKepler_2_Equinoctial      = JacECI_2_Equinoctial*inv(JacECI_2_Kepler);
            %JacKepler_2_Equinoctial
            %this.JacKepler_2_Equinoctial = JacKepler_2_Equinoctial;
            
            %this.JacECI_2_Kepler*this.JacKepler_2_ECI(1:6,2:7)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            epsilon     = 0.000001*ones(6,1);
            DummySensor = [this.StateVector7(1:3) + epsilon(1:3); zeros(6,1)];
            DummyLos    = [this.StateVector7(1:3) - epsilon(1:3)];
            InvCovUV    = eye(2);
            [Rinit,Vinit,ParamList,JacKepler_2_ECI,Hessian,GravityCan] = this.KeplerToECI.OrbitDerivatives(tExtrap, DummySensor, DummyLos, InvCovUV);
            %[Rinit,Vinit,ParamList] = OrbitAtTime(this.KeplerToECI, tExtrap, DummySensor, DummyLos, InvCovUV);
            % OrbitalCov   = Jacobian*Cov*Jacobian';  % transformed Covariances
            %JacKepler_2_ECI      = JacKepler_2_ECI(1:6,2:7)
            JacKepler_2_ECI = this.KeplerToECI.JacKepler_2_ECI;
            JacKepler_2_ECI;
            this.JacKepler_2_ECI = JacKepler_2_ECI;
            
            %JacKepler_2_ECI*JacECI_2_Kepler
            
            
            %JacKepler_2_Equinoctial
            JacKepler_2_Equinoctial      = this.KeplerToECI.JacKepler_2_Equinoctial;
            JacKepler_2_Equinoctial;
            this.JacKepler_2_Equinoctial = JacKepler_2_Equinoctial;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % % 'Debug Point'
            % tExtrap
            %Jac1 = this.JacKepler_2_ECI;
            %Jac2 = this.KeplerToECI.JacKepler_2_ECI;
            %Jac3 = this.JacECI_2_Kepler;
            %Jac1
            %Jac2
            %Jac3
            % Jac1*Jac2
            Okay = 1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        function F = backdiff(this, f, LineCheck)
            % eliminate some repetitions:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Canonical Coordinates assumed throughout -- length in DU, time in TU
            % Modifications can be made to the choice of units in wgs84Constants!
            % Orbital elements e, p, I, omega, Omega, PTime
            % or lines [17, 24, 64, 50, 47, 62]
            % and the parameters rx, ry, rz, vx, vy, vz
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Constants
            %units = wgs84Constants;
            twopi       = this.twopi;
            % TU          = units.TU;  % Canonical Time UnitTime Unit in Seconds
            % DU          = units.DU;  % Canonical Distance Unit
            % VU          = units.VU;  % Canonical Velocity Unit
            % AU          = units.AU;        % Canonical Acceleration Unit
            mu          = this.mu;  % Canonical Gravitational Constant
            sqrtmu      = this.sqrtmu;
            hmag        = this.hmag;
            SignOrbit   = this.SignOrbit; 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % TU          = 1.0;
            % DU          = 1.0;
            % VU          = 1.0;
            % AU          = 1.0;
            % this.TU     = TU;
            % this.DU     = DU;
            % this.VU     = VU;
            % this.AU     = AU;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Eprime      = this.Eprime;
            eDenom      = this.eDenom;
            dE_dM       = this.dE_dM;
            dE_de       = this.dE_de;
            d2E_dMdM    = this.d2E_dMdM;
            d2E_dMde    = this.d2E_dMde;
            d2E_dedM    = this.d2E_dedM;
            d2E_dede    = this.d2E_dede;
            
            sinO    = this.sinO;
            cosO    = this.cosO;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tEpoch      = this.DataList( 1); % = t;               %  1) t
            rx          = this.DataList( 2); % = rx;              %  2) Rx
            ry          = this.DataList( 3); % = ry;              %  3) Ry
            rz          = this.DataList( 4); % = rz;              %  4) Rz
            vx          = this.DataList( 5); % = vx;              %  5) Vx
            vy          = this.DataList( 6); % = vy;              %  6) Vy
            vz          = this.DataList( 7); % = vz;              %  7) Vz
            tExtrap     = this.DataList( 8); % = Extrapolation time.
            
            rmag        = this.DataList(10); % = rmag;            % 10) rmag
            vsq         = this.DataList(11); % = vsq;             % 11) vsq
            er          = this.DataList(12); % = er;              % 12) er
            ev          = this.DataList(13); % = ev;              % 13) ev
            ex          = this.DataList(14); % = ex;              % 14) Evec(1)
            ey          = this.DataList(15); % = ey;              % 15) Evec(2)
            ez          = this.DataList(16); % = ez;              % 16) Evec(3)
            e           = this.DataList(17); % = e;               % 17) e;
            Px          = this.DataList(18); % = Px;              % 18) Pvec(1)
            Py          = this.DataList(19); % = Py;              % 19) Pvec(2)
            Pz          = this.DataList(20); % = Pz;              % 20) Pvec(3)
            hx          = this.DataList(21); % = hx;              % 21) hVec(1)
            hy          = this.DataList(22); % = hy;              % 22) hVec(2)
            hz          = this.DataList(23); % = hz;              % 23) hVec(3)
            p           = this.DataList(24); % = p;;              % 24) p = hVec^2/mu
            Wx          = this.DataList(25); % = Wx;              % 25) Wvec(1)
            Wy          = this.DataList(26); % = Wy;              % 26) Wvec(2)
            Wz          = this.DataList(27); % = Wz;              % 27) Wvec(3)
            Qx          = this.DataList(28); % = Qx;              % 28) Qvec(1)
            Qy          = this.DataList(29); % = Qy;              % 29) Qvec(2)
            Qz          = this.DataList(30); % = Qz;              % 30) Qvec(3)
            rUnitx      = this.DataList(31); % = rUnitx;          % 31) rUnit(1)
            rUnity      = this.DataList(32); % = rUnity;          % 32) rUnit(2)
            rUnitz      = this.DataList(33); % = rUnitz;          % 33) rUnit(3)
            Ax          = this.DataList(34); % = Ax;              % 34) Ax
            Ay          = this.DataList(35); % = Ay;              % 35) Ay
            normN       = this.DataList(36); % = normN;           % 36) normN
            Nx          = this.DataList(37); % = Nx;              % 37) Nvec(1)
            Ny          = this.DataList(38); % = Ny;              % 38) Nvec(2)
            onePe       = this.DataList(39); % = onePe;           % 39)
            oneMe       = this.DataList(40); % = oneMe;           % 40)
            fac         = this.DataList(41); % = fac;             % 41)
            rootfac     = this.DataList(42); % = rootfac;         % 42)
            a           = this.DataList(43); % = a;               % 43)
            PeriodTU    = this.DataList(44); % = PeriodTU;        % 44) Period in TU's
            Period      = this.DataList(45); % = Period;          % 45) Work in TU's for now
            meanMotion  = this.DataList(46); % = meanMotion;      % 46) Radians per TU.
            Omega       = this.DataList(47); % = Omega;           % 47) Omega = atan(Ny/Nx);
            cosP        = this.DataList(48); % = cosP;            % 48) cos(omega)
            sinP        = this.DataList(49); % = sinP;            % 49) sin(omega)
            omega       = this.DataList(50); % = omega;           % 50) omega = atan(sinP/cosP)
            cosnu       = this.DataList(51); % = cosnu;           % 51) cos(nu)
            sinnu       = this.DataList(52); % = sinnu;           % 52) sin(nu)
            nuEpoch     = this.DataList(53); % = nuEpoch;         % 53) nu = atan(sinnu/cosnu)
            cosT        = this.DataList(54); % = cosT;            % 54) ~cos(E)
            sinT        = this.DataList(55); % = sinT;            % 55) ~sin(E)
            EccentricAnomalyEpoch = this.DataList(56); %= EccentricAnomalyEpoch;  % 56) E = atan(sinE/cosE)
            cosE        = this.DataList(57); % = cosE;                   % 57  (old 54)
            sinE        = this.DataList(58); % = sinE;                   % 58  (old 55)
            %MeanAnomalyEpoch   = this.DataList(59); % = MeanAnomalyEpoch;       % 59)
            %TimeSincePeriapsis = this.DataList(60); % = TimeSincePeriapsis;     % 60)
            %DeltaTime          = this.DataList(61); % = DeltaTime;              % 61)
            PTime       = this.DataList(62); % = PTime;  % DEFAULT% 62) in TUs
            Mp          = this.DataList(63); % = Mp;              % 63 in radians
            Inclination = this.DataList(64); % = Inclination;     % 64)
            M           = this.DataList(65); % = M;               % 65)
            EM          = this.DataList(66); % = EM;              % 66) EM = E(e,M)
            cosK        = this.DataList(67); % = cosK;            % 67)
            sinK        = this.DataList(68); % = sinK;            % 68)
            rtpinv      = this.DataList(69); %;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tanX        = this.DataList(70); %= tanX;             % 70)
            tanY        = this.DataList(71); %= tanY;             % 71)
            nu          = this.DataList(72); %= nu;               % 72)
            coss        = this.DataList(73); %= coss;             % 73)
            sins        = this.DataList(74); %= sins;             % 74)
            rorbit      = this.DataList(75); %= rorbit;           % 75)
            rorbitx     = this.DataList(76); %= rorbitx;          % 76)
            rorbity     = this.DataList(77); %= rorbity;          % 77)
            rorbitz     = this.DataList(78); %= rorbitz;          % 78)
            vorbitx     = this.DataList(79); %= vorbitx;          % 79)
            vorbity     = this.DataList(80); %= vorbity;          % 80)
            vorbitz     = this.DataList(81); %= vorbitz;          % 81)
            Aeq         = this.DataList(82); %= Aeq;              % line 82
            cosAeq      = this.DataList(83); %= cosAeq;           % line 83
            sinAeq      = this.DataList(84); %= sinAeq;           % line 84
            tanHalf     = this.DataList(85); %= tanHalf;          % line 85
            Feq         = this.DataList(86); %= Feq;              % line 86
            Geq         = this.DataList(87); %= Geq;              % line 87
            Heq         = this.DataList(88); %                    % line 88
            Keq         = this.DataList(89); %                    % line 89
            Leq         = this.DataList(90); %                    % line 90
            CosL        = this.DataList(91); %                    % line 91 CosL
            SinL        = this.DataList(92); %                    % line 92 SinL
            alphaSq     = this.DataList(93); %                    % line 93 alphaSq
            Seq         = this.DataList(94); %                    % line 94
            Weq         = this.DataList(95); %                    % line 95
            Req         = this.DataList(96); %                    % line 96
            RovS        = this.DataList(97); %                    % line 97
            srtpinv     = this.DataList(98); %                    % line 98
            HK          = this.DataList(99); %                    % line 99
            OnePalphaSqCosL  = this.DataList(100);                % line 100
            OneMalphaSqCosL  = this.DataList(101);                % line 101
            OnePalphaSqSinL  = this.DataList(102);                % line 102
            OneMalphaSqSinL  = this.DataList(103);                % line 103
            Xfac             = this.DataList(104);                % line 104
            Yfac             = this.DataList(105);                % line 105
            Zfac             = this.DataList(106);                % line 106
            VXfac            = this.DataList(107);                % line 107
            VYfac            = this.DataList(108);                % line 108
            VZfac            = this.DataList(109);                % line 109
            Xeq              = this.DataList(110);                % line 110
            Yeq              = this.DataList(111);                % line 111
            Zeq              = this.DataList(112);                % line 112
            VXeq             = this.DataList(113);                % line 113
            VYeq             = this.DataList(114);                % line 114
            VZeq             = this.DataList(115);                % line 115
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %//////////////////////////////////////////////////////////////////////
            if LineCheck > 64
                if LineCheck > 81
                    % from line 115: VZeq    = 2.0*srtpinv*VZfac; % line 115
                    f( 98) = f( 98) + 2.0*f(115)*VZfac;
                    f(109) = f(109) + 2.0*f(115)*srtpinv;
                    % from line 114: VYeq    =    -srtpinv*VYfac; % line 114
                    f( 98) = f( 98) - f(114)*VYfac;
                    f(108) = f(108) - f(114)*srtpinv;
                    % from line 113: VXeq    =    -srtpinv*VXfac; % line 113
                    f( 98) = f( 98) - f(113)*VXfac;
                    f(107) = f(107) - f(113)*srtpinv;
                    %from line 112: Zeq      = 2.0*RovS*Zfac;     % line 112
                    f( 97) = f( 97) + 2.0*f(112)*Zfac;
                    f(106) = f(106) + 2.0*f(112)*RovS;
                    % from line 111: Yeq     =     RovS*Yfac;     % line 111
                    f( 97) = f( 97) + f(111)*Yfac;
                    f(105) = f(105) + f(111)*RovS;
                    % from line 110: Xeq     =     RovS*Xfac;     % line 110
                    f( 97) = f( 97) + f(110)*Xfac;
                    f(104) = f(104) + f(110)*RovS;
                    % from line: VZfac   =  Heq*(Feq + CosL) + Keq*(Geq + SinL) % line 109
                    f( 86) = f( 86) + f(109)*Heq;               % Feq
                    f( 87) = f( 87) + f(109)*Keq;               % Geq
                    f( 88) = f( 88) + f(109)*(Feq+CosL);        % Heq
                    f( 89) = f( 89) + f(109)*(Geq+SinL);        % Keq
                    f( 91) = f( 91) + f(109)*Heq;               % CosL
                    f( 92) = f( 92) + f(109)*Keq;               % SinL
                    
                    % from line 108: VYfac = -OneMalphaSqCosL + 2.0*HK*(Geq + SinL) + Feq*(alphaSq - 1.0);% line 108
                    f( 86) = f( 86) +     f(108)*(alphaSq - 1.0); % Feq
                    f( 87) = f( 87) + 2.0*f(108)*HK;              % Geq
                    f( 92) = f( 92) + 2.0*f(108)*HK;              % SinL
                    f( 93) = f( 93) +     f(108)*Feq;             % alphaSq
                    f( 99) = f( 99) + 2.0*f(108)*(Geq+SinL);      % HK
                    f(101) = f(101) -     f(108);                 % OneMAlphaSqSinL
                    
                    % from line 107: VXfac =  OnePalphaSqSinL - 2.0*HK*(Feq + CosL) + Geq*(1.0 + alphaSq);% line 107
                    f( 86) = f( 86) - 2.0*f(107)*HK;              % Feq
                    f( 87) = f( 87) +     f(107)*(1 + alphaSq);   % Geq
                    f( 91) = f( 91) - 2.0*f(107)*HK;              % CosL
                    f( 93) = f( 93) +     f(107)*Geq;             % alphaSq
                    f( 99) = f( 99) - 2.0*f(107)*(Feq+CosL);      % HK
                    f(102) = f(102) +     f(107);                 % OnePAlphaSqSinL
                    
                    % from line 106: Zfac    = Heq*SinL - Keq*CosL;          % line 106
                    f(88) = f(88) + f(106)*SinL;
                    f(89) = f(89) - f(106)*CosL;
                    f(91) = f(91) - f(106)*Keq;
                    f(92) = f(92) + f(106)*Heq;
                    % from line 105: Yfac    = OneMalphaSqSinL + 2.0*HK*CosL; % line 105
                    f(91) = f(91) + 2.0*f(105)*HK;
                    f(99) = f(99) + 2.0*f(105)*CosL;
                    f(103) = f(103) + f(105);
                    % from line 104: Xfac    = OnePalphaSqCosL + 2.0*HK*SinL; % line 104
                    f(92) = f(92) + 2.0*f(104)*HK;
                    f(99) = f(99) + 2.0*f(104)*SinL;
                    f(100) = f(100) + f(104);
                    % from line 103: OneMalphaSqSinL = (1.-alphaSq)*SinL;% line 103
                    f(92) = f(92) + f(103)*(1.0 - alphaSq);
                    f(93) = f(93) - f(103)*SinL;
                    % from line 102: OnePalphaSqSinL = (1.+alphaSq)*SinL;% line 102
                    f(92) = f(92) + f(102)*(1.0 + alphaSq);
                    f(93) = f(93) + f(102)*SinL;
                    % from line 101: OneMalphaSqCosL = (1.-alphaSq)*CosL;% line 101
                    f(91) = f(91) + f(101)*(1.0 - alphaSq);
                    f(93) = f(93) - f(101)*CosL;
                    % from line 100: OnePalphaSqCosL = (1.+alphaSq)*CosL;% line 100
                    f(91) = f(91) + f(100)*(1.0 + alphaSq);
                    f(93) = f(93) + f(100)*CosL;
                    % from line 99: HK      = Heq*Keq;            % line 99
                    f(88) = f(88) + f(99)*Keq;
                    f(89) = f(89) + f(99)*Heq;
                    % from line 98: srtpinv = rtpinv/Seq;         % line 98
                    f(69) = f(69) + f(98)/Seq;
                    f(94) = f(94) - f(98)*srtpinv/Seq;
                    % from line 97: RovS    = Req/Seq;            % line 97
                    f(96) = f(96) + f(97)/Seq;
                    f(94) = f(94) - f(97)*RovS/Seq;
                    % from line 96: Req     = p/Weq;              % line 96
                    f(24) = f(24) + f(96)/Weq;
                    f(95) = f(95) - f(96)*Req/Weq;
                    % from line 95: Weq = 1.0 + Feq*CosL + Geq*SinL; % line 95
                    f(86) = f(86) + f(95)*CosL;
                    f(87) = f(87) + f(95)*SinL;
                    f(91) = f(91) + f(95)*Feq;
                    f(92) = f(92) + f(95)*Geq;
                    % from line 94: Seq     = 1.0 + Heq^2 + Keq^2;   % line 94
                    f(88) = f(88) + 2.0*f(94)*Heq;
                    f(89) = f(89) + 2.0*f(94)*Keq;
                    % from line 93: alphaSq = Heq^2 - Keq^2;         % line 93
                    f(88) = f(88) + 2.0*f(93)*Heq;
                    f(89) = f(89) - 2.0*f(93)*Keq;
                    % from line 92: SinL     = sin(Leq);          % line 92
                    f(90) = f(90) + f(92)*CosL;
                    % from line 91: CosL     = cos(Leq);          % line 91
                    f(90) = f(90) - f(91)*SinL;
                    % from line 90: Leq      = Aeq + nu;          % line 90
                    f(72) = f(72) + f(90);
                    f(82) = f(82) + f(90);
                    % from line 89: Keq      = tanHalf*sin(Omega);    % line 89
                    f(85) = f(85) + f(89)*sinO;
                    f(47) = f(47) + f(89)*tanHalf*cosO;
                    % from line 88: Heq      = tanHalf*cos(Omega);    % line 88
                    f(85) = f(85) + f(88)*cosO;
                    f(47) = f(47) - f(88)*tanHalf*sinO;
                    
                    % from line 87: Geq      = e*sinAeq;              % line 87
                    f(17) = f(17) + f(87)*sinAeq;
                    f(84) = f(84) + f(87)*e;
                    %f(82) = f(82) + f(87)*e*cosAeq;
                    % from line 86: Feq      = e*cosAeq;              % line 86
                    f(17) = f(17) + f(86)*cosAeq;
                    f(83) = f(83) + f(86)*e;
                    %f(82) = f(82) - f(86)*e*sinAeq;
                    
                    % from line 85: tanHalf  = tan(Inclination/2.0);  % line 85
                    f(64) = f(64) + 0.5*f(85)*(1.0 + tanHalf^2);
                    
                    % from line 84: sinAeq   = sin(Aeq);              % line 84
                    f(82) = f(82) + f(84)*cosAeq;
                    % from line 83: cosAeq   = cos(Aeq);              % line 83
                    f(82) = f(82) - f(83)*sinAeq;
                    % from line 82: Aeq      = omega + Omega;         % line 82
                    f(47) = f(47) + f(82);
                    f(50) = f(50) + f(82);
                end
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
                % from line 69: rtpinv   = sqrt(mu/p);           % line 69
                f(24) = f(24) - 0.5*f(69)*rtpinv/p;
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
                %//////////////////////////////////////////////////////////////////////
                % % from line 75: double rorbit   = p/(1.0 + e*coss);
                % denom = 1.0/(1.0 + e*coss);
                % f(73) = f(73) - f(75)*p*e*denom*denom;     % coss  #73
                % f(24) = f(24) + f(75)*denom;               % p     #24
                % f(17) = f(17) - f(75)*p*coss*denom*denom;  % e     #17
                %//////////////////////////////////////////////////////////////////////
                % % from line 75: rorbit   = a*(1.0 - e*cosK);    % line 75
                denom = 1.0 - e*cosK;
                f(67) = f(67) - f(75)*a*e;                      % cosK   42
                %f(66) = f(66) + f(75)*a*e*sinK;                % cos(E) 41
                f(43) = f(43) + f(75)*denom;                    %  a      3
                f(17) = f(17) - f(75)*a*cosK;                   %  e      2
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
                % % from line 65: M = meanMotion*t + MeanAnomalyEpoch;
                %f(59) = f(59) + f(65);                                % MeanAnomalyEpoch #59
                % %f(46) = f(46) + f(65)*t;                              % meanMotion       #46
                % %f( 1) = f( 1) + f(65)*meanMotion;                     % t       # 1
                % from line 65: M = meanMotion*tExtrap + Mp;
                f(63) = f(63) + f(65);              % Mp                #63
                f(46) = f(46) + f(65)*tExtrap;      % meanMotion        #46
                f( 8) = f( 8) + f(65)*meanMotion;   %tExtrap            # 8
            end
            %/////////////////////Block 64-54//////////////////////////////////////
            % from line 64: Inclination = acos(Wz);
            f(27) = f(27) - f(64)/sqrt(1.0 - Wz*Wz);              % Wz      #27
            % from line 63: Tnext     = PTime*time_unit;
            % f(62) = f(62) + f(63)*time_unit;                      % PTime   # 62
            % from line 63: Mp = MeanAnomalyEpoch - meanMotion*tEpoch;
            f(59) = f(59) + f(63);                   % MeanAnomalyEpoch #59
            f(46) = f(46) - f(63)*tEpoch;            % meanMotion       #46
            f( 1) = f( 1) - f(63)*meanMotion;        % tEpoch           # 1
            %//////////////////////////////////////////////////////////////////////
            % from line 62 : if DeltaTime > 0 PTime = DeltaTime;
            %                if DeltaTime < 0 PTime = DeltaTime + Period
            %  PTime = DeltaTime  % DEFAULT                       % 62)
            %  if DeltaTime < 0
            %    Then OVERWRITE PTime HERE:
            %    PTime = DeltaTime + Period;
            %  end
            % f(61) = f(61) + f(62);                                % DeltaTime  #61
            % if DeltaTime < 0
            %     f(45) = f(45) + f(62);                            % Period  #45
            % end
            % % from line 61: DeltaTime = tEpoch - TimeSincePeriapsis;
            % f(60) = f(60) - f(61);                                % TimeSincePeriapsis #60
            % f( 1) = f( 1) + f(61);                                % tEpoch        # 1
            % % from line 60: TimeSincePeriapsis = MeanAnomalyEpoch/meanMotion;
            % f(59) = f(59) + f(60)/meanMotion;                     % MeanAnomalyEpoch #59
            % f(46) = f(46) - f(60)*TimeSincePeriapsis/meanMotion;  % meanMotion #46
            % from line 59: MeanAnomalyEpoch   = EccentricAnomalyEpoch - e*sinE;
            %//////////////////////////////////////////////////////////////////////
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
            % from line: 44 PeriodTU = twopi*(SignOrbit*a)^(1.5)/sqrtmu;        % Period in TU's
            f(43) = f(43) + f(44)*PeriodTU*1.5*SignOrbit/a;                 % a       #43
            % Semi-Major Axis:
            % from line 43: a = p/fac;
            f(41) = f(41) - f(43)*a/fac;                          % fac     #41
            f(24) = f(24) + f(43)/fac;                            % p       #24
            %if (fac > 0.0)
                % from line 42: rootfac = sqrt(SignOrbit*fac);
                f(41) = f(41) + f(42)*(0.5*SignOrbit)/rootfac;               % fac     #41
            %else
            %    rootfac = 0;
            %end
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
            % Test Section
            % f( 1) = f( 1) + f( 8);  % for tExtrap == tEpoch
            %         vz   = Rdot(3);                                   % Vz      # 7
            %         vy   = Rdot(2);                                   % Vy      # 6
            %         vx   = Rdot(1);                                   % Vx      # 5
            %         rz   = Rpos(3);                                   % Rz      # 4
            %         ry   = Rpos(2);                                   % Ry      # 3
            %         rx   = Rpos(1);                                   % Rx      # 2
            %         t    = time;                                      % t       # 1
            F = f;
        end
        
        function Q = fordiff(this, q, LineCheck)
            %     eliminate some repetitions:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Canonical Coordinates assumed throughout -- length in DU, time in TU
            % Modifications can be made to the choice of units in wgs84Constants!
            % Orbital elements e, p, I, omega, Omega, PTime
            % or lines [17, 24, 64, 50, 47, 62]
            % and the parameters rx, ry, rz, vx, vy, vz
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Constants
            %units = wgs84Constants;
            twopi       = this.twopi;
            % TU          = units.TU;  % Canonical Time UnitTime Unit in Seconds
            % DU          = units.DU;  % Canonical Distance Unit
            % VU          = units.VU;  % Canonical Velocity Unit
            % AU          = units.AU;        % Canonical Acceleration Unit
            mu          = this.mu;  % Canonical Gravitational Constant
            sqrtmu      = this.sqrtmu;
            hmag        = this.hmag;
            SignOrbit   = this.SignOrbit; 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % TU          = 1.0;
            % DU          = 1.0;
            % VU          = 1.0;
            % AU          = 1.0;
            % this.TU     = TU;
            % this.DU     = DU;
            % this.VU     = VU;
            % this.AU     = AU;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Eprime      = this.Eprime;
            eDenom      = this.eDenom;
            dE_dM       = this.dE_dM;
            dE_de       = this.dE_de;
            d2E_dMdM    = this.d2E_dMdM;
            d2E_dMde    = this.d2E_dMde;
            d2E_dedM    = this.d2E_dedM;
            d2E_dede    = this.d2E_dede;
            
            sinO    = this.sinO;
            cosO    = this.cosO;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tEpoch      = this.DataList( 1); % = t;               %  1) t
            rx          = this.DataList( 2); % = rx;              %  2) Rx
            ry          = this.DataList( 3); % = ry;              %  3) Ry
            rz          = this.DataList( 4); % = rz;              %  4) Rz
            vx          = this.DataList( 5); % = vx;              %  5) Vx
            vy          = this.DataList( 6); % = vy;              %  6) Vy
            vz          = this.DataList( 7); % = vz;              %  7) Vz
            tExtrap     = this.DataList( 8); % = Extrapolation time.
            
            rmag        = this.DataList(10); % = rmag;            % 10) rmag
            vsq         = this.DataList(11); % = vsq;             % 11) vsq
            er          = this.DataList(12); % = er;              % 12) er
            ev          = this.DataList(13); % = ev;              % 13) ev
            ex          = this.DataList(14); % = ex;              % 14) Evec(1)
            ey          = this.DataList(15); % = ey;              % 15) Evec(2)
            ez          = this.DataList(16); % = ez;              % 16) Evec(3)
            e           = this.DataList(17); % = e;               % 17) e;
            Px          = this.DataList(18); % = Px;              % 18) Pvec(1)
            Py          = this.DataList(19); % = Py;              % 19) Pvec(2)
            Pz          = this.DataList(20); % = Pz;              % 20) Pvec(3)
            hx          = this.DataList(21); % = hx;              % 21) hVec(1)
            hy          = this.DataList(22); % = hy;              % 22) hVec(2)
            hz          = this.DataList(23); % = hz;              % 23) hVec(3)
            p           = this.DataList(24); % = p;;              % 24) p = hVec^2/mu
            Wx          = this.DataList(25); % = Wx;              % 25) Wvec(1)
            Wy          = this.DataList(26); % = Wy;              % 26) Wvec(2)
            Wz          = this.DataList(27); % = Wz;              % 27) Wvec(3)
            Qx          = this.DataList(28); % = Qx;              % 28) Qvec(1)
            Qy          = this.DataList(29); % = Qy;              % 29) Qvec(2)
            Qz          = this.DataList(30); % = Qz;              % 30) Qvec(3)
            rUnitx      = this.DataList(31); % = rUnitx;          % 31) rUnit(1)
            rUnity      = this.DataList(32); % = rUnity;          % 32) rUnit(2)
            rUnitz      = this.DataList(33); % = rUnitz;          % 33) rUnit(3)
            Ax          = this.DataList(34); % = Ax;              % 34) Ax
            Ay          = this.DataList(35); % = Ay;              % 35) Ay
            normN       = this.DataList(36); % = normN;           % 36) normN
            Nx          = this.DataList(37); % = Nx;              % 37) Nvec(1)
            Ny          = this.DataList(38); % = Ny;              % 38) Nvec(2)
            onePe       = this.DataList(39); % = onePe;           % 39)
            oneMe       = this.DataList(40); % = oneMe;           % 40)
            fac         = this.DataList(41); % = fac;             % 41)
            rootfac     = this.DataList(42); % = rootfac;         % 42)
            a           = this.DataList(43); % = a;               % 43)
            PeriodTU    = this.DataList(44); % = PeriodTU;        % 44) Period in TU's
            Period      = this.DataList(45); % = Period;          % 45) Work in TU's for now
            meanMotion  = this.DataList(46); % = meanMotion;      % 46) Radians per TU.
            Omega       = this.DataList(47); % = Omega;           % 47) Omega = atan(Ny/Nx);
            cosP        = this.DataList(48); % = cosP;            % 48) cos(omega)
            sinP        = this.DataList(49); % = sinP;            % 49) sin(omega)
            omega       = this.DataList(50); % = omega;           % 50) omega = atan(sinP/cosP)
            cosnu       = this.DataList(51); % = cosnu;           % 51) cos(nu)
            sinnu       = this.DataList(52); % = sinnu;           % 52) sin(nu)
            nuEpoch     = this.DataList(53); % = nuEpoch;         % 53) nu = atan(sinnu/cosnu)
            cosT        = this.DataList(54); % = cosT;            % 54) ~cos(E)
            sinT        = this.DataList(55); % = sinT;            % 55) ~sin(E)
            EccentricAnomalyEpoch = this.DataList(56); %= EccentricAnomalyEpoch;  % 56) E = atan(sinE/cosE)
            cosE        = this.DataList(57); % = cosE;                   % 57  (old 54)
            sinE        = this.DataList(58); % = sinE;                   % 58  (old 55)
            %MeanAnomalyEpoch   = this.DataList(59); % = MeanAnomalyEpoch;       % 59)
            %TimeSincePeriapsis = this.DataList(60); % = TimeSincePeriapsis;     % 60)
            %DeltaTime          = this.DataList(61); % = DeltaTime;              % 61)
            PTime       = this.DataList(62); % = PTime;  % DEFAULT% 62) in TUs
            Mp          = this.DataList(63); % = Mp;              % 63 in radians
            Inclination = this.DataList(64); % = Inclination;     % 64)
            M           = this.DataList(65); % = M;               % 65)
            EM          = this.DataList(66); % = EM;              % 66) EM = E(e,M)
            cosK        = this.DataList(67); % = cosK;            % 67)
            sinK        = this.DataList(68); % = sinK;            % 68)
            rtpinv      = this.DataList(69); %;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tanX        = this.DataList(70); %= tanX;             % 70)
            tanY        = this.DataList(71); %= tanY;             % 71)
            nu          = this.DataList(72); %= nu;               % 72)
            coss        = this.DataList(73); %= coss;             % 73)
            sins        = this.DataList(74); %= sins;             % 74)
            rorbit      = this.DataList(75); %= rorbit;           % 75)
            rorbitx     = this.DataList(76); %= rorbitx;          % 76)
            rorbity     = this.DataList(77); %= rorbity;          % 77)
            rorbitz     = this.DataList(78); %= rorbitz;          % 78)
            vorbitx     = this.DataList(79); %= vorbitx;          % 79)
            vorbity     = this.DataList(80); %= vorbity;          % 80)
            vorbitz     = this.DataList(81); %= vorbitz;          % 81)
            Aeq         = this.DataList(82); %= Aeq;              % line 82
            cosAeq      = this.DataList(83); %= cosAeq;           % line 83
            sinAeq      = this.DataList(84); %= sinAeq;           % line 84
            tanHalf     = this.DataList(85); %= tanHalf;          % line 85
            Feq         = this.DataList(86); %= Feq;              % line 86
            Geq         = this.DataList(87); %= Geq;              % line 87
            Heq         = this.DataList(88); %                    % line 88
            Keq         = this.DataList(89); %                    % line 89
            Leq         = this.DataList(90); %                    % line 90
            CosL        = this.DataList(91); %                    % line 91 CosL
            SinL        = this.DataList(92); %                    % line 92 SinL
            alphaSq     = this.DataList(93); %                    % line 93 alphaSq
            Seq         = this.DataList(94); %                    % line 94
            Weq         = this.DataList(95); %                    % line 95
            Req         = this.DataList(96); %                    % line 96
            RovS        = this.DataList(97); %                    % line 97
            srtpinv     = this.DataList(98); %                    % line 98
            HK          = this.DataList(99); %                    % line 99
            OnePalphaSqCosL  = this.DataList(100);                % line 100
            OneMalphaSqCosL  = this.DataList(101);                % line 101
            OnePalphaSqSinL  = this.DataList(102);                % line 102
            OneMalphaSqSinL  = this.DataList(103);                % line 103
            Xfac             = this.DataList(104);                % line 104
            Yfac             = this.DataList(105);                % line 105
            Zfac             = this.DataList(106);                % line 106
            VXfac            = this.DataList(107);                % line 107
            VYfac            = this.DataList(108);                % line 108
            VZfac            = this.DataList(109);                % line 109
            Xeq              = this.DataList(110);                % line 110
            Yeq              = this.DataList(111);                % line 111
            Zeq              = this.DataList(112);                % line 112
            VXeq             = this.DataList(113);                % line 113
            VYeq             = this.DataList(114);                % line 114
            VZeq             = this.DataList(115);                % line 115
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %     rtpinv   = rtpinv;
            %     Parameters 1-7:
            %     t    = time;                        %  1) t
            %     rx   = Rpos(1);                     %  2) Rx
            %     ry   = Rpos(2);                     %  3) Ry
            %     rz   = Rpos(3);                     %  4) Rz
            %     vx   = Rdot(1);                     %  5) Vx
            %     vy   = Rdot(2);                     %  6) Vy
            %     vz   = Rdot(3);                     %  7) Vz
            % Test Section
            % q( 8) = q( 8) + q( 1);
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
            % 	// from line 42: if(fac > 0) rootfac = sqrt(SignOrbit*fac);
            %if (fac > 0.0)
                q(42) = q(42) + q(41)*0.5*SignOrbit/rootfac;
            %end
            % 	// from line 43: a = p/fac;
            q(43) = q(43) + q(24)/fac - q(41)*a/fac;
            %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %     % Period
            % 	// from line 44: PeriodTU = twopi*(SignOrbit*a)^(1.5)/sqrt(mu);
            q(44) = q(44) + q(43)*1.5*SignOrbit*PeriodTU/a;
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
            %//////////////////////////////////////////////////////////////
            % % // from line 60: TimeSincePeriapsis = MeanAnomalyEpoch/meanMotion;
            % q(60) = q(60) + q(59)/meanMotion - q(46)*TimeSincePeriapsis/meanMotion;
            % % PTime                  = TimeSincePeriapsis;
            % % // from line 61: DeltaTime = tEpoch - TimeSincePeriapsis;
            % q(61) = q(61) + q( 1) - q(60);
            % % // from line 62: PTime:
            % if DeltaTime > 0
            %     q(62) = q(62) + q(61);
            % else
            %     q(62) = q(62) + q(61) + q(45);
            % end
            %//////////////////////////////////////////////////////////////
            % from line 63: Tnext     = PTime*time_unit;
            %  q(63) = q(63) + q(62)*time_unit;                     % PTime   # 62
            % from line 63: Mp = MeanAnomalyEpoch - meanMotion*tEpoch;
            q(63) = q(63) + q(59);                   % MeanAnomalyEpoch #59
            q(63) = q(63) - q(46)*tEpoch;            %      meanMotion  #46
            q(63) = q(63) - q( 1)*meanMotion;        %       tEpoch     # 1
            % Kepler 6: Inclination
            % // from line 64: Inclination = acos(Wz);
            q(64) = q(64) - q(27)/sqrt(1.0 - Wz*Wz);
            if LineCheck > 64
                % Comment this out for the moment!
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % // from line 65: M = MeanAnomalyEpoch + meanMotion*tExtrap;
                % %q(65) = q(65) + q(59) + q(46)*tExtrap + q( 8)*meanMotion;
                %q(65) = q(65) +  q(59);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % // from line 65: M = Mp + meanMotion*tExtrap;
                q(65) = q(65) + q(63) + q(46)*tExtrap + q( 8)*meanMotion;
                %q(65) = q(65) +  q(63);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                %denom = 1.0/(1.0 + e*coss);
                %      q(75) = q(75) + q(72)*p*e*sins*denom*denom;        % line 72 # nu
                %q(75) = q(75) - q(73)*p*e*denom*denom;              % line 73 # coss
                %q(75) = q(75) + q(24)*denom;
                %q(75) = q(75) - q(17)*p*coss*denom*denom;
                % // from line 49: rorbit   = a*(1.0 - e*cosK);
                q(75) = q(75) + q(43)*(1.0 - e*cosK) - q(67)*a*e - q(17)*a*cosK;
                %q(49) = q(49) + q( 3)*(1.0 - e*cosK) - q(42)*a*e - q(2)*a*cosK;
                % %q(49) = q(49) + q( 3)*(1.0 - e*cosK) + q(41)*a*e*sinK - q(2)*a*cosK;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                % from line 69: rtpinv = sqrt(mu/p);
                q(69) = q(69) - 0.5*q(24)*rtpinv/p;
                % // from line 79: vorbitx  = rtpinv*(-sins*Px + (e + coss)*Qx);
                %      q(79) = q(79) - q(72)*rtpinv*(coss*Px + sins*Qx);%line 72 # nu
                q(79) = q(79) - q(74)*rtpinv*Px;                  %line 74 # sins
                q(79) = q(79) + q(73)*rtpinv*Qx;                  %line 73 # coss
                q(79) = q(79) + q(28)*rtpinv*(e + coss);
                q(79) = q(79) - q(24)*(0.5/(p*sqrt(p)))*(-sins*Px + (e + coss)*Qx);
                q(79) = q(79) - q(18)*rtpinv*sins;
                q(79) = q(79) + q(17)*rtpinv*Qx;
                % // from line 80: vorbity  = rtpinv*(-sins*Py + (e + coss)*Qy);
                %      q(80) = q(80) - q(72)*rtpinv*(coss*Py + sins*Qy);%line 72 # nu
                q(80) = q(80) - q(74)*rtpinv*Py;                  %line 74 # sins
                q(80) = q(80) + q(73)*rtpinv*Qy;                  %line 73 # coss
                q(80) = q(80) + q(29)*rtpinv*(e + coss);
                q(80) = q(80) - q(24)*(0.5/(p*sqrt(p)))*(-sins*Py + (e + coss)*Qy);
                q(80) = q(80) - q(19)*rtpinv*sins;
                q(80) = q(80) + q(17)*rtpinv*Qy;
                % // from line 81: vorbitz  = rtpinv*(-sins*Pz + (e + coss)*Qz);
                %      q(81) = q(81) - q(72)*rtpinv*(coss*Pz + sins*Qz);%line 72 # nu
                q(81) = q(81) - q(74)*rtpinv*Pz;                  %line 74 # sins
                q(81) = q(81) + q(73)*rtpinv*Qz;                  %line 73 # coss
                q(81) = q(81) + q(30)*rtpinv*(e + coss);
                q(81) = q(81) - q(24)*(0.5/(p*sqrt(p)))*(-sins*Pz + (e + coss)*Qz);
                q(81) = q(81) - q(20)*rtpinv*sins;
                q(81) = q(81) + q(17)*rtpinv*Qz;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if LineCheck > 81
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
                    q(88) = q(88) + q(85)*cosO - q(47)*tanHalf*sinO;
                    % from line 89: Keq      = tanHalf*sin(Omega);    % line 89
                    q(89) = q(89) + q(85)*sinO + q(47)*tanHalf*cosO;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % from line 90: Leq      = Aeq + nu;         % line 90
                    q(90) = q(90) + q(72) + q(82);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % from line 91: CosL     = cos(Leq);          % line 91
                    q(91) = q(91) - q(90)*SinL;
                    % from line 92: SinL     = sin(Leq);          % line 92
                    q(92) = q(92) + q(90)*CosL;
                    % from line 93: alphaSq = Heq^2 - Keq^2;         % line 93
                    q(93) = q(93) + 2.0*(q(88)*Heq - q(89)*Keq);
                    % from line 94: Seq     = 1.0 + Heq^2 + Keq^2;   % line 94
                    q(94) = q(94) + 2.0*(q(88)*Heq + q(89)*Keq);
                    % from line 95: Weq = 1.0 + Feq*CosL + Geq*SinL; % line 95
                    q(95) = q(95) + q(86)*CosL + q(87)*SinL + q(91)*Feq + q(92)*Geq;
                    % from line 96: Req     = p/Weq;              % line 96
                    q(96) = q(96) + q(24)/Weq - q(95)*Req/Weq;
                    % from line 97: RovS    = Req/Seq;            % line 97
                    q(97) = q(97) + q(96)/Seq - q(94)*RovS/Seq;
                    % from line 98: srtpinv = rtpinv/Seq;         % line 98
                    q(98) = q(98) + q(69)/Seq - q(94)*srtpinv/Seq;
                    % from line 99: HK  = Heq*Keq;            % line 99
                    q(99) = q(99) + q(88)*Keq + q(89)*Heq;
                    % from line 100: OnePalphaSqCosL = (1.+alphaSq)*CosL;% line 100
                    q(100) = q(100) + q(91)*(1.0 + alphaSq) + q(93)*CosL;
                    % from line 101: OneMalphaSqCosL = (1.-alphaSq)*CosL;% line 101
                    q(101) = q(101) + q(91)*(1.0 - alphaSq) - q(93)*CosL;
                    % from line 102: OnePalphaSqSinL = (1.+alphaSq)*SinL;% line 102
                    q(102) = q(102) + q(92)*(1.0 + alphaSq) + q(93)*SinL;
                    % from line 103: OneMalphaSqSinL = (1.-alphaSq)*SinL;% line 103
                    q(103) = q(103) + q(92)*(1.0 - alphaSq) - q(93)*SinL;
                    % from line 104: Xfac    = OnePalphaSqCosL + 2.0*HK*SinL; % line 104
                    q(104) = q(104) + 2.0*q(92)*HK + 2.0*q(99)*SinL + q(100);
                    % from line 105: Yfac    = OneMalphaSqSinL + 2.0*HK*CosL; % line 105
                    q(105) = q(105) + 2.0*q(91)*HK + 2.0*q(99)*CosL + q(103);
                    
                    % from line 106: Zfac    = Heq*SinL - Keq*CosL;          % line 106
                    q(106) = q(106) + q(88)*SinL - q(89)*CosL  - q(91)*Keq  + q(92)*Heq;
                    
                    % from line 107: VXfac =  OnePalphaSqSinL - 2.0*HK*(Feq+CosL) + Geq*(1+alphaSq);% line 107
                    q(107) = q(107) -2.0*HK*(q(86)+q(91)) +q(87)*(1.0+alphaSq) +q(93)*Geq -2.0*q(99)*(Feq+CosL) +q(102);
                    
                    % from line 108: VYfac = -OneMalphaSqCosL + 2.0*HK*(Geq+SinL) + Feq*(alphaSq-1);% line 108
                    q(108) = q(108) + q(86)*(alphaSq-1.0) +2.0*HK*(q(87)+q(92))+q(93)*Feq +2.0*q(99)*(Geq+SinL) -q(101);
                    
                    % from line: VZfac   =  Heq*(Feq + CosL) + Keq*(Geq + SinL) % line 109
                    q(109) = q(109) +q(86)*Heq +q(87)*Keq +q(88)*(Feq+CosL) +q(89)*(Geq+SinL) +q(91)*Heq +q(92)*Keq;
                    
                    % from line 110: Xeq     =     RovS*Xfac;     % line 110
                    q(110) = q(110) + q(97)*Xfac + q(104)*RovS;
                    
                    % from line 111: Yeq     =     RovS*Yfac;     % line 111
                    q(111) = q(111) + q(97)*Yfac + q(105)*RovS;
                    
                    %from line 112: Zeq      = 2.0*RovS*Zfac;     % line 112
                    q(112) = q(112) + 2.0*(q(97)*Zfac + q(106)*RovS);
                    
                    % from line 113: VXeq    =    -srtpinv*VXfac; % line 113
                    q(113) = q(113) - q(98)*VXfac - q(107)*srtpinv;
                    
                    % from line 114: VYeq    =    -srtpinv*VYfac; % line 114
                    q(114) = q(114) - q(98)*VYfac - q(108)*srtpinv;
                    
                    % from line 115: VZeq    = 2.0*srtpinv*VZfac; % line 115
                    q(115) = q(115) + 2.0*q(98)*VZfac + 2.0*q(109)*srtpinv;
                end
            end
            
            Q = q;
        end
        
        function S = secdiff(this, f, q, s, LineCheck, varargin)
            %function S = secdiff(f, q, s, LineCheck)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Canonical Coordinates assumed throughout -- length in DU, time in TU
            % Modifications can be made to the choice of units in wgs84Constants!
            % Orbital elements e, p, I, omega, Omega, PTime
            % or lines [17, 24, 64, 50, 47, 62]
            % and the parameters rx, ry, rz, vx, vy, vz
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Constants
            %units = wgs84Constants;
            twopi       = this.twopi;
            % TU          = units.TU;  % Canonical Time UnitTime Unit in Seconds
            % DU          = units.DU;  % Canonical Distance Unit
            % VU          = units.VU;  % Canonical Velocity Unit
            % AU          = units.AU;        % Canonical Acceleration Unit
            mu          = this.mu;  % Canonical Gravitational Constant
            sqrtmu      = this.sqrtmu;
            hmag        = this.hmag;
            SignOrbit   = this.SignOrbit; 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % TU          = 1.0;
            % DU          = 1.0;
            % VU          = 1.0;
            % AU          = 1.0;
            % this.TU     = TU;
            % this.DU     = DU;
            % this.VU     = VU;
            % this.AU     = AU;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Eprime      = this.Eprime;
            eDenom      = this.eDenom;
            dE_dM       = this.dE_dM;
            dE_de       = this.dE_de;
            d2E_dMdM    = this.d2E_dMdM;
            d2E_dMde    = this.d2E_dMde;
            d2E_dedM    = this.d2E_dedM;
            d2E_dede    = this.d2E_dede;
            
            sinO    = this.sinO;
            cosO    = this.cosO;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tEpoch      = this.DataList( 1); % = t;               %  1) t
            rx          = this.DataList( 2); % = rx;              %  2) Rx
            ry          = this.DataList( 3); % = ry;              %  3) Ry
            rz          = this.DataList( 4); % = rz;              %  4) Rz
            vx          = this.DataList( 5); % = vx;              %  5) Vx
            vy          = this.DataList( 6); % = vy;              %  6) Vy
            vz          = this.DataList( 7); % = vz;              %  7) Vz
            tExtrap     = this.DataList( 8); % = Extrapolation time.
            
            rmag        = this.DataList(10); % = rmag;            % 10) rmag
            vsq         = this.DataList(11); % = vsq;             % 11) vsq
            er          = this.DataList(12); % = er;              % 12) er
            ev          = this.DataList(13); % = ev;              % 13) ev
            ex          = this.DataList(14); % = ex;              % 14) Evec(1)
            ey          = this.DataList(15); % = ey;              % 15) Evec(2)
            ez          = this.DataList(16); % = ez;              % 16) Evec(3)
            e           = this.DataList(17); % = e;               % 17) e;
            Px          = this.DataList(18); % = Px;              % 18) Pvec(1)
            Py          = this.DataList(19); % = Py;              % 19) Pvec(2)
            Pz          = this.DataList(20); % = Pz;              % 20) Pvec(3)
            hx          = this.DataList(21); % = hx;              % 21) hVec(1)
            hy          = this.DataList(22); % = hy;              % 22) hVec(2)
            hz          = this.DataList(23); % = hz;              % 23) hVec(3)
            p           = this.DataList(24); % = p;;              % 24) p = hVec^2/mu
            Wx          = this.DataList(25); % = Wx;              % 25) Wvec(1)
            Wy          = this.DataList(26); % = Wy;              % 26) Wvec(2)
            Wz          = this.DataList(27); % = Wz;              % 27) Wvec(3)
            Qx          = this.DataList(28); % = Qx;              % 28) Qvec(1)
            Qy          = this.DataList(29); % = Qy;              % 29) Qvec(2)
            Qz          = this.DataList(30); % = Qz;              % 30) Qvec(3)
            rUnitx      = this.DataList(31); % = rUnitx;          % 31) rUnit(1)
            rUnity      = this.DataList(32); % = rUnity;          % 32) rUnit(2)
            rUnitz      = this.DataList(33); % = rUnitz;          % 33) rUnit(3)
            Ax          = this.DataList(34); % = Ax;              % 34) Ax
            Ay          = this.DataList(35); % = Ay;              % 35) Ay
            normN       = this.DataList(36); % = normN;           % 36) normN
            Nx          = this.DataList(37); % = Nx;              % 37) Nvec(1)
            Ny          = this.DataList(38); % = Ny;              % 38) Nvec(2)
            onePe       = this.DataList(39); % = onePe;           % 39)
            oneMe       = this.DataList(40); % = oneMe;           % 40)
            fac         = this.DataList(41); % = fac;             % 41)
            rootfac     = this.DataList(42); % = rootfac;         % 42)
            a           = this.DataList(43); % = a;               % 43)
            PeriodTU    = this.DataList(44); % = PeriodTU;        % 44) Period in TU's
            Period      = this.DataList(45); % = Period;          % 45) Work in TU's for now
            meanMotion  = this.DataList(46); % = meanMotion;      % 46) Radians per TU.
            Omega       = this.DataList(47); % = Omega;           % 47) Omega = atan(Ny/Nx);
            cosP        = this.DataList(48); % = cosP;            % 48) cos(omega)
            sinP        = this.DataList(49); % = sinP;            % 49) sin(omega)
            omega       = this.DataList(50); % = omega;           % 50) omega = atan(sinP/cosP)
            cosnu       = this.DataList(51); % = cosnu;           % 51) cos(nu)
            sinnu       = this.DataList(52); % = sinnu;           % 52) sin(nu)
            nuEpoch     = this.DataList(53); % = nuEpoch;         % 53) nu = atan(sinnu/cosnu)
            cosT        = this.DataList(54); % = cosT;            % 54) ~cos(E)
            sinT        = this.DataList(55); % = sinT;            % 55) ~sin(E)
            EccentricAnomalyEpoch = this.DataList(56); %= EccentricAnomalyEpoch;  % 56) E = atan(sinE/cosE)
            cosE        = this.DataList(57); % = cosE;                   % 57  (old 54)
            sinE        = this.DataList(58); % = sinE;                   % 58  (old 55)
            MeanAnomalyEpoch   = this.DataList(59); % = MeanAnomalyEpoch;       % 59)
            %TimeSincePeriapsis = this.DataList(60); % = TimeSincePeriapsis;     % 60)
            %DeltaTime          = this.DataList(61); % = DeltaTime;              % 61)
            %PTime       = this.DataList(62); % = PTime;  % DEFAULT% 62) in TUs
            Mp          = this.DataList(63); % = Mp;              % 63 in radians
            Inclination = this.DataList(64); % = Inclination;     % 64)
            M           = this.DataList(65); % = M;               % 65)
            EM          = this.DataList(66); % = EM;              % 66) EM = E(e,M)
            cosK        = this.DataList(67); % = cosK;            % 67)
            sinK        = this.DataList(68); % = sinK;            % 68)
            rtpinv      = this.DataList(69); %;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tanX        = this.DataList(70); %= tanX;             % 70)
            tanY        = this.DataList(71); %= tanY;             % 71)
            nu          = this.DataList(72); %= nu;               % 72)
            coss        = this.DataList(73); %= coss;             % 73)
            sins        = this.DataList(74); %= sins;             % 74)
            rorbit      = this.DataList(75); %= rorbit;           % 75)
            rorbitx     = this.DataList(76); %= rorbitx;          % 76)
            rorbity     = this.DataList(77); %= rorbity;          % 77)
            rorbitz     = this.DataList(78); %= rorbitz;          % 78)
            vorbitx     = this.DataList(79); %= vorbitx;          % 79)
            vorbity     = this.DataList(80); %= vorbity;          % 80)
            vorbitz     = this.DataList(81); %= vorbitz;          % 81)
            Aeq         = this.DataList(82); %= Aeq;              % line 82
            cosAeq      = this.DataList(83); %= cosAeq;           % line 83
            sinAeq      = this.DataList(84); %= sinAeq;           % line 84
            tanHalf     = this.DataList(85); %= tanHalf;          % line 85
            Feq         = this.DataList(86); %= Feq;              % line 86
            Geq         = this.DataList(87); %= Geq;              % line 87
            Heq         = this.DataList(88); %                    % line 88
            Keq         = this.DataList(89); %                    % line 89
            Leq         = this.DataList(90); %                    % line 90
            CosL        = this.DataList(91); %                    % line 91 CosL
            SinL        = this.DataList(92); %                    % line 92 SinL
            alphaSq     = this.DataList(93); %                    % line 93 alphaSq
            Seq         = this.DataList(94); %                    % line 94
            Weq         = this.DataList(95); %                    % line 95
            Req         = this.DataList(96); %                    % line 96
            RovS        = this.DataList(97); %                    % line 97
            srtpinv     = this.DataList(98); %                    % line 98
            HK          = this.DataList(99); %                    % line 99
            OnePalphaSqCosL  = this.DataList(100);                % line 100
            OneMalphaSqCosL  = this.DataList(101);                % line 101
            OnePalphaSqSinL  = this.DataList(102);                % line 102
            OneMalphaSqSinL  = this.DataList(103);                % line 103
            Xfac             = this.DataList(104);                % line 104
            Yfac             = this.DataList(105);                % line 105
            Zfac             = this.DataList(106);                % line 106
            VXfac            = this.DataList(107);                % line 107
            VYfac            = this.DataList(108);                % line 108
            VZfac            = this.DataList(109);                % line 109
            Xeq              = this.DataList(110);                % line 110
            Yeq              = this.DataList(111);                % line 111
            Zeq              = this.DataList(112);                % line 112
            VXeq             = this.DataList(113);                % line 113
            VYeq             = this.DataList(114);                % line 114
            VZeq             = this.DataList(115);                % line 115
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  mu    = 1.0;
            % if nargin >= 5
            %     TU = this.TU;
            %     DU = this.DU;
            %     VU = this.VU;
            %     AU = this.AU;
            % end
            %disp('Hello World');
            denom    = rmag*rmag*rmag;
            edenom   = e*e*e;
            cesqr    = (coss*e+1.0)*(coss*e+1.0);
            cecube   = cesqr*(coss*e+1.0);
            psq      = p*p;
            pcubed   = p*p*p;
            pdenom   = sqrt(pcubed);
            cosDenom = 0.0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % t    = time;                                   % line   1 t
            % rx   = Rpos(1);                                % line   2 Rx
            % ry   = Rpos(2);                                % line   3 Ry
            % rz   = Rpos(3);                                % line   4 Rz
            % vx   = Rdot(1);                                % line   5 Vx
            % vy   = Rdot(2);                                % line   6 Vy
            % vz   = Rdot(3);                                % line   7 Vz
            
            % from line 10: rmag = sqrt(rx*rx + ry*ry + rz*rz);
            %// rx J = 2, KJ = rx/rmag I = 2,3,4
            s( 2) = s( 2) + f(10)*q( 2)*(1.0 - (rx/rmag)^2)/rmag;
            s( 2) = s( 2) - f(10)*q( 3)*(rx*ry)/rmag^3;
            s( 2) = s( 2) - f(10)*q( 4)*(rx*rz)/rmag^3;
            %// ry J = 3, KJ = ry/rmag I = 2,3,4
            s( 3) = s( 3) - f(10)*q( 2)*(ry*rx)/rmag^3;
            s( 3) = s( 3) + f(10)*q( 3)*(1.0 - (ry/rmag)^2)/rmag;
            s( 3) = s( 3) - f(10)*q( 4)*(ry*rz)/rmag^3;
            %// rz J = 4, KJ = rz/rmag I = 2,3,4
            s( 4) = s( 4) - f(10)*q( 2)*(rz*rx)/rmag^3;
            s( 4) = s( 4) - f(10)*q( 3)*(rz*ry)/rmag^3;
            s( 4) = s( 4) + f(10)*q( 4)*(1.0 - (rz/rmag)^2)/rmag;
            % from line 11: vsq  = vx*vx + vy*vy + vz*vz;
            %// vx J = 5, KJ = 2*vx I = 5
            s( 5) = s( 5) + 2.0*f(11)*q( 5);
            %// vy J = 6, KJ = 2*vy I = 6
            s( 6) = s( 6) + 2.0*f(11)*q( 6);
            %// vy J = 7, KJ = 2*vz I = 7
            s( 7) = s( 7) + 2.0*f(11)*q( 7);
            % from line 12: er   = vsq - mu/rmag;
            %// vsq J = 11, KJ = 1 all second derivatives vanish
            %// rmag J = 10: KJ = mu*rmag^(-2), I = 10
            s(10) = s(10) -2.0*f(12)*q(10)*mu/rmag^3;
            % from line 13: ev   = rx*vx + ry*vy + rz*vz;
            %// rx J = 2: KJ = vx I = 5
            s( 2) = s( 2) + f(13)*q( 5);
            %// ry J = 3: KJ = vy I = 6
            s( 3) = s( 3) + f(13)*q( 6);
            %// rz J = 4: KJ = vz I = 7
            s( 4) = s( 4) + f(13)*q( 7);
            %// vx J = 5: KJ = rx I = 2
            s( 5) = s( 5) + f(13)*q( 2);
            %// vy J = 6: KJ = ry I = 3
            s( 6) = s( 6) + f(13)*q( 3);
            %// vz J = 7: KJ = rz I = 4
            s( 7) = s( 7) + f(13)*q( 4);
            % Evec = er*Rpos - ev*Vpos;
            % from line 14: ex   = (er*rx - ev*vx)/mu;
            %// rx J = 2: KJ = er/mu I = 12
            s( 2) = s( 2) + f(14)*q(12)/mu;
            %// er J=12 : KJ = rx/mu I = 2
            s(12) = s(12) + f(14)*q( 2)/mu;
            %// vx J = 5: KJ = -ev/mu I = 13
            s( 5) = s( 5) - f(14)*q(13)/mu;
            %// ev J=13 : KJ = -vx/mu I = 5
            s(13) = s(13) - f(14)*q( 5)/mu;
            % from line 15: ey   = (er*ry - ev*vy)/mu;
            %// ry J = 3: KJ = er/mu I = 12
            s( 3) = s( 3) + f(15)*q(12)/mu;
            %// er J=12 : KJ = ry/mu I = 3
            s(12) = s(12) + f(15)*q( 3)/mu;
            %// vy J = 6: KJ = -ev/mu I = 13
            s( 6) = s( 6) - f(15)*q(13)/mu;
            %// ev J=13 : KJ = -vy/mu I = 6
            s(13) = s(13) - f(15)*q( 6)/mu;
            % from line 16: ez   = (er*rz - ev*vz)/mu;
            %// rz J = 4: KJ = er/mu I = 12
            s( 4) = s( 4) + f(16)*q(12)/mu;
            %// er J=12 : KJ = rz/mu I = 4
            s(12) = s(12) + f(16)*q( 4)/mu;
            %// vy J = 7: KJ = -ev/mu I = 13
            s( 7) = s( 7) - f(16)*q(13)/mu;
            %// ev J=13 : KJ = -vz/mu I = 7
            s(13) = s(13) - f(16)*q( 7)/mu;
            %Evec = [ex ey ez];
            % from line 17: e    = sqrt(ex*ex + ey*ey + ez*ez);
            %// ex J = 14, KJ = ex/e I = 14, 15, 16
            s(14) = s(14) + f(17)*( q(14)*(ey*ey + ez*ez) - q(15)*ex*ey           - q(16)*ex*ez          )/e^3;
            %// ey J = 15, KJ = ey/e I = 14,15,16
            s(15) = s(15) + f(17)*(-q(14)*(ey*ex)         + q(15)*(ex*ex + ez*ez) - q(16)*ey*ez          )/e^3;
            %// ez J = 16, KJ = ez/e I = 14,15,16
            s(16) = s(16) + f(17)*(-q(14)*ez*ex           - q(15)*ez*ey           + q(16)*(ex*ex + ey*ey))/e^3;
            % from line 18: Px   = ex/e;
            %// ex J = 14: KJ = 1/e I = 17
            s(14) = s(14) - f(18)*q(17)/e^2;
            %// e J = 17: KJ = -ex/e^2, I = 14, 17
            s(17) = s(17) + f(18)*(2.0*q(17)*Px - q(14))/e^2;
            % from line 19: Py   = ey/e;
            %// ey J = 15: KJ = 1/e  I = 17
            s(15) = s(15) - f(19)*q(17)/e^2;
            %// e J = 17: KJ = -ey/e^2, I = 15, 17
            s(17) = s(17) + f(19)*(2.0*q(17)*Py - q(15))/e^2;
            % from line 20: Pz   = ez/e;
            %// ez J = 16: KJ = 1/e  I = 17
            s(16) = s(16) - f(20)*q(17)/e^2;
            %// e J = 17: KJ = -ez/e^2, I = 16, 17
            s(17) = s(17) + f(20)*(2.0*q(17)*Pz - q(16))/e^2;
            % from line 21: hx   = ry*vz - rz*vy;            // line 21
            %// ry J = 3: KJ = vz I = 7
            s( 3) = s( 3) + f(21)*q( 7);
            s( 7) = s( 7) + f(21)*q( 3);
            s( 4) = s( 4) - f(21)*q( 6);
            s( 6) = s( 6) - f(21)*q( 4);
            % from line 22: hy   = rz*vx - rx*vz;            // line 22
            s( 4) = s( 4) + f(22)*q( 5);
            s( 5) = s( 5) + f(22)*q( 4);
            s( 2) = s( 2) - f(22)*q( 7);
            s( 7) = s( 7) - f(22)*q( 2);
            % from line 23: hz   = rx*vy - ry*vx;            // line 23
            s( 2) = s( 2) + f(23)*q( 6);
            s( 6) = s( 6) + f(23)*q( 2);
            s( 3) = s( 3) - f(23)*q( 5);
            s( 5) = s( 5) - f(23)*q( 3);
            % from line 24: p    = hx*hx + hy*hy + hz*hz;    // line 24
            s( 21) = s(21) + f(24)*q(21)*2.0;
            s( 22) = s(22) + f(24)*q(22)*2.0;
            s( 23) = s(23) + f(24)*q(23)*2.0;
            %
            % hmag = sqrt(mu*p) -- "p" scaled by mu
            % from line 25: Wx   = hx/hmag;                  // line 25
            %// hx J = 21: KJ = 1/sqrt(mu*p) I = p
            s(21) = s(21) - 0.5*f(25)*q(24)/(hmag*p);
            %// p J = 24: KJ = -0.5*hx/( sqrt(mu)*p^(3/2) ), I = 21, 24
            s(24) = s(24) + f(25)*(0.75*Wx*q(24)/p - 0.5*q(21)/hmag)/p;
            % from line 26: Wy   = hy/hmag;                  // line 26
            %// hy J = 22: KJ = 1/sqrt(mu*p) I = p
            s(22) = s(22) - f(26)*0.5*q(24)/(hmag*p);
            %// p J = 24: KJ = -0.5*hy/( sqrt(mu)*p^(3/2) ), I = 22, 24
            s(24) = s(24) + f(26)*(0.75*Wy*q(24)/p - 0.5*q(22)/hmag)/p;
            % from line 27: Wz   = hz/hmag;                  // line 27
            %// hz J = 23: KJ = 1/sqrt(mu*p) I = p
            s(23) = s(23) - f(27)*0.5*q(24)/(hmag*p);
            %// p J = 24: KJ = -0.5*hy/( sqrt(mu)*p^(3/2) ), I = 22, 24
            s(24) = s(24) + f(27)*(0.75*Wz*q(24)/p - 0.5*q(23)/hmag)/p;
            % from line 28: Qx   = Wy*Pz - Wz*Py;            // line 28
            s(20) = s(20) + f(28)*q(26);
            s(26) = s(26) + f(28)*q(20);
            s(19) = s(19) - f(28)*q(27);
            s(27) = s(27) - f(28)*q(19);
            % from line 29: Qy   = Wz*Px - Wx*Pz;            // line 29
            s(18) = s(18) + f(29)*q(27);
            s(27) = s(27) + f(29)*q(18);
            s(20) = s(20) - f(29)*q(25);
            s(25) = s(25) - f(29)*q(20);
            % from line 30: Qz   = Wx*Py - Wy*Px;            // line 30
            s(19) = s(19) + f(30)*q(25);
            s(25) = s(25) + f(30)*q(19);
            s(18) = s(18) - f(30)*q(26);
            s(26) = s(26) - f(30)*q(18);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % from line 31: rUnitx = rx/rmag;                // line 31
            %//rx J = 2, KJ = 1/rmag I = 10
            s( 2) = s( 2) - f(31)*q(10)/rmag^2;
            %// rmag J = 10: KJ = -rx/rmag^2, I = 2, 10
            s(10) = s(10) + f(31)*(2.0*q(10)*rUnitx - q( 2))/rmag^2;
            % from line 32: rUnity = ry/rmag;                // line 32
            %//ry J = 3, KJ = 1/rmag I = 10
            s( 3) = s( 3) - f(32)*q(10)/rmag^2;
            %// rmag J = 10: KJ = -ry/rmag^2, I = 3, 10
            s(10) = s(10) + f(32)*(2.0*q(10)*rUnity - q( 3))/rmag^2;
            % from line 33: rUnitz = rz/rmag;                // line 33
            %//rz J = 4, KJ = 1/rmag I = 10
            s( 4) = s( 4) - f(33)*q(10)/rmag^2;
            %// rmag J = 10: KJ = -rz/rmag^2, I = 4, 10
            s(10) = s(10) + f(33)*(2.0*q(10)*rUnitz - q( 4))/rmag^2;
            % from line 34: Ax     = -Wy;                    // line 34
            % All second derivatives vanish.
            % from line 35: Ay     =  Wx;                    // line 35
            % All second derivatives vanish
            % from line 36: normN  = sqrt(Ax*Ax + Ay*Ay);    // line 36
            %//Ax J = 34: KJ = Ax/sqrt(Ax^2 + Ay^2), I = 34, 35
            s(34) = s(34) + f(36)*(q(34)*Ay*Ay - q(35)*Ax*Ay)/normN^3;
            %//Ay J = 35: KJ = Ay/sqrt(Ax^2 + Ay^2), I = 34, 35
            s(35) = s(35) + f(36)*(q(35)*Ax*Ax - q(34)*Ax*Ay)/normN^3;
            % from line 37: Nx     = Ax/normN;               // line 37
            %// Ax J = 34: KJ = 1/normN I = 36
            s(34) = s(34) - f(37)*q(36)/normN^2;
            %// normN J = 36: KJ = -Ax/normN^2 I = 34, 36
            s(36) = s(36) + f(37)*(2.0*Nx*q(36) - q(34))/normN^2;
            % from line 38: Ny     = Ay/normN;               // line 38
            %// Ay J = 35: KJ = 1/normN I = 36
            s(35) = s(35) - f(38)*q(36)/normN^2;
            %// normN J = 36: KJ = -Ay/normN^2 I = 35, 36
            s(36) = s(36) + f(38)*(2.0*Ny*q(36) - q(35))/normN^2;
            % from line 39 onePe = 1.0 + e;                  // line 39
            % All second derivatives vanish
            % from line 40: oneMe = 1.0 - e;                 // line 40
            % All second derivatives vanish
            % from line 41: fac   = onePe*oneMe;             // line 41
            %//onePe J = 39: KJ = oneMe I = 40
            s(39) = s(39) + f(41)*q(40);
            %//oneMe J = 40: Kj = onePe I = 39
            s(40) = s(40) + f(41)*q(39);
            % from line 42: rootfac = sqrt(SignOrbit*fac);             // line 42
            %//fac J = 41: KJ = SignOrbit/sqrt(SignOrbit*fac) I = 41
            s(41) = s(41) - 0.25*f(42)*q(41)/(SignOrbit*fac*rootfac);
            % from line 43: a = p/fac;                       // line 43
            %//p J = 24: KJ = 1/fac
            s(24) = s(24) - f(43)*q(41)/fac^2;
            %//fac J = 41: KJ = -p/fac^2 I = 24, 41
            s(41) = s(41) + f(43)*(2.0*a*q(41) - q(24))/fac^2;
            % from line 44: PeriodTU   = twopi*(SignOrbit*a)^(1.5)/sqrtmu;// line 44
            %//a J = 43: KJ = 1.5*twopi*(SignOrbit*a)^(1/2)/sqrtmu
            s(43) = s(43) + 0.75*f(44)*q(43)*PeriodTU/a^2;
            % Period     = PeriodTU*time_unit;    % 45) Period in Seconds
            % from line 45: Period     = PeriodTU;           // line 45
            % All second derivatives vanish
            % from line 46: meanMotion = twopi/Period;       // line 46
            %//Period J = 45: KJ = -twopi/Period^2
            s(45) = s(45) + 2.0*f(46)*q(45)*meanMotion/Period^2;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Kepler 3) Omega
            % Longitude of the Ascending Node
            % from line 47: Omega = atan2(Ny, Nx);           // line 47
            cosDenom = (Nx*Nx + Ny*Ny)^2;
            % line K = 47 dependencies (  Ny 38,   Nx 37)
            % q(47) = q(47) - q(37)*Ny/(Nx*Nx + Ny*Ny) + q(38)*Nx/(Nx*Nx + Ny*Ny);
            % Nx J = 37 KJ = -Ny/(Nx*Nx + Ny*Ny);
            m00 = 2.0*Nx*Ny/cosDenom;
            m11 = (Ny*Ny - Nx*Nx)/cosDenom;
            s(37) = s(37) + f(47)*(q(37)*m00 + q(38)*m11);
            % Ny J = 38: KJ = Nx/(Nx*Nx + Ny*Ny);
            s(38) = s(38) + f(47)*(q(37)*m11 - q(38)*m00);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Kepler 4) omega
            % Argument of Periapsis
            %cosPeriapsis = dot(Nvec, Pvec);
            % from line: 48: cosP = Nx*Px + Ny*Py;           // line 48
            %//Px J = 18: KJ = Nx, I = 37
            s(18) = s(18) + f(48)*q(37);
            s(19) = s(19) + f(48)*q(38);
            s(37) = s(37) + f(48)*q(18);
            s(38) = s(38) + f(48)*q(19);
            %sinPeriapsis = dot(Wvec, cross(Nvec, Pvec));
            % from line 49: sinP =(Wx*Ny - Wy*Nx)*Pz + (Nx*Py - Ny*Px)*Wz; // line 49
            %//Nx J = 37: KJ = Py*Wz - Pz*Wy I = 19, 20, 26, 27
            s(37) = s(37) + f(49)*(q(19)*Wz + q(27)*Py - q(20)*Wy - q(26)*Pz);
            %//Ny J = 38: KJ = Pz*Wx - Px*Wz I = 18, 20, 25, 27
            s(38) = s(38) + f(49)*(q(20)*Wx + q(25)*Pz - q(18)*Wz - q(27)*Px);
            %//Wx J = 25: KJ = Pz*Ny, I = 20, 38
            s(25) = s(25) + f(49)*(q(20)*Ny + q(38)*Pz);
            %//Wy J = 26: KJ =-Pz*Nx, I = 20, 37
            s(26) = s(26) - f(49)*(q(20)*Nx + q(37)*Pz);
            %//Wz J = 27: KJ = (Py*Nx - Px*Ny);
            s(27) = s(27) + f(49)*(q(19)*Nx + q(37)*Py - q(18)*Ny - q(38)*Px);
            %//Px J = 18: KJ = -Ny*Wz
            s(18) = s(18) - f(49)*(q(27)*Ny +q(38)*Wz);
            %//Py J = 19: KJ = Nx*Wz;
            s(19) = s(19) + f(49)*(q(27)*Nx + q(37)*Wz);
            %//Pz J = 20: KJ = (Wx*Ny - Wy*Nx); I = 25, 26, 37, 38
            s(20) = s(20) + f(49)*(q(25)*Ny + q(38)*Wx - q(26)*Nx - q(37)*Wy);
            % from line 50: omega=atan2(sinP, cosP);                       % 50) omega = atan(sinP/cosP)
            % from line 47: Omega = atan2(Ny, Nx);           // line 47
            cosDenom = (cosP*cosP + sinP*sinP)^2;
            % line K = 47 dependencies (  sinP 49,   cosP 48)
            % line K = 47 dependencies (    Ny 38,     Nx 37)
            % q(47) = q(47) - q(37)*Ny/(Nx*Nx + Ny*Ny) + q(38)*Nx/(Nx*Nx + Ny*Ny);
            % cosP J = 48 KJ = -sinP/(cosP*cosP + sinP*sinP);
            m00 = 2.0*cosP*sinP/cosDenom;
            m11 = (sinP*sinP - cosP*cosP)/cosDenom;
            s(48) = s(48) + f(50)*(q(48)*m00 + q(49)*m11);
            % sinP J = 49: KJ = cosP/(cosP*cosP + sinP*sinP);
            s(49) = s(49) + f(50)*(q(48)*m11 - q(49)*m00);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % from line 51: cosnu    = rUnitx*Px + rUnity*Py + rUnitz*Pz;
            %//Px J = 18: KJ = rUnitx I = 31
            s(18) = s(18) + f(51)*q(31);
            %//Py J = 19: KJ = rUnity I = 32
            s(19) = s(19) + f(51)*q(32);
            %//Pz J = 20: KJ = rUnitz I = 33
            s(20) = s(20) + f(51)*q(33);
            %//rUnitx J = 31: KJ = Px I = 18
            s(31) = s(31) + f(51)*q(18);
            %//rUnity J = 32: KJ = Py I = 19
            s(32) = s(32) + f(51)*q(19);
            %//rUnitz J = 33: KJ = Pz I = 20
            s(33) = s(33) + f(51)*q(20);
            % from line: 52: sinnu    = Wx*(Py*rUnitz - Pz*rUnity) + ...
            %                           Wy*(Pz*rUnitx - Px*rUnitz) + ...
            %                           Wz*(Px*rUnity - Py*rUnitx);    % 52) sin(nu)
            %//Px J = 18: KJ = Wz*rUnity - Wy*rUnitz), I = 27, 32, 26, 33
            s(18) = s(18) + f(52)*(q(27)*rUnity + q(32)*Wz - q(26)*rUnitz - q(33)*Wy);
            %//Py J = 19: KJ = Wx*rUnitz - Wz*rUnitx), I = 25, 33, 27, 31
            s(19) = s(19) + f(52)*(q(25)*rUnitz + q(33)*Wx - q(27)*rUnitx - q(31)*Wz);
            %//Pz J = 20: KJ = Wy*rUnitx - Wx*rUnity), I = 26, 31, 25, 32
            s(20) = s(20) + f(52)*(q(26)*rUnitx + q(31)*Wy - q(25)*rUnity - q(32)*Wx);
            %//Wx J = 25: KJ = Py*rUnitz - Pz*rUnity), I = 19, 33, 20,  32
            s(25) = s(25) + f(52)*(q(19)*rUnitz + q(33)*Py - q(20)*rUnity - q(32)*Pz);
            %//Wy J = 26: KJ = Pz*rUnitx - Px*rUnitz), I = 20, 31, 18,  33
            s(26) = s(26) + f(52)*(q(20)*rUnitx + q(31)*Pz - q(18)*rUnitz - q(33)*Px);
            %//Wz J = 27: KJ = Px*rUnity - Py*rUnitx), I = 18, 32, 19   31
            s(27) = s(27) + f(52)*(q(18)*rUnity + q(32)*Px - q(19)*rUnitx - q(31)*Py);
            %//rUnitx J = 31: KJ = Pz*Wy - Py*Wz), I = 20, 26, 19,  27
            s(31) = s(31) + f(52)*(q(20)*Wy + q(26)*Pz - q(19)*Wz - q(27)*Py);
            %//rUnity J = 32: KJ = Px*Wz - Pz*Wx), I = 18, 27, 20,  25
            s(32) = s(32) + f(52)*(q(18)*Wz + q(27)*Px - q(20)*Wx - q(25)*Pz);
            %//rUnitz J = 33: KJ = Py*Wx- Px*Wy), I = 19, 25, 18,  26
            s(33) = s(33) + f(52)*(q(19)*Wx + q(25)*Py - q(18)*Wy - q(26)*Px);
            % from line 53: nuEpoch = atan2(sinnu, cosnu);            % 53)
            % from line 47: Omega = atan2(Ny, Nx);           // line 47
            cosDenom = (cosnu*cosnu + sinnu*sinnu)^2;
            % line K = 53 dependencies (   sinnu 52,   cosnu 51)
            % line K = 47 dependencies (    sinP 49,    cosP 48)
            % cosnu J = 51 KJ = -sinnu/(cosnu*cosnu + sinnu*sinnu);
            m00 = 2.0*cosnu*sinnu/cosDenom;
            m11 = (sinnu*sinnu - cosnu*cosnu)/cosDenom;
            s(51) = s(51) + f(53)*(q(51)*m00 + q(52)*m11);
            % sinP J = 52: KJ = cosnu/(cosnu*cosnu + sinnu*sinnu);
            s(52) = s(52) + f(53)*(q(51)*m11 - q(52)*m00);
            % from line 54: cosT    = e + cosnu;
            % All second derivatives vanish
            % from line 55: sinT    = sinnu*rootfac;                 55
            %//rootfac J = 42: KJ = sinnu I = 52
            s(42) = s(42) + f(55)*q(52);
            %//sinnu J = 52: KJ = rootfac I = 42
            s(52) = s(52) + f(55)*q(42);
            % from line 56: EccentricAnomalyEpoch = atan2(sinT, cosT); % 56)
            % from line 53:       nuEpoch = atan2(sinnu, cosnu);            % 53)
            cosDenom = (cosT*cosT + sinT*sinT)^2;
            % line K = 56 dependencies (    sinT 55,    cosT 54)
            % line K = 53 dependencies (   sinnu 52,   cosnu 51)
            % cosT J = 54 KJ = -sinT/(cosT*cosT + sinT*sinT);
            m00 = 2.0*cosT*sinT/cosDenom;
            m11 = (sinT*sinT - cosT*cosT)/cosDenom;
            s(54) = s(54) + f(56)*(q(54)*m00 + q(55)*m11);
            % sinT J = 55: KJ = cosnu/(cosnu*cosnu + sinnu*sinnu);
            s(55) = s(55) + f(56)*(q(54)*m11 - q(55)*m00);
            % from line 57: cosE    = cos(EccentricAnomalyEpoch);         % 57  (old 54)
            %//E J = 56: KJ = -sin(E), I = 56
            s(56) = s(56) - f(57)*q(56)*cosE;
            % from line 58: sinE    = sin(EccentricAnomalyEpoch);         % 58  (old 55)
            %//E J = 56: KJ = cos(E), I = 56
            s(56) = s(56) - f(58)*q(56)*sinE;
            % from line 59: MeanAnomalyEpoch = EccentricAnomalyEpoch - e*sinE;% 59)
            %//e J = 17: KJ = -sinE
            s(17) = s(17) - f(59)*q(58);
            %//sinE J = 58: KJ = -e
            s(58) = s(58) - f(59)*q(17);
            % %//e J = 17: KJ = -sinE
            % s(17) = s(17) - f(59)*q(56)*cosE;
            % %//sinE J = 58: KJ = -e*cosE
            % s(58) = s(58) - f(59)*q(17)*cosE + f(59)*q(56)*e*sinE;
            %///////////////////////////////////////////////////////////////
            % % from line 60: TimeSincePeriapsis = MeanAnomalyEpoch/meanMotion;   % 60)
            % %//MeanAnomalyEpoch J = 59: KJ = 1/meanMotion
            % s(59) = s(59) - f(60)*q(46)/meanMotion^2;
            % %//meanMotion J = 46: KJ = -MeanAnomalyEpoch/meanMotion^2 I = 46, 59
            % s(46) = s(46) + f(60)*(2.0*q(46)*MeanAnomalyEpoch - q(59))/meanMotion^2;
            % % from line 61: DeltaTime = tEpoch - TimeSincePeriapsis;     % 61)
            % % All second derivatives vanish
            % % from line 62: PTime     = DeltaTime;  % DEFAULT       % 62) in TUs
            % %if DeltaTime < 0
            % %   % The OVERWRITE PTime HERE:
            % %   PTime = DeltaTime + Period;
            % % end
            % % All second derivatives vanish
            %///////////////////////////////////////////////////////////////
            % from line 63: Mp = MeanAnomalyEpoch - meanMotion*tExtrap;      % 63 in radians
            %//tExtrap J = 8: KJ = meanMotion
            s( 8) = s( 8) - f(63)*q(46);
            %/meanMotion J = 46: KJ = tExtrap
            s(46) = s(46) - f(63)*q( 8);
            % Kepler 6: Inclination
            % from line 64: Inclination = acos(Wz);               % 64)
            %//Wz J = 27: KJ = -(1 - Wz^2)^(-1/2)
            s(27) = s(27) - f(64)*q(27)*Wz/((1.0 - Wz^2)*sqrt((1.0 - Wz^2)));
            if LineCheck > 64
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %from line 65: M = MeanAnomalyEpoch + meanMotion*tExtrap;    % 65)
                % K = 65 dependencies (MeanAnomalyEpoch 59, meanMotion 46, time 1)
                %//MeanAnomalyEpoch J = 59: KJ = 1; All second derivatives vanish
                %// meanMotion J = 46: KJ = tExtrap
                %s(46) = s(46) + f(65)*q( 8);
                % Mp J = 7: KJ = 1 (2nd derivative vanishes)!
                %// tExtrap J = 1: KJ = meanMotion
                %s( 8) = s( 8) + f(65)*q(46);
                % // from line 65: M = Mp + meanMotion*tExtrap;
                %Mp J = 63: KJ = 1 (All second derivative vanishes)!
                %meanMotion J = 46: KJ = tExtrap, I = 8
                % s(46) = s(46) + f(65)*q( 8);
                % %// tExtrap J = 8: KJ = meanMotion
                % s( 8) = s( 8) + f(65)*q(46);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line 66: E:  M = EM - e*sin(EM); % 66) EM = E(e,M)
                % Solve Kepler's Equation for E:  M = E - e*sin(E);
                % % from line 66: EM = E(e,M)
                % % K = 66 dependencies (M 65, e 2)
                % % M J = 65 KJ = dE_dM (i = 65, 2)
                % s(65) = s(65) + f(66)*(q( 2)*d2E_dMde + q(65)*d2E_dMdM);
                % % e J = 2 KJ = dE_de (i = 40,2)
                % s( 2) = s( 2) + f(66)*(q( 2)*d2E_dede + q(65)*d2E_dedM);
                % from line 66: EM = E(e,M)
                % K = 66 dependencies (M 65, e 17)
                % M J = 65 KJ = dE_dM (i = 65, 17)
                s(65) = s(65) + f(66)*(q(17)*d2E_dMde + q(65)*d2E_dMdM);
                % e J = 2 KJ = dE_de (i = 40,2)
                s(17) = s(17) + f(66)*(q(17)*d2E_dede + q(65)*d2E_dedM);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line 67: cosK  = cos(EM);                  % line 67
                %//EM J = 66: KJ = -sin(EM), I = 66
                s(66) = s(66) - f(67)*q(66)*cosK;
                % from line 68: sinK  = sin(EM);                  % line 68
                %//EM J = 66: KJ = cos(EM), I = 66
                s(66) = s(66) - f(68)*q(66)*sinK;
                % from line 70: tanX  = cosK - e;                 % line 70
                %// All second derivatives vanish
                % from line 71: tanY  = rootfac*sinK;             % line 71
                %//rootfac J = 42: KJ = sinK, I = 68
                s(42) = s(42) + f(71)*q(68);
                %//sinK J = 68: KJ = rootfac, I = 42
                s(68) = s(68) + f(71)*q(42);
                % from line 72: nu     = atan2(tanY, tanX);       % line 72
                % from line 56: EccentricAnomalyEpoch = atan2(sinT, cosT); % 56)
                cosDenom = (tanX*tanX + tanY*tanY)^2;
                % line K = 72 dependencies (    tanY 71,    tanX 70)
                % line K = 56 dependencies (    sinT 55,    cosT 54)
                m00 = 2.0*tanX*tanY/cosDenom;
                m11 = (tanY*tanY - tanX*tanX)/cosDenom;
                % tanX J = 70 KJ = -tanY/(tanX*tanX + tanY*tanY);
                s(70) = s(70) + f(72)*(q(70)*m00 + q(71)*m11);
                % tanY J = 71: KJ = tanX/(tanX*tanX + tanY*tanY);
                s(71) = s(71) + f(72)*(q(70)*m11 - q(71)*m00);
                % from line 73: coss  = cos(nu);                  % line 73
                %//nu J = 72: KJ = -sin(nu), I = 72
                s(72) = s(72) - f(73)*q(72)*coss;
                % from line 74: sins  = sin(nu);                  % line 74
                %//nu J = 72: KJ = cos(nu);
                s(72) = s(72) - f(74)*q(72)*sins;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % % from line 75: rorbit   = p/(1.0 + e*coss)       % line 75
                % %// from line 28: double rorbit   = p/(1.0 + e*coss);                          // line 28
                % %// nu: J = 72
                % m00  = ((coss*e+1.0)*coss + 2.0*sins*sins*e)*e*p/cecube;
                % m01  = -(coss*e-1.0)*sins*p/cecube;
                % m02  = sins*e/cesqr;
                % s(72) = s(72) + f(75)*(q(72)*m00 + q(17)*m01 + q(24)*m02);
                %
                % %// e: J = 17
                % m10  = -(coss*e-1.0)*sins*p/cecube;
                % m11  = 2.0*coss*coss*p/cecube;
                % m12  = -coss/cesqr;
                % s(17) = s(17) + f(75)*(q(72)*m10 + q(17)*m11 + q(24)*m12);
                %
                % %// p: J = 24
                % m20  =  sins*e/cesqr;
                % m21  = -coss/cesqr;
                % m22  =  0.0;
                % s(24) = s(24) + f(75)*( q(72)*m20 + q(17)*m21 + q(24)*m22);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % // from line 75: rorbit   = a*(1.0 - e*cosK);
                % cosK: J=67 KJ = -a*e
                s(67) = s(67) - f(75)*(          q(43)*e + q(17)*a);
                % a: J = 43  KJ = 1 - e*cosK
                s(43) = s(43) - f(75)*(q(67)*e +           q(17)*cosK);
                % e: J = 7 KJ = -a*cosK
                s(17) = s(17) - f(75)*(q(67)*a + q(43)*cosK);
                % E: J = 41 KJ = +a*e*sin(E)
                % m15 = -e;
                % m02 = -a;
                % s(42) = s(42) + f(49)*(q(15)*m15 + q( 2)*m02);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line 76: rorbitx  = rorbit*(coss*Px + sins*Qx);           % line 76
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line 29: double rorbitx  = rorbit*(coss*Px + sins*Qx);              // line 29
                %//nu: J = 72
                m00 = -(coss*Px + sins*Qx)*rorbit;
                m01 =  coss*Qx - sins*Px;
                m02 = -sins*rorbit;
                m03 =  coss*rorbit;
                s(72) = s(72) +  f(76)*(q(72)*m00 + q(75)*m01 + q(18)*m02 + q(28)*m03);
                %// rorbit: J = 75
                m10 =  coss*Qx - sins*Px;
                m11 =  0.0;
                m12 =  coss;
                m13 =  sins;
                s(75) = s(75) + f(76)*(q(72)*m10 + q(75)*m11 + q(18)*m12 + q(28)*m13);
                %// Px: J = 18
                m20 = -sins*rorbit;
                m21 =  coss;
                m22 =  0.0;
                m23 =  0.0;
                s(18) = s(18) + f(76)*(q(72)*m20 + q(75)*m21 + q(18)*m22 + q(28)*m23);
                %// Qx: J = 28
                m30 =  coss*rorbit;
                m31 =  sins;
                m32 =  0.0;
                m33 =  0.0;
                s(28) = s(28) + f(76)*(q(72)*m30 + q(75)*m31 + q(18)*m32 + q(28)*m33);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %rorbity  = rorbit*(coss*Py + sins*Qy);           % line 77
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %// from line 30: double rorbity  = rorbit*(coss*Py + sins*Qy);              // line 30
                %// nu: J = 72
                m00 = -(coss*Py + sins*Qy)*rorbit;
                m01 =  coss*Qy - sins*Py;
                m02 = -sins*rorbit;
                m03 =  coss*rorbit;
                s(72) = s(72) + f(77)*(q(72)*m00 + q(75)*m01 + q(19)*m02 + q(29)*m03);
                %// rorbit: J = 75
                m10 =  coss*Qy - sins*Py;
                m11 =  0.0;
                m12 =  coss;
                m13 =  sins;
                s(75) = s(75) + f(77)*(q(72)*m10 + q(75)*m11 + q(19)*m12 + q(29)*m13);
                %// Py: J = 19
                m20 = -sins*rorbit;
                m21 =  coss;
                m22 =  0.0;
                m23 =  0.0;
                s(19) = s(19) + f(77)*(q(72)*m20 + q(75)*m21 + q(19)*m22 + q(29)*m23);
                %// Qy: J = 29
                m30 =  coss*rorbit;
                m31 =  sins;
                m32 =  0.0;
                m33 =  0.0;
                s(29) = s(29) + f(77)*(q(72)*m30 + q(75)*m31 + q(19)*m32 + q(29)*m33);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %rorbitz  = rorbit*(coss*Pz + sins*Qz);           % line 78
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %// from line 31: double rorbitz  = rorbit*(coss*Pz + sins*Qz);              // line 31
                %// nu: J = 72
                m00 = -(coss*Pz + sins*Qz)*rorbit;
                m01 =  coss*Qz - sins*Pz;
                m02 = -sins*rorbit;
                m03 =  coss*rorbit;
                s(72) = s(72) + f(78)*(q(72)*m00 + q(75)*m01 + q(20)*m02 + q(30)*m03);
                %// rorbit: J = 75
                m10 =  coss*Qz - sins*Pz;
                m11 =  0.0;
                m12 =  coss;
                m13 =  sins;
                s(75) = s(75) + f(78)*(q(72)*m10 + q(75)*m11 + q(20)*m12 + q(30)*m13);
                %// Pz: J = 20
                m20 = -sins*rorbit;
                m21 =  coss;
                m22 =  0.0;
                m23 =  0.0;
                s(20) = s(20) + f(78)*(q(72)*m20 + q(75)*m21 + q(20)*m22 + q(30)*m23);
                %// Qz: J = 30
                m30 =  coss*rorbit;
                m31 =  sins;
                m32 =  0.0;
                m33 =  0.0;
                s(30) = s(30) + f(78)*(q(72)*m30 + q(75)*m31 + q(20)*m32 + q(30)*m33);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line 69: rtpinv   = sqrt(mu/p);
                %p J = 24: KJ = -sqrt(mu)*p^(-1.5);
                s(24) = s(24) + 0.75*f(69)*q(24)*rtpinv/psq;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %vorbitx  = rtpinv*(-sins*Px + (e + coss)*Qx);    % line 79
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %// from line 32: double vorbitx  = rtpinv*(-sin(s)*Px + (e + cos(s))*Qx);  // line 32
                %// nu: J = 72
                m00  = -rtpinv*(coss*Qx-(sins*Px));
                m01  =  0.0;
                m02  =  (coss*Px+sins*Qx)/(2.0*rtpinv*psq);
                m03  = -rtpinv*coss;
                m04  = -rtpinv*sins;
                s(72) = s(72) + f(79)*(q(72)*m00 + q(17)*m01 + q(24)*m02 + q(18)*m03 + q(28)*m04);
                %// e: J = 17
                m10  =  0.0;
                m11  =  0.0;
                m12  = -Qx/(2.0*rtpinv*psq);
                m13  =  0.0;
                m14  =  rtpinv;
                s(17) = s(17) + f(79)*(q(72)*m10 + q(17)*m11 + q(24)*m12 + q(18)*m13 + q(28)*m14);
                %// p: J = 24
                m20  =  (coss*Px+sins*Qx)/(2.0*rtpinv*psq);
                m21  = -Qx/(2.0*rtpinv*psq);
                m22  =  (3.0*((coss+e)*Qx-(sins*Px)))/(4.0*rtpinv*pcubed);
                m23  =  sins/(2.0*rtpinv*psq);
                m24  = -(coss+e)/(2.0*rtpinv*psq);
                s(24) = s(24) + f(79)*(q(72)*m20 + q(17)*m21 + q(24)*m22 + q(18)*m23 + q(28)*m24);
                %// Px: J = 18
                m30  = -rtpinv*coss;
                m31  =  0.0;
                m32  =  sins/(2.0*rtpinv*psq);
                m33  =  0.0;
                m34  =  0.0;
                s(18) = s(18) + f(79)*(q(72)*m30 + q(17)*m31 + q(24)*m32 + q(18)*m33 + q(28)*m34);
                %// Qx: J = 28
                m40  = -rtpinv*sins;
                m41  =  rtpinv;
                m42  = -(coss+e)/(2.0*rtpinv*psq);
                m43  =  0.0;
                m44  =  0.0;
                s(28) = s(28) + f(79)*(q(72)*m40 + q(17)*m41 + q(24)*m42 + q(18)*m43 + q(28)*m44);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %vorbity  = rtpinv*(-sins*Py + (e + coss)*Qy);    % line 80
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %// from line 33: double vorbity  = rtpinv*(-sins*Py + (e + coss)*Qy);  // line 33
                %// nu: J = 72
                m00  = -rtpinv*(coss*Qy-(sins*Py));
                m01  =  0.0;
                m02  =  (coss*Py+sins*Qy)/(2.0*rtpinv*psq);
                m03  = -rtpinv*coss;
                m04  = -rtpinv*sins;
                s(72) = s(72) + f(80)*(q(72)*m00 + q(17)*m01 + q(24)*m02 + q(19)*m03 + q(29)*m04);
                %// e: J = 17
                m10  =  0.0;
                m11  =  0.0;
                m12  = -Qy/(2.0*rtpinv*psq);
                m13  =  0.0;
                m14  =  rtpinv;
                s(17) = s(17) + f(80)*(q(72)*m10 + q(17)*m11 + q(24)*m12 + q(19)*m13 + q(29)*m14);
                %// p: J = 24
                m20  =  (coss*Py+sins*Qy)/(2.0*rtpinv*psq);
                m21  = -Qy/(2.0*rtpinv*psq);
                m22  =  (3.0*((coss+e)*Qy-(sins*Py)))/(4.0*rtpinv*pcubed);
                m23  =  sins/(2.0*rtpinv*psq);
                m24  = -(coss+e)/(2.0*rtpinv*psq);
                s(24) = s(24) + f(80)*(q(72)*m20 + q(17)*m21 + q(24)*m22 + q(19)*m23 + q(29)*m24);
                %// Py: J = 19
                m30  = -rtpinv*coss;
                m31  =  0.0;
                m32  =  sins/(2.0*rtpinv*psq);
                m33  =  0.0;
                m34  =  0.0;
                s(19) = s(19) + f(80)*(q(72)*m30 + q(17)*m31 + q(24)*m32 + q(19)*m33 + q(29)*m34);
                %// Qy: J = 29
                m40  = -rtpinv*sins;
                m41  =  rtpinv;
                m42  = -(coss+e)/(2.0*rtpinv*psq);
                m43  =  0.0;
                m44  =  0.0;
                s(29) = s(29) + f(80)*(q(72)*m40 + q(17)*m41 + q(24)*m42 + q(19)*m43 + q(29)*m44);                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %vorbitz  = rtpinv*(-sins*Pz + (e + coss)*Qz);    % line 81
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %// from line 34:  double vorbitz  = rtpinv*(-sins*Pz + (e + coss)*Qz);  // line 34
                %// nu: J = 72
                m00  = -rtpinv*(coss*Qz-(sins*Pz));
                m01  =  0.0;
                m02  =  (coss*Pz+sins*Qz)/(2.0*rtpinv*psq);
                m03  = -rtpinv*coss;
                m04  = -rtpinv*sins;
                s(72) = s(72) + f(81)*(q(72)*m00 + q(17)*m01 + q(24)*m02 + q(20)*m03 + q(30)*m04);
                %// e: J = 17
                m10  =  0.0;
                m11  =  0.0;
                m12  = -Qz/(2.0*rtpinv*psq);
                m13  =  0.0;
                m14  =  rtpinv;
                s(17) = s(17) + f(81)*(q(72)*m10 + q(17)*m11 + q(24)*m12 + q(20)*m13 + q(30)*m14);
                %// p: J = 24
                m20  =  (coss*Pz+sins*Qz)/(2.0*rtpinv*psq);
                m21  = -Qz/(2.0*rtpinv*psq);
                m22  =  (3.0*((coss+e)*Qz-(sins*Pz)))/(4.0*rtpinv*pcubed);
                m23  =  sins/(2.0*rtpinv*psq);
                m24  = -(coss+e)/(2.0*rtpinv*psq);
                s(24) = s(24) + f(81)*(q(72)*m20 + q(17)*m21 + q(24)*m22 + q(20)*m23 + q(30)*m24);
                %// Pz: J = 20
                m30  = -rtpinv*coss;
                m31  =  0.0;
                m32  =  sins/(2.0*rtpinv*psq);
                m33  =  0.0;
                m34  =  0.0;
                s(20) = s(20) + f(81)*(q(72)*m30 + q(17)*m31 + q(24)*m32 + q(20)*m33 + q(30)*m34);
                %// Qz: J = 30
                m40  = -rtpinv*sins;
                m41  =  rtpinv;
                m42  = -(coss+e)/(2.0*rtpinv*psq);
                m43  =  0.0;
                m44  =  0.0;
                s(30) = s(30) + f(81)*(q(72)*m40 + q(17)*m41 + q(24)*m42 + q(20)*m43 + q(30)*m44);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if LineCheck > 81
                    % from line 82: Aeq      = omega + Omega;        % line 82
                    % All second derivatives vanish
                    % from line 83: cosAeq   = cos(Aeq);
                    %//Aeq J = 82: KJ = -sin(Aeq)
                    s(82) = s(82) - f(83)*q(82)*cosAeq;
                    % from line 84: sinAeq   = sin(Aeq);             % line 84
                    %//Aeq J = 82: KJ = cos(Aeq)
                    s(82) = s(82) - f(84)*q(82)*sinAeq;
                    % from line 85: tanHalf  = tan(Inclination/2.0); % line 85
                    %//Inclination J = 64: KJ = (1 + tanHalf^2)/2
                    s(64) = s(64) + 0.5*f(85)*q(64)*tanHalf*(1.0 + tanHalf^2);
                    % from line 86: Feq      = e*cosAeq;             % line 86
                    %// e J =  17: KJ = cosAeq I = 83
                    s(17) = s(17) + f(86)*q(83);
                    %//cosAeq J = 83: KJ = e
                    s(83) = s(83) + f(86)*q(17);
                    % from line 87: Geq      = e*sinAeq;             % line 87
                    %//e J = 17: KJ = sinAeq
                    s(17) = s(17) + f(87)*q(84);
                    %//sinAeq J = 84: KJ = e
                    s(84) = s(84) + f(87)*q(17);
                    
                    % from line 88: Heq      = tanHalf*cos(Omega);    % line 88
                    %q(88) = q(88) + q(85)*cosO - q(47)*tanHalf*sinO;
                    % from line 89: Keq      = tanHalf*sin(Omega);    % line 89
                    %q(89) = q(89) + q(85)*sinO + q(47)*tanHalf*cosO;
                    
                    % from line 88: Heq      = tanHalf*cosO;     % line 88
                    %//Omega J = 47: KJ = -tanHalf*sinO
                    s(47) = s(47) - f(88)*(q(47)*Heq + q(85)*sinO);
                    %//tanHalf J = 85: KJ = cosO
                    s(85) = s(85) - f(88)*q(47)*sinO;
                    % from line 89: Keq      = tanHalf*sinO;     % line 89
                    %//Omega J = 47: KJ = tanHalf*cosO
                    s(47) = s(47) + f(89)*(-q(47)*Keq + q(85)*cosO);
                    %//tanHalf J = 85: KJ = sinO
                    s(85) = s(85) + f(89)*q(47)*cosO;
                    % from line 90: Leq    = Aeq + nu;               % line 90
                    % All second derivatives vanish
                    % from line 91: CosL     = cos(Leq);             % line 91
                    %//Leq J = 90: KJ = -sin(Leg) I = 90
                    s(90) = s(90) - f(91)*q(90)*CosL;
                    % from line 92: SinL     = sin(Leq);             % line 92
                    %//LEQ J = 90: KJ = cos(Leq), I = 90
                    s(90) = s(90) - f(92)*q(90)*SinL;
                    % from line 93: alphaSq = Heq^2 - Keq^2;         % line 93
                    %//Heq J = 88: KJ = 2.0*Heq, I = 88
                    s(88) = s(88) + 2.0*f(93)*q(88);
                    %//Keq J = 89: KJ = -2.0*Keq, I = 89
                    s(89) = s(89) - 2.0*f(93)*q(89);
                    % from line 94: Seq     = 1.0 + Heq^2 + Keq^2;   % line 94
                    %//Heq J = 88: KJ = 2.0*Heq, I = 88
                    s(88) = s(88) + 2.0*f(94)*q(88);
                    %//Keq J = 89: KJ = 2.0*Keq, I = 89
                    s(89) = s(89) + 2.0*f(94)*q(89);
                    % from line 95: Weq = 1.0 + Feq*CosL + Geq*SinL; % line 95
                    %//Feq J = 86: KJ = CosL, I = 91
                    s(86) = s(86) + f(95)*q(91);
                    %//Geq J = 87: KJ = SinL, I = 92
                    s(87) = s(87) + f(95)*q(92);
                    %//CosL J = 91: KJ = Feq, I = 86
                    s(91) = s(91) + f(95)*q(86);
                    %//SinL J = 92: KJ = Geq, I = 87
                    s(92) = s(92) + f(95)*q(87);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % from line 96: Req     = p/Weq;              % line 96
                    %//p J = 24: KJ = 1/Weq
                    s(24) = s(24) -f(96)*q(95)/Weq^2;
                    %//Weq J = 95: KJ = -p/Weq^2, I = 24, 95
                    s(95) = s(95) + f(96)*(2.0*q(95)*Req - q(24))/Weq^2;
                    
                    % from line 97: RovS    = Req/Seq;            % line 97
                    %//Req J = 96: KJ = 1/Seq
                    s(96) = s(96) -f(97)*q(94)/Seq^2;
                    %//Seq J = 94: KJ = -Req/Seq^2, I = 96, 94
                    s(94) = s(94) + f(97)*(2.0*q(94)*RovS - q(96))/Seq^2;
                    %    q(97) = q(97) + q(96)/Seq - q(94)*RovS/Seq;
                    
                    % from line 98: srtpinv = rtpinv/Seq;         % line 98
                    %//rtpinv J = 69: KJ = 1/Seq, I = 94
                    s(69) = s(69) - f(98)*q(94)/Seq^2;
                    %//Seq J = 94: KJ = -rtpinv/Seq^2, I = 69, 94
                    s(94) = s(94) + f(98)*(2.0*q(94)*srtpinv - q(69))/Seq^2;
                    % q(98) = q(98) + q(69)/Seq - q(94)*srtpinv/Seq;
                    
                    % from line 99: HK  = Heq*Keq;            % line 99
                    %//Heq J = 88: KJ = Keq
                    s(88) = s(88) + f(99)*q(89);
                    s(89) = s(89) + f(99)*q(88);
                    %  q(99) = q(99) + q(88)*Keq + q(89)*Heq;
                    
                    % from line 100: OnePalphaSqCosL = (1.+alphaSq)*CosL;% line 100
                    %//alphaSq J = 93: KJ = CosL, I = 91
                    s(93) = s(93) + f(100)*q(91);
                    %//CosL J = 91: KJ = (1 + alphaSq), I = 93
                    s(91) = s(91) + f(100)*q(93);
                    % q(100) = q(100) + q(91)*(1.0 + alphaSq) + q(93)*CosL;
                    
                    % from line 101: OneMalphaSqCosL = (1.-alphaSq)*CosL;% line 101
                    %//alphaSq J = 93: KJ =-CosL, I = 91
                    s(93) = s(93) - f(101)*q(91);
                    %//CosL J = 91: KJ = (1 - alphaSq), I = 93
                    s(91) = s(91) - f(101)*q(93);
                    % q(101) = q(101) + q(91)*(1.0 - alphaSq) - q(93)*CosL;
                    
                    % from line 102: OnePalphaSqSinL = (1.+alphaSq)*SinL;% line 102
                    %//alphaSq J = 93: KJ = SinL, I = 92
                    s(93) = s(93) + f(102)*q(92);
                    %//SinL J = 92: KJ = (1 + alphaSq), I = 93
                    s(92) = s(92) + f(102)*q(93);
                    % q(102) = q(102) + q(92)*(1.0 + alphaSq) + q(93)*SinL;
                    
                    % from line 103: OneMalphaSqSinL = (1.-alphaSq)*SinL;% line 103
                    %//alphaSq J = 93: KJ =-SinL, I = 92
                    s(93) = s(93) - f(103)*q(92);
                    %//SinL J = 92: KJ = (1 - alphaSq), I = 93
                    s(92) = s(92) - f(103)*q(93);
                    % q(103) = q(103) + q(92)*(1.0 - alphaSq) - q(93)*SinL;
                    
                    % from line 104: Xfac    = OnePalphaSqCosL + 2.0*HK*SinL; % line 104
                    %//HK J = 99: KJ = 2.0*SinL
                    s(99) = s(99) + f(104)*2.0*q(92);
                    %//SinL J = 92: KJ = 2.0*HK
                    s(92) = s(92) + f(104)*2.0*q(99);
                    % OnePalphaSqCosL J = 100 -- All second derivatives vanish
                    % q(104) = q(104) + 2.0*q(92)*HK + 2.0*q(99)*SinL + q(100);
                    
                    % from line 105: Yfac    = OneMalphaSqSinL + 2.0*HK*CosL; % line 105
                    %//HK J = 99: KJ = 2.0*CosL
                    s(99) = s(99) + f(105)*2.0*q(91);
                    %//CosL J = 91: KJ = 2.0*HK
                    s(91) = s(91) + f(105)*2.0*q(99);
                    % OnePalphaSqCosL J = 100 -- All second derivatives vanish
                    % q(105) = q(105) + 2.0*q(91)*HK + 2.0*q(99)*CosL + q(103);
                    
                    % from line 106: Zfac    = Heq*SinL - Keq*CosL;          % line 106
                    s(88) = s(88) + f(106)*q(92);
                    s(89) = s(89) - f(106)*q(91);
                    s(91) = s(91) - f(106)*q(89);
                    s(92) = s(92) + f(106)*q(88);
                    % q(106) = q(88) + q(88)*SinL - q(89)*CosL  - q(91)*Keq  + q(92)*Heq;
                    
                    % from line 107: VXfac =  OnePalphaSqSinL - 2.0*HK*(Feq+CosL) + Geq*(1+alphaSq);% line 107
                    %//Feq J = 86: KJ = -2.0*HK, I = 99
                    s(86) = s(86) - 2.0*f(107)*q(99);
                    %//Geq J = 87: KJ = (1+alphaSq), I = 93
                    s(87) = s(87) + f(107)*q(93);
                    %//CosL J = 91: KJ = -2.0*HK, I = 99
                    s(91) = s(91) - 2.0*f(107)*q(99);
                    %//alphaSq J = 93: KJ = Geq, I =87
                    s(93) = s(93) + f(107)*q(87);
                    %//HK J = 99: KJ = -2.0*(Feq+CosL)
                    s(99) = s(99) -2.0*f(107)*(q(86) + q(91));
                    %//OnePalphaSqSinL J=102 -- All second derivatives vanish
                    % q(107) = q(107) -2.0*q(86)*HK +q(87)*(1.0 + alphaSq) -2.0*q(91)*HK +q(93)*Geq -2.0*q(99)*(Feq+CosL) +q(102);
                    
                    % from line 108: VYfac = -OneMalphaSqCosL + 2.0*HK*(Geq+SinL) + Feq*(alphaSq-1) ;% line 108
                    %//Feq J = 86: KJ = (alphaSq-1), I = 93
                    s(86) = s(86) + f(108)*q(93);
                    %//Geq J = 87: KJ = 2.0*HK, I = 99
                    s(87) = s(87) + 2.0*f(108)*q(99);
                    %//SinL J = 92: KJ = 2.0*HK, I = 99
                    s(92) = s(92) + 2.0*f(108)*q(99);
                    %//alphaSq J = 93: KJ = Feq, I =86
                    s(93) = s(93) + f(108)*q(86);
                    %//HK J = 99: KJ = 2.0*(Geq+SinL)
                    s(99) = s(99) +2.0*f(108)*(q(87) + q(92));
                    %//OneMalphaSqCosL J=101 -- All second derivatives vanish
                    %q(108) = q(108) +2.0*q(87)*HK +q(86)*(alphaSq - 1.0) +2.0*q(92)*HK +q(93)*Feq +2.0*q(99)*(Geq+SinL) -q(101);
                    
                    % from line: VZfac   =  Heq*(Feq + CosL) + Keq*(Geq + SinL) % line 109
                    %//Feq J = 86: KJ = Heq, I = 88
                    s(86) = s(86) + f(109)*q(88);
                    %//Geq J = 87: KJ = Keq, I = 89
                    s(87) = s(87) + f(109)*q(89);
                    %//Heq J = 88: KJ = (Feq + CosL)
                    s(88) = s(88) + f(109)*(q(86) + q(91));
                    %//Keq J = 89: KJ = (Geq + SinL)
                    s(89) = s(89) + f(109)*(q(87) + q(92));
                    %//CosL J = 91: KJ = Heq, I = 88
                    s(91) = s(91) + f(109)*q(88);
                    %//SinL J = 92: KJ = Keq, I = 89
                    s(92) = s(92) + f(109)*q(89);
                    %q(109) = q(109) +q(86)*Heq +q(87)*Keq +q(88)*(Feq+CosL) +q(89)*(Geq+SinL) +q(91)*Heq +q(92)*Keq;
                    
                    % from line 110: Xeq     =     RovS*Xfac;     % line 110
                    %//RovS J = 97: KJ = Xfac, I = 104
                    s( 97) = s( 97) + f(110)*q(104);
                    %//Xfac J = 104: KJ = RovS, I = 97
                    s(104) = s(104) + f(110)*q(97);
                    %q(110) = q(110) + q(97)*Xfac + q(104)*RovS;
                    
                    % from line 111: Yeq     =     RovS*Yfac;     % line 111
                    %//RovS J = 97: KJ = Yfac, I = 105
                    s( 97) = s( 97) + f(111)*q(105);
                    %//Yfac J = 105: KJ = RovS, I = 97
                    s(105) = s(105) + f(111)*q(97);
                    %q(111) = q(111) + q(97)*Yfac + q(105)*RovS;
                    
                    %from line 112: Zeq      = 2.0*RovS*Zfac;     % line 112
                    %//RovS J = 97: KJ = 2.0*Zfac, I = 106
                    s( 97) = s( 97) + 2.0*f(112)*q(106);
                    %//Zfac J = 106: KJ = 2.0*RovS, I = 97
                    s(106) = s(106) + 2.0*f(112)*q(97);
                    %q(112) = q(112) + 2.0*(q(97)*Zfac + q(106)*RovS);
                    
                    % from line 113: VXeq    =    -srtpinv*VXfac; % line 113
                    %//srtpinv J = 98: KJ = -VXfac, I = 107
                    s( 98) = s( 98) - f(113)*q(107);
                    %//VXfac J = 107: KJ = -srtpinv, I = 98
                    s(107) = s(107) - f(113)*q( 98);
                    %q(113) = q(113) - q(98)*VXfac - q(107)*srtpinv;
                    
                    % from line 114: VYeq    =    -srtpinv*VYfac; % line 114
                    %//srtpinv J = 98: KJ = -VYfac, I = 108
                    s( 98) = s( 98) - f(114)*q(108);
                    %//VYfac J = 108: KJ = -srtpinv, I = 98
                    s(108) = s(108) - f(114)*q( 98);
                    %q(114) = q(114) - q(98)*VYfac - q(108)*srtpinv;
                    
                    % from line 115: VZeq    = 2.0*srtpinv*VZfac; % line 115
                    %//srtpinv J = 98: KJ = -VZfac, I = 109
                    s( 98) = s( 98) - f(115)*q(109);
                    %//VZfac J = 109: KJ = -srtpinv, I = 98
                    s(109) = s(109) - f(115)*q( 98);
                    %q(115) = q(115) + 2.0*q(98)*VZfac + 2.0*q(109)*srtpinv;
                end
            end
            
            %
            %                 t       = this.time;                                 % line   1
            %                 e       = this.e;                                    % line   2
            %                s a       = this.a;                                    % line   3
            %                 %p       = this.p;                                    % line   3
            %                 I       = this.Inclination;                          % line   4
            %                 omega   = this.omega;                                % line   5
            %                 Omega   = this.Omega;                                % line   6
            %                 Mp      = this.Mp;                                   % line   7
            %
            %                 onePe   = this.onePe;                                % line  11
            %                 oneMe   = this.oneMe;                                % line  12
            %                 fac     = this.fac;                                  % line  13
            %                 rootfac = this.rootfac;                              % line  14
            %                 %a       = this.a;                                   % line  15
            %                 p       = this.p;                                    % line  15
            %                 meanMotion = this.meanMotion;                        % line  16
            %
            %                 cosI    = this.cosI;                                 % line  17
            %                 sinI    = this.sinI;                                 % line  18
            %                 cosom   = this.cosom;                                % line  19
            %                 sinom   = this.sinom;                                % line  20
            %                 cosO    = this.cosO;                                 % line  21
            %                 sinO    = this.sinO;                                 % line  22
            %
            %                 Px      = this.Px;                                   % line  23
            %                 Py      = this.Py;                                   % line  24
            %                 Pz      = this.Pz;                                   % line  25
            %                 Qx      = this.Qx;                                   % line  26
            %                 Qy      = this.Qy;                                   % line  27
            %                 Qz      = this.Qz;                                   % line  28
            %
            %                 M       = this.M;                                    % line  40
            %                 E       = this.E;                                    % line  41
            %                 cosK    = this.cosK;                                 % line  42
            %                 sinK    = this.sinK;                                 % line  43
            %                 tanX    = this.tanX;                                 % line  44
            %                 tanY    = this.tanY;                                 % line  45
            %                 nu      = this.nu;                                   % line  46
            %                 coss    = this.coss;                                 % line  47
            %                 sins    = this.sins;                                 % line  48
            %                 rorbit  = this.rorbit;                               % line  49
            %                 rorbitx = this.rorbitx;                              % line  50
            %                 rorbity = this.rorbity;                              % line  51
            %                 rorbitz = this.rorbitz;                              % line  52
            %                 rtpinv  = this.rtpinv;                               % line  53
            %                 vorbitx = this.vorbitx;                              % line  54
            %                 vorbity = this.vorbity;                              % line  55
            %                 vorbitz = this.vorbitz;                              % line  56
            %
            %                 dE_dM   = this.dE_dM;
            %                 dE_de   = this.dE_de;
            %                 d2E_dMdM = this.d2E_dMdM;
            %                 d2E_dMde = this.d2E_dMde;
            %                 d2E_dedM = this.d2E_dedM;
            %                 d2E_dede = this.d2E_dede;
            %                 rmag     = this.rmag;
            %                 denom    = rmag*rmag*rmag;
            %                 edenom   = e*e*e;
            %                 cesqr    = (coss*e+1.0)*(coss*e+1.0);
            %                 cecube   = cesqr*(coss*e+1.0);
            %                 psq      = p*p;
            %                 pcubed   = p*p*p;
            %                 pdenom   = sqrt(pcubed);
            %                 cosDenom = 0.0;
            %
            %                 CK = this.CK;
            %                 C2 = this.C2;
            %                 C3 = this.C3;
            %                 C4 = this.C4;
            %                 C5 = this.C5;
            %                 C6 = this.C6;
            %
            %                 % Escobal Gravity Terms
            %                 X        = this.X;
            %                 Y        = this.Y;
            %                 Z        = this.Z;
            %                 R2       = this.R2;
            %                 R1       = this.R1;
            %                 R3       = this.R3;
            %                 R4       = this.R4;
            %                 R5       = this.R5;
            %                 R6       = this.R6;
            %                 %sinomnu  = this.sinomnu;
            %                 %cosomnu  = this.cosomnu;
            %                 %sd1chk   = this.sd1chk;
            %                 %cd1chk   = this.cd1chk;
            %                 sd1      = this.sd1;
            %                 sd2      = this.sd2;
            %                 sd3      = this.sd3;
            %                 sd4      = this.sd4;
            %                 sd5      = this.sd5;
            %                 sd6      = this.sd6;
            %                 F2       = this.F2;
            %                 F3       = this.F3;
            %                 F4       = this.F4;
            %                 F5       = this.F5;
            %                 F6       = this.F6;
            %                 V1       = this.V1;
            %                 V2       = this.V2;
            %                 V3       = this.V3;
            %                 V4       = this.V4;
            %                 V5       = this.V5;
            %                 V6       = this.V6;
            %                 V        = this.V;
            %                 D        = this.D;
            %                 Phi      = this.Phi;
            %
            %                 fac1      = this.fac1;
            %                 fac2      = this.fac2;
            %                 fac3      = this.fac3;
            %                 fac4      = this.fac4;
            %                 fac5      = this.fac5;
            %                 Mwe       = this.Mwe;
            %                 MeM       = this.MeM;
            %                 MaM       = this.MaM;
            %                 MWi       = this.MWi;
            %                 Miw       = this.Miw;
            %                 %this.MLagrange = MLagrange;
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 los       = this.los;   % Rextrap - Rsensor
            %                 lx        = this.lx;                                 % line 101
            %                 ly        = this.ly;                                 % line 102
            %                 lz        = this.lz;                                 % line 103
            %                 lperpSq   = this.lperpSq;                            % line 104
            %                 rangeSq   = this.rangeSq;                            % line 105
            %                 range     = this.range;                              % line 106
            %                 zLos      = this.zLos;                               % line 107
            %                 thetaLos  = this.thetaLos;                           % line 108
            %                 phiLos    = this.phiLos;                             % line 109
            %                 vlos      = this.vlos;                          % lines 110-112
            %                 vx        = this.vx;                                 % line 110
            %                 vy        = this.vy;                                 % line 111
            %                 vz        = this.vz;                                 % line 112
            %                 ldotv     = this.ldotv;                              % line 113
            %                 ldot      = this.ldot;                               % line 114
            %                 lperp     = this.lperp;                              % line 115
            %                 thdotfac1 = this.thdotfac1;                          % line 116
            %                 thdotfac2 = this.thdotfac2;                          % line 117
            %                 thdot     = this.thdot;                              % line 118
            %                 phiDotNum = this.phiDotNum;                          % line 119
            %                 phidot    = this.phidot;                             % line 120
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 thetaM      = this.thetaM;
            %                 cosThetaM   = this.cosThetaM;
            %                 sinThetaM   = this.sinThetaM;
            %                 phiM        = this.phiM;
            %                 cosPhiM     = this.cosPhiM;
            %                 sinPhiM     = this.sinPhiM;
            %                 cosThetaLos = this.cosThetaLos;                      % line 121
            %                 sinThetaLos = this.sinThetaLos;                      % line 122
            %                 phiDiff     = this.phiDiff;                          % line 123
            %                 cosPhiDiff  = this.cosPhiDiff;                       % line 124
            %                 sinPhiDiff  = this.sinPhiDiff;                       % line 125
            %                 delUV       = this.delUV;                            % line 126
            %                 Uview       = this.Uview;                            % line 127
            %                 Vview       = this.Vview;                            % line 128
            %                 xLos        = this.xLos;                             % line 129
            %                 yLos        = this.yLos;                             % line 130
            %                 xLosM       = this.xLosM;
            %                 yLosM       = this.yLosM;
            %                 zLosM       = this.zLosM;
            %                 Collinearity  = this.Collinearity;                   % line 131
            %                 Mdelta      = this.Mdelta;
            %                 InvCovM     = this.InvCovM;
            %                 losDiff     = this.losDiff;                          % line 132
            %                 xUnit       = this.xUnit;
            %                 yUnit       = this.yUnit;
            %                 zUnit       = this.zUnit;
            %                 ChiSqLos    = this.ChiSqLos;                         % line 133
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 % % Semi-Major Axis:
            %                 % % 	// from line 11: onePe = 1.0 + e;
            %                 % % e: J=2  KJ = 1 All 2nd derivatives vanish
            %                 % % 	// from line 12: oneMe = 1.0 - e;
            %                 % % e: J=2 KJ = 1 All 2nd derivatives vanish
            %                 % % 	// from line 13: fac   = onePe*oneMe;
            %                 % % onePe: J = 11 KJ = oneMe
            %                 % s(11) = s(11) + f(13)*q(12);
            %                 % % oneMe: J = 12 KJ = onePe
            %                 % s(12) = s(12) + f(13)*q(11);
            %                 % // from line 13: fac   = 1 - e^2;
            %                 % J = 2: KJ = -2.*e
            %                 s(2) = s(2) - f(13)*q(2)*2.0;
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 % 	// from line 14: if(fac > 0) rootfac = sqrt(fac);
            %                 % K = 14
            %                 if fac>0
            %                     % fac: J=13  KJ = 0.5/rootfac = 0.5*fac^(1/2), I = 13
            %                     s(13) = s(13) - f(14)*q(13)*0.25/(rootfac^3);
            %                     % try I = 14
            %                     %s(14) = s(14) - f(14)*q(14)*0.5/(rootfac^2);
            %                 end
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 % 	// from line 15: a = p/fac;
            %                 % fac: J = 13 KJ = -p*fac^(-2)
            %                 % fac: I = 13 & p I = 3
            %                 %s(13) = s(13) + f(15)*(q(13)*2.0*a/fac^2 - q(3)/fac^2);
            %                 % p: J = 3 KJ = fac(-1);
            %                 %s(3) = s(3) - f(15)*q(13)/fac^2;
            %                 % 	// from line 15: p = a*fac;
            %                 % fac: J = 13 KJ = a I = 3
            %                 s(13) = s(13) + f(15)*q( 3);
            %                 % a: J = 3 KJ = fac, I =13;
            %                 s( 3) = s( 3) + f(15)*q(13);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 % 	// from line 16: meanMotion = sqrt(mu)*a^(-1.5);
            %                 % a: J = 15 KJ = -1.5*a(-2.5)
            %                 % a I = 15
            %                 %s(15) = s(15) + f(16)*q(15)*3.75*meanMotion/(a^2);
            %                 % a: J = 3 KJ = -1.5*a(-2.5)
            %                 % a I = 3
            %                 s( 3) = s( 3) + f(16)*q( 3)*3.75*meanMotion/(a^2);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 %   // from line 17: cosI             = cos(this.Inclination);
            %                 % I: J = 4 KJ = -sinI;  (i = 18)
            %                 %s(4) = s(4) - f(17)*q(18);
            %                 % Alternate method
            %                 s(4)  = s(4) - f(17)*q(4)*cosI;
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 %   // from line 18: sinI             = sin(this.Inclination);
            %                 % I  J = 4: KJ = cosI (i = 17)
            %                 %s(4) = s(4) + f(18)*q(17);
            %                 % Alternate method
            %                 s(4)  = s(4) - f(18)*q(4)*sinI;
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 %   // from line 19: cosom            = cos(this.omega);
            %                 % omega J = 5: KJ = -sinom (i = 20)
            %                 %s(5) = s(5) - f(19)*q(20);
            %                 s(5) = s(5) - f(19)*q(5)*cosom;
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 %   // from line 20: sinom            = sin(this.omega);
            %                 % omega J = 5: KJ = cosom (i = 19)
            %                 %s(5) = s(5) + f(20)*q(19);
            %                 s(5) = s(5) - f(20)*q(5)*sinom;
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 %   // from line 21: cosO             = cos(this.Omega);
            %                 % Omega J = 6: KJ = -sinO (i = 22)
            %                 %s(6) = s(6) - f(21)*q(22);
            %                 s(6) = s(6) - f(21)*q(6)*cosO;
            %                 %   // from line 22: sinO             = sin(this.Omega);
            %                 % Omega J = 6: KJ = cosO (i = 21)
            %                 %s(6) = s(6) + f(22)*q(21);
            %                 s(6) = s(6) - f(22)*q(6)*sinO;
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 %   // from line K = 23: Px   =  cosO*cosom - sinO*sinom*cosI;
            %                 % K = 23 dependencies sinO 22, cosO 21, sinom 20, cosom 19, cosI 17
            %                 % sinO J = 22: KJ = -sinom*cosI (i = 20, 17)
            %                 s(22) = s(22) - f(23)*(q(20)*cosI + q(17)*sinom);
            %                 % cosO J = 21: KJ = cosom (i = 19)
            %                 s(21) = s(21) + f(23)*q(19);
            %                 % sinom J = 20: KJ = -sinO*cosI (i = 22, 17)
            %                 s(20) = s(20) - f(23)*(q(22)*cosI + q(17)*sinO);
            %                 % cosom J = 19: KJ = cosO (i=21)
            %                 s(19) = s(19) + f(23)*q(21);
            %                 % cosI J = 17: KJ = - sinO*sinom (i=22, 20)
            %                 s(17) = s(17) - f(23)*(q(22)*sinom + q(20)*sinO);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 %   // from line K = 24: Py   =  sinO*cosom + cosO*sinom*cosI;
            %                 % K = 24 dependencies sinO 22, cosO 21, sinom 20, cosom 19, cosI 17
            %                 % sinO J = 22: KJ = cosom (i = 19)
            %                 s(22) = s(22) + f(24)*q(19);
            %                 % cosO J = 21: KJ = sinom*cosI (i = 20, 17)
            %                 s(21) = s(21) + f(24)*(q(20)*cosI + q(17)*sinom);
            %                 % sinom J = 20: KJ = cosO*cosI (i = 21, 17)
            %                 s(20) = s(20) + f(24)*(q(21)*cosI + q(17)*cosO);
            %                 % cosom J = 19: KJ = sinO (i=22)
            %                 s(19) = s(19) + f(24)*q(22);
            %                 % cosI J = 17: KJ = cosI*sinom (i=21, 20)
            %                 s(17) = s(17) + f(24)*(q(21)*sinom + q(20)*cosO);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 %   // from line K = 25: Pz   =  sinom*sinI;
            %                 % K = 25 dependencies sinom 20, sinI 18
            %                 % sinom J = 20: KJ = sinI (i=18)
            %                 s(20) = s(20) + f(25)*q(18);
            %                 % sinI J = 18: KJ = sinom (i=20)
            %                 s(18) = s(18) + f(25)*q(20);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 %   // from line 26: Qx   = -cosO*sinom - sinO*cosom*cosI;
            %                 % K = 26 dependencies (sinO 22, cosO 21,  sinom 20, cosom 19, cosI 17
            %                 % sinO J = 22: KJ = -cosom*cosI (i = 19, 17)
            %                 s(22) = s(22) - f(26)*(q(19)*cosI + q(17)*cosom);
            %                 % cosO J = 21: KJ = -sinom (i=20)
            %                 s(21) = s(21) - f(26)*q(20);
            %                 % sinom J = 20: KJ = -cosO (i=21)
            %                 s(20) = s(20) - f(26)*q(21);
            %                 % cosom J = 19: KJ = -sinO*cosI (i=22, 17)
            %                 s(19) = s(19) - f(26)*(q(22)*cosI + q(17)*sinO);
            %                 % cosI J = 17: KJ = -sinO*cosom (i=22, 19)
            %                 s(17) = s(17) - f(26)*(q(22)*cosom + q(19)*sinO);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 %   // from line 27: Qy   = -sinO*sinom + cosO*cosom*cosI;
            %                 % K = 27 dependencies (sinO 22, cosO 21, sinom 20, cosom 19, cosI 17)
            %                 % sinO J = 22: KJ = -sinom
            %                 s(22) = s(22) - f(27)*q(20);
            %                 % cosO J = 21: KJ = cosom*cosI (i=19, 17)
            %                 s(21) = s(21) + f(27)*(q(19)*cosI + q(17)*cosom);
            %                 % sinom J = 20: KJ = -sinO (i=22);
            %                 s(20) = s(20) - f(27)*q(22);
            %                 % cosom J = 19: KJ = cosO*cosI (i=21, 17)
            %                 s(19) = s(19) + f(27)*(q(21)*cosI + q(17)*cosO);
            %                 % cosI J = 17: KJ = cosO*cosom (i=21, 19)
            %                 s(17) = s(17) + f(27)*(q(21)*cosom + q(19)*cosO);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 %   // from line 28: Qz   =  cosom*sinI;
            %                 % K = 28 dependencies (cosom 19, sinI 18)
            %                 % cosom J = 19: KJ = sinI (i=18)
            %                 s(19) = s(19) + f(28)*q(18);
            %                 % sinI J = 18: KJ = cosom (i=19)
            %                 s(18) = s(18) + f(28)*q(19);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 %   // from line 29: Wx   =  sinO*sinI;
            %                 % K = 29 dependencies (sinO 22, sinI 18)
            %                 % sinO J = 22: KJ = sinI (i=18)
            %                 s(22) = s(22) + f(29)*q(18);
            %                 % sinI J = 18: KJ = sinO (i=22)
            %                 s(18) = s(18) + f(29)*q(22);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 %   // from line 30: Wy   = -cosO*sinI;
            %                 % K = 30 dependencies (cosI 21, sinI 18)
            %                 % cosI J = 21: KJ = -sinI (i=18)
            %                 s(21) = s(21) - f(30)*q(18);
            %                 % sinI J = 18: KJ = -cosI (i=21)
            %                 s(18) = s(18) - f(30)*q(21);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 %   // from line 31: Wz   =  cosI;
            %                 % K = 31 dependencies cosI 17
            %                 % J = 17 KJ = -sinI (i=18)
            %                 %s(17) = s(17) - f(31)*q(18);
            %                 s(4) = s(4) - f(31)*q(4)*cosI;
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 % // from line 40: M = meanMotion*t + Mp;
            %                 % K = 40 dependencies (meanMotion 16, Mp 7, time 1)
            %                 % meanMotion J = 16: KJ = t
            %                 s(16) = s(16) + f(40)*q(1);
            %                 % Mp J = 7: KJ = 1 (2nd derivative vanishes)!
            %                 % t J = 1: KJ = meanMotion
            %                 s(1) = s(1) + f(40)*q(16);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 % Solve Kepler's Equation for E:  M = E - e*sin(E);
            %                 % from line 41: E = E(e,M)
            %                 % K = 41 dependencies (M 40, e 2)
            %                 % M J = 40 KJ = dE_dM (i = 40, 2)
            %                 s(40) = s(40) + f(41)*(q(2)*d2E_dMde + q(40)*d2E_dMdM);
            %                 % e J = 2 KJ = dE_de (i = 40,2)
            %                 s( 2) = s( 2) + f(41)*(q(2)*d2E_dede + q(40)*d2E_dedM);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 % from line 42: cosK  = cos(E);
            %                 % line K = 42 dependencies (E 41)
            %                 % // E J = 41: KJ = -sinK (i = 43)
            %                 %s(41) = s(41) - f(42)*q(43);
            %                 s(41) = s(41) - f(42)*q(41)*cosK;
            %                 % from line 43: sinK  = sin(E);
            %                 % line K = 43 dependencies (E 41);
            %                 % E J = 41: KJ = cosK (i = 42)
            %                 %s(41) = s(41) + f(43)*q(42);
            %                 s(41) = s(41) - f(43)*q(41)*sinK;
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 % from line 44: tanX  = cosK - e;
            %                 % line K = 44 dependencies (cosK 42, e 2)
            %                 % cosK J = 42: KJ = 1                              (2nd derivatives vanish)
            %                 % e    J =  2: KJ =-1                              (2nd derivatives vanish)
            %                 % // from line 45: tanY  = rootfac*sinK;
            %                 % line K = 45 dependencies (sinK 43, rootfac 14)
            %                 % sinK J = 43: KJ = rootfac (i = 14)
            %                 s(43) = s(43) + f(45)*q(14);
            %                 % rootfac J = 14: KJ = sinK (i = 43)
            %                 s(14) = s(14) + f(45)*q(43);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 % //from line 46 double nu;// nu = atan2(Ytan, Xtan);                 // line 77
            %                 % // from line 46: nu     = atan2(tanY, tanX);
            %                 cosDenom = (tanX*tanX + tanY*tanY)^2;
            %                 % line K = 46 dependencies (tanY 45, tanX 44)
            %                 % %q(46) = q(46) + (q(45)*tanX - q(44)*tanY)/(tanX*tanX + tanY*tanY);
            %                 % tanX J = 44 KJ = -tanY/(tanX*tanX + tanY*tanY);
            %                 m00 = 2.0*tanX*tanY/cosDenom;
            %                 m11 = (tanY*tanY - tanX*tanX)/cosDenom;
            %                 s(44) = s(44) + f(46)*(q(44)*m00 + q(45)*m11);
            %                 % tanY J = 45: KJ = tanX/(tanX*tanX + tanY*tanY);
            %                 s(45) = s(45) + f(46)*(q(44)*m11 - q(45)*m00);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 % // from line 47: coss  = cos(nu);
            %                 % line K = 47 dependencies (nu 46)
            %                 % nu J = 46: KJ = -sin(nu) = -sins (i = 48);
            %                 %s(46) = s(46) - f(47)*q(48);
            %                 s(46) = s(46) - f(47)*q(46)*coss;
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 % // from line 48: sins  = sin(nu);
            %                 % line K = 48 dependencies (nu 46)
            %                 % nu J = 46: KJ = cos(nu) = coss (i= 47)
            %                 % s(46) = s(46) + f(48)*q(47);
            %                 s(46) = s(46) - f(48)*q(46)*sins;
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 % // from line 49: rorbit   = a*(1.0 - e*cosK);
            %                 % cosK: J=42 KJ = -a*e
            %                 m03 = -e;
            %                 m02 = -a;
            %                 s(42) = s(42) + f(49)*(q( 3)*m03 + q( 2)*m02);
            %                 % E: J = 41 KJ = +a*e*sin(E)
            %                 % m15 = -e;
            %                 % m02 = -a;
            %                 % s(42) = s(42) + f(49)*(q(15)*m15 + q( 2)*m02);
            %                 % a: J = 3  KJ = 1 - e*cosK
            %                 m42 = -e;
            %                 m02 = -cosK;
            %                 s( 3) = s( 3) + f(49)*(q(42)*m42 + q( 2)*m02);
            %                 % e: J = 2 KJ = -a*cosK
            %                 m42 = -a;
            %                 m03 = -cosK;
            %                 s( 2) = s( 2) + f(49)*(q(42)*m42 + q(3)*m03);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 % // from line 50: rorbity  = rorbit*(coss*Px + sins*Qx);
            %                 % K = 50 Dependencies rorbit 49, sins 48, coss 47, Qx 26, Px 23
            %                 % rorbit: J=49 KJ = coss*Px + sins*Qx (i = 48, 47, 26, 23)
            %                 s(49) = s(49) + f(50)*(q(48)*Qx + q(47)*Px + q(26)*sins + q(23)*coss);
            %                 % sins J=48 KJ = rorbit*Qx  (i = 49, 26)
            %                 s(48) = s(48) + f(50)*(q(49)*Qx + q(26)*rorbit);
            %                 % coss J=47 KJ = rorbit*Px (i = 49, 23)
            %                 s(47) = s(47) + f(50)*(q(49)*Px + q(23)*rorbit);
            %                 % Qx J=26 KJ = rorbit*sins (i = 49, 48)
            %                 s(26) = s(26) + f(50)*(q(49)*sins + q(48)*rorbit);
            %                 % Px J=23 KJ = rorbit*coss (i = 49, 47)
            %                 s(23) = s(23) + f(50)*(q(49)*coss + q(47)*rorbit);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 % // from line 51: rorbity  = rorbit*(coss*Py + sins*Qy);
            %                 % rorbit: J=49 KJ = coss*Py + sins*Qy
            %                 m48 = Qy;
            %                 m47 = Py;
            %                 m27 = sins;
            %                 m24 = coss;
            %                 s(49) = s(49) + f(51)*(q(48)*m48 + q(47)*m47 + q(27)*m27 + q(24)*m24);
            %                 % sins 48 KJ = rorbit*Qy
            %                 m49 = Qy;
            %                 m27 = rorbit;
            %                 s(48) = s(48) + f(51)*(q(49)*m49 + q(27)*m27);
            %                 % sins 47 KJ = rorbit*Py
            %                 m49 = Py;
            %                 m24 = rorbit;
            %                 s(47) = s(47) + f(51)*(q(49)*m49 + q(24)*m24);
            %                 % sins 27 KJ = rorbit*sins
            %                 m49 = sins;
            %                 m48 = rorbit;
            %                 s(27) = s(27) + f(51)*(q(49)*m49 + q(48)*m48);
            %                 % sins 24 KJ = rorbit*coss
            %                 m49 = coss;
            %                 m47 = rorbit;
            %                 s(24) = s(24) + f(51)*(q(49)*m49 + q(47)*m47);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 % // from line 52: rorbitz  = rorbit*(coss*Pz + sins*Qz);
            %                 % rorbit: 49 KJ = coss*Pz + sins*Qz
            %                 m48 = Qz;
            %                 m47 = Pz;
            %                 m28 = sins;
            %                 m25 = coss;
            %                 s(49) = s(49) + f(52)*(q(48)*m48 + q(47)*m47 + q(28)*m28 + q(25)*m25);
            %                 % sins 48 KJ = rorbit*Qz
            %                 m49 = Qz;
            %                 m28 = rorbit;
            %                 s(48) = s(48) + f(52)*(q(49)*m49 + q(28)*m28);
            %                 % sins 47 KJ = rorbit*Pz
            %                 m49 = Pz;
            %                 m25 = rorbit;
            %                 s(47) = s(47) + f(52)*(q(49)*m49 + q(25)*m25);
            %                 % sins 28 KJ = rorbit*sins
            %                 m49 = sins;
            %                 m48 = rorbit;
            %                 s(28) = s(28) + f(52)*(q(49)*m49 + q(48)*m48);
            %                 % sins 25 KJ = rorbit*coss
            %                 m49 = coss;
            %                 m47 = rorbit;
            %                 s(25) = s(25) + f(52)*(q(49)*m49 + q(47)*m47);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 % from line 53: rtpinv   = rtpinv;
            %                 % p: J = 3 KJ = -0.5*rtpinv/p
            %                 %s(3)  = s(3) + f(53)*q(3)*0.75*rtpinv/psq;
            %                 % p: J = 15 KJ = -0.5*rtpinv/p
            %                 s(15)  = s(15) + f(53)*q(15)*0.75*rtpinv/psq;
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 %from line 54:  double vorbitx  = rtpinv*(-sins*Px + (e + coss)*Qx);  // line 54
            %                 % K = 54 Dependencies: rtpinv = 53, sins 48, coss 47, Qx 26, Px 23, e 2
            %                 %rtpinv: J = 53 KJ = -sins*Px + (e + coss)*Qx (i = 48, 47, 26, 23, 2)
            %                 s(53) = s(53) + f(54)*(-q(48)*Px + q(47)*Qx + q(26)*(e + coss) -q(23)*sins + q(2)*Qx);
            %                 % m48   = -Px;
            %                 % m47   =  Qx;
            %                 % m26   =  e + coss;
            %                 % m23   = -sins;
            %                 % m02   =  Qx;
            %                 % s(53) =  s(53) + f(54)*(q(48)*m48 + q(47)*m47 + q(26)*m26 + q(23)*m23 + q( 2)*m02);
            %                 % // sins: J = 48   KJ = -rtpinv*Px (i = 49, 23)
            %                 s(48) = s(48) - f(54)*(q(53)*Px + q(23)*rtpinv);
            %                 % m53   = -Px;
            %                 % m23   = -rtpinv;
            %                 % s(48) = s(48) + f(54)*(q(53)*m53 + q(23)*m23);
            %                 % // coss: J = 47   KJ = rtpinv*Qx (i = 49, 26)
            %                 s(47) = s(47) + f(54)*(q(53)*Qx + q(26)*rtpinv);
            %                 % m53   = Qx;
            %                 % m26   = rtpinv;
            %                 % s(47) = s(47) + f(54)*(q(53)*m53 + q(26)*m26);
            %                 %// Qx: J = 28   KJ = rtpinv*(e+coss) (i = 53, 47, 2)
            %                 s(26) = s(26) + f(54)*(q(53)*(e+coss) + q(47)*rtpinv + q( 2)*rtpinv);
            %                 % m53   =  e + coss;
            %                 % m47   =  rtpinv;
            %                 % m02   =  rtpinv;
            %                 % s(26) = s(26) + f(54)*(q(53)*m53 + q(47)*m47 + q( 2)*m02);
            %                 %// Px: J = 23  KJ = -rtpinv*sins (i = 53, 48)
            %                 s(23) = s(23) - f(54)*(q(53)*sins + q(48)*rtpinv);
            %                 % m53   = -sins;
            %                 % m48   = -rtpinv;
            %                 % s(24) = s(24) + f(54)*(q(53)*m53 + q(48)*m48);
            %                 %// e: J = 2  KJ = rtpinv*Qx (i = 53, 26)
            %                 s( 2) = s( 2) + f(54)*(q(53)*Qx + q(26)*rtpinv);
            %                 % m53   = Qx;
            %                 % m26   = rtpinv;
            %                 % s( 2) = s( 2) + f(54)*(q(53)*m53 + q(26)*m26);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 %// from line 55:  double vorbity  = rtpinv*(-sins*Py + (e + coss)*Qy);  // line 55
            %                 %//rtpinv: J = 53 KJ = -sins*Py + (e + coss)*Qy
            %                 m48   = -Py;
            %                 m47   =  Qy;
            %                 m27   =  e + coss;
            %                 m24   = -sins;
            %                 m02   =  Qy;
            %                 s(53) =  s(53) + f(55)*(q(48)*m48 + q(47)*m47 + q(27)*m27 + q(24)*m24 + q( 2)*m02);
            %                 % // sins: J = 48   KJ = -rtpinv*Py
            %                 m53   = -Py;
            %                 m24   = -rtpinv;
            %                 s(48) = s(48) + f(55)*(q(53)*m53 + q(24)*m24);
            %                 % // coss: J = 47   KJ = rtpinv*Qy
            %                 m53   = Qy;
            %                 m27   = rtpinv;
            %                 s(47) = s(47) + f(55)*(q(53)*m53 + q(27)*m27);
            %                 %// Qy: J = 27   KJ = rtpinv*(e+coss)
            %                 m53   =  e + coss;
            %                 m47   =  rtpinv;
            %                 m02   =  rtpinv;
            %                 s(27) = s(27) + f(55)*(q(53)*m53 + q(47)*m47 + q( 2)*m02);
            %                 %// Py: J = 24  KJ = -rtpinv*sins
            %                 m53   = -sins;
            %                 m48   = -rtpinv;
            %                 s(24) = s(24) + f(55)*(q(53)*m53 + q(48)*m48);
            %                 %// e: J = 2  KJ = rtpinv*Qy
            %                 m53   = Qy;
            %                 m27   = rtpinv;
            %                 s( 2) = s( 2) + f(55)*(q(53)*m53 + q(27)*m27);
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 %// from line 56:  double vorbitz  = rtpinv*(-sins*Pz + (e + coss)*Qz);  // line 56
            %                 %//rtpinv: J = 53 KJ = -sins*Pz + (e + coss)*Qz
            %                 m48   = -Pz;
            %                 m47   =  Qz;
            %                 m28   =  e + coss;
            %                 m25   = -sins;
            %                 m02   =  Qz;
            %                 s(53) =  s(53) + f(56)*(q(48)*m48 + q(47)*m47 + q(28)*m28 + q(25)*m25 + q( 2)*m02);
            %                 % // sins: J = 48   KJ = -rtpinv*Pz
            %                 m53   = -Pz;
            %                 m25   = -rtpinv;
            %                 s(48) = s(48) + f(56)*(q(53)*m53 + q(25)*m25);
            %                 % // coss: J = 47   KJ = rtpinv*Qz
            %                 m53   = Qz;
            %                 m28   = rtpinv;
            %                 s(47) = s(47) + f(56)*(q(53)*m53 + q(28)*m28);
            %                 %// Qz: J = 28   KJ = rtpinv*(e+coss)
            %                 m53   =  e + coss;
            %                 m47   =  rtpinv;
            %                 m02   =  rtpinv;
            %                 s(28) = s(28) + f(56)*(q(53)*m53 + q(47)*m47 + q( 2)*m02);
            %                 %// Pz: J = 25  KJ = -rtpinv*sins
            %                 m53   = -sins;
            %                 m48   = -rtpinv;
            %                 s(25) = s(25) + f(56)*(q(53)*m53 + q(48)*m48);
            %                 %// e: J = 2  KJ = rtpinv*Qz
            %                 m53   = Qz;
            %                 m28   = rtpinv;
            %                 s( 2) = s( 2) + f(56)*(q(53)*m53 + q(28)*m28);
            %                 if LineCheck > 60
            %                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                     % Escobal Potential starts here
            %                     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                     % from line 61: X    = rorbitx = Kvec(1);                     % line  61
            %                     % All 2nd Derivatives vanish
            %                     % from line 62: Y    = rorbity = Kvec(2);                     % line  62
            %                     % All 2nd Derivatives vanish
            %                     % from line 63: Z    = rorbitz = Kvec(3);                     % line  63
            %                     % All 2nd Derivatives vanish
            %                     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                     % from line 64: R2   = X^2 + Y^2 + Z^2;                         % line  64
            %                     % K = 64
            %                     %// X: J = 61  KJ = 2.0*X
            %                     s(61) = s(61) + 2.0*f(64)*q(61);
            %                     %// Y: J = 62  KJ = 2.0*Y
            %                     s(62) = s(62) + 2.0*f(64)*q(62);
            %                     %// Z: J = 63  KJ = 2.0*Z
            %                     s(63) = s(63) + 2.0*f(64)*q(63);
            %                     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                     % from line 65: R1   = rorbit;                       % line  64
            %                     % All 2nd Derivatives vanish
            %                     % from line 65: R1   = sqrt(R2);                                % line  65
            %                     %// R2: J = 64  KJ = 0.5/R1   check this
            %                     s(64) = s(64) - f(65)*q(64)*0.25/R3;
            %                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                     % from line 64: R2   = R1^2;                                % line  65
            %                     %// R1: J = 65  KJ = 2.0*R1   I = 65
            %                     %s(65) = s(65) + 2.0*f(64)*q(65);
            %
            %                     % from line 66: R3   = R1*R2;
            %                     %// R1: J = 65 KJ = R2  I = 64
            %                     s(65) = s(65) + f(66)*q(64);
            %                     %// R2: J = 64 KJ = R1 I = 65
            %                     s(64) = s(64) + f(66)*q(65);
            %
            %                     % from line 67: R4   = R1*R3;                                   % line  67
            %                     %// R1: J = 65 KJ = R3  I = 66
            %                     s(65) = s(65) + f(67)*q(66);
            %                     %// R3: J = 66 KJ = R1
            %                     s(66) = s(66) + f(67)*q(65);
            %
            %                     % from line 68: R5   = R1*R4;                                   % line  68
            %                     % K = 68
            %                     %// R1: J = 65  KJ = R4 I = 67
            %                     s(65) = s(65) + f(68)*q(67);
            %                     %// R4: J = 67  KJ = RI I = 65
            %                     s(67) = s(67) + f(68)*q(65);
            %
            %                     % from line 69: R6   = R1*R5;                                   % line  69
            %                     %// R1: J = 65  KJ = R5 I = 68
            %                     s(65) = s(65) + f(69)*q(68);
            %                     %// R5: J = 68  KJ = R1 I = 65
            %                     s(68) = s(68) + f(69)*q(65);
            %                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                     % from line 70: sd1  = rorbitz/R1;
            %                     %//rorbitz: J = 52 KJ = 1/R1 I = 65
            %                     s(52) = s(52) - f(70)*q(65)/R2;
            %                     %// R1: J = 65 KJ = -Z*R1^(-2) I = 52, 65
            %                     s(65) = s(65) + f(70)*(2.0*q(65)*sd1 - q(52))/R2;
            %
            %                     % from line 70: sd1  = Z/R1;
            %                     %// Z: J = 63 KJ = 1/R1 I = 65
            %                     %s(63) = s(63) - f(70)*q(65)/R2;
            %                     %// R1: J = 65 KJ = -Z*R1^(-2) I = 63, 65
            %                     %s(65) = s(65) + f(70)*(2.0*q(65)*sd1 - q(63))/R2;
            %                     % from line 70: sd1  = sin(I)*sin(omega + nu);
            %                     %q(70)  = q(70) + q(18)*sinomnu + sinI*cosomnu*(q(5) + q(46));
            %                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                     % from line 71: sd2  = sd1^2;
            %                     %//sd1: J = 70 KJ = 2.0*sd1 I = 70
            %                     s(70) = s(70) + f(71)*q(70)*2.0;
            %
            %                     % from line 72: sd3  = sd1*sd2;
            %                     %//sd1: J = 70 KJ = sd2 I = 71
            %                     s(70) = s(70) + f(72)*q(71);
            %                     %//sd2: J = 71 KJ = sd1 I = 70
            %                     s(71) = s(71) + f(72)*q(70);
            %
            %                     % from line 73: sd4  = sd1*sd3;
            %                     %//sd1: J = 70 KJ = sd3 I = 72
            %                     s(70) = s(70) + f(73)*q(72);
            %                     %//sd3: J = 72 KJ = sd1 I = 70
            %                     s(72) = s(72) + f(73)*q(70);
            %
            %                     % from line 74: sd5  = sd1*sd4;
            %                     %//sd1: J = 70 KJ = sd4 I = 73
            %                     s(70) = s(70) + f(74)*q(73);
            %                     %//sd4: J = 73 KJ = sd1 I = 70
            %                     s(73) = s(73) + f(74)*q(70);
            %
            %                     % from line 75: sd6  = sd1*sd5;
            %                     %//sd1: J = 70 KJ = sd5 I = 74
            %                     s(70) = s(70) + f(75)*q(74);
            %                     %//sd5: J = 73 KJ = sd1 I = 70
            %                     s(74) = s(74) + f(75)*q(70);
            %
            %                     % All Second Derivatives Vanish here
            %                     % from line 76: F2   =  1.0      -  3.0*sd2;
            %                     % from line 77: F3   =  3.0*sd1  -  5.0*sd3;
            %                     % from line 78: F4   =  3.0      -  30.*sd2  + 35*sd4;
            %                     % from line 79: F5   =  15.*sd1  -  70.*sd3  + 63*sd5;
            %                     % from line 80: F6   =  5.       - 105.*sd2  + 315*sd4 - 231*sd6;
            %
            %                     % from line 81: V1    = 1/R1;          Kepler Potential                 % line 21
            %                     %//R1: J = 65 KJ = -1/R1^2 I = 65
            %                     s(65) = s(65) + 2.0*f(81)*q(65)/R3;
            %
            %                     % from line 82: V2      = C2*F2/R2;
            %                     %//F2 J = 76 KJ = C2/R2 I = 64
            %                     s(76) = s(76) - f(82)*q(64)*C2/R4;
            %                     %//R2 J = 64 KJ = -C2*F2/R2^2 I = 76, 67
            %                     s(64) = s(64) + f(82)*(2.0*q(64)*V2 - C2*q(76))/R4;
            %
            %                     % from line 83: V3      = C3*F3/R3;
            %                     %//F3 J = 77 KJ = C3/R3 I = 66
            %                     s(77) = s(77) - f(83)*q(66)*C3/R6;
            %                     %//R3 J = 66 KJ = -C3*F3/R3^2 I = 77, 66
            %                     s(66) = s(66) + f(83)*(2.0*q(66)*V3 - C3*q(77))/R6;
            %
            %                     % from line 84: V4      = C4*F4/R4;
            %                     %//F4 J = 78 KJ = C4/R4 I = 67
            %                     s(78) = s(78) - f(84)*q(67)*C4/R4^2;
            %                     %//R4 J = 67 KJ = -C4*F4/R4^-2 I = 78, 67
            %                     s(67) = s(67) + f(84)*(2.0*q(67)*V4 - C4*q(78))/R4^2;
            %
            %                     % from line 85: V5      = C5*F5/R5;
            %                     %//F5 J = 79 KJ = C5/R5 I = 68
            %                     s(79) = s(79) - f(85)*q(68)*C5/R5^2;
            %                     %//R5 J = 68 KJ = -C5*F5/R5^-2 I = 79, 68
            %                     s(68) = s(68) + f(85)*(2.0*q(68)*V5 - C5*q(79))/R5^2;
            %
            %                     % from line 86: V6      = C6*F6/R6;
            %                     %//F6 J = 80 KJ = C5/R6 I = 69
            %                     s(80) = s(80) - f(86)*q(69)*C6/R6^2;
            %                     %//R6 J = 69 KJ = -C6*F6/R6^-2 I = 80, 69
            %                     s(69) = s(69) + f(86)*(2.0*q(69)*V6 - C6*q(80))/R6^2;
            %
            %                     % All Second Derivative Vanish
            %                     % from line 87: V       = mu*V1;
            %                     % from line 88: D       = V2 + V3 + V4 + V5 + V6;
            %
            %                     % %from line 89: R       = mu*D;
            %                     % %q(89) = q(89) + q(88)*mu;
            %                     % % from line 90: Phi     = V + R;
            %                     % %q(90) = q(90) + q(87) + q(89);
            %
            %                     % from line 90: Phi     = V*(CK + D);
            %                     %//V J = 87  KJ = (CK + D)  I = 88
            %                     s(87) = s(87) + f(90)*q(88);
            %                     %//D J = 88  KJ = V  I = 87
            %                     s(88) = s(88) + f(90)*q(87);
            %                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                     if LineCheck > 90
            %                         fac1 = this.fac1;
            %                         fac2 = this.fac2;
            %                         fac3 = this.fac3;
            %                         fac4 = this.fac4;
            %                         fac5 = this.fac5;
            %                         Mwe  = this.Mwe;
            %                         MeM  = this.MeM;
            %                         MaM  = this.MaM;
            %                         MWi  = this.MWi;
            %                         Miw  = this.Miw;
            %                         MLagrange = this.MLagrange;
            %                         % from line 91: fac1 = n*a
            %                         % K = 91
            %                         %//n J = 16, KJ = a I = 3
            %                         s(16) = s(16) + f(91)*q( 3);
            %                         %//a J = 3, KJ = n I = 16
            %                         s( 3) = s( 3) + f(91)*q(16);
            %
            %                         % from line 92: fac2 = fac1*a % line 92
            %                         % K = 92
            %                         %//a J = 3, KJ = fac1, I = 91)
            %                         s( 3) = s( 3) + f(92)*q(91);
            %                         %// fac1 J = 91 KJ = a I = 3
            %                         s(91) = s(91) + f(92)*q( 3);
            %
            %                         % from line 93: fac3 = fac2*e % line 93
            %                         % K = 93
            %                         %//e J = 2, KJ = fac2 I = 92
            %                         s( 2) = s( 2) + f(93)*q(92);
            %                         %//fac2 J = 92, KJ = e I = 2
            %                         s(92) = s(92) + f(93)*q( 2);
            %
            %                         % from line 94: fac4 = fac2*rootfac  % line 94
            %                         % K = 94
            %                         %//rootfac J = 14, KJ = fac2 I = 92
            %                         s(14) = s(14) + f(94)*q(92);
            %                         %//fac2 J = 92, KJ = rootfac I = 14
            %                         s(92) = s(92) + f(94)*q(14);
            %
            %                         % from line 95: fac5 = fac4*sinI        % line 95
            %                         % K = 95
            %                         %//I J = 4, KJ = fac4*cosI, I = 4, 94
            %                         s( 4) = s( 4) + f(95)*(q(94)*cosI - q( 4)*fac4*sinI);
            %                         %// fac4 J = 94, KJ = sinI, I = 4
            %                         s(94) = s(94) + f(95)*q( 4)*cosI;
            %
            %                         %from line 96: Mwe  = rootfac/fac3      % line 96
            %                         %K = 96
            %                         %//rootfac J = 14, KJ = 1/fac3 I = 93
            %                         s(14) = s(14) - f(96)*q(93)/fac3^2;
            %                         %//fac3 J = 93, KJ = -rootfac/fac3^2, I = 14, 93
            %                         s(93) = s(93)  + f(96)*(2.0*q(93)*Mwe - q(14))/fac3^2;
            %
            %                         % from line 97: MeM  = fac/fac3    % line 97
            %                         %// K = 97
            %                         %//fac J = 13, KJ = 1/fac3
            %                         s(13) = s(13) - f(97)*q(93)/fac3^2;
            %                         %//fac3 J = 93, KJ = -fac/fac3^2 I = 13, 93
            %                         s(93) = s(93) + f(97)*(2.0*q(93)*MeM - q(13))/fac3^2;
            %
            %                         % from line 98: MaM  = 2.0/fac1         % line 98
            %                         %//fac1 J = 91, KJ = -2.0/fac1^2 I = 91
            %                         s(91) = s(91) + f(98)*2.0*q(91)*MaM/fac1^2;
            %
            %                         % from line 99: MWi  = 1/fac5       % line 99
            %                         %//fac5 J = 95, KJ = -1/fac5^2 I = 95
            %                         s(95) = s(95) + 2.0*f(99)*q(95)/fac5^3;
            %
            %                         % from line 100: Miw  = cosI*MWi    % line 100
            %                         %//I J = 4, KJ = -sinI*MWi I = 4, 99
            %                         s( 4) = s( 4) - f(100)*(cosI*MWi*q( 4) + sinI*q(99));
            %                         %//MWi J = 99, KJ = cosI, I = 4
            %                         s(99) = s(99) - f(100)*sinI*q( 4);
            %                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                         % All Second Derivative Vanish
            %                         if LineCheck > 100
            %                             %los       = Rextrap - state_vector_sat(1:3);    % lines 101-103
            %                             %lx = los(1);                                    % line 101
            %                             %ly = los(2);                                    % line 102
            %                             %lz = los(3);                                    % line 103
            %                             % from line K = 85: V5      = C5*F5/R5;
            %                             %//F5 J = 79 KJ = C5/R5 I = 68
            %                             %s(J) = s(J) - f(K)*q(I)*C5/R5^2;
            %
            %                             % from line 104: lperpSq  = lx^2 + ly^2;         % line 104
            %                             % K = 104:
            %                             %//lx J = 101, KJ = 2.0*lx, I = 101
            %                             s(101)  = s(101) + 2.0*f(104)*q(101);
            %                             %//ly J = 102, KJ = 2.0*ly, I = 102
            %                             s(102)  = s(102) + 2.0*f(104)*q(102);
            %
            %                             % from line: 105: rangeSq  = lperpSq + lz^2;     % line 105
            %                             % K = 105
            %                             %//lz J = 103, KJ = 2.0*lz, I = 103
            %                             s(103)  = s(103) + 2.0*f(105)*q(103);
            %                             %//lperpSq J = 104, KJ = 1
            %                             % Second Derivative vanishes
            %
            %                             % from line 106: range    = sqrt(rangeSq);       % line 106
            %                             % K = 106
            %                             %// rangeSq J = 105, KJ = 0.5/range, I = 106
            %                             s(105)    = s(105) - 0.25*f(106)*q(105)/range^3;
            %
            %                             % from line 107: zLos     = lz/range;
            %                             % K = 107
            %                             %//lz J = 103, KJ = 1/range I = 106
            %                             s(103) = s(103) - f(107)*q(106)/rangeSq;
            %                             %//range J = 106, KJ = -lz/range^2 I = 103, 106
            %                             s(106) = s(106) + f(107)*(2.0*q(106)*zLos - q(103))/rangeSq;
            %
            %                             % from line 108: thetaLos = acos(zLos);          % line 108
            %                             % K = 108
            %                             %//zLos J = 107, KJ = -1/sqrt(1-zLos^2) I = 107
            %                             denom = sqrt(1.0-zLos^2);
            %                             s(107) = s(107) - f(108)*q(107)*zLos/denom^3;
            %
            %                             %phiLos   = atan2(ly,lx);                        % line 109
            %                             % K = 109
            %                             %//lx J = 101, KJ = -ly/(lx^2 + ly^2), I = 101, 102
            %                             s(101) = s(101) + f(109)*(2.0*q(101)*lx*ly + q(102)*(ly^2 - lx^2))/lperpSq^2;
            %                             %//ly J = 102, KJ = lx/(lx^2 + ly^2), I = 101, 102
            %                             s(102) = s(102) + f(109)*(q(101)*(ly^2 - lx^2) - 2.0*q(102)*lx*ly)/lperpSq^2;
            %
            %                             % All Second Derivative Vanish
            %                             %vlos     = Vextrap - state_vector_sat(4:6);       % lines 110-112
            %                             %vx       = vlos(1);                               % line 110
            %                             %vy       = vlos(2);                               % line 111
            %                             %vz       = vlos(3);                               % line 112
            %
            %                             % from line 113: ldotv    = lx*vx + ly*vy + lz*vz; % line 113
            %                             % K = 113
            %                             %//lx J = 101, KJ = vx I = 110
            %                             s(101)    = s(101) + f(113)*q(110);
            %                             %//vx J = 110, KJ = lx, I = 101
            %                             s(110)    = s(110) + f(113)*q(101);
            %
            %                             %//ly J = 102, KJ = vy I = 111
            %                             s(102)    = s(102) + f(113)*q(111);
            %                             %//vy J = 111, KJ = ly, I = 102
            %                             s(111)    = s(111) + f(113)*q(102);
            %
            %                             %//lz J = 103, KJ = vz I = 112
            %                             s(103)    = s(103) + f(113)*q(112);
            %                             %//vz J = 112, KJ = lz, I = 102
            %                             s(112)    = s(112) + f(113)*q(103);
            %
            %                             % from line 114: ldot     = ldotv/range;         % line 114
            %                             % K = 114
            %                             %//ldotv J = 113, KJ = 1/range I = 106
            %                             s(113) = s(113) - f(114)*q(106)/rangeSq;
            %                             %//range J = 106, KJ = -ldotv/rangeSq I = 106, 113
            %                             s(106) = s(106) + f(114)*(2.0*q(106)*ldot - q(113))/rangeSq;
            %
            %                             % from line 115: lperp    = sqrt(lperpSq);       % line 115
            %                             % K = 115
            %                             %//lperpSq J = 104, KJ = 0.5/sqrt(lperpSq) I = 104
            %                             s(104) = s(104) - 0.25*f(115)*q(104)/(lperpSq*lperp);
            %
            %                             % from line 116: thdotfac1 = range*vz - lz*ldot; % line 116
            %                             % K = 116
            %                             %//range J = 106, KJ = vz, I = 103
            %                             s(106) = s(106) + f(116)*q(112);
            %                             %//vz J = 112, KJ = range, I = 106
            %                             s(112) = s(112) + f(116)*q(106);
            %
            %                             %//lz J = 103, KJ = -ldot, I = 114
            %                             s(103) = s(103) - f(116)*q(114);
            %                             %//ldot J = 114, KJ = -vz, I = 103
            %                             s(114) = s(114) - f(116)*q(103);
            %
            %                             % from line 117: thdotfac2 = thdotfac1/range;    % line 117
            %                             %//thdotfac1 J = 116, KJ = 1/range I = 106
            %                             s(116) = s(116) - f(117)*q(106)/rangeSq;
            %                             %//range J = 106, KJ = -thdotfac1/rangeSq I = 106, 116
            %                             s(106) = s(106) + f(117)*(2.0*q(106)*thdotfac2 - q(116))/rangeSq;
            %
            %                             % from line 118: thdot     = -thdotfac2/lperp;   % line 118
            %                             %//thdotfac2 J = 117, KJ = -1/lperp I = 115
            %                             s(117) = s(117) + f(118)*q(115)/lperpSq;
            %                             %//lperp J = 115, KJ = thdotfac2/lperpSq I = 115, 117
            %                             s(115) = s(115) + f(118)*(q(117) + 2.0*q(115)*thdot)/lperpSq;
            %
            %                             % from line 119: phiDotNum = lx*vy - ly*vx;      % line 119
            %                             %K = 119
            %                             %//lx J = 101, KJ = vy I = 111
            %                             s(101) = s(101) + f(119)*q(111);
            %                             %//vy J = 111, KJ = lx I = 101
            %                             s(111) = s(111) + f(119)*q(101);
            %
            %                             %//ly J = 102, KJ = -vx I = 110
            %                             s(102) = s(102) - f(119)*q(110);
            %                             %//vx J = 110, KJ = -ly I = 102
            %                             s(110) = s(110) - f(119)*q(102);
            %
            %                             %phidot    = phiDotNum/lperpSq;                  % line 120
            %                             %K = 120
            %                             %//phiDotNum J = 119, KJ = 1/lperpSq I = 104
            %                             s(119) = s(119) - f(120)*q(104)/lperpSq^2;
            %                             %//lperpSq J = 104, KJ = -phiDotNum/lperpSq^2 I = 104, 119
            %                             s(104) = s(104) + f(120)*(2.0*q(104)*phidot - q(119))/lperpSq^2;
            %                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                             %cosThetaLos = cos(thetaLos);                        % line 121
            %                             % K = 121
            %                             %//thetaLos J = 108, KJ = -sin(thetaLos), I = 108
            %                             s(108) = s(108) - f(121)*q(108)*cosThetaLos;
            %                             %sinThetaLos = sin(thetaLos);                        % line 122
            %                             % K = 122
            %                             %//thetaLos J = 108, KJ = cos(thetaLos), I = 108
            %                             s(108) = s(108) - f(122)*q(108)*sinThetaLos;
            %                             %phiDiff     = phiLos - phiM;                        % line 123
            %                             % K = 123
            %                             %//phiLos J = 109, KJ = 1 -- All Second Derivatives vanish
            %                             %cosPhiDiff  = cos(phiDiff);                         % line 124
            %                             % K = 124
            %                             %//phiDiff J = 123, KJ = -sin(phiDiff) I = 123
            %                             s(123) = s(123) - f(124)*q(123)*cosPhiDiff;
            %                             %sinPhiDiff  = sin(phiDiff);                         % line 125
            %                             % K = 125
            %                             %//phiDiff J = 123, KJ = cos(phiDiff), I = 123
            %                             s(123) = s(123) - f(125)*q(123)*sinPhiDiff;
            %                             %delUV       = cosThetaLos*sinThetaM*cosPhiDiff;    % line 126
            %                             % K = 126
            %                             %cosThetaLos J = 121 KJ = sinThetaM*cosPhiDiff I = 124
            %                             s(121) = s(121) + f(126)*q(124)*sinThetaM;
            %                             %//cosPhiDiff J = 124 KJ = cosThetaLos*sinThetaM I = 121
            %                             s(124) = s(124) + f(126)*q(121)*sinThetaM;
            %                             %Uview       = sinThetaLos*cosThetaM - delUV;        % line 127
            %                             % K = 127
            %                             %//delUV K = 126, KJ = 1 -- All Second Derivatives vanish
            %                             %//sinThetaLos J = 122 KJ = cosThetaM -- All second derivatives vanish
            %                             %Vview       = sinThetaM*sinPhiDiff;                 % line 128
            %                             % K = 128
            %                             %//sinPhiDiff J = 123, KJ = sinThetaM -- All second derivatives vanish
            %
            %                             % from line 129: xLos     = lx/range;
            %                             % K = 129
            %                             %//lx J = 101, KJ = 1/range I = 106
            %                             s(101) = s(101) - f(129)*q(106)/rangeSq;
            %                             %//range J = 106, KJ = -lz/range^2 I = 101, 106
            %                             s(106) = s(106) + f(129)*(2.0*q(106)*xLos - q(101))/rangeSq;
            %
            %                             % from line 130: yLos     = ly/range;
            %                             % K = 130
            %                             %//ly J = 102, KJ = 1/range I = 106
            %                             s(102) = s(102) - f(130)*q(106)/rangeSq;
            %                             %//range J = 106, KJ = -ly/range^2 I = 102, 106
            %                             s(106) = s(106) + f(130)*(2.0*q(106)*yLos - q(102))/rangeSq;
            %
            %                             % from line 131: Collinearity = 1 - xLos*xLosM - yLos*yLosM - zLos*zLosM;
            %                             %  all second derivatives vanish
            %                             % from line 133: ChiSqLos = losDiff'*InvCovM*losDiff;
            %                             % K = 133
            %                             % xLos  J = 129 KJ = xUnit'*InvCovM*losDiff
            %                             % I = xLos 129, yLos 130, zLos 107
            %                             s(129)  = s(129) + f(133)*(  q(129)*(xUnit'*InvCovM*xUnit) ...
            %                                 + q(130)*(xUnit'*InvCovM*yUnit) ...
            %                                 + q(107)*(xUnit'*InvCovM*zUnit));
            %
            %                             % yLos  J = 130 KJ = yUnit'*InvCovM*losDiff
            %                             % I = xLos 129, yLos 130, zLos 107
            %                             s(130)  = s(130) + f(133)*(  q(129)*(yUnit'*InvCovM*xUnit) ...
            %                                 + q(130)*(yUnit'*InvCovM*yUnit) ...
            %                                 + q(107)*(yUnit'*InvCovM*zUnit));
            %
            %                             % zLos  J = 107 KJ = xUnit'*InvCovM*losDiff
            %                             % I = xLos 129, yLos 130, zLos 107
            %                             s(107)  = s(107) + f(133)*(  q(129)*(zUnit'*InvCovM*xUnit) ...
            %                                 + q(130)*(zUnit'*InvCovM*yUnit) ...
            %                                 + q(107)*(zUnit'*InvCovM*zUnit));
            %                         end
            %                     end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            S = s;
        end
        
        
        function [Jacob, JacobFinite, Hess, HessFinite, ForDerivs, DerivsDataList] = DigitalJustice(this,mu,I)
            
            Nvar          = 7;
            time          = this.time;
            Rpos          = this.Rpos;
            Rdot          = this.Rdot;
            Kepler        = [Rpos(1), Rpos(2), Rpos(3), Rdot(1), Rdot(2), Rdot(3)]';
            
            Kdelta        = KeplerFromECI(time, Rpos, Rdot, mu);
            InMotion      = Extrapolator(Kdelta, time);
            [H0, Jacob, Hess] = WorkOrder(Kdelta, I);
            DataList0     = Kdelta.DataList';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Checking First Derivatives
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            format short
            iStart = 1;
            Nlines = 115;
            CombinedDifferences = 0.0;
            ForDerivs  = [];
            for LineCheck = 1:Nlines
                %BackDerivs = [];
                %Derivs     = [];
                %Combined   = [];
                
                %F            = zeros(Nlines,1);
                %F(LineCheck) = 1;
                %F            = backdiff(this, F, LineCheck);
                % BackDerivs   = [F(1) F(2) F(3) F(4) F(5) F(6) F(7)];
                % Combined     = [Index F];
                % Kate's Combined'
                
                forward      = [];
                %for jk = 1:Nlines % check range of derivatives wrt forwards differentiation
                for jk = 1:Nvar % check range of derivatives wrt forwards differentiation
                    %jk
                    Q        = zeros(Nlines,1);
                    Q(jk)     = 1;
                    Q        = fordiff(this, Q, LineCheck);
                    % Q(LineCheck) = dK_LineCheck/dLj
                    %forward  = [forward; j Q(LineCheck)];
                    %forward  = [forward; jk+1-iStart Q(LineCheck)];
                    forward  = [forward; Q(LineCheck)];
                    %Derivs   = [Derivs; Q(LineCheck)];
                    %Combined = [Combined Q];
                end
                ForDerivs    = [ForDerivs, forward];
                %Derivs       = [ForDerivs F (F-ForDerivs(:,2))];
                %Derivs(iStart:LineCheck,:)
                %Difference   = sum(F-ForDerivs(:,2));
                %CombinedDifferences = CombinedDifferences + Difference;
                %disp([' Line = ', num2str(LineCheck),'  Sum of Differences = ',num2str(Difference)]);
                % Kate's ForDerivs'
                % Derivs      = [F Derivs (F-Derivs)];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %     ForDerivs    = [ForDerivs; forward];
                %     Derivs = [ForDerivs F (F-ForDerivs(:,2))];
                %     Derivs = Derivs(1:LineCheck,:)
                %     %Derivs = [F Derivs (F-Derivs)];
                % 'Debug Point'
            end
            %disp([' Combined Sum of Differences = ',num2str(CombinedDifferences)]);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Do time Increment Separately
            
            epsilon     = 0.000000001;
            KeplerDelta = epsilon*eye(6);
            Kdelta = KeplerFromECI(time+epsilon, Kepler(1:3), Kepler(4:6), mu);
            InMotion      = Extrapolator(Kdelta, time+epsilon);
            [HT, JacobT, HessT] = Kdelta.WorkOrder(I);
            DataListT    = Kdelta.DataList';
            
            Kepler1     = Kepler + KeplerDelta(:,1);
            Kdelta1 = KeplerFromECI(time, Kepler1(1:3), Kepler1(4:6), mu);
            InMotion       = Extrapolator(Kdelta1, time);
            [H1, Jacob1, Hess1] = Kdelta1.WorkOrder(I);
            DataList1    = Kdelta1.DataList';
            
            Kepler2     = Kepler + KeplerDelta(:,2);
            Kdelta2 = KeplerFromECI(time, Kepler2(1:3), Kepler2(4:6), mu);
            InMotion       = Extrapolator(Kdelta2, time);
            [H2, Jacob2, Hess2] = Kdelta2.WorkOrder(I);
            DataList2    = Kdelta2.DataList';
            
            Kepler3     = Kepler + KeplerDelta(:,3);
            Kdelta3 = KeplerFromECI(time, Kepler3(1:3), Kepler3(4:6), mu);
            InMotion       = Extrapolator(Kdelta3, time);
            [H3, Jacob3, Hess3] = Kdelta3.WorkOrder(I);
            DataList3    = Kdelta3.DataList';
            
            Kepler4     = Kepler + KeplerDelta(:,4);
            Kdelta4 = KeplerFromECI(time, Kepler4(1:3), Kepler4(4:6), mu);
            InMotion       = Extrapolator(Kdelta4, time);
            [H4, Jacob4, Hess4] = Kdelta4.WorkOrder(I);
            DataList4    = Kdelta4.DataList';
            
            Kepler5     = Kepler + KeplerDelta(:,5);
            Kdelta5 = KeplerFromECI(time, Kepler5(1:3), Kepler5(4:6), mu);
            InMotion       = Extrapolator(Kdelta5, time);
            [H5, Jacob5, Hess5] = Kdelta5.WorkOrder(I);
            DataList5    = Kdelta5.DataList';
            
            Kepler6     = Kepler + KeplerDelta(:,6);
            Kdelta6 = KeplerFromECI(time, Kepler6(1:3), Kepler6(4:6), mu);
            InMotion       = Extrapolator(Kdelta6, time);
            [H6, Jacob6, Hess6] = Kdelta6.WorkOrder(I);
            DataList6    = Kdelta6.DataList';
            
            JacobFinite = [HT-H0,H1-H0, H2-H0, H3-H0, H4-H0, H5-H0, H6-H0]/epsilon;
            
            HessFinite = [(JacobT-Jacob);...
                (Jacob1-Jacob);...
                (Jacob2-Jacob);...
                (Jacob3-Jacob);...
                (Jacob4-Jacob);...
                (Jacob5-Jacob);...
                (Jacob6-Jacob)]/epsilon;
            
            DerivsDataList = [DataListT - DataList0;...
                DataList1 - DataList0;...
                DataList2 - DataList0;...
                DataList3 - DataList0;...
                DataList4 - DataList0;...
                DataList5 - DataList0;...
                DataList6 - DataList0]/epsilon;
            
            % NumericalDerivs = [];
            % for J = 1:Nlines
            %     DataListTemp = DataList0;
            %     %DataListTemp(J) = DataListTemp(J) + epsilon;
            %     tTemp = DataListTemp(1);
            %     rTemp = DataListTemp(2:4);
            %     vTemp = DataListTemp(4:7);
            %     KdeltaTemp = KeplerFromECI(tTemp, rTemp, vTemp, mu);
            %     InMotion   = Extrapolator(KdeltaTemp, time);
            %     DataListTemp = KdeltaTemp.DataList';
            %     Delta = (DataListTemp - DataList0)/epsilon
            %     NumericalDerivs = [NumericalDerivs; Delta];
            % end
            % ForDerivs - NumericalDerivs
            
        end
        
        function [H_I, Jacob, Hess] = WorkOrder(this,I)
            
            H_I          = this.DataList(I);
            ipoint = [1, 2, 3, 4, 5, 6, 7];
            %ipoint = [8, 2, 3, 4, 5, 6, 7];
            
            Nvar         = 7;
            LineCheck    = I;
            F            = zeros(115,1);
            F(LineCheck) = 1;
            F = backdiff(this, F, LineCheck);
            Jacob        = [F(1) F(2) F(3) F(4) F(5) F(6) F(7)];
            % e,a,...         t    e    a    I    w    W    Mp
            %Jacob        = [F(ipoint)]';
            % e,a,...   tExtrap    e    a    I    w    W    Mp
            Hess   = [];
            for j = 1: Nvar
                Q = zeros(115,1);
                S = zeros(115,1);
                %Q(ipoint(j)) = 1.0;
                Q(j) = 1.0;
                Q = fordiff(this, Q, LineCheck);
                S = secdiff(this, F, Q, S, LineCheck);
                S = backdiff(this, S, LineCheck);
                row = [];
                for k = 1 : Nvar
                    row = [row, S(k)];
                    %row = [row, S(ipoint(k))];
                end
                %row = [S(ipoint)]';
                Hess  = [Hess; row];
            end
            
        end
        
        function [Nlines] = CheckDerivatives(time, varargin)
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Checking First Derivatives
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %Index = [1:56]';
            % %    for LineCheck = 76:81 %LineCheck  = 90
            % format short
            % iStart = 1;
            % for LineCheck = iStart:133
            %     ForDerivs  = [];
            %     %BackDerivs = [];
            %     %Derivs     = [];
            %     %Combined   = [];
            %
            %     F            = zeros(133,1);
            %     F(LineCheck) = 1;
            %     F            = backdiff(this, F, LineCheck);
            %     % BackDerivs   = [F(1) F(2) F(3) F(4) F(5) F(6) F(7)];
            %     % Combined     = [Index F];
            %     % Kate's Combined'
            %
            %     forward      = [];
            %     for jk = 1:133 % check range of derivatives wrt forwards differentiation
            %         %jk
            %         Q        = zeros(133,1);
            %         Q(jk)     = 1;
            %         Q        = fordiff(this, Q, LineCheck);
            %         % Q(LineCheck) = dK_LineCheck/dLj
            %         %forward  = [forward; j Q(LineCheck)];
            %         forward  = [forward; jk+1-iStart Q(LineCheck)];
            %         %Derivs   = [Derivs; Q(LineCheck)];
            %         %Combined = [Combined Q];
            %     end
            %     ForDerivs    = [ForDerivs; forward];
            %     Derivs       = [ForDerivs F (F-ForDerivs(:,2))];
            %     Derivs(iStart:LineCheck,:)
            %     disp([' Line = ', num2str(LineCheck),'  Sum of Differences = ',num2str(sum(F-ForDerivs(:,2)))]);
            %     % Kate's ForDerivs'
            %     % Derivs      = [F Derivs (F-Derivs)];
            %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %     %     ForDerivs    = [ForDerivs; forward];
            %     %     Derivs = [ForDerivs F (F-ForDerivs(:,2))];
            %     %     Derivs = Derivs(1:LineCheck,:)
            %     %     %Derivs = [F Derivs (F-Derivs)];
            %     % 'Debug Point'
            % end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %%% TEST SECOND DERIVATIVES:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Debug Second Derivatives by checking symmetry:
            %Begin Checking Out Second Derivatives:
            Nvar     = 7;
            % LineCheck    = 56;  % Kepler Part
            % LineCheck    = 90;  % Escobal Gravitation to J6
            % F            = zeros(115,1);
            % F(LineCheck) = 1;
            % F = backdiff(this, F, LineCheck);
            Nlines = 115;
            format short
            for LineCheck = 1:Nlines
                ListMax   = LineCheck + 1;
                %Nvar = LineCheck;
                Hessian      = [];
                F            = zeros(Nlines,1);
                F(LineCheck) = 1;
                F = backdiff(this, F, LineCheck);
                %F(ListMax:Nlines) = 0;
                % for j = 1:Nvar
                for j = 1:LineCheck
                    Q = zeros(Nlines,1);
                    S = zeros(Nlines,1);
                    Q(j) = 1.0;
                    Q = fordiff(this, Q, LineCheck);
                    %Q(ListMax:Nlines) = 0;
                    S = secdiff(this, F,Q,S, LineCheck);
                    %S(ListMax:Nlines) = 0;
                    S = backdiff(this, S, LineCheck);
                    %S(ListMax:Nlines) = 0;
                    row = [];
                    %for k = 1:Nvar
                    for k = 1:LineCheck
                        row = [row S(k)];
                    end
                    Hessian  = [Hessian; row];
                end
                %disp([' Line = ', num2str(Nvar)]);
                maxElement = max(max(abs(Hessian)));
                HessianTest = Hessian;
                if maxElement > 0
                    HessianTest = Hessian/maxElement
                end
                Symmetry = HessianTest - HessianTest'
                max(max(Symmetry))
                disp([' Line Check = ', num2str(LineCheck)]);
                %disp('Debug Point')
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end% end of member methods
    
end
