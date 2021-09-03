        function [Inhomogeneous, this] = LagrangePlanetary_odefun(TrackParams, time, Y, units)
            % Initialize the "this" structure.
            this = InitialValues();

            % For Choice of units
            TU = units.TU;
            DU = units.DU;
            VU = units.VU;
            AU = units.AU;
            mu = units.mu;
            this.mu     = mu;
            this.sqrtmu = sqrt(mu);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            wgs84       = wgs84Constants();
            this.J2     = wgs84.J2;        % J2 perturbation constant
            this.J3     = wgs84.J3;
            this.J4     = wgs84.J4;
            this.J5     = wgs84.J5;
            this.J6     = wgs84.J6;
            % Truncate Gravity Model to Kepler only
            %this.J2     = 0.0;
            %this.J3     = 0.0;
            %this.J4     = 0.0;
            %this.J5     = 0.0;
            %this.J6     = 0.0;
            %this.CK     = 0.0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % To Compute the Disturbing Potential Only set CK = 0
            %this.CK     = 1.0;   % Kepler contribution to Gravity Potential
            this.CK     = 0.0;
            this.C2     = 0.5000*this.J2;
            this.C3     = 0.5000*this.J3;
            this.C4     =-0.1250*this.J4;
            this.C5     =-0.1250*this.J5;
            this.C6     = 0.0625*this.J6;
            % J4 only for now: 
            %this.C5     = 0.0;
            %this.C6     = 0.0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            TrackParams = TrackParams + Y';
            this.time        = time; 
            this.t           = time;                             % line   1
            this.e           = TrackParams(1);                   % line   2
            %this.p           = TrackParams(2);                  % line   3
            this.a           = TrackParams(2);                   % line   3
            this.Inclination = TrackParams(3);                   % line   4
            this.omega       = TrackParams(4);                   % line   5
            this.Omega       = TrackParams(5);                   % line   6
            this.Mp          = TrackParams(6);                   % line   7
            this.ClassicalElements = [this.e; this.a; this.Inclination; this.omega; this.Omega; this.Mp; this.t];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Computations based on the above input variables begins ---
            % start with line 11 to leave room for more input parameters in
            % lines 1-10:
            
            this.onePe = 1.0 + this.e;                           % line  11
            this.oneMe = 1.0 - this.e;                           % line  12
            this.fac   = this.onePe*this.oneMe;                  % line  13
            if this.fac > 0.0
                this.rootfac = sqrt(this.fac);                   % line  14
            else
                this.rootfac = 0;                                % line  14
            end
            %this.a           = this.p/this.fac;                 % line  15
            this.p           = this.a*this.fac;                  % line  15
            this.meanMotion  = this.sqrtmu*(this.a)^(-1.5);      % line  16
            
            cosI             = cos(this.Inclination);            % line  17
            sinI             = sin(this.Inclination);            % line  18
            cosom            = cos(this.omega);                  % line  19
            sinom            = sin(this.omega);                  % line  20
            cosO             = cos(this.Omega);                  % line  21
            sinO             = sin(this.Omega);                  % line  22
            
            this.Px   =  cosO*cosom - sinO*sinom*cosI;           % line  23
            this.Py   =  sinO*cosom + cosO*sinom*cosI;           % line  24
            this.Pz   =  sinom*sinI;                             % line  25
            
            this.Qx   = -cosO*sinom - sinO*cosom*cosI;           % line  26
            this.Qy   = -sinO*sinom + cosO*cosom*cosI;           % line  27
            this.Qz   =  cosom*sinI;                             % line  28
            
            % For Pure "Kepler" orbit we don't need Wvec, but it might play
            % a role later when we include Perturbations, so leave it in
            this.Wx   =  sinO*sinI;                              % line  29
            this.Wy   = -cosO*sinI;                              % line  30
            this.Wz   =  cosI;                                   % line  31
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Class Variable Storage Assignments -- not "computations"
            this.cosI   = cosI;
            this.sinI   = sinI;
            this.cosom  = cosom;
            this.sinom  = sinom;
            this.cosO   = cosO;
            this.sinO   = sinO;
            this.Pvec = [this.Px; this.Py; this.Pz];
            this.Qvec = [this.Qx; this.Qy; this.Qz];
            this.Wvec = [this.Wx; this.Wy; this.Wz];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % variables for use in Derivatives based on Implicit Function Theorem
            % on M = E - e*sin(E)
            this.eDenom   = 0.;
            this.Eprime   = 0.;
            this.dE_dM    = 0.;
            this.dE_de    = 0.;
            this.d2E_dMdM = 0.;
            this.d2E_dMde = 0.;
            this.d2E_dedM = 0.;
            this.d2E_dede = 0.;
            
            % invert Kepler's Equation
            % Brute Force iteration (good way to seed Newton's method which follows)
            % M = MeanAnomalyEpoch;
            %  meanAnomaly M to eccentric Anomaly E
            M = this.meanMotion*this.time + this.Mp;
            E = M + this.e;
            if (M > pi)
                E = M - this.e;
            elseif (-pi < M)
                if (M < 0)
                    E = M - this.e;
                end
            end
            %         //cout << " Input meanAnomaly = " << meanAnomaly << endl;
            for k=1:1:10
                E  = M + this.e*sin(E);
                %//cout << " Mean Anomaly Solution " << E << endl;
            end
            % //      10 rounds of Newton's root finding method based on the above "seed".
            for k=1:1:10
                Eprime      = 1.0 - this.e*cos(E);
                E           = E + (M - E + this.e*sin(E))./Eprime;
            end
            %KeplerInv = E;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cosK  = cos(E);                                   % line 42
            sinK  = sin(E);                                   % line 43
            %E     = atan2(sinK,cosK)
            %%%%%%%%%%%%%%%%%%% Implicit Function Theorem Derivatives
            eDenom   = Eprime*Eprime*Eprime;
            this.dE_dM    = 1.0/Eprime;
            this.dE_de    = sinK/Eprime;
            this.d2E_dMdM = -this.e*sinK/eDenom;
            %// Mike Cain Corrections!
            this.d2E_dMde = (cosK - this.e)/eDenom;
            this.d2E_dedM = (cosK - this.e)/eDenom;
            this.d2E_dede = ((2.0 - this.e*cosK)*cosK*sinK - this.e*sinK)/eDenom;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   % for convenience -- you can eliminate kfac;
            tanX  = cosK - this.e;                                 % line 44
            tanY  = this.rootfac*sinK;                             % line 45
            nu     = atan2(tanY, tanX);                       % line 46
            coss  = cos(nu);                                  % line 47
            sins  = sin(nu);                                  % line 48
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
            %rorbit   = p/(1.0 + e*coss);                        % line  49
            rorbit   = this.a*(1.0 - this.e*cosK);                    % line  49
            rorbitx  = rorbit*(coss*this.Px + sins*this.Qx);               % line  50
            rorbity  = rorbit*(coss*this.Py + sins*this.Qy);               % line  51
            rorbitz  = rorbit*(coss*this.Pz + sins*this.Qz);               % line  52
            rtpinv   = sqrt(mu/this.p);                               % line  53
            vorbitx  = rtpinv*(-sins*this.Px + (this.e + coss)*this.Qx);        % line  54
            vorbity  = rtpinv*(-sins*this.Py + (this.e + coss)*this.Qy);        % line  55
            vorbitz  = rtpinv*(-sins*this.Pz + (this.e + coss)*this.Qz);        % line  56
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Rextrap  = [rorbitx; rorbity; rorbitz];
            Vextrap  = [vorbitx; vorbity; vorbitz];
            rmag     = norm(Rextrap);
            this.rmag = rmag;
            Rpos  = Rextrap;
            Rdot  = [vorbitx' vorbity' vorbitz'];
            %Rextrap  = rorbit*( coss*this.Pvec + sins*this.Qvec);
            %Vextrap  = rtpinv*(-sins*this.Pvec + (e + coss)*this.Qvec);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            this.M       = M;
            this.E       = E;
            this.cosK    = cosK;
            this.sinK    = sinK;
            this.tanX    = tanX;
            this.tanY    = tanY;
            this.nu      = nu;
            this.coss    = coss;
            this.sins    = sins;
            this.rorbit  = rorbit;
            this.rorbitx = rorbitx;
            this.rorbity = rorbity;
            this.rorbitz = rorbitz;
            this.rtpinv  = rtpinv;
            this.vorbitx = vorbitx;
            this.vorbity = vorbity;
            this.vorbitz = vorbitz;
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % this.DataList(  1) = this.time;                      % line   1
            % this.DataList(  2) = this.e;                         % line   2
            % this.DataList(  3) = this.a;                         % line   3
            % this.DataList(  4) = this.Inclination;               % line   4
            % this.DataList(  5) = this.omega;                     % line   5
            % this.DataList(  6) = this.Omega;                     % line   6
            % this.DataList(  7) = this.Mp;                        % line   7
            % this.DataList( 11) = this.onePe;                     % line  11
            % this.DataList( 12) = this.oneMe;                     % line  12
            % this.DataList( 13) = this.fac;                       % line  13
            % if this.fac > 0.0
            %     this.DataList( 14) = this.rootfac;               % line  14
            % else
            %     this.DataList( 14) = 0;                          % line  14
            % end
            % this.DataList( 15) = this.p;                         % line  15
            % this.DataList( 16) = this.meanMotion;                % line  16
            % this.DataList( 17) = this.cosI;                      % line  17
            % this.DataList( 18) = this.sinI;                      % line  18
            % this.DataList( 19) = this.cosom;                     % line  19
            % this.DataList( 20) = this.sinom;                     % line  20
            % this.DataList( 21) = this.cosO;                      % line  21
            % this.DataList( 22) = this.sinO;                      % line  22
            %
            % this.DataList( 23) = this.Px;                        % line  23
            % this.DataList( 24) = this.Py;                        % line  24
            % this.DataList( 25) = this.Pz;                        % line  25
            %
            % this.DataList( 26) = this.Qx;                        % line  26
            % this.DataList( 27) = this.Qy;                        % line  27
            % this.DataList( 28) = this.Qz;                        % line  28
            %
            % this.DataList( 29) = this.Wx;                        % line  29
            % this.DataList( 30) = this.Wy;                        % line  30
            % this.DataList( 31) = this.Wz;                        % line  31
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Class Variable Storage Assignments -- not "computations"
            % %this.DataList() = this.Pvec;
            % %this.DataList() = this.Qvec;
            % %this.DataList() = this.Wvec;
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % variables for use in Derivatives based on Implicit Function Theorem
            % % on M = E - e*sin(E)
            % %this.DataList() = this.eDenom;
            % %this.DataList() = this.Eprime;
            % %this.DataList() = this.dE_dM;
            % %this.DataList() = this.dE_de;
            % %this.DataList() = this.d2E_dMdM;
            % %this.DataList() = this.d2E_dMde;
            % %this.DataList() = this.d2E_dedM;
            % %this.DataList() = this.d2E_dede;
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % this.DataList( 40) = this.M;                         % line  40
            % this.DataList( 41) = this.E;                         % line  41
            % this.DataList( 42) = this.cosK;                      % line  42
            % this.DataList( 43) = this.sinK;                      % line  43
            % %%%%%%%%%%%%%%%%%%% Implicit Function Theorem Derivatives
            % this.DataList( 44) = this.tanX;                      % line  44
            % this.DataList( 45) = this.tanY;                      % line  45
            % this.DataList( 46) = this.nu;                        % line  46
            % this.DataList( 47) = this.coss;                      % line  47
            % this.DataList( 49) = this.sins;                      % line  48
            % % OLD: rorbit   = p/(1.0 + e*coss);                  % line  49
            % this.DataList( 49) = this.rorbit;                    % line  49
            % this.DataList( 50) = this.rorbitx;                   % line  50
            % this.DataList( 51) = this.rorbity;                   % line  51
            % this.DataList( 52) = this.rorbitz;                   % line  52
            % this.DataList( 53) = this.rtpinv;                    % line  53
            % this.DataList( 54) = this.vorbitx;                   % line  54
            % this.DataList( 55) = this.vorbity;                   % line  55
            % this.DataList( 56) = this.vorbitz;                   % line  56
            % %  rmag
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %Kvec   = Rvec(1:3)/DU;
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %X    = Kvec(1);                                     % line  61
            % %Y    = Kvec(2);                                     % line  62
            % %Z    = Kvec(3);                                     % line  63
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Change Order of R1 and R2 and skip X, Y, Z
            X    = rorbitx;                                      % line  61
            Y    = rorbity;                                      % line  62
            Z    = rorbitz;                                      % line  63
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            R2   = X^2 + Y^2 + Z^2;                              % line  64
            R1   = sqrt(R2);                                     % line  65
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %R1   = rorbit;                                      % line  64
            %R2   = R1^2;                                        % line  65
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            R3   = R1*R2;                                        % line  66
            R4   = R1*R3;                                        % line  67
            R5   = R1*R4;                                        % line  68
            R6   = R1*R5;                                        % line  69
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sd1  = rorbitz/R1;                                   % line  70
            % Compute Perturbing Function:
            %sinomnu = sin(omega + nu);
            %cosomnu = cos(omega + nu);
            %sd1chk  = sinI*sinomnu;
            %cd1chk  = cosI*sinomnu;
            %testDiff  = sd1chk - sd1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sd2  = sd1^2;                                        % line  71
            sd3  = sd1*sd2;                                      % line  72
            sd4  = sd1*sd3;                                      % line  73
            sd5  = sd1*sd4;                                      % line  74
            sd6  = sd1*sd5;                                      % line  75
            
            F2   =  1.0      -  3.0*sd2;                         % line  76
            F3   =  3.0*sd1  -  5.0*sd3;                         % line  77
            F4   =  3.0      -  30.*sd2  + 35*sd4;               % line  78
            F5   =  15.*sd1  -  70.*sd3  + 63*sd5;               % line  79
            F6   =  5.       - 105.*sd2  + 315*sd4 - 231*sd6;    % line  80
            
            V1      = 1.0/R1;                                    % line  81
            V2      = this.C2*F2/R2;                                  % line  82
            V3      = this.C3*F3/R3;                                  % line  83
            V4      = this.C4*F4/R4;                                  % line  84
            V5      = this.C5*F5/R5;                                  % line  85
            V6      = this.C6*F6/R6;                                  % line  86
            
            V       = mu*V1;           % Kepler Potential        % line  87
            % J6:
            D       = V2 + V3 + V4 + V5 + V6;                    % line  88
            % R is the disturbing potential
            % R       = mu*D;                                    % line  89
            Phi     = V*(this.CK + D);                                % line  90
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            this.X   = X;
            this.Y   = Y;
            this.Z   = Z;
            this.R2  = R2;
            this.R1  = R1;
            this.R3  = R3;
            this.R4  = R4;
            this.R5  = R5;
            this.R6  = R6;
            %this.sinomnu = sinomnu;
            %this.cosomnu = cosomnu;
            %this.sd1chk  = sd1chk;
            %this.cd1chk  = cd1chk;
            this.sd1 = sd1;
            this.sd2 = sd2;
            this.sd3 = sd3;
            this.sd4 = sd4;
            this.sd5 = sd5;
            this.sd6 = sd6;
            this.F2  = F2;
            this.F3  = F3;
            this.F4  = F4;
            this.F5  = F5;
            this.F6  = F6;
            this.V1  = V1;
            this.V2  = V2;
            this.V3  = V3;
            this.V4  = V4;
            this.V5  = V5;
            this.V6  = V6;
            this.V   = V;
            this.D   = D;
            this.Phi = Phi;
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Treat Equinoctial Variables:
            % Aeq      = this.omega + this.Omega;                            % line 142
            % cosAeq   = cos(Aeq);                                 % line 143
            % sinAeq   = sin(Aeq);                                 % line 144
            % tanHalf  = tan(this.Inclination/2.0);                               % line 145
            % Feq      = this.e*cosAeq;                                 % line 146
            % Geq      = this.e*sinAeq;                                 % line 147
            % Heq      = tanHalf*cosO;                             % line 148
            % Keq      = tanHalf*sinO;                             % line 149
            % Leq      = Aeq + nu;                                 % line 150
            % CosL     = cos(Leq);                                 % line 151
            % SinL     = sin(Leq);                                 % line 152
            % alphaSq = Heq^2 - Keq^2;                             % line 153
            % Seq     = 1.0 + Heq^2 + Keq^2;                       % line 154
            % Weq     = 1.0 + Feq*CosL + Geq*SinL;                 % line 155
            % Req     = this.p/Weq;                                     % line 156
            % RovS    = Req/Seq;                                   % line 157
            % srtpinv = rtpinv/Seq;                                % line 158
            % HK      = Heq*Keq;                                   % line 159
            % OnePalphaSqCosL = (1.+alphaSq)*CosL;                 % line 160
            % OneMalphaSqCosL = (1.-alphaSq)*CosL;                 % line 161
            % OnePalphaSqSinL = (1.+alphaSq)*SinL;                 % line 162
            % OneMalphaSqSinL = (1.-alphaSq)*SinL;                 % line 163
            % Xfac    = OnePalphaSqCosL + 2.0*HK*SinL;             % line 164
            % Yfac    = OneMalphaSqSinL + 2.0*HK*CosL;             % line 165
            % Zfac    = Heq*SinL - Keq*CosL;                       % line 166
            % VXfac   =  OnePalphaSqSinL - 2.0*HK*(Feq + CosL) + Geq*(1. + alphaSq);% line 167
            % VYfac   = -OneMalphaSqCosL + 2.0*HK*(Geq + SinL) + Feq*(alphaSq - 1.);% line 168
            % VZfac   =  Heq*(Feq+CosL) + Keq*(Geq+SinL);          % line 169
            % Xeq     =     RovS*Xfac;                             % line 170
            % Yeq     =     RovS*Yfac;                             % line 171
            % Zeq     = 2.0*RovS*Zfac;                             % line 172
            % VXeq    =    -srtpinv*VXfac;                         % line 173
            % VYeq    =    -srtpinv*VYfac;                         % line 174
            % VZeq    = 2.0*srtpinv*VZfac;                         % line 175
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % this.Aeq             = Aeq;                      % line 142 Aeq
            % this.cosAeq          = cosAeq;                   % line 143 cosAeq
            % this.sinAeq          = sinAeq;                   % line 144 sinAeq
            % this.tanHalf         = tanHalf;                  % line 145 tanHalf
            % this.Feq             = Feq;                      % line 146 Feq
            % this.Geq             = Geq;                      % line 147 Geq
            % this.Heq             = Heq;                      % line 148 Heq
            % this.Keq             = Keq;                      % line 149 Keq
            % this.Leq             = Leq;                      % line 150 Leq
            % this.CosL            = CosL;                     % line 151 CosL
            % this.SinL            = SinL;                     % line 152 SinL
            % this.alphaSq         = alphaSq;                  % line 153 alphaSq
            % this.Seq             = Seq;                      % line 154
            % this.Weq             = Weq;                      % line 155
            % this.Req             = Req;                      % line 156
            % this.RovS            = RovS;                     % line 157
            % this.srtpinv         = srtpinv;                  % line 158
            % this.HK              = HK;                       % line 159
            % this.OnePalphaSqCosL = OnePalphaSqCosL;          % line 160
            % this.OneMalphaSqCosL = OneMalphaSqCosL;          % line 161
            % this.OnePalphaSqSinL = OnePalphaSqSinL;          % line 162
            % this.OneMalphaSqSinL = OneMalphaSqSinL;          % line 163
            % this.Xfac            = Xfac;                     % line 164
            % this.Yfac            = Yfac;                     % line 165
            % this.Zfac            = Zfac;                     % line 166
            % this.VXfac           = VXfac;                    % line 167
            % this.VYfac           = VYfac;                    % line 168
            % this.VZfac           = VZfac;                    % line 169
            % this.Xeq             = Xeq;                      % line 170
            % this.Yeq             = Yeq;                      % line 171
            % this.Zeq             = Zeq;                      % line 172
            % this.VXeq            = VXeq;                     % line 173
            % this.VYeq            = VYeq;                     % line 174
            % this.VZeq            = VZeq;                     % line 175
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %B = Runit';
            % %B = B .^2;
            % %sum(B);
            % Jacobian = [];
            % %Hessian  = cell(6,6);
            %
            % [H_I, Jacob, Hess] = WorkOrder(this, 50);
            % %Hessian{1}  = Hess;
            % Jacobian      = [Jacobian; Jacob];   % row 1
            %
            % %LineCheck    = 51;                % Ry
            % [H_I, Jacob, Hess] = WorkOrder(this, 51);
            % %Hessian{2}  = Hess;
            % Jacobian      = [Jacobian; Jacob];   % row 1
            %
            % %LineCheck    = 52;                % Rz
            % [H_I, Jacob, Hess] = WorkOrder(this, 52);
            % %Hessian{3}  = Hess;
            % Jacobian      = [Jacobian; Jacob];   % row 1
            %
            % %LineCheck    = 54;                % Vx
            % [H_I, Jacob, Hess] = WorkOrder(this, 54);
            % %Hessian{4}  = Hess;
            % Jacobian      = [Jacobian; Jacob];   % row 1
            %
            % %LineCheck    = 55;                % Vy
            % [H_I, Jacob, Hess] = WorkOrder(this, 55);
            % %Hessian{5}  = Hess;
            % Jacobian      = [Jacobian; Jacob];   % row 1
            %
            % %LineCheck    = 56;                % Vz
            % [H_I, Jacob, Hess] = WorkOrder(this, 56);
            % %Hessian{6}  = Hess;
            % Jacobian      = [Jacobian; Jacob];   % row 1
            % %this.JacKepler_2_ECI = Jacobian(:,2:7);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            this.StateVector7 = [rorbitx; rorbity; rorbitz; vorbitx; vorbity; vorbitz; time];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Change Order of R1 and R2 and skip X, Y, Z
            % X    = rorbitx;                                      % line  61
            % Y    = rorbity;                                      % line  62
            % Z    = rorbitz;                                      % line  63
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % R2   = X^2 + Y^2 + Z^2;                              % line  64
            % R1   = sqrt(R2);                                     % line  65
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %R1   = rorbit;                                      % line  64
            % %R2   = R1^2;                                        % line  65
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % R3   = R1*R2;                                        % line  66
            % R4   = R1*R3;                                        % line  67
            % R5   = R1*R4;                                        % line  68
            % R6   = R1*R5;                                        % line  69
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % sd1  = rorbitz/R1;                                   % line  70
            % % Compute Perturbing Function:
            % %sinomnu = sin(omega + nu);
            % %cosomnu = cos(omega + nu);
            % %sd1chk  = sinI*sinomnu;
            % %cd1chk  = cosI*sinomnu;
            % %testDiff  = sd1chk - sd1
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % sd2  = sd1^2;                                        % line  71
            % sd3  = sd1*sd2;                                      % line  72
            % sd4  = sd1*sd3;                                      % line  73
            % sd5  = sd1*sd4;                                      % line  74
            % sd6  = sd1*sd5;                                      % line  75
            %
            % F2   =  1.0      -  3.0*sd2;                         % line  76
            % F3   =  3.0*sd1  -  5.0*sd3;                         % line  77
            % F4   =  3.0      -  30.*sd2  + 35*sd4;               % line  78
            % F5   =  15.*sd1  -  70.*sd3  + 63*sd5;               % line  79
            % F6   =  5.       - 105.*sd2  + 315*sd4 - 231*sd6;    % line  80
            %
            % V1      = 1.0/R1;                                    % line  81
            % V2      = C2*F2/R2;                                  % line  82
            % V3      = C3*F3/R3;                                  % line  83
            % V4      = C4*F4/R4;                                  % line  84
            % V5      = C5*F5/R5;                                  % line  85
            % V6      = C6*F6/R6;                                  % line  86
            %
            % V       = mu*V1;           % Kepler Potential        % line  87
            % % J6:
            % D       = V2 + V3 + V4 + V5 + V6;                    % line  88
            % % R is the disturbing potential
            % % R       = mu*D;                                    % line  89
            % Phi     = V*(CK + D);                                % line  90
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % this.X   = X;
            % this.Y   = Y;
            % this.Z   = Z;
            % this.R2  = R2;
            % this.R1  = R1;
            % this.R3  = R3;
            % this.R4  = R4;
            % this.R5  = R5;
            % this.R6  = R6;
            % %this.sinomnu = sinomnu;
            % %this.cosomnu = cosomnu;
            % %this.sd1chk  = sd1chk;
            % %this.cd1chk  = cd1chk;
            % this.sd1 = sd1;
            % this.sd2 = sd2;
            % this.sd3 = sd3;
            % this.sd4 = sd4;
            % this.sd5 = sd5;
            % this.sd6 = sd6;
            % this.F2  = F2;
            % this.F3  = F3;
            % this.F4  = F4;
            % this.F5  = F5;
            % this.F6  = F6;
            % this.V1  = V1;
            % this.V2  = V2;
            % this.V3  = V3;
            % this.V4  = V4;
            % this.V5  = V5;
            % this.V6  = V6;
            % this.V   = V;
            % this.D   = D;
            % this.Phi = Phi;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Coefficient Matrix in Lagrange Planetary Equations
            n    = this.meanMotion;                              % line  16
            fac1 = n*this.a         ;% n*a                            % line  91
            fac2 = fac1*this.a      ;% n*a^2                          % line  92
            fac3 = fac2*this.e      ;% n*a^2*e                        % line  93
            fac4 = fac2*this.rootfac;% n*a^2*sqrt(1-e^2)              % line  94
            fac5 = fac4*sinI   ;% n*a^2*sqrt(1-e^2)*sinI         % line  95
            Mwe  = this.rootfac/fac3;% sqrt(1-e^2)/(n*a^2*e)          % line  96
            MeM  = this.fac/fac3    ;% (1-e^2)/(n*a^2*e)              % line  97
            MaM  = 2.0/fac1    ;% 2/(n*a)                        % line  98
            MWi  = 1/fac5      ;% 1/(n*a^2*sqrt(1-e^2)*sinI)     % line  99
            Miw  = cosI*MWi    ;% cosI/(n*a^2*sqrt(1-e^2)*sinI)  % line 100
            
            MLagrange = [ 0,    0,    0,  -Mwe,     0,  MeM;...
                0,    0,    0,     0,     0,  MaM;...
                0,    0,    0,   Miw,  -MWi,    0;...
                Mwe,    0,  -Miw,    0,     0,    0;...
                0,    0,   MWi,    0,     0,    0;...
                -MeM, -MaM,    0,     0,     0,    0];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            this.fac1 = fac1;
            this.fac2 = fac2;
            this.fac3 = fac3;
            this.fac4 = fac4;
            this.fac5 = fac5;
            this.Mwe  = Mwe;
            this.MeM  = MeM;
            this.MaM  = MaM;
            this.MWi  = MWi;
            this.Miw  = Miw;
            this.MLagrange = MLagrange;
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % los       = Rextrap - state_vector_sat(1:3);     %lines 101-103
            % %los       = los/norm(los);    % JRS EXPERIMENTAL
            % % On the fly Backwards Differentiation:
            % % Prediction:
            % lx = los(1);                                         % line 101
            % ly = los(2);                                         % line 102
            % lz = los(3);                                         % line 103
            % %range   = norm(los)
            % %theta   = acos(los(3)/range)
            % %phi     = atan2(los(2),los(1))
            % lperpSq  = lx^2 + ly^2;                              % line 104
            % rangeSq  = lperpSq + lz^2;                           % line 105
            % range    = sqrt(abs(rangeSq));                       % line 106
            % zLos     = lz/range;                                 % line 107
            % thetaLos = acos(zLos);                               % line 108
            % phiLos   = atan2(ly,lx);                             % line 109
            % vlos     = Vextrap - state_vector_sat(4:6);     % lines 110-112
            % vx       = vlos(1);                                  % line 110
            % vy       = vlos(2);                                  % line 111
            % vz       = vlos(3);                                  % line 112
            % ldotv    = lx*vx + ly*vy + lz*vz;                    % line 113
            % ldot     = ldotv/range;                              % line 114
            % lperp    = sqrt(lperpSq);                            % line 115
            % thdotfac1 = range*vz - lz*ldot;                      % line 116
            % thdotfac2 = thdotfac1/range;                         % line 117
            % thdot     = -thdotfac2/lperp;                        % line 118
            % phiDotNum = lx*vy - ly*vx;                           % line 119
            % phidot    = phiDotNum/lperpSq;                       % line 120
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % rangeM      = norm(losM);
            % losM        = losM/rangeM;
            % xLosM       = losM(1);
            % yLosM       = losM(2);
            % zLosM       = losM(3);
            % thetaM      = acos(losM(3));
            % phiM        = atan2(losM(2), losM(1));
            % cosThetaM   = cos(thetaM);
            % sinThetaM   = sin(thetaM);
            % cosPhiM     = cos(phiM);
            % sinPhiM     = sin(phiM);
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % cosThetaLos = cos(thetaLos);                         % line 121
            % sinThetaLos = sin(thetaLos);                         % line 122
            % phiDiff     = phiLos - phiM;                         % line 123
            % cosPhiDiff  = cos(phiDiff);                          % line 124
            % sinPhiDiff  = sin(phiDiff);                          % line 125
            % delUV       = cosThetaLos*sinThetaM*cosPhiDiff;      % line 126
            % Uview       = sinThetaLos*cosThetaM - delUV;         % line 127
            % Vview       = sinThetaM*sinPhiDiff;                  % line 128
            % xLos        = lx/range;                              % line 129
            % yLos        = ly/range;                              % line 130
            % Collinearity = 1-xLos*xLosM-yLos*yLosM-zLos*zLosM;   % line 131
            %
            % Mdelta      =[cosThetaM*cosPhiM, cosThetaM*sinPhiM, -sinThetaM;...
            %     -sinPhiM,           cosPhiM,       0.0 ];
            % InvCovM  = Mdelta'*InvCovUV*Mdelta;
            % xUnit    = [1.0; 0.0; 0.0];
            % yUnit    = [0.0; 1.0; 0.0];
            % zUnit    = [0.0; 0.0; 1.0];
            %
            % losDiff  = [xLos - xLosM; yLos - yLosM; zLos - zLosM];% line 132
            % ChiSqLos = 0.5*losDiff'*InvCovM*losDiff;              % line 133
            %
            % ParamList   = [range, thetaLos, phiLos, ldot, thdot, phidot, lx,ly,lz, Uview, Vview, Collinearity, ChiSqLos]';
            % this.ParamList = ParamList;
            % %ParamList = [range, thetaLos, phiLos, ldot, thdot, phidot, Uview, Vview]';
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % this.los       = los;                            % line 101-103
            % this.lx        = lx;                                 % line 101
            % this.ly        = ly;                                 % line 102
            % this.lz        = lz;                                 % line 103
            % this.lperpSq   = lperpSq;                            % line 104
            % this.rangeSq   = rangeSq;                            % line 105
            % this.range     = range;                              % line 106
            % this.zLos      = zLos;                               % line 107
            % this.thetaLos  = thetaLos;                           % line 108
            % this.phiLos    = phiLos;                             % line 109
            % this.vlos      = vlos;                           % line 110-112
            % this.vx        = vx;                                 % line 110
            % this.vy        = vy;                                 % line 111
            % this.vz        = vz;                                 % line 112
            % this.ldotv     = ldotv;                              % line 113
            % this.ldot      = ldot;                               % line 114
            % this.lperp     = lperp;                              % line 115
            % this.thdotfac1 = thdotfac1;                          % line 116
            % this.thdotfac2 = thdotfac2;                          % line 117
            % this.thdot     = thdot;                              % line 118
            % this.phiDotNum = phiDotNum;                          % line 119
            % this.phidot    = phidot;                             % line 120
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % this.thetaM      = thetaM;
            % this.cosThetaM   = cosThetaM;
            % this.sinThetaM   = sinThetaM;
            % this.phiM        = phiM;
            % this.cosPhiM     = cosPhiM;
            % this.sinPhiM     = sinPhiM;
            % this.cosThetaLos = cosThetaLos;                      % line 121
            % this.sinThetaLos = sinThetaLos;                      % line 122
            % this.phiDiff     = phiDiff;                          % line 123
            % this.cosPhiDiff  = cosPhiDiff;                       % line 124
            % this.sinPhiDiff  = sinPhiDiff;                       % line 125
            % this.delUV       = delUV;                            % line 126
            % this.Uview       = Uview;                            % line 127
            % this.Vview       = Vview;                            % line 128
            % this.xLos        = xLos;                             % line 129
            % this.yLos        = yLos;                             % line 130
            % this.xLosM       = xLosM;
            % this.yLosM       = yLosM;
            % this.zLosM       = zLosM;
            % this.Collinearity = Collinearity;                    % line 131
            % this.Mdelta      = Mdelta;
            % this.InvCovM     = InvCovM;
            % this.losDiff     = losDiff;                          % line 132
            % this.xUnit       = xUnit;
            % this.yUnit       = yUnit;
            % this.zUnit       = zUnit;
            % this.ChiSqLos    = ChiSqLos;                         % line 133
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Treat Equinoctial Variables:
            % Aeq      = omega + Omega;                            % line 142
            % cosAeq   = cos(Aeq);                                 % line 143
            % sinAeq   = sin(Aeq);                                 % line 144
            % tanHalf  = tan(I/2.0);                               % line 145
            % Feq      = e*cosAeq;                                 % line 146
            % Geq      = e*sinAeq;                                 % line 147
            % Heq      = tanHalf*cosO;                             % line 148
            % Keq      = tanHalf*sinO;                             % line 149
            % Leq      = Aeq + nu;                                 % line 150
            % CosL     = cos(Leq);                                 % line 151
            % SinL     = sin(Leq);                                 % line 152
            % alphaSq = Heq^2 - Keq^2;                             % line 153
            % Seq     = 1.0 + Heq^2 + Keq^2;                       % line 154
            % Weq     = 1.0 + Feq*CosL + Geq*SinL;                 % line 155
            % Req     = p/Weq;                                     % line 156
            % RovS    = Req/Seq;                                   % line 157
            % srtpinv = rtpinv/Seq;                                % line 158
            % HK      = Heq*Keq;                                   % line 159
            % OnePalphaSqCosL = (1.+alphaSq)*CosL;                 % line 160
            % OneMalphaSqCosL = (1.-alphaSq)*CosL;                 % line 161
            % OnePalphaSqSinL = (1.+alphaSq)*SinL;                 % line 162
            % OneMalphaSqSinL = (1.-alphaSq)*SinL;                 % line 163
            % Xfac    = OnePalphaSqCosL + 2.0*HK*SinL;             % line 164
            % Yfac    = OneMalphaSqSinL + 2.0*HK*CosL;             % line 165
            % Zfac    = Heq*SinL - Keq*CosL;                       % line 166
            % VXfac   =  OnePalphaSqSinL - 2.0*HK*(Feq + CosL) + Geq*(1. + alphaSq);% line 167
            % VYfac   = -OneMalphaSqCosL + 2.0*HK*(Geq + SinL) + Feq*(alphaSq - 1.);% line 168
            % VZfac   =  Heq*(Feq+CosL) + Keq*(Geq+SinL);          % line 169
            % Xeq     =     RovS*Xfac;                             % line 170
            % Yeq     =     RovS*Yfac;                             % line 171
            % Zeq     = 2.0*RovS*Zfac;                             % line 172
            % VXeq    =    -srtpinv*VXfac;                         % line 173
            % VYeq    =    -srtpinv*VYfac;                         % line 174
            % VZeq    = 2.0*srtpinv*VZfac;                         % line 175
            % this.EquinoctialElements = [p; Feq; Geq; Heq; Keq; Leq; t];
            %
            % RVeq    = [Xeq; Yeq; Zeq];
            % %Rpos'
            % VVeq     = [VXeq; VYeq; VZeq];
            % %Rdot'
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % this.Aeq             = Aeq;                      % line 142 Aeq
            % this.cosAeq          = cosAeq;                   % line 143 cosAeq
            % this.sinAeq          = sinAeq;                   % line 144 sinAeq
            % this.tanHalf         = tanHalf;                  % line 145 tanHalf
            % this.Feq             = Feq;                      % line 146 Feq
            % this.Geq             = Geq;                      % line 147 Geq
            % this.Heq             = Heq;                      % line 148 Heq
            % this.Keq             = Keq;                      % line 149 Keq
            % this.Leq             = Leq;                      % line 150 Leq
            % this.CosL            = CosL;                     % line 151 CosL
            % this.SinL            = SinL;                     % line 152 SinL
            % this.alphaSq         = alphaSq;                  % line 153 alphaSq
            % this.Seq             = Seq;                      % line 154
            % this.Weq             = Weq;                      % line 155
            % this.Req             = Req;                      % line 156
            % this.RovS            = RovS;                     % line 157
            % this.srtpinv         = srtpinv;                  % line 158
            % this.HK              = HK;                       % line 159
            % this.OnePalphaSqCosL = OnePalphaSqCosL;          % line 160
            % this.OneMalphaSqCosL = OneMalphaSqCosL;          % line 161
            % this.OnePalphaSqSinL = OnePalphaSqSinL;          % line 162
            % this.OneMalphaSqSinL = OneMalphaSqSinL;          % line 163
            % this.Xfac            = Xfac;                     % line 164
            % this.Yfac            = Yfac;                     % line 165
            % this.Zfac            = Zfac;                     % line 166
            % this.VXfac           = VXfac;                    % line 167
            % this.VYfac           = VYfac;                    % line 168
            % this.VZfac           = VZfac;                    % line 169
            % this.Xeq             = Xeq;                      % line 170
            % this.Yeq             = Yeq;                      % line 171
            % this.Zeq             = Zeq;                      % line 172
            % this.VXeq            = VXeq;                     % line 173
            % this.VYeq            = VYeq;                     % line 174
            % this.VZeq            = VZeq;                     % line 175
                       
            % DataSetList = SetDataSet(this);
            % this = DataSetList;
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Checking First Derivatives
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %Index = [1:56]';
            % %    for LineCheck = 76:81 %LineCheck  = 90
            % CombinedDifferences = 0.0;
            % format short
            % iStart = 1;
            % for LineCheck = iStart:175
            %     ForDerivs  = [];
            %     %BackDerivs = [];
            %     %Derivs     = [];
            %     %Combined   = [];
            %
            %     F            = zeros(175,1);
            %     F(LineCheck) = 1;
            %     F            = this.backdiff(F, LineCheck);
            %     % BackDerivs   = [F(1) F(2) F(3) F(4) F(5) F(6) F(7)];
            %     % Combined     = [Index F];
            %     % Kate's Combined'
            %
            %     forward      = [];
            %     for jk = 1:175 % check range of derivatives wrt forwards differentiation
            %         %jk
            %         Q        = zeros(175,1);
            %         Q(jk)     = 1;
            %         Q        = this.fordiff(Q, LineCheck);
            %         % Q(LineCheck) = dK_LineCheck/dLj
            %         %forward  = [forward; j Q(LineCheck)];
            %         forward  = [forward; jk+1-iStart Q(LineCheck)];
            %         %Derivs   = [Derivs; Q(LineCheck)];
            %         %Combined = [Combined Q];
            %     end
            %     ForDerivs    = [ForDerivs; forward];
            %     Derivs       = [ForDerivs F (F-ForDerivs(:,2))];
            %     Difference   = sum(F-ForDerivs(:,2));
            %     CombinedDifferences = CombinedDifferences + abs(Difference);
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
            % disp([' Combined Sum of Differences = ',num2str(CombinedDifferences)]);
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % %%% TEST SECOND DERIVATIVES:
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %Debug Second Derivatives by checking symmetry:
            % %Begin Checking Out Second Derivatives:
            % HessianDifferences = 0.0;
            % Nvar         = 7;
            % % LineCheck    = 56;  % Kepler Part
            % % LineCheck    = 90;  % Escobal Gravitation to J6
            % % F            = zeros(90,1);
            % % F(LineCheck) = 1;
            % % F = this.backdiff(F, LineCheck);
            % Nlines = 175;
            % format short
            % for LineCheck = 1:Nlines
            %     ListMax   = LineCheck + 1;
            %     %Nvar = LineCheck;
            %     Hessian      = [];
            %     F            = zeros(Nlines,1);
            %     F(LineCheck) = 1;
            %     F = this.backdiff(F, LineCheck);
            %     %F(ListMax:Nlines) = 0;
            %     % for j = 1:Nvar
            %     for j = 1:LineCheck
            %         Q = zeros(Nlines,1);
            %         S = zeros(Nlines,1);
            %         Q(j) = 1.0;
            %         Q = this.fordiff(Q, LineCheck);
            %         %Q(ListMax:Nlines) = 0;
            %         S = this.secdiff(F,Q,S, LineCheck);
            %         %S(ListMax:Nlines) = 0;
            %         S = this.backdiff(S, LineCheck);
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
            %         HessianTest = Hessian/maxElement
            %     end
            %     Symmetry = HessianTest - HessianTest'
            %     Symmetry = HessianTest - HessianTest';
            %     Differences = max(max(abs(Symmetry)));
            %     HessianDifferences = HessianDifferences + Differences;
            % % disp(['  Total Difference with Hessians = ',num2str(HessianDifferences)]);
            %
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
            %     max(max(Symmetry))
            %     disp([' Line Check = ', num2str(LineCheck),' max Hessian = ', num2str(maxElement),'  Summed Differences = ', num2str(Differences)]);
            %     %disp('Debug Point')
            % end
            % disp(['  Total Difference with Hessians = ',num2str(HessianDifferences)]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %// independent variable locater: Index: t = 1, rx = 2, ry = 3, rz = 4, vx = 5, vy = 6, vz = 7
            % ipoint    = zeros(7,1);
            % ipoint(1) = 1;
            % ipoint(2) = 2;
            % ipoint(3) = 3;
            % ipoint(4) = 4;
            % ipoint(5) = 5;
            % ipoint(6) = 6;
            % ipoint(7) = 40;
            % %Important Notation for F array indexes shown below.
            % %grads        = [F(1) F(2) F(3) F(4) F(5) F(6) F(40)];
            % e,a,...         t    e    a    I    w    W    M
            % Important Notation for F array indexes shown below.
            % %grads        = [F(1) F(2) F(3) F(4) F(5) F(6) F(7)];
            % e,a,...         t    e    a    I    w    W    Mp
            % %grads        = [F(1) F(2) F(15) F(4) F(5) F(6) F(7)];
            % e,p,...         t     e     p    I    w    W    Mp
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Get Jacobians J(Rx,Ry,Rz,Vx,Vy, Vz)
            %                (e,  a, I, w, W, Mp)
            
            % Hessian = cell(6, 6);
            % Jacobian = [];
            
            % %LineCheck    = 50;                % Rx
            % [H_I, Jacob, Hess] = WorkOrder(this, 50);
            % Hessian{1}  = Hess;
            % Jacobian      = [Jacobian; Jacob];   % row 1
            %
            % %LineCheck    = 51;                % Ry
            % [H_I, Jacob, Hess] = WorkOrder(this, 51);
            % Hessian{2}  = Hess;
            % Jacobian      = [Jacobian; Jacob];   % row 1
            %
            % %LineCheck    = 52;                % Rz
            % [H_I, Jacob, Hess] = WorkOrder(this, 52);
            % Hessian{3}  = Hess;
            % Jacobian      = [Jacobian; Jacob];   % row 1
            %
            % %LineCheck    = 54;                % Vx
            % [H_I, Jacob, Hess] = WorkOrder(this, 54);
            % Hessian{4}  = Hess;
            % Jacobian      = [Jacobian; Jacob];   % row 1
            %
            % %LineCheck    = 55;                % Vy
            % [H_I, Jacob, Hess] = WorkOrder(this, 55);
            % Hessian{5}  = Hess;
            % Jacobian      = [Jacobian; Jacob];   % row 1
            %
            % %LineCheck    = 56;                % Vz
            % [H_I, Jacob, Hess] = WorkOrder(this, 56);
            % Hessian{6}  = Hess;
            % Jacobian      = [Jacobian; Jacob];   % row 1
            %
            % this.JacKepler_2_ECI = Jacobian(:,2:7);
            %
            %this.Jacobian = Jacobian;
            %this.RVHessian = RVHessian;
            %Initialization of WorkVector for Backwards Differentiation:
                                    
            % Get Jacobians for J(los,theta_los,phi_los,ldot, theta_dot, phi_dot)
            %                          (e,  a, I, w, W, Mp)
            %JacobianLos = [];
            %HessianLos = cell(1, 6);
            
            % %LineCheck    = 106;  % Los Range
            % [H_I, Jacob, Hess] = WorkOrder(this, 106);
            % Jacobian   = [Jacobian; Jacob];   % row 1
            % Hessian{8} = Hess;
            %
            % %LineCheck    = 108;  % Los Theta
            % [H_I, Jacob, Hess] = WorkOrder(this, 108);
            % Jacobian   = [Jacobian; Jacob];   % row 2
            % Hessian{9} = Hess;
            %
            % %LineCheck    = 109;  % Los Phi
            % [H_I, Jacob, Hess] = WorkOrder(this, 109);
            % Jacobian   = [Jacobian; Jacob];   % row 3
            % Hessian{10} = Hess;
            %
            % %LineCheck    = 114;  % Los range-dot
            % [H_I, Jacob, Hess] = WorkOrder(this, 114);
            % Jacobian   = [Jacobian; Jacob];   % row 4
            % Hessian{11} = Hess;
            %
            % %LineCheck    = 118;  % Los Theta-dot
            % [H_I, Jacob, Hess] = WorkOrder(this, 118);
            % Jacobian   = [Jacobian; Jacob];   % row 5
            % Hessian{12} = Hess;
            %
            % %LineCheck    = 120;  % Los Phi-dot
            % [H_I, Jacob, Hess] = WorkOrder(this, 120);
            % Jacobian   = [Jacobian; Jacob];   % row 6
            % Hessian{13} = Hess;
            %
            % LineCheck    = 129;  % xLos
            % [H_I, Jacob, Hess] = WorkOrder(this, LineCheck);
            % Jacobian   = [Jacobian; Jacob];   % row 7
            % Hessian{14} = Hess;
            %
            % LineCheck    = 130;  % yLos Vview
            % [H_I, Jacob, Hess] = WorkOrder(this, LineCheck);
            % Jacobian   = [Jacobian; Jacob];   % row 8
            % Hessian{15} = Hess;
            %
            % LineCheck    = 107;  % zLos Vview
            % [H_I, Jacob, Hess] = WorkOrder(this, LineCheck);
            % Jacobian   = [Jacobian; Jacob];   % row 9
            % Hessian{16} = Hess;
            %
            % LineCheck    = 127;  % Los Uview
            % [H_I, Jacob, Hess] = WorkOrder(this, LineCheck);
            % Jacobian   = [Jacobian; Jacob];   % row 10
            % Hessian{17} = Hess;
            %
            % LineCheck    = 128;  % Los Vview
            % [H_I, Jacob, Hess] = WorkOrder(this, LineCheck);
            % Jacobian   = [Jacobian; Jacob];   % row 11
            % Hessian{18} = Hess;
            %
            % % The Equinoctials:
            %
            % LineCheck    = 15;  % p
            % [H_I, Jacob, Hess] = WorkOrder(this, LineCheck);
            % Jacobian   = [Jacobian; Jacob];   % row 11
            % Hessian{19} = Hess;
            %
            % LineCheck    = 146;  % f
            % [H_I, Jacob, Hess] = WorkOrder(this, LineCheck);
            % Jacobian   = [Jacobian; Jacob];   % row 11
            % Hessian{20} = Hess;
            %
            % LineCheck    = 147;  % g
            % [H_I, Jacob, Hess] = WorkOrder(this, LineCheck);
            % Jacobian   = [Jacobian; Jacob];   % row 11
            % Hessian{21} = Hess;
            %
            % LineCheck    = 148;  % h
            % [H_I, Jacob, Hess] = WorkOrder(this, LineCheck);
            % Jacobian   = [Jacobian; Jacob];   % row 11
            % Hessian{22} = Hess;
            %
            % LineCheck    = 149;  % k
            % [H_I, Jacob, Hess] = WorkOrder(this, LineCheck);
            % Jacobian   = [Jacobian; Jacob];   % row 11
            % Hessian{23} = Hess;
            %
            % LineCheck    = 150;  % L
            % [H_I, Jacob, Hess] = WorkOrder(this, LineCheck);
            % Jacobian   = [Jacobian; Jacob];   % row 11
            % Hessian{24} = Hess;
            %
            % JacKepler_2_Equinoctial = Jacobian(19:24,2:7);
            % this.JacKepler_2_Equinoctial = JacKepler_2_Equinoctial;
            
            % Aternative Jacobian Definition
            % Jac = [p, f, g, h, k, L]
            %        [e, a, I, omega, Omega, nu]
                    
            % JacKepler_2_Equinoctial = [];
            %
            % LineCheck    = 15;  %                   p
            % F            = zeros(175,1);
            % F(LineCheck) = 1;
            % F            = this.backdiff(F, LineCheck);
            % JacKepler_2_Equinoctial = [JacKepler_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(46)];
            %
            % LineCheck    = 146;  %                   f
            % F            = zeros(175,1);
            % F(LineCheck) = 1;
            % F            = this.backdiff(F, LineCheck);
            % JacKepler_2_Equinoctial = [JacKepler_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(46)];
            %
            % LineCheck    = 147;  %                   g
            % F            = zeros(175,1);
            % F(LineCheck) = 1;
            % F            = this.backdiff(F, LineCheck);
            % JacKepler_2_Equinoctial = [JacKepler_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(46)];
            %
            % LineCheck    = 148;  %                   h
            % F            = zeros(175,1);
            % F(LineCheck) = 1;
            % F            = this.backdiff(F, LineCheck);
            % JacKepler_2_Equinoctial = [JacKepler_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(46)];
            %
            % LineCheck    = 149;  %                   k
            % F            = zeros(175,1);
            % F(LineCheck) = 1;
            % F            = this.backdiff(F, LineCheck);
            % JacKepler_2_Equinoctial = [JacKepler_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(46)];
            %
            % LineCheck    = 150;  %                   L
            % F            = zeros(175,1);
            % F(LineCheck) = 1;
            % F            = this.backdiff(F, LineCheck);
            % JacKepler_2_Equinoctial = [JacKepler_2_Equinoctial; F(2) F(3) F(4) F(5) F(6) F(46)];
            % this.JacKepler_2_Equinoctial = JacKepler_2_Equinoctial;
            
            %this.Gradient    = Gradient;
            %this.JacobianLos = JacobianLos;
            %this.HessianLos  = HessianLos;
                        
  %      function [Zdot, Inhomogeneous, Perturbation, DMatrix, LagrangeTerms, LagrangeGrads] = EquationsOfMotion(this,Z)
            %Nvar      = 7;
            %meanMotion = this.meanMotion;
            %ipoint    = zeros(Nvar,1);
            %ipoint(1) = 1;
            %ipoint(2) = 2;
            %ipoint(3) = 3;
            %ipoint(4) = 4;
            %ipoint(5) = 5;
            %ipoint(6) = 6;
            %ipoint(7) = 40;
            
            LineCheck    = 90;  % Escobal Gravitation to J6
            F            = zeros(175,1);
            F(LineCheck) = 1;
            F            = backdiff(this, F, LineCheck);
            Jacob        = [F(1) F(2) F(3) F(4) F(5) F(6) F(7)];
            % e,a,...         t    e    a    I    w    W    Mp
            % this.GravityCan    = [F(61); F(62); F(63)];
            % GravityMKS  = Gravity
            %[J6 Jperturb Jhom] = J6Gravity(Rextrap);
            %J6Compare = Jperturb'; 
            %J6Compare - GravityCan
            %[H_I, Jacob, Hess] = WorkOrder(this, LineCheck);
            %Jacobian      = [Jacobian; Jacob];   % row 1
            %hess1     = Hess(2:7,2:7);
            %Hessian   = (hess1+hess1')/2.0;
            %Hessian{1}  = Hess;
            %this.Gradient  = Jacob(2:7)';
            %Gradient = this.Gradient;
            %MLagrange     = this.MLagrange;
            Inhomogeneous = this.MLagrange*Jacob(2:7)';   % The Zeroth Order Term
            %Gradient       = Jacob(2:7)';
            %Derivatives    = (Hess + Hess')/2.0;
            %Hessian        = Derivatives(2:7,2:7);
            %MLagrange*Gradient;   % Zeroth Order term

            %this.Jacobian    = Jacobian;
            %this.Hessian     = Hessian;            
            
            % % Important Notation for F array indexes shown below.
            % LineCheck = 96;   % Mwe derivatives
            % [H_I, Jacob, Hess] = WorkOrder(this, LineCheck);
            % gradsMwe  = Jacob(2:7)';
            %
            % LineCheck = 97;   % MeM derivatives
            % [H_I, Jacob, Hess] = WorkOrder(this, LineCheck);
            % gradsMeM  = Jacob(2:7)';
            %
            % LineCheck    = 98;   % MaM derivatives
            % [H_I, Jacob, Hess] = WorkOrder(this, LineCheck);
            % gradsMaM     = Jacob(2:7)';
            %
            % LineCheck    = 99;   % MWi derivatives
            % [H_I, Jacob, Hess] = WorkOrder(this, LineCheck);
            % gradsMWi     = Jacob(2:7)';
            %
            % LineCheck    = 100;   % Miw derivatives
            % [H_I, Jacob, ~] = WorkOrder(this, LineCheck);
            % gradsMiw     = Jacob(2:7)';
            %
            % LineCheck    = 16;   % meanMotion derivatives
            % [H_I, Jacob, Hess] = WorkOrder(this, LineCheck);
            % nDerivs      = Jacob(2:7)';
            % nMatrix      = zeros(6);
            % %nMatrix(6,:) = nDerivs;    % turn this gradient of n wrt elements off for now.
            %
            % DMatrix = cell(6,6);
            % for k = 1:6
            %     DMatrix{k}  = [       0,         0,             0, -gradsMwe(k),  0,        gradsMeM(k);...
            %         0,         0,             0,            0,  0,        gradsMaM(k);...
            %         0,         0,             0,  gradsMiw(k), -gradsMWi(k),        0;...
            %         gradsMwe(k),      0,  -gradsMiw(k),            0,  0,                  0;...
            %         0,         0,   gradsMWi(k),            0,  0,                  0;...
            %         -gradsMeM(k), -gradsMaM(k),       0,            0,  0,                  0];
            %
            % end
            % %DMatrix{6} = DMatrix{6} + nMatrix;
            % this.DMatrix = DMatrix;
            %
            % %Pertubation = zeros(6,1);
            % %for k = 1:6
            % %    Perturbation = Pertubation + DMatrix{k}*Gradient + MLagrange*Hessian(:,k);
            % %end
            %
            % PerturbationGrad = [];
            % for k = 1:6
            %     PerturbationGrad = [PerturbationGrad, DMatrix{k}*Gradient];
            % end
            % PerturbationHess = MLagrange*Hessian;
            % Perturbation     = PerturbationGrad + PerturbationHess;
            % Zdot             = Inhomogeneous + Perturbation*Z;
            %
            % Mwe = this.Mwe;
            % MeM = this.MeM;
            % MaM = this.MaM;
            % MWi = this.MWi;
            % Miw = this.Miw;
            % LagrangeTerms  = [Mwe; MeM; MWi; MaM; Miw];
            %
            % this.gradsMwe  = gradsMwe;
            % this.gradsMeM  = gradsMeM;
            % this.gradsMWi  = gradsMWi;
            % this.gradsMaM  = gradsMaM;
            % this.gradsMiw  = gradsMiw;
            % LagrangeGrads  = [gradsMwe, gradsMeM, gradsMWi, gradsMaM, gradsMiw];
            
            % this.Inhomogeneous = Inhomogeneous;
            % this.Perturbation  = Perturbation;
            
            %[V,D] = eig(Perturbation)
            %invPert = inv(Perturbation)
            %
            % Inhomogeneous
            % Perturbation
            %
            % Phi = expm(Perturbation*2)
            %
            % [V,D] = eig(Perturbation)
            % tProp = 1
            % PhiFundamental = [];
            % D = diag(D)
            % for k = 1:6
            %    PhiFundamental = [PhiFundamental, V(:,k)]
            % end
            % invPhi = inv(PhiFundamental)
            %
            % PhiFundamental = [];
            % for k = 1:6
            %    PhiFundamental = [PhiFundamental, V(:,k)*exp(tProp*D(k))]
            % end
            %
            % Phi = PhiFundamental*invPhi
            %
            % Initial = -invPhi*invPert*Inhomogeneous
            
%        end
        
        %
        % function [H_I, Jacob, Hess] = GravityDerivs(this, I)
        %
     end
        
        % function [ChiSq, DerivChiSq, HessianChiSq] = ChiSquare(this, DataList, Covariance, PhiFundamental, varargin)
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     DU       = 1.0;
        %     VU       = 1.0;
        %     AU       = 1.0;
        %     TU       = 1.0;
        %     mu       = this.mu;
        %
        %     range    = this.range;      % line 106
        %     thetaLos = this.thetaLos;   % line 108
        %     phiLos   = this.phiLos;     % line 109
        %
        %     ldot     = this.ldot;       % line 114
        %     thdot    = this.thdot;      % line 118
        %     phidot   = this.phidot      % line 120
        %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     DMatrix = cell(6,6);
        %     for k = 1:6
        %       DMatrix{k}  = [       0,         0,             0, -gradsMwe(k),  0,        gradsMeM(k);...
        %                             0,         0,             0,            0,  0,        gradsMaM(k);...
        %                             0,         0,             0,  gradsMiw(k), -gradsMWi(k),        0;...
        %                      gradsMwe(k),      0,  -gradsMiw(k),            0,  0,                  0;...
        %                             0,         0,   gradsMWi(k),            0,  0,                  0;...
        %                     -gradsMeM(k), -gradsMaM(k),       0,            0,  0,                  0];
        %
        %     end
        %     Prediction = this.
        %     PerturbationGrad = [];
        %     for k = 1:6
        %        PerturbationGrad = [PerturbationGrad, DMatrix{k}*Gradient];
        %     end
        % end
        
        
        %function F = backdiff(this, f, varargin)
        function F = backdiff(this, f, LineCheck)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % DU    = 1.0;
            % VU    = 1.0;
            % TU    = 1.0;
            mu    = this.mu;
            % if nargin >= 3
            %     TU = this.TU;
            %     DU = this.DU;
            %     VU = this.VU;
            %     AU = this.AU;
            % end
            %disp('Hello World');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            t       = this.time;                              % line  1
            e       = this.e;                                 % line  2
            %p       = this.p;                                 % line  3
            a       = this.a;                                 % line  3
            I       = this.Inclination;                       % line  4
            omega   = this.omega;                             % line  5
            Omega   = this.Omega;                             % line  6
            Mp      = this.Mp;                                % line  7
            
            onePe   = this.onePe;                             % line 11
            oneMe   = this.oneMe;                             % line 12
            fac     = this.fac;                               % line 13
            rootfac = this.rootfac;                           % line 14
            %a       = this.a;                                 % line 15
            p       = this.p;                                 % line 15
            meanMotion = this.meanMotion;                     % line 16
            
            cosI    = this.cosI;                              % line 17
            sinI    = this.sinI;                              % line 18
            cosom   = this.cosom;                             % line 19
            sinom   = this.sinom;                             % line 20
            cosO    = this.cosO;                              % line 21
            sinO    = this.sinO;                              % line 22
            
            Px      = this.Px;                                % line 23
            Py      = this.Py;                                % line 24
            Pz      = this.Pz;                                % line 25
            Qx      = this.Qx;                                % line 26
            Qy      = this.Qy;                                % line 27
            Qz      = this.Qz;                                % line 28
            
            M       = this.M;                                 % line 40
            E       = this.E;                                 % line 41
            cosK    = this.cosK;                              % line 42
            sinK    = this.sinK;                              % line 43
            tanX    = this.tanX;                              % line 44
            tanY    = this.tanY;                              % line 45
            nu      = this.nu;                                % line 46
            coss    = this.coss;                              % line 47
            sins    = this.sins;                              % line 48
            rorbit  = this.rorbit;                            % line 49
            rorbitx = this.rorbitx;                           % line 50
            rorbity = this.rorbity;                           % line 51
            rorbitz = this.rorbitz;                           % line 52
            rtpinv  = this.rtpinv;                            % line 53
            vorbitx = this.vorbitx;                           % line 54
            vorbity = this.vorbity;                           % line 55
            vorbitz = this.vorbitz;                           % line 56
            
            dE_dM   = this.dE_dM;
            dE_de   = this.dE_de;
            d2E_dMdM = this.d2E_dMdM;
            d2E_dMde = this.d2E_dMde;
            d2E_dedM = this.d2E_dedM;
            d2E_dede = this.d2E_dede;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if LineCheck > 60
                %J2 = this.J2;        % J2 perturbation constant
                %J3 = this.J3;
                %J4 = this.J4;
                %J5 = this.J5;
                %J6 = this.J6;
                % Truncate Gravity Model
                % J2     = 0.0;
                % J3     = 0.0;
                % J4     = 0.0;
                % J5     = 0.0;
                % J6     = 0.0;
                CK = this.CK;
                C2 = this.C2;
                C3 = this.C3;
                C4 = this.C4;
                C5 = this.C5;
                C6 = this.C6;
                
                % Escobal Gravity Terms
                X        = this.X;
                Y        = this.Y;
                Z        = this.Z;
                R2       = this.R2;
                R1       = this.R1;
                R3       = this.R3;
                R4       = this.R4;
                R5       = this.R5;
                R6       = this.R6;
                %sinomnu  = this.sinomnu;
                %cosomnu  = this.cosomnu;
                %sd1chk   = this.sd1chk;
                %cd1chk   = this.cd1chk;
                sd1      = this.sd1;
                sd2      = this.sd2;
                sd3      = this.sd3;
                sd4      = this.sd4;
                sd5      = this.sd5;
                sd6      = this.sd6;
                F2       = this.F2;
                F3       = this.F3;
                F4       = this.F4;
                F5       = this.F5;
                F6       = this.F6;
                V1       = this.V1;
                V2       = this.V2;
                V3       = this.V3;
                V4       = this.V4;
                V5       = this.V5;
                V6       = this.V6;
                V        = this.V;
                D        = this.D;
                Phi      = this.Phi;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fac1      = this.fac1;
                fac2      = this.fac2;
                fac3      = this.fac3;
                fac4      = this.fac4;
                fac5      = this.fac5;
                Mwe       = this.Mwe;
                MeM       = this.MeM;
                MaM       = this.MaM;
                MWi       = this.MWi;
                Miw       = this.Miw;
                %this.MLagrange = MLagrange;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if LineCheck > 90
                    if LineCheck > 100
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        los       = this.los;   % Rextrap - Rsensor
                        lx        = this.lx;                                 % line 101
                        ly        = this.ly;                                 % line 102
                        lz        = this.lz;                                 % line 103
                        lperpSq   = this.lperpSq;                            % line 104
                        rangeSq   = this.rangeSq;                            % line 105
                        range     = this.range;                              % line 106
                        zLos      = this.zLos;                               % line 107
                        thetaLos  = this.thetaLos;                           % line 108
                        phiLos    = this.phiLos;                             % line 109
                        vlos      = this.vlos;                          % lines 110-112
                        vx        = this.vx;                                 % line 110
                        vy        = this.vy;                                 % line 111
                        vz        = this.vz;                                 % line 112
                        ldotv     = this.ldotv;                              % line 113
                        ldot      = this.ldot;                               % line 114
                        lperp     = this.lperp;                              % line 115
                        thdotfac1 = this.thdotfac1;                          % line 116
                        thdotfac2 = this.thdotfac2;                          % line 117
                        thdot     = this.thdot;                              % line 118
                        phiDotNum = this.phiDotNum;                          % line 119
                        phidot    = this.phidot;                             % line 120
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        thetaM      = this.thetaM;
                        cosThetaM   = this.cosThetaM;
                        sinThetaM   = this.sinThetaM;
                        phiM        = this.phiM;
                        cosPhiM     = this.cosPhiM;
                        sinPhiM     = this.sinPhiM;
                        cosThetaLos = this.cosThetaLos;                      % line 121
                        sinThetaLos = this.sinThetaLos;                      % line 122
                        phiDiff     = this.phiDiff;                          % line 123
                        cosPhiDiff  = this.cosPhiDiff;                       % line 124
                        sinPhiDiff  = this.sinPhiDiff;                       % line 125
                        delUV       = this.delUV;                            % line 126
                        Uview       = this.Uview;                            % line 127
                        Vview       = this.Vview;                            % line 128
                        xLos        = this.xLos;                             % line 129
                        yLos        = this.yLos;                             % line 130
                        xLosM       = this.xLosM;
                        yLosM       = this.yLosM;
                        zLosM       = this.zLosM;
                        Collinearity  = this.Collinearity;                   % line 131
                        Mdelta      = this.Mdelta;
                        InvCovM     = this.InvCovM;
                        losDiff     = this.losDiff;                          % line 132
                        xUnit       = this.xUnit;
                        yUnit       = this.yUnit;
                        zUnit       = this.zUnit;
                        ChiSqLos    = this.ChiSqLos;                         % line 133
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if LineCheck > 133
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            Aeq             = this.Aeq;                      % line 142 Aeq
                            cosAeq          = this.cosAeq;                   % line 143 cosAeq
                            sinAeq          = this.sinAeq;                   % line 144 sinAeq
                            tanHalf         = this.tanHalf;                  % line 145 tanHalf
                            Feq             = this.Feq;                      % line 146 Feq
                            Geq             = this.Geq;                      % line 147 Geq
                            Heq             = this.Heq;                      % line 148 Heq
                            Keq             = this.Keq;                      % line 149 Keq
                            Leq             = this.Leq;                      % line 150 Leq
                            CosL            = this.CosL;                     % line 151 CosL
                            SinL            = this.SinL;                     % line 152 SinL
                            alphaSq         = this.alphaSq;                  % line 153 alphaSq
                            Seq             = this.Seq;                      % line 154
                            Weq             = this.Weq;                      % line 155
                            Req             = this.Req;                      % line 156
                            RovS            = this.RovS;                     % line 157
                            srtpinv         = this.srtpinv;                  % line 158
                            HK              = this.HK;                       % line 159
                            OnePalphaSqCosL = this.OnePalphaSqCosL;          % line 160
                            OneMalphaSqCosL = this.OneMalphaSqCosL;          % line 161
                            OnePalphaSqSinL = this.OnePalphaSqSinL;          % line 162
                            OneMalphaSqSinL = this.OneMalphaSqSinL;          % line 163
                            Xfac            = this.Xfac;                     % line 164
                            Yfac            = this.Yfac;                     % line 165
                            Zfac            = this.Zfac;                     % line 166
                            VXfac           = this.VXfac;                    % line 167
                            VYfac           = this.VYfac;                    % line 168
                            VZfac           = this.VZfac;                    % line 169
                            Xeq             = this.Xeq;                      % line 170
                            Yeq             = this.Yeq;                      % line 171
                            Zeq             = this.Zeq;                      % line 172
                            VXeq            = this.VXeq;                     % line 173
                            VYeq            = this.VYeq;                     % line 174
                            VZeq            = this.VZeq;                     % line 175
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % from line 175: VZeq    = 2.0*srtpinv*VZfac; % line 175
                            f(158) = f(158) + 2.0*f(175)*VZfac;
                            f(169) = f(169) + 2.0*f(175)*srtpinv;
                            % from line 174: VYeq    =    -srtpinv*VYfac; % line 174
                            f(158) = f(158) - f(174)*VYfac;
                            f(168) = f(168) - f(174)*srtpinv;
                            %from line 173: VXeq    =    -srtpinv*VXfac; % line 173
                            f(158) = f(158) - f(173)*VXfac;
                            f(167) = f(167) - f(173)*srtpinv;
                            %from line 172: Zeq      = 2.0*RovS*Zfac;     % line 172
                            f(157) = f(157) + 2.0*f(172)*Zfac;
                            f(166) = f(166) + 2.0*f(172)*RovS;
                            % from line 171: Yeq     =     RovS*Yfac;     % line 171
                            f(157) = f(157) + f(171)*Yfac;
                            f(165) = f(165) + f(171)*RovS;
                            % from line 170: Xeq     =     RovS*Xfac;     % line 170
                            f(157) = f(157) + f(170)*Xfac;
                            f(164) = f(164) + f(170)*RovS;
                            % from line 169: VZfac   =  Heq*(Feq + CosL) + Keq*(Geq + SinL) % line 169
                            f(146) = f(146) + f(169)*Heq;               % Feq
                            f(147) = f(147) + f(169)*Keq;               % Geq
                            f(148) = f(148) + f(169)*(Feq+CosL);        % Heq
                            f(149) = f(149) + f(169)*(Geq+SinL);        % Keq
                            f(151) = f(151) + f(169)*Heq;               % CosL
                            f(152) = f(152) + f(169)*Keq;               % SinL
                            
                            % from line 168: VYfac = -OneMalphaSqCosL + 2.0*HK*(Geq + SinL) + Feq*(alphaSq - 1.0);% line 168
                            f(146) = f(146) +     f(168)*(alphaSq - 1.0); % Feq
                            f(147) = f(147) + 2.0*f(168)*HK;              % Geq
                            f(152) = f(152) + 2.0*f(168)*HK;              % SinL
                            f(153) = f(153) +     f(168)*Feq;             % alphaSq
                            f(159) = f(159) + 2.0*f(168)*(Geq+SinL);      % HK
                            f(161) = f(161) -     f(168);                 % OneMAlphaSqSinL
                            
                            % from line 167: VXfac =  OnePalphaSqSinL - 2.0*HK*(Feq + CosL) + Geq*(1.0 + alphaSq);% line 107
                            f(146) = f(146) - 2.0*f(167)*HK;              % Feq
                            f(147) = f(147) +     f(167)*(1 + alphaSq);   % Geq
                            f(151) = f(151) - 2.0*f(167)*HK;              % CosL
                            f(153) = f(153) +     f(167)*Geq;             % alphaSq
                            f(159) = f(159) - 2.0*f(167)*(Feq+CosL);      % HK
                            f(162) = f(162) +     f(167);                 % OnePAlphaSqSinL
                            
                            % from line 166: Zfac    = Heq*SinL - Keq*CosL;          % line 106
                            f(148) = f(148) + f(166)*SinL;
                            f(149) = f(149) - f(166)*CosL;
                            f(151) = f(151) - f(166)*Keq;
                            f(152) = f(152) + f(166)*Heq;
                            % from line 165: Yfac    = OneMalphaSqSinL + 2.0*HK*CosL; % line 105
                            f(151) = f(151) + 2.0*f(165)*HK;
                            f(159) = f(159) + 2.0*f(165)*CosL;
                            f(163) = f(163) + f(165);
                            % from line 164: Xfac    = OnePalphaSqCosL + 2.0*HK*SinL; % line 104
                            f(152) = f(152) + 2.0*f(164)*HK;
                            f(159) = f(159) + 2.0*f(164)*SinL;
                            f(160) = f(160) + f(164);
                            % from line 163: OneMalphaSqSinL = (1.-alphaSq)*SinL;% line 103
                            f(152) = f(152) + f(163)*(1.0 - alphaSq);
                            f(153) = f(153) - f(163)*SinL;
                            % from line 162: OnePalphaSqSinL = (1.+alphaSq)*SinL;% line 102
                            f(152) = f(152) + f(162)*(1.0 + alphaSq);
                            f(153) = f(153) + f(162)*SinL;
                            % from line 161: OneMalphaSqCosL = (1.-alphaSq)*CosL;% line 101
                            f(151) = f(151) + f(161)*(1.0 - alphaSq);
                            f(153) = f(153) - f(161)*CosL;
                            % from line 160: OnePalphaSqCosL = (1.+alphaSq)*CosL;% line 100
                            f(151) = f(151) + f(160)*(1.0 + alphaSq);
                            f(153) = f(153) + f(160)*CosL;
                            % from line 159: HK      = Heq*Keq;            % line 99
                            f(148) = f(148) + f(159)*Keq;
                            f(149) = f(149) + f(159)*Heq;
                            % from line 158: srtpinv = rtpinv/Seq;         % line 98
                            f( 53) = f( 53) + f(158)/Seq;
                            f(154) = f(154) - f(158)*srtpinv/Seq;
                            % from line 157: RovS    = Req/Seq;            % line 97
                            f(156) = f(156) + f(157)/Seq;
                            f(154) = f(154) - f(157)*RovS/Seq;
                            % from line 156: Req     = p/Weq;              % line 96
                            f( 15) = f( 15) + f(156)/Weq;
                            f(155) = f(155) - f(156)*Req/Weq;
                            % from line 155: Weq = 1.0 + Feq*CosL + Geq*SinL; % line 95
                            f(146) = f(146) + f(155)*CosL;
                            f(147) = f(147) + f(155)*SinL;
                            f(151) = f(151) + f(155)*Feq;
                            f(152) = f(152) + f(155)*Geq;
                            % from line 154: Seq     = 1.0 + Heq^2 + Keq^2;   % line 94
                            f(148) = f(148) + 2.0*f(154)*Heq;
                            f(149) = f(149) + 2.0*f(154)*Keq;
                            % from line 153: alphaSq = Heq^2 - Keq^2;         % line 93
                            f(148) = f(148) + 2.0*f(153)*Heq;
                            f(149) = f(149) - 2.0*f(153)*Keq;
                            % from line 152: SinL     = sin(Leq);          % line 92
                            f(150) = f(150) + f(152)*CosL;
                            % from line 151: CosL     = cos(Leq);          % line 91
                            f(150) = f(150) - f(151)*SinL;
                            % from line 150: Leq      = Aeq + nu;          % line 90
                            f( 46) = f( 46) + f(150);
                            f(142) = f(142) + f(150);
                            % from line 149: Keq      = tanHalf*sin(Omega);    % line 89
                            f(145) = f(145) + f(149)*sinO;
                            f(  6) = f(  6) + f(149)*tanHalf*cosO;
                            % from line 148: Heq      = tanHalf*cos(Omega);    % line 88
                            f(145) = f(145) + f(148)*cosO;
                            f( 6) = f( 6) - f(148)*tanHalf*sinO;
                            % from line 147: Geq      = e*sinAeq;              % line 87
                            f( 2) = f( 2) + f(147)*sinAeq;
                            f(144) = f(144) + f(147)*e;
                            % from line 146: Feq      = e*cosAeq;              % line 86
                            f( 2) = f( 2) + f(146)*cosAeq;
                            f(143) = f(143) + f(146)*e;
                            % from line 145: tanHalf  = tan(Inclination/2.0);  % line 85
                            f( 4) = f( 4) + 0.5*f(145)*(1.0 + tanHalf^2);
                            % from line 144: sinAeq   = sin(Aeq);              % line 84
                            f(142) = f(142) + f(144)*cosAeq;
                            % from line 143: cosAeq   = cos(Aeq);              % line 83
                            f(142) = f(142) - f(143)*sinAeq;
                            % from line 142: Aeq      = omega + Omega;         % line 142
                            f( 5) = f( 5) + f(142);
                            f( 6) = f( 6) + f(142);
                        end
                        % Computation of backwards derivative starting at the bottom of the algorithm
                        %////////////  /////Block 81-75////////////////////////////////////////
                        % from line 133: ChiSqLos = 0.5*losDiff'*InvCovM*losDiff;
                        f(129)  = f(129) + f(133)*(xUnit'*InvCovM*losDiff);
                        f(130)  = f(130) + f(133)*(yUnit'*InvCovM*losDiff);
                        f(107)  = f(107) + f(133)*(zUnit'*InvCovM*losDiff);
                        
                        %losDiff  = [xLos - xLosM; yLos - yLosM; zLos - zLosM];% line 132
                        
                        %from line 131: Collinearity = 1-xLos*xLosM-yLos*yLosM-zLos*zLosM;
                        f(107)  = f(107) - f(131)*zLosM;
                        f(130)  = f(130) - f(131)*yLosM;
                        f(129)  = f(129) - f(131)*xLosM;
                        % from line 130: yLos     = ly/range;
                        f(102)  = f(102) + f(130)/range;
                        f(106)  = f(106) - f(130)*yLos/range;
                        % from line 129: xLos     = lx/range;
                        f(101)  = f(101) + f(129)/range;
                        f(106)  = f(106) - f(129)*xLos/range;
                        % from line 128: Vview       = sinThetaM*sinPhiDiff;
                        f(125) = f(125) + f(128)*sinThetaM;
                        % from line 127: Uview       = sinThetaLos*cosThetaM - delUV;
                        f(122)  = f(122) + f(127)*cosThetaM;
                        f(126)  = f(126) - f(127);
                        % from line 126: delUV       = cosThetaLos*sinThetaM*cosPhiDiff;
                        f(121)  = f(121) + f(126)*sinThetaM*cosPhiDiff;
                        f(124)  = f(124) + f(126)*sinThetaM*cosThetaLos;
                        % from line 125: sinPhiDiff  = sin(phiDiff);
                        f(123)  = f(123) + f(125)*cosPhiDiff;
                        % from line 124: cosPhiDiff  = cos(phiDiff);
                        f(123)  = f(123) - f(124)*sinPhiDiff;
                        % from line 123: phiDiff     = phiLos - phiM;
                        f(109)  = f(109) + f(123);
                        % from line 122: sinThetaLos = sin(thetaLos);
                        f(108)  = f(108) + f(122)*cosThetaLos;
                        f(108)  = f(108) - f(121)*sinThetaLos;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % from line 120: phidot    = phiDotNum/lperpSq;
                        f(119)  = f(119) + f(120)/lperpSq;
                        f(104)  = f(104) - f(120)*phidot/lperpSq;
                        % from line 119: phiDotNum = lx*vy - ly*vx;
                        f(101)  = f(101) + f(119)*vy;
                        f(102)  = f(102) - f(119)*vx;
                        f(110)  = f(110) - f(119)*ly;
                        f(111)  = f(111) + f(119)*lx;
                        % from line 118: thdot     = -thdotfac2/lperp;
                        f(115)  = f(115) - f(118)*thdot/lperp;
                        f(117)  = f(117) - f(118)/lperp;
                        % from line 117: thdotfac2 = thdotfac1/range;
                        f(116)  = f(116) + f(117)/range;
                        f(106)  = f(106) - f(117)*thdotfac2/range;
                        % from line 116: thdotfac1 = range*vz - lz*ldot;
                        f(103)  = f(103) - f(116)*ldot;    % lz
                        f(106)  = f(106) + f(116)*vz;      % range
                        f(112)  = f(112) + f(116)*range;   % vz
                        f(114)  = f(114) - f(116)*lz;      % ldot
                        % from line 115: lperp    = sqrt(lperpSq);
                        f(104)  = f(104) + 0.5*f(115)/lperp;
                        % from line 114: ldot     = ldotv/range;
                        f(106)  = f(106) - f(114)*ldot/range;
                        f(113)  = f(113) + f(114)/range;
                        % from line 113: ldotv    = lx*vx + ly*vy + lz*vz;
                        f(101)  = f(101) + f(113)*vx;
                        f(102)  = f(102) + f(113)*vy;
                        f(103)  = f(103) + f(113)*vz;
                        f(110)  = f(110) + f(113)*lx;
                        f(111)  = f(111) + f(113)*ly;
                        f(112)  = f(112) + f(113)*lz;
                        %vlos     = Vextrap - state_vector_sat(4:6);     % lines 110-112
                        % from line 112: vz  = vlos(3);
                        f( 56)  = f( 56) + f(112);
                        % from line 111: vy       = vlos(2);
                        f( 55)  = f( 55) + f(111);
                        % from line 110: vx       = vlos(1);
                        f( 54)  = f( 54) + f(110);
                        % from line 109: phiLos   = atan2(ly,lx);
                        f(101)   = f(101) - f(109)*ly/lperpSq;
                        f(102)   = f(102) + f(109)*lx/lperpSq;
                        % from line 108: thetaLos = acos(zLos);
                        f(107)   = f(107) - f(108)/sqrt(1.0 - zLos^2);
                        % from line 107: zLos     = lz/range;
                        f(103)  = f(103) + f(107)/range;
                        f(106)  = f(106) - f(107)*zLos/range;
                        % from line 106: range    = sqrt(rangeSq);
                        f(105)  = f(105) + 0.5*f(106)/range;
                        % from line 105: rangeSq  = lperpSq + lz^2;
                        f(104)  = f(104) + f(105);
                        f(103)  = f(103) + 2.0*f(105)*lz;
                        % from line 104: lperpSq  = lx^2 + ly^2;
                        f(101) = f(101) + 2.0*f(104)*lx;
                        f(102) = f(102) + 2.0*f(104)*ly;
                        % from line 103: lz        = this.lz;
                        f( 52) = f( 52) + f(103);
                        % from line 102: ly        = this.ly;
                        f( 51) = f( 51) + f(102);
                        % from line 101: lx        = this.lx;
                        f( 50) = f( 50) + f(101);
                    end
                    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % from line 100: Miw  = cosI*MWi      % cosI/(n*a^2*sqrt(1-e^2)*sinI) % line 100
                    f( 4) = f( 4) - f(100)*sinI*MWi;
                    f(99) = f(99) + f(100)*cosI;
                    % from line  99: MWi  = 1/fac5        % 1/(n*a^2*sqrt(1-e^2)*sinI)    % line  99
                    f(95) = f(95) - f(99)/fac5^2;
                    %from line   98: MaM  = 2.0/fac1      % 2/(n*a)                       % line  98
                    f(91) = f(91) - 2.0*f(98)/fac1^2;
                    % from line  97: MeM  = fac/fac3      % (1-e^2)/(n*a^2*e)             % line  97
                    f(13) = f(13) + f(97)/fac3;
                    f(93) = f(93) - f(97)*fac/fac3^2;
                    % from line  96: Mwe  = rootfac/fac3  % sqrt(1-e^2)/(n*a^2*e)         % line  96
                    f(14) = f(14) + f(96)/fac3;
                    f(93) = f(93) - f(96)*rootfac/fac3^2;
                    % from line  95: fac5 = fac4*sinI     % n*a^2*sqrt(1-e^2)*sinI        % line  95
                    f( 4) = f( 4) + f(95)*fac4*cosI;
                    f(94) = f(94) + f(95)*sinI;
                    % from line  94: fac4 = fac2*rootfac  % n*a^2*sqrt(1-e^2)             % line  94
                    f(14) = f(14) + f(94)*fac2;
                    f(92) = f(92) + f(94)*rootfac;
                    % from line  93: fac3 = fac2*e        % n*a^2*e                       % line  93
                    f( 2) = f( 2) + f(93)*fac2;
                    f(92) = f(92) + f(93)*e;
                    % from line 92: fac2 = fac1*a        % n*a^2                         % line  92
                    f( 3) = f( 3) + f(92)*fac1;
                    f(91) = f(91) + f(92)*a;
                    % from line 91: fac1 = meanMotion*a           % n*a                           % line  91
                    f( 3) = f( 3) + f(91)*meanMotion;
                    f(16) = f(16) + f(91)*a;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line 90: Phi     = V*(CK + D);
                f(87) = f(87) + f(90)*(CK + D);
                f(88) = f(88) + f(90)*V;
                % R is the disturbing potential
                % from line 89: R       = mu*D;
                % f(88) = f(88) + f(89)*mu;
                % Up to J6:
                % from line 88: D       = V2 + V3 + V4 + V5 + V6;
                f(82) = f(82) + f(88);
                f(83) = f(83) + f(88);
                f(84) = f(84) + f(88);
                f(85) = f(85) + f(88);
                f(86) = f(86) + f(88);
                % from line 87: V       = mu*V1;
                f(81) = f(81) + f(87)*mu;
                % from line 86: V6      = C6*F6/R6;
                f(80) = f(80) + f(86)*C6/R6;
                f(69) = f(69) - f(86)*V6/R6;
                % from line 85: V5      = C5*F5/R5;
                f(79) = f(79) + f(85)*C5/R5;
                f(68) = f(68) - f(85)*V5/R5;
                % from line 84: V4      = C4*F4/R4;
                f(78) = f(78) + f(84)*C4/R4;
                f(67) = f(67) - f(84)*V4/R4;
                % from line 83: V3      = C3*F3/R3;
                f(77) = f(77) + f(83)*C3/R3;
                f(66) = f(66) - f(83)*V3/R3;
                % from line 82: V2      = C2*F2/R2;
                f(76) = f(76) + f(82)*C2/R2;
                f(64) = f(64) - f(82)*V2/R2;
                % from line 81: V1      = 1/R1;
                f(65) = f(65) - f(81)*V1/R1;
                % from line 80: F6   =  5.       - 105.*sd2  + 315*sd4 - 231*sd6;
                f(71) = f(71) - f(80)*105.;
                f(73) = f(73) + f(80)*315.;
                f(75) = f(75) - f(80)*231.;
                % from line 79: F5   =  15.*sd1  -  70.*sd3  + 63*sd5;
                f(70) = f(70) + f(79)*15.;
                f(72) = f(72) - f(79)*70.;
                f(74) = f(74) + f(79)*63.;
                % from line 78: F4   =  3.0      -  30.*sd2  + 35*sd4;
                f(71) = f(71) - f(78)*30.;
                f(73) = f(73) + f(78)*35.;
                % from line 77: F3   =  3.0*sd1  -  5.0*sd3;
                f(70) = f(70) + f(77)*3.;
                f(72) = f(72) - f(77)*5.;
                % from line 76: F2   =  1.0      -  3.0*sd2;
                f(71) = f(71) - f(76)*3.;
                % from line 75: sd6  = sd1*sd5;
                f(70) = f(70) + f(75)*sd5;
                f(74) = f(74) + f(75)*sd1;
                % from line 74: sd5  = sd1*sd4;
                f(70) = f(70) + f(74)*sd4;
                f(73) = f(73) + f(74)*sd1;
                % from line 73: sd4  = sd1*sd3;
                f(70) = f(70) + f(73)*sd3;
                f(72) = f(72) + f(73)*sd1;
                % from line 72: sd3  = sd1*sd2;
                f(70) = f(70) + f(72)*sd2;
                f(71) = f(71) + f(72)*sd1;
                % from line 71: sd2  = sd1^2;
                f(70) = f(70) + f(71)*2.0*sd1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line 70: sd1  = Z/R1; (63 + 65)
                f(63) = f(63) + f(70)/R1;
                f(65) = f(65) - f(70)*sd1/R1;
                % from line 70: sd1  = sinI*sin(omega + nu);
                %f(18) = f(18) + f(70)*sinomnu;
                %f( 5) = f( 5) + f(70)*sinI*cosomnu;
                %f(46) = f(46) + f(70)*sinI*cosomnu;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line  69: R6   = R1*R5;
                f(65) = f(65) + f(69)*R5;
                f(68) = f(68) + f(69)*R1;
                % from line  68: R5   = R1*R4;
                f(65) = f(65) + f(68)*R4;
                f(67) = f(67) + f(68)*R1;
                % from line  67: R4   = R1*R3;
                f(65) = f(65) + f(67)*R3;
                f(66) = f(66) + f(67)*R1;
                % from line  66: R3   = R1*R2;
                f(65) = f(65) + f(66)*R2;
                f(64) = f(64) + f(66)*R1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line  64: R2   = R1^2;     --  keep original line numbers
                %f(65) = f(65) + 2.0*f(64)*R1;
                % from line  65: R1   = rorbit;    --      keep the line numbers
                %f(49) = f(49) + f(65);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line  65: R1   = sqrt(R2);
                f(64) = f(64) + f(65)*0.5/R1;
                % from line  64: R2   = X^2 + Y^2 + Z^2;
                f(61) = f(61) + f(64)*2.0*X;
                f(62) = f(62) + f(64)*2.0*Y;
                f(63) = f(63) + f(64)*2.0*Z;
                %
                % from line 63: Z    = rorbitz;
                f(52) = f(52) + f(63);
                % from line 62: Y    = rorbity;
                f(51) = f(51) + f(62);
                % from line 61: X    = rorbitx;
                f(50) = f(50) + f(61);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            %//////////////////////////////////////////////////////////////////////
            % Computation of backwards derivative starting at the bottom of the algorithm
            %////////////  /////Block 56-50////////////////////////////////////////
            % from line 56: double vorbitz  = sqrt(1.0/p)*(-sins*Pz + (e + coss)*Qz);
            f(53) = f(53) + f(56)*(-sins*Pz + (e+coss)*Qz);  % rtpinv #53
            f(48) = f(48) - f(56)*rtpinv*Pz;                 % sins 48
            f(47) = f(47) + f(56)*rtpinv*Qz;                 % coss 47
            f(28) = f(28) + f(56)*rtpinv*(e + coss);         % Qz   28
            % skip over line 53 and go all the way back to p -- line 3
            %f( 3) = f( 3) - f(56)*rtpinv*(0.5/p)*(-sins*Pz + (e+coss)*Qz); % p #3
            f(25) = f(25) - f(56)*rtpinv*sins;               % Pz   25
            f( 2) = f( 2) + f(56)*rtpinv*Qz;                 % e     2
            % from line 55: double vorbity  = sqrt(1.0/p)*(-sins*Py + (e +coss)*Qy);
            f(53) = f(53) + f(55)*(-sins*Py + (e + coss)*Qy);% rtpinv 53
            f(48) = f(48) - f(55)*rtpinv*Py;                 % sins  48
            f(47) = f(47) + f(55)*rtpinv*Qy;                 % coss  47
            f(27) = f(27) + f(55)*rtpinv*(e + coss);         % Qy    27
            % skip over line 53 and go all the way back to p -- line 3
            %f( 3) = f( 3) - f(55)*rtpinv*(0.5/p)*(-sins*Py + (e + coss)*Qy); % p 3
            f(24) = f(24) - f(55)*rtpinv*sins;               % Py    24
            f( 2) = f( 2) + f(55)*rtpinv*Qy;                 % e      2
            % from line 54: double vorbitx  = sqrt(1.0/p)*(-sins*Px + (e + coss)*Qx);
            f(53) = f(53) + f(54)*(-sins*Px + (e + coss)*Qx); % rtpinv 53
            f(48) = f(48) - f(54)*rtpinv*Px;                 % sins  48
            f(47) = f(47) + f(54)*rtpinv*Qx;                 % coss  47
            f(26) = f(26) + f(54)*rtpinv*(e + coss);         % Qx    26
            % skip over line 53 and go all the way back to p -- line 3
            %f( 3) = f( 3) - f(54)*rtpinv*(0.5/p)*(-sins*Px + (e + coss)*Qx); % p 3
            f(23) = f(23) - f(54)*rtpinv*sins;               % Px    23
            f( 2) = f( 2) + f(54)*rtpinv*Qx;                 % e      2
            % Don't skip over line 53:
            % from line 53: rtpinv   = sqrt(1.0/p);
            %f( 3) = f(3) - f(53)*rtpinv*(0.5/p);             % p      3
            f(15) = f(15) - f(53)*rtpinv*(0.5/p);             % p      15
            % skip over line 53 and go all the way back to p -- line 3
            % from line 52: double rorbitz  = rorbit*(coss*Pz + sins*Qz);
            f(49) = f(49) + f(52)*(coss*Pz + sins*Qz);       % rorbit  49
            f(48) = f(48) + f(52)*rorbit*Qz;                 % sins    48
            f(47) = f(47) + f(52)*rorbit*Pz;                 % coss    47
            f(28) = f(28) + f(52)*rorbit*sins;               % Qz      28
            f(25) = f(25) + f(52)*rorbit*coss;               % Pz      25
            % from line 51: double rorbity  = rorbit*(coss*Py + sins*Qy);
            f(49) = f(49) + f(51)*(coss*Py + sins*Qy);       % rorbit  49
            f(48) = f(48) + f(51)*rorbit*Qy;                 % sins    48
            f(47) = f(47) + f(51)*rorbit*Py;                 % coss    47
            f(27) = f(27) + f(51)*rorbit*sins;               % Qy      27
            f(24) = f(24) + f(51)*rorbit*coss;               % Py      24
            % from line 50: double rorbitx  = rorbit*(coss*Px + sins*Qx);
            f(49) = f(49) + f(50)*(coss*Px + sins*Qx);       % rorbit  49
            f(48) = f(48) + f(50)*rorbit*Qx;                 % sins    48
            f(47) = f(47) + f(50)*rorbit*Px;                 % coss    47
            f(26) = f(26) + f(50)*rorbit*sins;               % Qx      26
            f(23) = f(23) + f(50)*rorbit*coss;               % Px      23
            %             % from line 49: double rorbit   = p/(1.0 + e*coss);
            %             denom = 1.0/(1.0 + e*coss);
            %             f(73) = f(73) - f(75)*p*e*denom*denom;         % coss  #73
            %             f(24) = f(24) + f(75)*denom;                   % p     #24
            %             f(17) = f(17) - f(75)*p*coss*denom*denom;      % e     #17
            % from line 49: rorbit   = a*(1.0 - e*cosK);     % line 49
            denom = 1.0 - e*cosK;
            f(42) = f(42) - f(49)*a*e;                         % cosK   42
            %f(41) = f(41) + f(49)*a*e*sinK;                     % cos(E) 41
            f( 3) = f( 3) + f(49)*denom;                        %  a      3
            f( 2) = f( 2) - f(49)*a*cosK;                       %  e      2
            %//////////////////////////////////////////////////////////////////////
            %     Block -48
            % from line 48: sins  = sin(nu);
            f(46) = f(46) + f(48)*coss;                      % nu     46
            % from line: 47: coss  = cos(nu);
            f(46) = f(46) - f(47)*sins;                      % nu     46
            %             %//////////////////////////////////////////////////////////////////////
            %             % from line: 46: nu     = atan2(tanY, tanX);
            %             % JRS: Comment out following two lines to take
            %             % (nu, rx, ry, rx, vx, vy, vz) as independent variables
            f(45) = f(45) + f(46)*tanX/(tanX*tanX + tanY*tanY);% tanY  45
            f(44) = f(44) - f(46)*tanY/(tanX*tanX + tanY*tanY);% tanX  44
            %/////////////////////Block 71-65//////////////////////////////////////
            % from line 71: sins  = rootfac*sinK/kfac;
            %  f(69) = f(69) - f(71)*rootfac*sinK/(kfac*kfac);       % kfac    #69
            %  f(68) = f(68) + f(71)*rootfac/kfac;                   % sinK    #68
            %  f(42) = f(42) + f(71)*sinK/kfac;                      % rootfac #42
            % from line 70: coss = (cosK - e)/kfac;
            %  f(69) = f(69) - f(70)*(cosK - e)/(kfac*kfac);         % kfac    #69
            %  f(67) = f(67) + f(70)/kfac;                           % cosK    #67
            %  f(17) = f(17) - f(70)/kfac;                           % e       #17
            % from line 45: tanY  = rootfac*sinK;
            f(43) = f(43) + f(45)*rootfac;                    % sinK    43
            f(14) = f(14) + f(45)*sinK;                       % rootfac 14
            % from line 44: tanX  = cosK - e;
            f(42) = f(42) + f(44);                            % cosK    42
            f( 2) = f( 2) - f(44);                            % e        2
            % from line 69: kfac  = 1.0 - e*cosK;
            %  f(67) = f(67) - f(69)*e;                              % cosK    #67
            %  f(17) = f(17) - f(69)*cosK;                           % e       #17
            % from line 43: sinK  = sin(E);
            f(41) = f(41) + f(43)*cosK;                       % E       41
            % from line 42: cosK  = cos(E);
            f(41) = f(41) - f(42)*sinK;                       % E       41
            % from line 41: Implicitly E = E(e, meanAnomaly)
            % eccentricity line 17, meanAnomaly line 65
            f(40) = f(40) + f(41)*dE_dM;                      % M       40
            f( 2) = f( 2) + f(41)*dE_de;                      % e        2
            % from line 40: M  = meanMotion*t + Mp;      % line 40
            f(16) = f(16) + f(40)*t;                     % meanMotion   16
            f( 7) = f( 7) + f(40);                       % Mp            7
            f( 1) = f( 1) + f(40)*meanMotion;            % t             1
            %/////////////////////Block 31-23//////////////////////////////////////
            % from line 31: Wz   =  cosI;
            %f(17) = f(17) + f(31);                       % cosI         17
            f( 4) = f( 4) - f(31)*sinI;                   %    I          4
            % from line 30: Wy   = -cosO*sinI;
            f(21) = f(21) - f(30)*sinI;                  % cosO         21
            f(18) = f(18) - f(30)*cosO;                  % sinI         18
            % from line 29: Wx   =  sinO*sinI;
            f(22) = f(22) + f(29)*sinI;                 % sinO         22
            f(18) = f(18) + f(29)*sinO;                 % sinI         18
            % from line 28: Qz   =  cosom*sinI;          % line 28
            f(19) = f(19) + f(28)*sinI;                 % cosom        19
            f(18) = f(18) + f(28)*cosom;                % sinI         18
            % from line 27: Qy   = -sinO*sinom + cosO*cosom*cosI;
            f(22) = f(22) - f(27)*sinom;                % sinO         22
            f(21) = f(21) + f(27)*cosom*cosI;           % cosO         21
            f(20) = f(20) - f(27)*sinO;                 % sinom        20
            f(19) = f(19) + f(27)*cosO*cosI;            % cosom        19
            f(17) = f(17) + f(27)*cosO*cosom;           % cosI         17
            % from line 26: Qx   = -cosO*sinom - sinO*cosom*cosI;
            f(22) = f(22) - f(26)*cosom*cosI;           % sinO         22
            f(21) = f(21) - f(26)*sinom;                % cosO         21
            f(20) = f(20) - f(26)*cosO;                 % sinom        20
            f(19) = f(19) - f(26)*sinO*cosI;            % cosom        19
            f(17) = f(17) - f(26)*sinO*cosom;           % cosI         17
            % from line 25: Pz   =  sinom*sinI;
            f(20) = f(20) + f(25)*sinI;                 % sinom        20
            f(18) = f(18) + f(25)*sinom;                % sinI         18
            % from line 24: Py   =  sinO*cosom + cosO*sinom*cosI;
            f(22) = f(22) + f(24)*cosom;                % sinO         22
            f(21) = f(21) + f(24)*sinom*cosI;           % cosO         21
            f(20) = f(20) + f(24)*cosO*cosI;            % sinom        20
            f(19) = f(19) + f(24)*sinO;                 % cosom        19
            f(17) = f(17) + f(24)*cosO*sinom;           % cosI         17
            % from line 23: Px   =  cosO*cosom - sinO*sinom*cosI;
            f(22) = f(22) - f(23)*sinom*cosI;           % sinO         22
            f(21) = f(21) + f(23)*cosom;                % cosO         21
            f(20) = f(20) - f(23)*sinO*cosI;            % sinom        20
            f(19) = f(19) + f(23)*cosO;                 % cosom        19
            f(17) = f(17) - f(23)*sinO*sinom;           % cosI         17
            %//////////////////////////////////////////////////////////////////////
            %/////////////////////Block 22-17//////////////////////////////////////
            % from line 22: sinO            = sin(this.Omega);
            f( 6) = f( 6) + f(22)*cosO;                 % Omega         6
            % from line 21: cosO            = cos(this.Omega);
            f( 6) = f( 6) - f(21)*sinO;                 % Omega         6
            %from line 20: sinom            = sin(this.omega);
            f( 5) = f( 5) + f(20)*cosom;                % omega         5
            %from line 19: cosom            = cos(this.omega);
            f( 5) = f( 5) - f(19)*sinom;                % omega         5
            % from line 18: sinI            = sin(this.Inclination);
            f( 4) = f( 4) + f(18)*cosI;                 % I             4
            % from line 17: cosI            = cos(this.Inclination);
            f( 4) = f( 4) - f(17)*sinI;                  % I            4
            %//////////////////////////////////////////////////////////////////////
            %/////////////////////Block 16-11//////////////////////////////////////
            % from line 16: meanMotion = sqrt(mu)*a^(-1.5);    % Radians per Second
            %f(15) = f(15) - f(16)*1.5*meanMotion/a;            % Period  #45
            f( 3) = f( 3) - f(16)*1.5*meanMotion/a;            % Period  #45
            % Semi-Major Axis:
            % % from line 15: a = p/fac;
            % f(13) = f(13) - f(15)*a/fac;                      % fac     13
            % f( 3) = f( 3) + f(15)/fac;                        % p        3
            % from line 15: p = a*fac;
            f(13) = f(13) + f(15)*a;                          % fac     13
            f( 3) = f( 3) + f(15)*fac;                        % a        3
            %disp(' backdiff fac')
            %fac
            % from line 14: rootfac = sqrt(fac);
            if (fac > 0.0)
                %    % from line 14: rootfac = sqrt(fac);
                f(13) = f(13) + f(14)*0.5/rootfac;             % fac     13
            else
                rootfac = 0;
            end
            %//////////////////////////////////////////////////////////////////////
            % % from line 13: fac   = onePe*oneMe;
            % f(12) = f(12) + f(13)*onePe;                      % oneMe   12
            % f(11) = f(11) + f(13)*oneMe;                      % onePe   11
            % % from line 12: oneMe = 1.0 - e;
            % f( 2) = f( 2) - f(12);                            % e        2
            % % from line 11: onePe = 1.0 + e;
            % f( 2) = f( 2) + f(11);                            % e        2
            % from line 13: fac   = 1 - e^2;
            f(2) = f(2) - f(13)*2.0*e;
            %//////////////////////////////////////////////////////////////////////
            % Mp    = TrackParams(6);                       % Mp       7
            % Omega = TrackParams(5);                       % Omega    6
            % omega = TrackParams(4);                       % omega    5
            % I     = TrackParams(3);                       % I        4
            % a     = TrackParams(2);                       % a        3
            % e     = TrackParams(1);                       % e        2
            % t     = time;                                 % t        1
            F = f;
        end
        
        %function Q = fordiff(this, q, varargin)
        function Q = fordiff(this, q, LineCheck)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % DU    = 1.0;
            % VU    = 1.0;
            % TU    = 1.0;
            mu    = this.mu;
            % if nargin >= 3
            %     TU = this.TU;
            %     DU = this.DU;
            %     VU = this.VU;
            %     AU = this.AU;
            % end
            %disp('Hello World');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            t       = this.time;                              % line  1
            e       = this.e;                                 % line  2
            a       = this.a;                                 % line 3
            %p       = this.p;                                % line  3
            I       = this.Inclination;                       % line  4
            omega   = this.omega;                             % line  5
            Omega   = this.Omega;                             % line  6
            Mp      = this.Mp;                                % line  7
            
            onePe   = this.onePe;                             % line 11
            oneMe   = this.oneMe;                             % line 12
            fac     = this.fac;                               % line 13
            rootfac = this.rootfac;                           % line 14
            %a       = this.a;                                % line 15
            p       = this.p;                                % line  15
            meanMotion = this.meanMotion;                     % line 16
            
            cosI    = this.cosI;                              % line 17
            sinI    = this.sinI;                              % line 18
            cosom   = this.cosom;                             % line 19
            sinom   = this.sinom;                             % line 20
            cosO    = this.cosO;                              % line 21
            sinO    = this.sinO;                              % line 22
            
            Px      = this.Px;                                % line 23
            Py      = this.Py;                                % line 24
            Pz      = this.Pz;                                % line 25
            Qx      = this.Qx;                                % line 26
            Qy      = this.Qy;                                % line 27
            Qz      = this.Qz;                                % line 28
            
            M       = this.M;                                    % line 40
            E       = this.E;                                    % line 41
            cosK    = this.cosK;                                 % line 42
            sinK    = this.sinK;                                 % line 43
            tanX    = this.tanX;                                 % line 44
            tanY    = this.tanY;                                 % line 45
            nu      = this.nu;                                   % line 46
            coss    = this.coss;                                 % line 47
            sins    = this.sins;                                 % line 48
            rorbit  = this.rorbit;                               % line 49
            rorbitx = this.rorbitx;                              % line 50
            rorbity = this.rorbity;                              % line 51
            rorbitz = this.rorbitz;                              % line 52
            rtpinv  = this.rtpinv;                               % line 53
            vorbitx = this.vorbitx;                              % line 54
            vorbity = this.vorbity;                              % line 55
            vorbitz = this.vorbitz;                              % line 56
            
            dE_dM   = this.dE_dM;
            dE_de   = this.dE_de;
            d2E_dMdM = this.d2E_dMdM;
            d2E_dMde = this.d2E_dMde;
            d2E_dedM = this.d2E_dedM;
            d2E_dede = this.d2E_dede;
            
            %J2 = this.J2;        % J2 perturbation constant
            %J3 = this.J3;
            %J4 = this.J4;
            %J5 = this.J5;
            %J6 = this.J6;
            % Truncate Gravity Model
            % J2     = 0.0;
            % J3     = 0.0;
            % J4     = 0.0;
            % J5     = 0.0;
            % J6     = 0.0;
            CK = this.CK;
            C2 = this.C2;
            C3 = this.C3;
            C4 = this.C4;
            C5 = this.C5;
            C6 = this.C6;
            
            % Escobal Gravity Terms
            X        = this.X;
            Y        = this.Y;
            Z        = this.Z;
            R2       = this.R2;
            R1       = this.R1;
            R3       = this.R3;
            R4       = this.R4;
            R5       = this.R5;
            R6       = this.R6;
            %sinomnu  = this.sinomnu;
            %cosomnu  = this.cosomnu;
            sd1      = this.sd1;
            sd2      = this.sd2;
            sd3      = this.sd3;
            sd4      = this.sd4;
            sd5      = this.sd5;
            sd6      = this.sd6;
            F2       = this.F2;
            F3       = this.F3;
            F4       = this.F4;
            F5       = this.F5;
            F6       = this.F6;
            V1       = this.V1;
            V2       = this.V2;
            V3       = this.V3;
            V4       = this.V4;
            V5       = this.V5;
            V6       = this.V6;
            V        = this.V;
            D        = this.D;
            Phi      = this.Phi;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fac1      = this.fac1;
            fac2      = this.fac2;
            fac3      = this.fac3;
            fac4      = this.fac4;
            fac5      = this.fac5;
            Mwe       = this.Mwe;
            MeM       = this.MeM;
            MaM       = this.MaM;
            MWi       = this.MWi;
            Miw       = this.Miw;
            %this.MLagrange = MLagrange;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            los       = this.los;   % Rextrap - Rsensor
            lx        = this.lx;                                 % line 101
            ly        = this.ly;                                 % line 102
            lz        = this.lz;                                 % line 103
            lperpSq   = this.lperpSq;                            % line 104
            rangeSq   = this.rangeSq;                            % line 105
            range     = this.range;                              % line 106
            zLos      = this.zLos;                               % line 107
            thetaLos  = this.thetaLos;                           % line 108
            phiLos    = this.phiLos;                             % line 109
            vlos      = this.vlos;                          % lines 110-112
            vx        = this.vx;                                 % line 110
            vy        = this.vy;                                 % line 111
            vz        = this.vz;                                 % line 112
            ldotv     = this.ldotv;                              % line 113
            ldot      = this.ldot;                               % line 114
            lperp     = this.lperp;                              % line 115
            thdotfac1 = this.thdotfac1;                          % line 116
            thdotfac2 = this.thdotfac2;                          % line 117
            thdot     = this.thdot;                              % line 118
            phiDotNum = this.phiDotNum;                          % line 119
            phidot    = this.phidot;                             % line 120
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            thetaM      = this.thetaM;
            cosThetaM   = this.cosThetaM;
            sinThetaM   = this.sinThetaM;
            phiM        = this.phiM;
            cosPhiM     = this.cosPhiM;
            sinPhiM     = this.sinPhiM;
            cosThetaLos = this.cosThetaLos;                      % line 121
            sinThetaLos = this.sinThetaLos;                      % line 122
            phiDiff     = this.phiDiff;                          % line 123
            cosPhiDiff  = this.cosPhiDiff;                       % line 124
            sinPhiDiff  = this.sinPhiDiff;                       % line 125
            delUV       = this.delUV;                            % line 126
            Uview       = this.Uview;                            % line 127
            Vview       = this.Vview;                            % line 128
            xLos        = this.xLos;                             % line 129
            yLos        = this.yLos;                             % line 130
            xLosM       = this.xLosM;
            yLosM       = this.yLosM;
            zLosM       = this.zLosM;
            Collinearity  = this.Collinearity;                   % line 131
            Mdelta      = this.Mdelta;
            InvCovM     = this.InvCovM;
            losDiff     = this.losDiff;                          % line 132
            xUnit       = this.xUnit;
            yUnit       = this.yUnit;
            zUnit       = this.zUnit;
            ChiSqLos    = this.ChiSqLos;                         % line 133
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Aeq             = this.Aeq;                      % line 142 Aeq
            cosAeq          = this.cosAeq;                   % line 143 cosAeq
            sinAeq          = this.sinAeq;                   % line 144 sinAeq
            tanHalf         = this.tanHalf;                  % line 145 tanHalf
            Feq             = this.Feq;                      % line 146 Feq
            Geq             = this.Geq;                      % line 147 Geq
            Heq             = this.Heq;                      % line 148 Heq
            Keq             = this.Keq;                      % line 149 Keq
            Leq             = this.Leq;                      % line 150 Leq
            CosL            = this.CosL;                     % line 151 CosL
            SinL            = this.SinL;                     % line 152 SinL
            alphaSq         = this.alphaSq;                  % line 153 alphaSq
            Seq             = this.Seq;                      % line 154
            Weq             = this.Weq;                      % line 155
            Req             = this.Req;                      % line 156
            RovS            = this.RovS;                     % line 157
            srtpinv         = this.srtpinv;                  % line 158
            HK              = this.HK;                       % line 159
            OnePalphaSqCosL = this.OnePalphaSqCosL;          % line 160
            OneMalphaSqCosL = this.OneMalphaSqCosL;          % line 161
            OnePalphaSqSinL = this.OnePalphaSqSinL;          % line 162
            OneMalphaSqSinL = this.OneMalphaSqSinL;          % line 163
            Xfac            = this.Xfac;                     % line 164
            Yfac            = this.Yfac;                     % line 165
            Zfac            = this.Zfac;                     % line 166
            VXfac           = this.VXfac;                    % line 167
            VYfac           = this.VYfac;                    % line 168
            VZfac           = this.VZfac;                    % line 169
            Xeq             = this.Xeq;                      % line 170
            Yeq             = this.Yeq;                      % line 171
            Zeq             = this.Zeq;                      % line 172
            VXeq            = this.VXeq;                     % line 173
            VYeq            = this.VYeq;                     % line 174
            VZeq            = this.VZeq;                     % line 175
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Semi-Major Axis:
            % % 	// from line 11: onePe = 1.0 + e;
            % q(11) = q(11) + q(2);
            % % 	// from line 12: oneMe = 1.0 - e;
            % q(12) = q(12) - q(2);
            % % 	// from line 13: fac   = onePe*oneMe;
            % q(13) = q(13) + q(11)*oneMe + q(12)*onePe;
            %  // from line 13: fac   = 1 - e^2;
            q(13) = q(13) - 2.0*e*q(2);
            % 	// from line 14: if(fac > 0) rootfac = sqrt(fac);
            %disp(' fordiff')
            %fac
            if (fac > 0.0)
                q(14) = q(14) + q(13)*0.5/rootfac;
            end
            % 	// from line 15: a = p/fac;
            % q(15) = q(15) + q(3)/fac - q(13)*a/fac;
            % 	// from line 15: p = a*fac;
            q(15) = q(15) + q(3)*fac + q(13)*a;
            %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %     % Period
            % 	// from line 16: meanMotion = sqrt(mu)*(a^(-1.5);
            %q(16) = q(16) - q(15)*1.5*meanMotion/a;
            q(16) = q(16) - q( 3)*1.5*meanMotion/a;
            %   // from line 17: cosI             = cos(this.Inclination);
            q(17) = q(17) - q(4)*sinI;
            %   // from line 18: sinI             = sin(this.Inclination);
            q(18) = q(18) + q(4)*cosI;
            %   // from line 19: cosom            = cos(this.omega);
            q(19) = q(19) - q(5)*sinom;
            %   // from line 20: sinom            = sin(this.omega);
            q(20) = q(20) + q(5)*cosom;
            %   // from line 21: cosO             = cos(this.Omega);
            q(21) = q(21) - q(6)*sinO;
            %   // from line 22: sinO             = sin(this.Omega);
            q(22) = q(22) + q(6)*cosO;
            %///////////////////////////////////////Block 07-27/////////////////////////////////////////
            %   // from line 23: Px   =  cosO*cosom - sinO*sinom*cosI;
            q(23) = q(23) + q(21)*cosom + q(19)*cosO;
            q(23) = q(23) - q(22)*sinom*cosI - q(20)*sinO*cosI - q(17)*sinO*sinom;
            %   // from line 24: Py   =  sinO*cosom + cosO*sinom*cosI;
            q(24) = q(24) + q(22)*cosom + q(19)*sinO;
            q(24) = q(24) + q(21)*sinom*cosI + q(20)*cosO*cosI + q(17)*cosO*sinom;
            %   // from line 25: Pz   =  sinom*sinI;
            q(25) = q(25) + q(20)*sinI + q(18)*sinom;
            %   // from line 26: Qx   = -cosO*sinom - sinO*cosom*cosI;
            q(26) = q(26) - q(21)*sinom - q(20)*cosO;
            q(26) = q(26) - q(22)*cosom*cosI - q(19)*sinO*cosI - q(17)*sinO*cosom;
            %   // from line 27: Qy   = -sinO*sinom + cosO*cosom*cosI;
            q(27) = q(27) - q(22)*sinom - q(20)*sinO;
            q(27) = q(27) + q(21)*cosom*cosI + q(19)*cosO*cosI + q(17)*cosO*cosom;
            %   // from line 28: Qz   =  cosom*sinI;
            q(28) = q(28) + q(19)*sinI + q(18)*cosom;
            %   // from line 29: Wx   =  sinO*sinI;
            q(29) = q(29) + q(22)*sinI + q(18)*sinO;
            %   // from line 30: Wy   = -cosO*sinI;
            q(30) = q(30) - q(21)*sinI - q(18)*cosO;
            %   // from line 31: Wz   =  cosI;
            %q(31) = q(31) + q(17);
            q(31) = q(31) - sinI*q( 4);
            %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % // from line 40: M = meanMotion*t + Mp;
            q(40) = q(40) + q( 1)*meanMotion + q(16)*t + q(7);
            % Solve Kepler's Equation for E:  M = E - e*sin(E);
            % // from line 41: E = E(e,M)
            q(41) = q(41) + q(2)*dE_de + q(40)*dE_dM;
            % from line 42: cosK  = cos(E);
            q(42) = q(42) - q(41)*sinK;
            % from line 43: sinK  = sin(E);
            q(43) = q(43) + q(41)*cosK;
            %///////////////////////////////////////Block 07-27/////////////////////////////////////////
            % from line 44: tanX  = cosK - e;
            q(44) = q(44) - q(2) + q(42);
            % // from line 45: tanY  = rootfac*sinK;
            q(45) = q(45) + q(14)*sinK + q(43)*rootfac;
            % // from line 46: nu     = atan2(tanY, tanX);
            q(46) = q(46) + (q(45)*tanX - q(44)*tanY)/(tanX*tanX + tanY*tanY);
            % // from line 47: coss  = cos(nu);
            q(47) = q(47) - q(46)*sins;
            % // from line 48: sins  = sin(nu);
            q(48) = q(48) + q(46)*coss;
            % // from line 49: rorbit   = a*(1.0 - e*cosK);
            %q(49) = q(49) + q(15)*(1.0 - e*cosK) - q(42)*a*e - q(2)*a*cosK;
            q(49) = q(49) + q( 3)*(1.0 - e*cosK) - q(42)*a*e - q(2)*a*cosK;
            %q(49) = q(49) + q( 3)*(1.0 - e*cosK) + q(41)*a*e*sinK - q(2)*a*cosK;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % // from line 50: rorbitx  = rorbit*(coss*Px + sins*Qx);
            %      q(50) = q(50)  + q(46)*rorbit*(-sins*Px + coss*Qx);% line 46 nu
            q(50) = q(50)  + q(49)*(coss*Px + sins*Qx);         % line 49 rorbit
            q(50) = q(50)  + q(48)*rorbit*Qx;                   % line 48 sins
            q(50) = q(50)  + q(47)*rorbit*Px;                   % line 47 coss
            q(50) = q(50)  + q(26)*rorbit*sins;                 % line 26 Qx
            q(50) = q(50)  + q(23)*rorbit*coss;                 % line 23 Px
            % // from line 51: rorbity  = rorbit*(coss*Py + sins*Qy);
            %      q(51) = q(51) + q(46)*rorbit*(-sins*Py + coss*Qy); % line 46 nu
            q(51) = q(51) + q(49)*(coss*Py + sins*Qy);          % line 49 rorbit
            q(51) = q(51) + q(48)*rorbit*Qy;                    % line 48 sins
            q(51) = q(51) + q(47)*rorbit*Py;                    % line 47 coss
            q(51) = q(51) + q(27)*rorbit*sins;                  % line 27 Qy
            q(51) = q(51) + q(24)*rorbit*coss;                  % line 24 Py
            % // from line 52: rorbitz  = rorbit*(coss*Pz + sins*Qz);
            %      q(52) = q(52) + q(46)*rorbit*(-sins*Pz + coss*Qz); % line 46 nu
            q(52) = q(52) + q(49)*(coss*Pz + sins*Qz);          % line 49 rorbit
            q(52) = q(52) + q(48)*rorbit*Qz;                    % line 48 sins
            q(52) = q(52) + q(47)*rorbit*Pz;                    % line 47 coss
            q(52) = q(52) + q(28)*rorbit*sins;                  % line 28 Qz
            q(52) = q(52) + q(25)*rorbit*coss;                  % line 25 Pz
            % from line 53: rtpinv   = sqrt(1.0/p);
            %q(53) = q(53) - q( 3)*rtpinv*(0.5/p);
            q(53) = q(53) - q(15)*rtpinv*(0.5/p);
            % // from line 54: vorbitx  = sqrt(1.0/p)*(-sins*Px + (e + coss)*Qx);
            %      q(54) = q(54) - q(46)*rtpinv*(coss*Px + sins*Qx);  % line 46 nu
            q(54) = q(54) - q(48)*rtpinv*Px;                    % line 48 sins
            q(54) = q(54) + q(47)*rtpinv*Qx;                    % line 47 coss
            q(54) = q(54) + q(26)*rtpinv*(e + coss);            % line 26 Qx
            %q(54) = q(54) - q( 3)*rtpinv*(0.5/p)*(-sins*Px + (e + coss)*Qx); % line 3 p
            q(54) = q(54) + q(53)*(-sins*Px + (e + coss)*Qx); % line 3 p
            q(54) = q(54) - q(23)*rtpinv*sins;                  % line 23 Px
            q(54) = q(54) + q( 2)*rtpinv*Qx;                    % line  2 e
            % // from line 55: vorbity  = sqrt(1.0/p)*(-sins*Py + (e + coss)*Qy);
            %      q(55) = q(55) - q(46)*rtpinv*(coss*Py + sins*Qy);  % line 46 nu
            q(55) = q(55) - q(48)*rtpinv*Py;                    % line 48 sins
            q(55) = q(55) + q(47)*rtpinv*Qy;                    % line 47 coss
            q(55) = q(55) + q(27)*rtpinv*(e + coss);            % line 27 Qy
            %q(55) = q(55) - q( 3)*rtpinv*(0.5/p)*(-sins*Py + (e + coss)*Qy); % line 3 p
            q(55) = q(55) + q(53)*(-sins*Py + (e + coss)*Qy); % line 3 p
            q(55) = q(55) - q(24)*rtpinv*sins;                  % line 24 Pq
            q(55) = q(55) + q( 2)*rtpinv*Qy;                    % line  2 e
            % // from line 56: vorbitz  = sqrt(1.0/p)*(-sins*Pz + (e + coss)*Qz);
            %      q(56) = q(56) - q(46)*rtpinv*(coss*Pz + sins*Qz);  % line 46 nu
            q(56) = q(56) - q(48)*rtpinv*Pz;                    % line 48 sins
            q(56) = q(56) + q(47)*rtpinv*Qz;                    % line 47 coss
            q(56) = q(56) + q(28)*rtpinv*(e + coss);            % line 28 Qz
            %q(56) = q(56) - q( 3)*rtpinv*(0.5/p)*(-sins*Pz + (e + coss)*Qz); % line 3 p
            q(56) = q(56) + q(53)*(-sins*Pz + (e + coss)*Qz);   % line 3 p
            q(56) = q(56) - q(25)*rtpinv*sins;                  % line 25 Pz
            q(56) = q(56) + q( 2)*rtpinv*Qz;                    % line  2 e
            if LineCheck > 60
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line 61: X    = rorbitx = Kvec(1);                     % line  61
                q(61)  = q(61) + q(50);
                % from line 62: Y    = rorbity = Kvec(2);                     % line  62
                q(62)  = q(62) + q(51);
                % from line 63: Z    = rorbitz = Kvec(3);                     % line  63
                q(63)  = q(63) + q(52);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line 65: R1   = rorbit;                       % line  64
                %q(65) = q(65) + q(49);
                % from line 64: R2   = R1^2;                         % line  65
                %q(64) = q(64) + 2.0*q(65)*R1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line 64: R2   = X^2 + Y^2 + Z^2;                         % line  64
                q(64) = q(64) + 2.0*(q(61)*X + q(62)*Y + q(63)*Z);
                % from line 65: R1   = sqrt(R2);                                % line  65
                q(65) = q(65) + q(64)*0.5/R1;
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line 66: R3   = R1*R2;
                q(66) = q(66) + q(65)*R2 + q(64)*R1;
                % from line 67: R4   = R1*R3;                                   % line  67
                q(67) = q(67) + q(65)*R3 + q(66)*R1;
                % from line 68: R5   = R1*R4;                                   % line  68
                q(68) = q(68) + q(65)*R4 + q(67)*R1;
                % from line 69: R6   = R1*R5;                                   % line  69
                q(69) = q(69) + q(65)*R5 + q(68)*R1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line 70: sd1  = Z/R1;
                q(70) = q(70) + q(63)/R1 - q(65)*sd1/R1;
                % from line 70: sd1  = sin(I)*sin(omega + nu);
                %q(70)  = q(70) + q(18)*sinomnu + sinI*cosomnu*(q(5) + q(46));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line 71: sd2  = sd1^2;
                q(71) = q(71) + q(70)*2.0*sd1;
                % from line 72: sd3  = sd1*sd2;
                q(72) = q(72) + q(70)*sd2 + q(71)*sd1;
                % from line 73: sd4  = sd1*sd3;
                q(73) = q(73) + q(70)*sd3 + q(72)*sd1;
                % from line 74: sd5  = sd1*sd4;
                q(74) = q(74) + q(70)*sd4 + q(73)*sd1;
                % from line 75: sd6  = sd1*sd5;
                q(75) = q(75) + q(70)*sd5 + q(74)*sd1;
                %
                % from line 76: F2   =  1.0      -  3.0*sd2;
                q(76) = q(76) - q(71)*3.;
                % from line 77: F3   =  3.0*sd1  -  5.0*sd3;
                q(77) = q(77) + q(70)*3.0 - q(72)*5.0;
                % from line 78: F4   =  3.0      -  30.*sd2  + 35*sd4;
                q(78) = q(78) - q(71)*30.0 + q(73)*35.0;
                % from line 79: F5   =  15.*sd1  -  70.*sd3  + 63*sd5;
                q(79) = q(79) + q(70)*15. - q(72)*70.  + q(74)*63.;
                % from line 80: F6   =  5.       - 105.*sd2  + 315*sd4 - 231*sd6;
                q(80) = q(80) - q(71)*105. + q(73)*315. - q(75)*231.;
                %
                % from line 81: V1    = 1/R1; Kepler Potential                 % line 21
                q(81) = q(81) - q(65)*V1/R1;
                % from line 82: V2      = C2*F2/R2;
                q(82) = q(82) + q(76)*C2/R2 - q(64)*V2/R2;
                % from line 83: V3      = C3*F3/R3;
                q(83) = q(83) + q(77)*C3/R3 - q(66)*V3/R3;
                % from line 84: V4      = C4*F4/R4;
                q(84) = q(84) + q(78)*C4/R4 - q(67)*V4/R4;
                % from line 85: V5      = C5*F5/R5;
                q(85) = q(85) + q(79)*C5/R5 - q(68)*V5/R5;
                % from line 86: V6      = C6*F6/R6;
                q(86) = q(86) + q(80)*C6/R6 - q(69)*V6/R6;
                % from line 87: V       = mu*V1;
                q(87) = q(87) + q(81)*mu;
                % Up to J6
                % from line 88: D       = V2 + V3 + V4 + V5 + V6;
                q(88) = q(88) + q(82) + q(83) + q(84) + q(85) + q(86);
                % from line 89: R       = mu*D;
                %q(89) = q(89) + q(88)*mu;
                % % from line 90: Phi     = V + R;
                %q(90) = q(90) + q(87) + q(89);
                % from line 90: Phi     = V*(CK + D);
                q(90) = q(90) + q(87)*(CK + D) + q(88)*V;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if LineCheck > 90
                    % from line 91: fac1 = meanMotion*a           % n*a               % line  91
                    q(91) = q(91) + q( 3)*meanMotion + q(16)*a;
                    % from line 92: fac2 = fac1*a        % n*a^2                         % line  92
                    q(92) = q(92) + q( 3)*fac1       + q(91)*a;
                    % from line  93: fac3 = fac2*e        % n*a^2*e                       % line  93
                    q(93) = q(93) + q( 2)*fac2       + q(92)*e;
                    % from line  94: fac4 = fac2*rootfac  % n*a^2*sqrt(1-e^2)             % line  94
                    q(94) = q(94) + q(14)*fac2       + q(92)*rootfac;
                    % from line  95: fac5 = fac4*sinI     % n*a^2*sqrt(1-e^2)*sinI        % line  95
                    q(95) = q(95) + q( 4)*fac4*cosI  + q(94)*sinI;
                    % from line  96: Mwe  = rootfac/fac3  % sqrt(1-e^2)/(n*a^2*e)         % line  96
                    q(96) = q(96) + q(14)/fac3       - q(93)*rootfac/fac3^2;
                    % from line  97: MeM  = fac/fac3      % (1-e^2)/(n*a^2*e)             % line  97
                    q(97) = q(97) + q(13)/fac3       - q(93)*fac/fac3^2;
                    %from line   98: MaM  = 2.0/fac1      % 2/(n*a)                       % line  98
                    q(98) = q(98) - 2.0*q(91)/fac1^2;
                    % from line  99: MWi  = 1/fac5        % 1/(n*a^2*sqrt(1-e^2)*sinI)    % line  99
                    q(99) = q(99) - q(95)/fac5^2;
                    % from line 100: Miw  = cosI*MWi      % cosI/(n*a^2*sqrt(1-e^2)*sinI) % line 100
                    q(100)= q(100)- q( 4)*sinI*MWi + q(99)*cosI;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if LineCheck > 100
                        % Forwards Differentiation:
                        % Prediction:
                        %los       = Rextrap - state_vector_sat(1:3);                  %lines 101-103
                        % from line 101: lx = los(1);
                        q(101) = q(101) + q(50);
                        % from line 102: ly = los(2);                                    %line 102
                        q(102) = q(102) + q(51);
                        % from line 103: lz = los(3);                                    %line 103
                        q(103) = q(103) + q(52);
                        %range   = norm(los)
                        %theta   = acos(los(3)/range)
                        %phi     = atan2(los(2),los(1))
                        % from line 104: lperpSq  = lx^2 + ly^2;                         % line 104
                        q(104) = q(104) + 2.0*(lx*q(101) + ly*q(102));
                        % from line 105: rangeSq  = lperpSq + lz^2;                      % line 105
                        q(105) = q(105) + q(104) + 2.0*lz*q(103);
                        % from line 106: range    = sqrt(rangeSq);                       % line 106
                        q(106) = q(106) + 0.5*q(105)/range;
                        % from line 107: zLos     = lz/range;                            % line 107
                        q(107) = q(107) + (q(103) - q(106)*zLos )/range;
                        % from line 108: thetaLos = acos(zLos);                          % line 108
                        q(108) = q(108) - q(107)/sqrt(1.0 - zLos^2);
                        % from line 109: phiLos   = atan2(ly,lx);                        % line 109
                        q(109) = q(109) + (q(102)*lx - q(101)*ly)/lperpSq;
                        %vlos     = Vextrap - state_vector_sat(4:6);     % lines 110-112
                        % from line 110: vx       = vlos(1);                             % line 110
                        q(110) = q(110) + q(54);
                        % from line 111: vy       = vlos(2);                             % line 111
                        q(111) = q(111) + q(55);
                        % from line 112: vz       = vlos(3);                             % line 112
                        q(112) = q(112) + q(56);
                        % from line 113: ldotv    = lx*vx + ly*vy + lz*vz;               % line 113
                        q(113) = q(113) + q(101)*vx + q(102)*vy + q(103)*vz + q(110)*lx + q(111)*ly + q(112)*lz;
                        % from line 114: ldot     = ldotv/range;                         % line 114
                        q(114) = q(114) + (q(113)- q(106)*ldot)/range;
                        % from line 115: lperp    = sqrt(lperpSq);                       % line 115
                        q(115) = q(115) + 0.5*q(104)/lperp;
                        % from line 116: thdotfac1 = range*vz - lz*ldot;                 % line 116
                        q(116) = q(116) - q(103)*ldot + q(106)*vz + q(112)*range - q(114)*lz;
                        % from line: thdotfac2 = thdotfac1/range;                    % line 117
                        q(117) = q(117) + (q(116) - q(106)*thdotfac2)/range;
                        % from line 118: thdot     = -thdotfac2/lperp;                   % line 118
                        q(118) = q(118) - (q(115)*thdot + q(117))/lperp;
                        % from line 119: phiDotNum = lx*vy - ly*vx;                      % line 119
                        q(119) = q(119) + q(101)*vy + q(111)*lx - q(102)*vx - q(110)*ly;
                        % from line 120: phidot    = phiDotNum/lperpSq;                  % line 120
                        q(120) = q(120) + (q(119) - q(104)*phidot)/lperpSq;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % from line 121: cosThetaLos = cos(thetaLos);
                        q(121) = q(121) - q(108)*sinThetaLos;
                        % from line 122: sinThetaLos = sin(thetaLos);
                        q(122) = q(122) + q(108)*cosThetaLos;
                        % from line 123: phiDiff     = phiLos - phiM;
                        q(123) = q(123) + q(109);
                        % from line 124: cosPhiDiff  = cos(phiDiff);
                        q(124) = q(124) - q(123)*sinPhiDiff;
                        % from line 125: sinPhiDiff  = sin(phiDiff);
                        q(125) = q(125) + q(123)*cosPhiDiff;
                        % from line 126: delUV       = cosThetaLos*sinThetaM*cosPhiDiff;
                        q(126) = q(126) + sinThetaM*(q(121)*cosPhiDiff + q(124)*cosThetaLos);
                        % from line 127: Uview       = sinThetaLos*cosThetaM - delUV;
                        q(127) = q(127) + q(122)*cosThetaM - q(126);
                        % from line 128: Vview       = sinThetaM*sinPhiDiff;
                        q(128) = q(128) + sinThetaM*q(125);
                        % from line 129: xLos     = lx/range;
                        q(129) = q(129) + (q(101) - q(106)*xLos )/range;
                        % from line 130: yLos     = ly/range;
                        q(130) = q(130) + (q(102) - q(106)*yLos )/range;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %from line 131: Collinearity = 1-xLos*xLosM-yLos*yLosM-zLos*zLosM;
                        q(131) = q(131) - q(129)*xLosM - q(130)*yLosM - q(107)*zLosM ;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % from line 133: ChiSqLos = losDiff'*InvCovM*losDiff;
                        q(133)  = q(133) + q(129)*(xUnit'*InvCovM*losDiff) ...
                            + q(130)*(yUnit'*InvCovM*losDiff) ...
                            + q(107)*(zUnit'*InvCovM*losDiff);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if LineCheck > 133
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % from line 142: Aeq      = omega + Omega;       % line 142
                            q(142) = q(142) + q( 5) + q( 6);
                            % from line 143: cosAeq   = cos(Aeq);              % line 143
                            q(143) = q(143) - q(142)*sinAeq;
                            % from line 84: sinAeq   = sin(Aeq);              % line 144
                            q(144) = q(144) + q(142)*cosAeq;
                            % from line 145: tanHalf  = tan(Inclination/2.0);  % line 145
                            q(145) = q(145) + 0.5*q( 4)*(1.0 + tanHalf^2);
                            % from line 146: Feq      = e*cosAeq;              % line 146
                            q(146) = q(146) + q( 2)*cosAeq + q(143)*e;
                            % from line 147: Geq      = e*sinAeq;              % line 147
                            q(147) = q(147) + q( 2)*sinAeq + q(144)*e;
                            % from line 148: Heq      = tanHalf*cos(Omega);    % line 148
                            q(148) = q(148) + q(145)*cosO - q( 6)*tanHalf*sinO;
                            % from line 149: Keq      = tanHalf*sin(Omega);    % line 149
                            q(149) = q(149) + q(145)*sinO + q( 6)*tanHalf*cosO;
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % from line 140: Leq      = Aeq + nu;         % line 150
                            q(150) = q(150) + q(46) + q(142);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % from line 151: CosL     = cos(Leq);          % line 151
                            q(151) = q(151) - q(150)*SinL;
                            % from line 152: SinL     = sin(Leq);          % line 152
                            q(152) = q(152) + q(150)*CosL;
                            % from line 153: alphaSq = Heq^2 - Keq^2;         % line 153
                            q(153) = q(153) + 2.0*(q(148)*Heq - q(149)*Keq);
                            % from line 154: Seq     = 1.0 + Heq^2 + Keq^2;   % line 154
                            q(154) = q(154) + 2.0*(q(148)*Heq + q(149)*Keq);
                            % from line 155: Weq = 1.0 + Feq*CosL + Geq*SinL; % line 155
                            q(155) = q(155) + q(146)*CosL + q(147)*SinL + q(151)*Feq + q(152)*Geq;
                            % from line 156: Req     = p/Weq;              % line 156
                            q(156) = q(156) + q(15)/Weq - q(155)*Req/Weq;
                            % from line 157: RovS    = Req/Seq;            % line 157
                            q(157) = q(157) + q(156)/Seq - q(154)*RovS/Seq;
                            % from line 158: srtpinv = rtpinv/Seq;         % line 158
                            q(158) = q(158) + q(53)/Seq - q(154)*srtpinv/Seq;
                            % from line 159: HK  = Heq*Keq;            % line 159
                            q(159) = q(159) + q(148)*Keq + q(149)*Heq;
                            % from line 160: OnePalphaSqCosL = (1.+alphaSq)*CosL;% line 160
                            q(160) = q(160) + q(151)*(1.0 + alphaSq) + q(153)*CosL;
                            % from line 161: OneMalphaSqCosL = (1.-alphaSq)*CosL;% line 161
                            q(161) = q(161) + q(151)*(1.0 - alphaSq) - q(153)*CosL;
                            % from line 162: OnePalphaSqSinL = (1.+alphaSq)*SinL;% line 162
                            q(162) = q(162) + q(152)*(1.0 + alphaSq) + q(153)*SinL;
                            % from line 163: OneMalphaSqSinL = (1.-alphaSq)*SinL;% line 163
                            q(163) = q(163) + q(152)*(1.0 - alphaSq) - q(153)*SinL;
                            % from line 164: Xfac    = OnePalphaSqCosL + 2.0*HK*SinL; % line 164
                            q(164) = q(164) + 2.0*q(152)*HK + 2.0*q(159)*SinL + q(160);
                            % from line 165: Yfac    = OneMalphaSqSinL + 2.0*HK*CosL; % line 165
                            q(165) = q(165) + 2.0*q(151)*HK + 2.0*q(159)*CosL + q(163);
                            
                            % from line 166: Zfac    = Heq*SinL - Keq*CosL;          % line 166
                            q(166) = q(166) + q(148)*SinL - q(149)*CosL  - q(151)*Keq  + q(152)*Heq;
                            
                            % from line 167: VXfac =  OnePalphaSqSinL - 2.0*HK*(Feq+CosL) + Geq*(1+alphaSq);% line 167
                            q(167) = q(167) -2.0*HK*(q(146)+q(151)) +q(147)*(1.0+alphaSq) +q(153)*Geq -2.0*q(159)*(Feq+CosL) +q(162);
                            
                            % from line 168: VYfac = -OneMalphaSqCosL + 2.0*HK*(Geq+SinL) + Feq*(alphaSq-1);% line 168
                            q(168) = q(168) + q(146)*(alphaSq-1.0) +2.0*HK*(q(147)+q(152))+q(153)*Feq +2.0*q(159)*(Geq+SinL) -q(161);
                            
                            % from line 169: VZfac   =  Heq*(Feq + CosL) + Keq*(Geq + SinL) % line 169
                            q(169) = q(169) +q(146)*Heq +q(147)*Keq +q(148)*(Feq+CosL) +q(149)*(Geq+SinL) +q(151)*Heq +q(152)*Keq;
                            
                            % from line 170: Xeq     =     RovS*Xfac;     % line 170
                            q(170) = q(170) + q(157)*Xfac + q(164)*RovS;
                            
                            % from line 171: Yeq     =     RovS*Yfac;     % line 171
                            q(171) = q(171) + q(157)*Yfac + q(165)*RovS;
                            
                            %from line 172: Zeq      = 2.0*RovS*Zfac;     % line 172
                            q(172) = q(172) + 2.0*(q(157)*Zfac + q(166)*RovS);
                            
                            % from line 173: VXeq    =    -srtpinv*VXfac; % line 173
                            q(173) = q(173) - q(158)*VXfac - q(167)*srtpinv;
                            
                            % from line 174: VYeq    =    -srtpinv*VYfac; % line 174
                            q(174) = q(174) - q(158)*VYfac - q(168)*srtpinv;
                            
                            % from line 175: VZeq    = 2.0*srtpinv*VZfac; % line 175
                            q(175) = q(175) + 2.0*q(158)*VZfac + 2.0*q(169)*srtpinv;
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                    end
                end
            end
            Q = q;
        end
        
        %function S = secdiff(this, f, q, s, varargin)
        function S = secdiff(this, f, q, s, LineCheck)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % DU    = 1.0;
            % VU    = 1.0;
            % TU    = 1.0;
            mu    = this.mu;
            % if nargin >= 5
            %     TU = this.TU;
            %     DU = this.DU;
            %     VU = this.VU;
            %     AU = this.AU;
            % end
            %disp('Hello World');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            t       = this.time;                                 % line   1
            e       = this.e;                                    % line   2
            a       = this.a;                                    % line   3
            %p       = this.p;                                    % line   3
            I       = this.Inclination;                          % line   4
            omega   = this.omega;                                % line   5
            Omega   = this.Omega;                                % line   6
            Mp      = this.Mp;                                   % line   7
            
            onePe   = this.onePe;                                % line  11
            oneMe   = this.oneMe;                                % line  12
            fac     = this.fac;                                  % line  13
            rootfac = this.rootfac;                              % line  14
            %a       = this.a;                                   % line  15
            p       = this.p;                                    % line  15
            meanMotion = this.meanMotion;                        % line  16
            
            cosI    = this.cosI;                                 % line  17
            sinI    = this.sinI;                                 % line  18
            cosom   = this.cosom;                                % line  19
            sinom   = this.sinom;                                % line  20
            cosO    = this.cosO;                                 % line  21
            sinO    = this.sinO;                                 % line  22
            
            Px      = this.Px;                                   % line  23
            Py      = this.Py;                                   % line  24
            Pz      = this.Pz;                                   % line  25
            Qx      = this.Qx;                                   % line  26
            Qy      = this.Qy;                                   % line  27
            Qz      = this.Qz;                                   % line  28
            
            M       = this.M;                                    % line  40
            E       = this.E;                                    % line  41
            cosK    = this.cosK;                                 % line  42
            sinK    = this.sinK;                                 % line  43
            tanX    = this.tanX;                                 % line  44
            tanY    = this.tanY;                                 % line  45
            nu      = this.nu;                                   % line  46
            coss    = this.coss;                                 % line  47
            sins    = this.sins;                                 % line  48
            rorbit  = this.rorbit;                               % line  49
            rorbitx = this.rorbitx;                              % line  50
            rorbity = this.rorbity;                              % line  51
            rorbitz = this.rorbitz;                              % line  52
            rtpinv  = this.rtpinv;                               % line  53
            vorbitx = this.vorbitx;                              % line  54
            vorbity = this.vorbity;                              % line  55
            vorbitz = this.vorbitz;                              % line  56
            
            dE_dM   = this.dE_dM;
            dE_de   = this.dE_de;
            d2E_dMdM = this.d2E_dMdM;
            d2E_dMde = this.d2E_dMde;
            d2E_dedM = this.d2E_dedM;
            d2E_dede = this.d2E_dede;
            rmag     = this.rmag;
            denom    = rmag*rmag*rmag;
            edenom   = e*e*e;
            cesqr    = (coss*e+1.0)*(coss*e+1.0);
            cecube   = cesqr*(coss*e+1.0);
            psq      = p*p;
            pcubed   = p*p*p;
            pdenom   = sqrt(pcubed);
            cosDenom = 0.0;
            
            CK = this.CK;
            C2 = this.C2;
            C3 = this.C3;
            C4 = this.C4;
            C5 = this.C5;
            C6 = this.C6;
            
            % Escobal Gravity Terms
            X        = this.X;
            Y        = this.Y;
            Z        = this.Z;
            R2       = this.R2;
            R1       = this.R1;
            R3       = this.R3;
            R4       = this.R4;
            R5       = this.R5;
            R6       = this.R6;
            %sinomnu  = this.sinomnu;
            %cosomnu  = this.cosomnu;
            %sd1chk   = this.sd1chk;
            %cd1chk   = this.cd1chk;
            sd1      = this.sd1;
            sd2      = this.sd2;
            sd3      = this.sd3;
            sd4      = this.sd4;
            sd5      = this.sd5;
            sd6      = this.sd6;
            F2       = this.F2;
            F3       = this.F3;
            F4       = this.F4;
            F5       = this.F5;
            F6       = this.F6;
            V1       = this.V1;
            V2       = this.V2;
            V3       = this.V3;
            V4       = this.V4;
            V5       = this.V5;
            V6       = this.V6;
            V        = this.V;
            D        = this.D;
            Phi      = this.Phi;
            
            fac1      = this.fac1;
            fac2      = this.fac2;
            fac3      = this.fac3;
            fac4      = this.fac4;
            fac5      = this.fac5;
            Mwe       = this.Mwe;
            MeM       = this.MeM;
            MaM       = this.MaM;
            MWi       = this.MWi;
            Miw       = this.Miw;
            %this.MLagrange = MLagrange;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            los       = this.los;   % Rextrap - Rsensor
            lx        = this.lx;                                 % line 101
            ly        = this.ly;                                 % line 102
            lz        = this.lz;                                 % line 103
            lperpSq   = this.lperpSq;                            % line 104
            rangeSq   = this.rangeSq;                            % line 105
            range     = this.range;                              % line 106
            zLos      = this.zLos;                               % line 107
            thetaLos  = this.thetaLos;                           % line 108
            phiLos    = this.phiLos;                             % line 109
            vlos      = this.vlos;                          % lines 110-112
            vx        = this.vx;                                 % line 110
            vy        = this.vy;                                 % line 111
            vz        = this.vz;                                 % line 112
            ldotv     = this.ldotv;                              % line 113
            ldot      = this.ldot;                               % line 114
            lperp     = this.lperp;                              % line 115
            thdotfac1 = this.thdotfac1;                          % line 116
            thdotfac2 = this.thdotfac2;                          % line 117
            thdot     = this.thdot;                              % line 118
            phiDotNum = this.phiDotNum;                          % line 119
            phidot    = this.phidot;                             % line 120
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            thetaM      = this.thetaM;
            cosThetaM   = this.cosThetaM;
            sinThetaM   = this.sinThetaM;
            phiM        = this.phiM;
            cosPhiM     = this.cosPhiM;
            sinPhiM     = this.sinPhiM;
            cosThetaLos = this.cosThetaLos;                      % line 121
            sinThetaLos = this.sinThetaLos;                      % line 122
            phiDiff     = this.phiDiff;                          % line 123
            cosPhiDiff  = this.cosPhiDiff;                       % line 124
            sinPhiDiff  = this.sinPhiDiff;                       % line 125
            delUV       = this.delUV;                            % line 126
            Uview       = this.Uview;                            % line 127
            Vview       = this.Vview;                            % line 128
            xLos        = this.xLos;                             % line 129
            yLos        = this.yLos;                             % line 130
            xLosM       = this.xLosM;
            yLosM       = this.yLosM;
            zLosM       = this.zLosM;
            Collinearity  = this.Collinearity;                   % line 131
            Mdelta      = this.Mdelta;
            InvCovM     = this.InvCovM;
            losDiff     = this.losDiff;                          % line 132
            xUnit       = this.xUnit;
            yUnit       = this.yUnit;
            zUnit       = this.zUnit;
            ChiSqLos    = this.ChiSqLos;                         % line 133
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Aeq             = this.Aeq;                      % line 142 Aeq
            cosAeq          = this.cosAeq;                   % line 143 cosAeq
            sinAeq          = this.sinAeq;                   % line 144 sinAeq
            tanHalf         = this.tanHalf;                  % line 145 tanHalf
            Feq             = this.Feq;                      % line 146 Feq
            Geq             = this.Geq;                      % line 147 Geq
            Heq             = this.Heq;                      % line 148 Heq
            Keq             = this.Keq;                      % line 149 Keq
            Leq             = this.Leq;                      % line 150 Leq
            CosL            = this.CosL;                     % line 151 CosL
            SinL            = this.SinL;                     % line 152 SinL
            alphaSq         = this.alphaSq;                  % line 153 alphaSq
            Seq             = this.Seq;                      % line 154
            Weq             = this.Weq;                      % line 155
            Req             = this.Req;                      % line 156
            RovS            = this.RovS;                     % line 157
            srtpinv         = this.srtpinv;                  % line 158
            HK              = this.HK;                       % line 159
            OnePalphaSqCosL = this.OnePalphaSqCosL;          % line 160
            OneMalphaSqCosL = this.OneMalphaSqCosL;          % line 161
            OnePalphaSqSinL = this.OnePalphaSqSinL;          % line 162
            OneMalphaSqSinL = this.OneMalphaSqSinL;          % line 163
            Xfac            = this.Xfac;                     % line 164
            Yfac            = this.Yfac;                     % line 165
            Zfac            = this.Zfac;                     % line 166
            VXfac           = this.VXfac;                    % line 167
            VYfac           = this.VYfac;                    % line 168
            VZfac           = this.VZfac;                    % line 169
            Xeq             = this.Xeq;                      % line 170
            Yeq             = this.Yeq;                      % line 171
            Zeq             = this.Zeq;                      % line 172
            VXeq            = this.VXeq;                     % line 173
            VYeq            = this.VYeq;                     % line 174
            VZeq            = this.VZeq;                     % line 175
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Semi-Major Axis:
            % % 	// from line 11: onePe = 1.0 + e;
            % % e: J=2  KJ = 1 All 2nd derivatives vanish
            % % 	// from line 12: oneMe = 1.0 - e;
            % % e: J=2 KJ = 1 All 2nd derivatives vanish
            % % 	// from line 13: fac   = onePe*oneMe;
            % % onePe: J = 11 KJ = oneMe
            % s(11) = s(11) + f(13)*q(12);
            % % oneMe: J = 12 KJ = onePe
            % s(12) = s(12) + f(13)*q(11);
            % // from line 13: fac   = 1 - e^2;
            % J = 2: KJ = -2.*e
            s(2) = s(2) - f(13)*q(2)*2.0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 	// from line 14: if(fac > 0) rootfac = sqrt(fac);
            % K = 14
            if fac>0
                % fac: J=13  KJ = 0.5/rootfac = 0.5*fac^(1/2), I = 13
                s(13) = s(13) - f(14)*q(13)*0.25/(rootfac^3);
                % try I = 14
                %s(14) = s(14) - f(14)*q(14)*0.5/(rootfac^2);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 	// from line 15: a = p/fac;
            % fac: J = 13 KJ = -p*fac^(-2)
            % fac: I = 13 & p I = 3
            %s(13) = s(13) + f(15)*(q(13)*2.0*a/fac^2 - q(3)/fac^2);
            % p: J = 3 KJ = fac(-1);
            %s(3) = s(3) - f(15)*q(13)/fac^2;
            % 	// from line 15: p = a*fac;
            % fac: J = 13 KJ = a I = 3
            s(13) = s(13) + f(15)*q( 3);
            % a: J = 3 KJ = fac, I =13;
            s( 3) = s( 3) + f(15)*q(13);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 	// from line 16: meanMotion = sqrt(mu)*a^(-1.5);
            % a: J = 15 KJ = -1.5*a(-2.5)
            % a I = 15
            %s(15) = s(15) + f(16)*q(15)*3.75*meanMotion/(a^2);
            % a: J = 3 KJ = -1.5*a(-2.5)
            % a I = 3
            s( 3) = s( 3) + f(16)*q( 3)*3.75*meanMotion/(a^2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   // from line 17: cosI             = cos(this.Inclination);
            % I: J = 4 KJ = -sinI;  (i = 18)
            %s(4) = s(4) - f(17)*q(18);
            % Alternate method
            s(4)  = s(4) - f(17)*q(4)*cosI;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   // from line 18: sinI             = sin(this.Inclination);
            % I  J = 4: KJ = cosI (i = 17)
            %s(4) = s(4) + f(18)*q(17);
            % Alternate method
            s(4)  = s(4) - f(18)*q(4)*sinI;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   // from line 19: cosom            = cos(this.omega);
            % omega J = 5: KJ = -sinom (i = 20)
            %s(5) = s(5) - f(19)*q(20);
            s(5) = s(5) - f(19)*q(5)*cosom;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   // from line 20: sinom            = sin(this.omega);
            % omega J = 5: KJ = cosom (i = 19)
            %s(5) = s(5) + f(20)*q(19);
            s(5) = s(5) - f(20)*q(5)*sinom;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   // from line 21: cosO             = cos(this.Omega);
            % Omega J = 6: KJ = -sinO (i = 22)
            %s(6) = s(6) - f(21)*q(22);
            s(6) = s(6) - f(21)*q(6)*cosO;
            %   // from line 22: sinO             = sin(this.Omega);
            % Omega J = 6: KJ = cosO (i = 21)
            %s(6) = s(6) + f(22)*q(21);
            s(6) = s(6) - f(22)*q(6)*sinO;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   // from line K = 23: Px   =  cosO*cosom - sinO*sinom*cosI;
            % K = 23 dependencies sinO 22, cosO 21, sinom 20, cosom 19, cosI 17
            % sinO J = 22: KJ = -sinom*cosI (i = 20, 17)
            s(22) = s(22) - f(23)*(q(20)*cosI + q(17)*sinom);
            % cosO J = 21: KJ = cosom (i = 19)
            s(21) = s(21) + f(23)*q(19);
            % sinom J = 20: KJ = -sinO*cosI (i = 22, 17)
            s(20) = s(20) - f(23)*(q(22)*cosI + q(17)*sinO);
            % cosom J = 19: KJ = cosO (i=21)
            s(19) = s(19) + f(23)*q(21);
            % cosI J = 17: KJ = - sinO*sinom (i=22, 20)
            s(17) = s(17) - f(23)*(q(22)*sinom + q(20)*sinO);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   // from line K = 24: Py   =  sinO*cosom + cosO*sinom*cosI;
            % K = 24 dependencies sinO 22, cosO 21, sinom 20, cosom 19, cosI 17
            % sinO J = 22: KJ = cosom (i = 19)
            s(22) = s(22) + f(24)*q(19);
            % cosO J = 21: KJ = sinom*cosI (i = 20, 17)
            s(21) = s(21) + f(24)*(q(20)*cosI + q(17)*sinom);
            % sinom J = 20: KJ = cosO*cosI (i = 21, 17)
            s(20) = s(20) + f(24)*(q(21)*cosI + q(17)*cosO);
            % cosom J = 19: KJ = sinO (i=22)
            s(19) = s(19) + f(24)*q(22);
            % cosI J = 17: KJ = cosI*sinom (i=21, 20)
            s(17) = s(17) + f(24)*(q(21)*sinom + q(20)*cosO);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   // from line K = 25: Pz   =  sinom*sinI;
            % K = 25 dependencies sinom 20, sinI 18
            % sinom J = 20: KJ = sinI (i=18)
            s(20) = s(20) + f(25)*q(18);
            % sinI J = 18: KJ = sinom (i=20)
            s(18) = s(18) + f(25)*q(20);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   // from line 26: Qx   = -cosO*sinom - sinO*cosom*cosI;
            % K = 26 dependencies (sinO 22, cosO 21,  sinom 20, cosom 19, cosI 17
            % sinO J = 22: KJ = -cosom*cosI (i = 19, 17)
            s(22) = s(22) - f(26)*(q(19)*cosI + q(17)*cosom);
            % cosO J = 21: KJ = -sinom (i=20)
            s(21) = s(21) - f(26)*q(20);
            % sinom J = 20: KJ = -cosO (i=21)
            s(20) = s(20) - f(26)*q(21);
            % cosom J = 19: KJ = -sinO*cosI (i=22, 17)
            s(19) = s(19) - f(26)*(q(22)*cosI + q(17)*sinO);
            % cosI J = 17: KJ = -sinO*cosom (i=22, 19)
            s(17) = s(17) - f(26)*(q(22)*cosom + q(19)*sinO);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   // from line 27: Qy   = -sinO*sinom + cosO*cosom*cosI;
            % K = 27 dependencies (sinO 22, cosO 21, sinom 20, cosom 19, cosI 17)
            % sinO J = 22: KJ = -sinom
            s(22) = s(22) - f(27)*q(20);
            % cosO J = 21: KJ = cosom*cosI (i=19, 17)
            s(21) = s(21) + f(27)*(q(19)*cosI + q(17)*cosom);
            % sinom J = 20: KJ = -sinO (i=22);
            s(20) = s(20) - f(27)*q(22);
            % cosom J = 19: KJ = cosO*cosI (i=21, 17)
            s(19) = s(19) + f(27)*(q(21)*cosI + q(17)*cosO);
            % cosI J = 17: KJ = cosO*cosom (i=21, 19)
            s(17) = s(17) + f(27)*(q(21)*cosom + q(19)*cosO);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   // from line 28: Qz   =  cosom*sinI;
            % K = 28 dependencies (cosom 19, sinI 18)
            % cosom J = 19: KJ = sinI (i=18)
            s(19) = s(19) + f(28)*q(18);
            % sinI J = 18: KJ = cosom (i=19)
            s(18) = s(18) + f(28)*q(19);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   // from line 29: Wx   =  sinO*sinI;
            % K = 29 dependencies (sinO 22, sinI 18)
            % sinO J = 22: KJ = sinI (i=18)
            s(22) = s(22) + f(29)*q(18);
            % sinI J = 18: KJ = sinO (i=22)
            s(18) = s(18) + f(29)*q(22);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   // from line 30: Wy   = -cosO*sinI;
            % K = 30 dependencies (cosI 21, sinI 18)
            % cosI J = 21: KJ = -sinI (i=18)
            s(21) = s(21) - f(30)*q(18);
            % sinI J = 18: KJ = -cosI (i=21)
            s(18) = s(18) - f(30)*q(21);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   // from line 31: Wz   =  cosI;
            % K = 31 dependencies cosI 17
            % J = 17 KJ = -sinI (i=18)
            %s(17) = s(17) - f(31)*q(18);
            s(4) = s(4) - f(31)*q(4)*cosI;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % // from line 40: M = meanMotion*t + Mp;
            % K = 40 dependencies (meanMotion 16, Mp 7, time 1)
            % meanMotion J = 16: KJ = t
            s(16) = s(16) + f(40)*q(1);
            % Mp J = 7: KJ = 1 (2nd derivative vanishes)!
            % t J = 1: KJ = meanMotion
            s(1) = s(1) + f(40)*q(16);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Solve Kepler's Equation for E:  M = E - e*sin(E);
            % from line 41: E = E(e,M)
            % K = 41 dependencies (M 40, e 2)
            % M J = 40 KJ = dE_dM (i = 40, 2)
            s(40) = s(40) + f(41)*(q(2)*d2E_dMde + q(40)*d2E_dMdM);
            % e J = 2 KJ = dE_de (i = 40,2)
            s( 2) = s( 2) + f(41)*(q(2)*d2E_dede + q(40)*d2E_dedM);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % from line 42: cosK  = cos(E);
            % line K = 42 dependencies (E 41)
            % // E J = 41: KJ = -sinK (i = 43)
            %s(41) = s(41) - f(42)*q(43);
            s(41) = s(41) - f(42)*q(41)*cosK;
            % from line 43: sinK  = sin(E);
            % line K = 43 dependencies (E 41);
            % E J = 41: KJ = cosK (i = 42)
            %s(41) = s(41) + f(43)*q(42);
            s(41) = s(41) - f(43)*q(41)*sinK;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % from line 44: tanX  = cosK - e;
            % line K = 44 dependencies (cosK 42, e 2)
            % cosK J = 42: KJ = 1                              (2nd derivatives vanish)
            % e    J =  2: KJ =-1                              (2nd derivatives vanish)
            % // from line 45: tanY  = rootfac*sinK;
            % line K = 45 dependencies (sinK 43, rootfac 14)
            % sinK J = 43: KJ = rootfac (i = 14)
            s(43) = s(43) + f(45)*q(14);
            % rootfac J = 14: KJ = sinK (i = 43)
            s(14) = s(14) + f(45)*q(43);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % //from line 46 double nu;// nu = atan2(Ytan, Xtan);                 // line 77
            % // from line 46: nu     = atan2(tanY, tanX);
            cosDenom = (tanX*tanX + tanY*tanY)^2;
            % line K = 46 dependencies (tanY 45, tanX 44)
            % %q(46) = q(46) + (q(45)*tanX - q(44)*tanY)/(tanX*tanX + tanY*tanY);
            % tanX J = 44 KJ = -tanY/(tanX*tanX + tanY*tanY);
            m00 = 2.0*tanX*tanY/cosDenom;
            m11 = (tanY*tanY - tanX*tanX)/cosDenom;
            s(44) = s(44) + f(46)*(q(44)*m00 + q(45)*m11);
            % tanY J = 45: KJ = tanX/(tanX*tanX + tanY*tanY);
            s(45) = s(45) + f(46)*(q(44)*m11 - q(45)*m00);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % // from line 47: coss  = cos(nu);
            % line K = 47 dependencies (nu 46)
            % nu J = 46: KJ = -sin(nu) = -sins (i = 48);
            %s(46) = s(46) - f(47)*q(48);
            s(46) = s(46) - f(47)*q(46)*coss;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % // from line 48: sins  = sin(nu);
            % line K = 48 dependencies (nu 46)
            % nu J = 46: KJ = cos(nu) = coss (i= 47)
            % s(46) = s(46) + f(48)*q(47);
            s(46) = s(46) - f(48)*q(46)*sins;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % // from line 49: rorbit   = a*(1.0 - e*cosK);
            % cosK: J=42 KJ = -a*e
            m03 = -e;
            m02 = -a;
            s(42) = s(42) + f(49)*(q( 3)*m03 + q( 2)*m02);
            % E: J = 41 KJ = +a*e*sin(E)
            % m15 = -e;
            % m02 = -a;
            % s(42) = s(42) + f(49)*(q(15)*m15 + q( 2)*m02);
            % a: J = 3  KJ = 1 - e*cosK
            m42 = -e;
            m02 = -cosK;
            s( 3) = s( 3) + f(49)*(q(42)*m42 + q( 2)*m02);
            % e: J = 2 KJ = -a*cosK
            m42 = -a;
            m03 = -cosK;
            s( 2) = s( 2) + f(49)*(q(42)*m42 + q(3)*m03);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % // from line 50: rorbity  = rorbit*(coss*Px + sins*Qx);
            % K = 50 Dependencies rorbit 49, sins 48, coss 47, Qx 26, Px 23
            % rorbit: J=49 KJ = coss*Px + sins*Qx (i = 48, 47, 26, 23)
            s(49) = s(49) + f(50)*(q(48)*Qx + q(47)*Px + q(26)*sins + q(23)*coss);
            % sins J=48 KJ = rorbit*Qx  (i = 49, 26)
            s(48) = s(48) + f(50)*(q(49)*Qx + q(26)*rorbit);
            % coss J=47 KJ = rorbit*Px (i = 49, 23)
            s(47) = s(47) + f(50)*(q(49)*Px + q(23)*rorbit);
            % Qx J=26 KJ = rorbit*sins (i = 49, 48)
            s(26) = s(26) + f(50)*(q(49)*sins + q(48)*rorbit);
            % Px J=23 KJ = rorbit*coss (i = 49, 47)
            s(23) = s(23) + f(50)*(q(49)*coss + q(47)*rorbit);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % // from line 51: rorbity  = rorbit*(coss*Py + sins*Qy);
            % rorbit: J=49 KJ = coss*Py + sins*Qy
            m48 = Qy;
            m47 = Py;
            m27 = sins;
            m24 = coss;
            s(49) = s(49) + f(51)*(q(48)*m48 + q(47)*m47 + q(27)*m27 + q(24)*m24);
            % sins 48 KJ = rorbit*Qy
            m49 = Qy;
            m27 = rorbit;
            s(48) = s(48) + f(51)*(q(49)*m49 + q(27)*m27);
            % sins 47 KJ = rorbit*Py
            m49 = Py;
            m24 = rorbit;
            s(47) = s(47) + f(51)*(q(49)*m49 + q(24)*m24);
            % sins 27 KJ = rorbit*sins
            m49 = sins;
            m48 = rorbit;
            s(27) = s(27) + f(51)*(q(49)*m49 + q(48)*m48);
            % sins 24 KJ = rorbit*coss
            m49 = coss;
            m47 = rorbit;
            s(24) = s(24) + f(51)*(q(49)*m49 + q(47)*m47);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % // from line 52: rorbitz  = rorbit*(coss*Pz + sins*Qz);
            % rorbit: 49 KJ = coss*Pz + sins*Qz
            m48 = Qz;
            m47 = Pz;
            m28 = sins;
            m25 = coss;
            s(49) = s(49) + f(52)*(q(48)*m48 + q(47)*m47 + q(28)*m28 + q(25)*m25);
            % sins 48 KJ = rorbit*Qz
            m49 = Qz;
            m28 = rorbit;
            s(48) = s(48) + f(52)*(q(49)*m49 + q(28)*m28);
            % sins 47 KJ = rorbit*Pz
            m49 = Pz;
            m25 = rorbit;
            s(47) = s(47) + f(52)*(q(49)*m49 + q(25)*m25);
            % sins 28 KJ = rorbit*sins
            m49 = sins;
            m48 = rorbit;
            s(28) = s(28) + f(52)*(q(49)*m49 + q(48)*m48);
            % sins 25 KJ = rorbit*coss
            m49 = coss;
            m47 = rorbit;
            s(25) = s(25) + f(52)*(q(49)*m49 + q(47)*m47);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % from line 53: rtpinv   = sqrt(1.0/p);
            % p: J = 3 KJ = -0.5*rtpinv/p
            %s(3)  = s(3) + f(53)*q(3)*0.75*rtpinv/psq;
            % p: J = 15 KJ = -0.5*rtpinv/p
            s(15)  = s(15) + f(53)*q(15)*0.75*rtpinv/psq;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %from line 54:  double vorbitx  = rtpinv*(-sins*Px + (e + coss)*Qx);  // line 54
            % K = 54 Dependencies: rtpinv = 53, sins 48, coss 47, Qx 26, Px 23, e 2
            %rtpinv: J = 53 KJ = -sins*Px + (e + coss)*Qx (i = 48, 47, 26, 23, 2)
            s(53) = s(53) + f(54)*(-q(48)*Px + q(47)*Qx + q(26)*(e + coss) -q(23)*sins + q(2)*Qx);
            % m48   = -Px;
            % m47   =  Qx;
            % m26   =  e + coss;
            % m23   = -sins;
            % m02   =  Qx;
            % s(53) =  s(53) + f(54)*(q(48)*m48 + q(47)*m47 + q(26)*m26 + q(23)*m23 + q( 2)*m02);
            % // sins: J = 48   KJ = -rtpinv*Px (i = 49, 23)
            s(48) = s(48) - f(54)*(q(53)*Px + q(23)*rtpinv);
            % m53   = -Px;
            % m23   = -rtpinv;
            % s(48) = s(48) + f(54)*(q(53)*m53 + q(23)*m23);
            % // coss: J = 47   KJ = rtpinv*Qx (i = 49, 26)
            s(47) = s(47) + f(54)*(q(53)*Qx + q(26)*rtpinv);
            % m53   = Qx;
            % m26   = rtpinv;
            % s(47) = s(47) + f(54)*(q(53)*m53 + q(26)*m26);
            %// Qx: J = 28   KJ = rtpinv*(e+coss) (i = 53, 47, 2)
            s(26) = s(26) + f(54)*(q(53)*(e+coss) + q(47)*rtpinv + q( 2)*rtpinv);
            % m53   =  e + coss;
            % m47   =  rtpinv;
            % m02   =  rtpinv;
            % s(26) = s(26) + f(54)*(q(53)*m53 + q(47)*m47 + q( 2)*m02);
            %// Px: J = 23  KJ = -rtpinv*sins (i = 53, 48)
            s(23) = s(23) - f(54)*(q(53)*sins + q(48)*rtpinv);
            % m53   = -sins;
            % m48   = -rtpinv;
            % s(24) = s(24) + f(54)*(q(53)*m53 + q(48)*m48);
            %// e: J = 2  KJ = rtpinv*Qx (i = 53, 26)
            s( 2) = s( 2) + f(54)*(q(53)*Qx + q(26)*rtpinv);
            % m53   = Qx;
            % m26   = rtpinv;
            % s( 2) = s( 2) + f(54)*(q(53)*m53 + q(26)*m26);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %// from line 55:  double vorbity  = rtpinv*(-sins*Py + (e + coss)*Qy);  // line 55
            %//rtpinv: J = 53 KJ = -sins*Py + (e + coss)*Qy
            m48   = -Py;
            m47   =  Qy;
            m27   =  e + coss;
            m24   = -sins;
            m02   =  Qy;
            s(53) =  s(53) + f(55)*(q(48)*m48 + q(47)*m47 + q(27)*m27 + q(24)*m24 + q( 2)*m02);
            % // sins: J = 48   KJ = -rtpinv*Py
            m53   = -Py;
            m24   = -rtpinv;
            s(48) = s(48) + f(55)*(q(53)*m53 + q(24)*m24);
            % // coss: J = 47   KJ = rtpinv*Qy
            m53   = Qy;
            m27   = rtpinv;
            s(47) = s(47) + f(55)*(q(53)*m53 + q(27)*m27);
            %// Qy: J = 27   KJ = rtpinv*(e+coss)
            m53   =  e + coss;
            m47   =  rtpinv;
            m02   =  rtpinv;
            s(27) = s(27) + f(55)*(q(53)*m53 + q(47)*m47 + q( 2)*m02);
            %// Py: J = 24  KJ = -rtpinv*sins
            m53   = -sins;
            m48   = -rtpinv;
            s(24) = s(24) + f(55)*(q(53)*m53 + q(48)*m48);
            %// e: J = 2  KJ = rtpinv*Qy
            m53   = Qy;
            m27   = rtpinv;
            s( 2) = s( 2) + f(55)*(q(53)*m53 + q(27)*m27);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %// from line 56:  double vorbitz  = rtpinv*(-sins*Pz + (e + coss)*Qz);  // line 56
            %//rtpinv: J = 53 KJ = -sins*Pz + (e + coss)*Qz
            m48   = -Pz;
            m47   =  Qz;
            m28   =  e + coss;
            m25   = -sins;
            m02   =  Qz;
            s(53) =  s(53) + f(56)*(q(48)*m48 + q(47)*m47 + q(28)*m28 + q(25)*m25 + q( 2)*m02);
            % // sins: J = 48   KJ = -rtpinv*Pz
            m53   = -Pz;
            m25   = -rtpinv;
            s(48) = s(48) + f(56)*(q(53)*m53 + q(25)*m25);
            % // coss: J = 47   KJ = rtpinv*Qz
            m53   = Qz;
            m28   = rtpinv;
            s(47) = s(47) + f(56)*(q(53)*m53 + q(28)*m28);
            %// Qz: J = 28   KJ = rtpinv*(e+coss)
            m53   =  e + coss;
            m47   =  rtpinv;
            m02   =  rtpinv;
            s(28) = s(28) + f(56)*(q(53)*m53 + q(47)*m47 + q( 2)*m02);
            %// Pz: J = 25  KJ = -rtpinv*sins
            m53   = -sins;
            m48   = -rtpinv;
            s(25) = s(25) + f(56)*(q(53)*m53 + q(48)*m48);
            %// e: J = 2  KJ = rtpinv*Qz
            m53   = Qz;
            m28   = rtpinv;
            s( 2) = s( 2) + f(56)*(q(53)*m53 + q(28)*m28);
            if LineCheck > 60
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Escobal Potential starts here
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line 61: X    = rorbitx = Kvec(1);                     % line  61
                % All 2nd Derivatives vanish
                % from line 62: Y    = rorbity = Kvec(2);                     % line  62
                % All 2nd Derivatives vanish
                % from line 63: Z    = rorbitz = Kvec(3);                     % line  63
                % All 2nd Derivatives vanish
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line 64: R2   = X^2 + Y^2 + Z^2;                         % line  64
                % K = 64
                %// X: J = 61  KJ = 2.0*X
                s(61) = s(61) + 2.0*f(64)*q(61);
                %// Y: J = 62  KJ = 2.0*Y
                s(62) = s(62) + 2.0*f(64)*q(62);
                %// Z: J = 63  KJ = 2.0*Z
                s(63) = s(63) + 2.0*f(64)*q(63);
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line 65: R1   = rorbit;                       % line  64
                % All 2nd Derivatives vanish
                % from line 65: R1   = sqrt(R2);                                % line  65
                %// R2: J = 64  KJ = 0.5/R1   check this
                s(64) = s(64) - f(65)*q(64)*0.25/R3;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line 64: R2   = R1^2;                                % line  65
                %// R1: J = 65  KJ = 2.0*R1   I = 65
                %s(65) = s(65) + 2.0*f(64)*q(65);
                
                % from line 66: R3   = R1*R2;
                %// R1: J = 65 KJ = R2  I = 64
                s(65) = s(65) + f(66)*q(64);
                %// R2: J = 64 KJ = R1 I = 65
                s(64) = s(64) + f(66)*q(65);
                
                % from line 67: R4   = R1*R3;                                   % line  67
                %// R1: J = 65 KJ = R3  I = 66
                s(65) = s(65) + f(67)*q(66);
                %// R3: J = 66 KJ = R1
                s(66) = s(66) + f(67)*q(65);
                
                % from line 68: R5   = R1*R4;                                   % line  68
                % K = 68
                %// R1: J = 65  KJ = R4 I = 67
                s(65) = s(65) + f(68)*q(67);
                %// R4: J = 67  KJ = RI I = 65
                s(67) = s(67) + f(68)*q(65);
                
                % from line 69: R6   = R1*R5;                                   % line  69
                %// R1: J = 65  KJ = R5 I = 68
                s(65) = s(65) + f(69)*q(68);
                %// R5: J = 68  KJ = R1 I = 65
                s(68) = s(68) + f(69)*q(65);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line 70: sd1  = rorbitz/R1;
                %//rorbitz: J = 52 KJ = 1/R1 I = 65
                s(52) = s(52) - f(70)*q(65)/R2;
                %// R1: J = 65 KJ = -Z*R1^(-2) I = 52, 65
                s(65) = s(65) + f(70)*(2.0*q(65)*sd1 - q(52))/R2;
                
                % from line 70: sd1  = Z/R1;
                %// Z: J = 63 KJ = 1/R1 I = 65
                %s(63) = s(63) - f(70)*q(65)/R2;
                %// R1: J = 65 KJ = -Z*R1^(-2) I = 63, 65
                %s(65) = s(65) + f(70)*(2.0*q(65)*sd1 - q(63))/R2;
                % from line 70: sd1  = sin(I)*sin(omega + nu);
                %q(70)  = q(70) + q(18)*sinomnu + sinI*cosomnu*(q(5) + q(46));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % from line 71: sd2  = sd1^2;
                %//sd1: J = 70 KJ = 2.0*sd1 I = 70
                s(70) = s(70) + f(71)*q(70)*2.0;
                
                % from line 72: sd3  = sd1*sd2;
                %//sd1: J = 70 KJ = sd2 I = 71
                s(70) = s(70) + f(72)*q(71);
                %//sd2: J = 71 KJ = sd1 I = 70
                s(71) = s(71) + f(72)*q(70);
                
                % from line 73: sd4  = sd1*sd3;
                %//sd1: J = 70 KJ = sd3 I = 72
                s(70) = s(70) + f(73)*q(72);
                %//sd3: J = 72 KJ = sd1 I = 70
                s(72) = s(72) + f(73)*q(70);
                
                % from line 74: sd5  = sd1*sd4;
                %//sd1: J = 70 KJ = sd4 I = 73
                s(70) = s(70) + f(74)*q(73);
                %//sd4: J = 73 KJ = sd1 I = 70
                s(73) = s(73) + f(74)*q(70);
                
                % from line 75: sd6  = sd1*sd5;
                %//sd1: J = 70 KJ = sd5 I = 74
                s(70) = s(70) + f(75)*q(74);
                %//sd5: J = 73 KJ = sd1 I = 70
                s(74) = s(74) + f(75)*q(70);
                
                % All Second Derivatives Vanish here
                % from line 76: F2   =  1.0      -  3.0*sd2;
                % from line 77: F3   =  3.0*sd1  -  5.0*sd3;
                % from line 78: F4   =  3.0      -  30.*sd2  + 35*sd4;
                % from line 79: F5   =  15.*sd1  -  70.*sd3  + 63*sd5;
                % from line 80: F6   =  5.       - 105.*sd2  + 315*sd4 - 231*sd6;
                
                % from line 81: V1    = 1/R1;          Kepler Potential                 % line 21
                %//R1: J = 65 KJ = -1/R1^2 I = 65
                s(65) = s(65) + 2.0*f(81)*q(65)/R3;
                
                % from line 82: V2      = C2*F2/R2;
                %//F2 J = 76 KJ = C2/R2 I = 64
                s(76) = s(76) - f(82)*q(64)*C2/R4;
                %//R2 J = 64 KJ = -C2*F2/R2^2 I = 76, 67
                s(64) = s(64) + f(82)*(2.0*q(64)*V2 - C2*q(76))/R4;
                
                % from line 83: V3      = C3*F3/R3;
                %//F3 J = 77 KJ = C3/R3 I = 66
                s(77) = s(77) - f(83)*q(66)*C3/R6;
                %//R3 J = 66 KJ = -C3*F3/R3^2 I = 77, 66
                s(66) = s(66) + f(83)*(2.0*q(66)*V3 - C3*q(77))/R6;
                
                % from line 84: V4      = C4*F4/R4;
                %//F4 J = 78 KJ = C4/R4 I = 67
                s(78) = s(78) - f(84)*q(67)*C4/R4^2;
                %//R4 J = 67 KJ = -C4*F4/R4^-2 I = 78, 67
                s(67) = s(67) + f(84)*(2.0*q(67)*V4 - C4*q(78))/R4^2;
                
                % from line 85: V5      = C5*F5/R5;
                %//F5 J = 79 KJ = C5/R5 I = 68
                s(79) = s(79) - f(85)*q(68)*C5/R5^2;
                %//R5 J = 68 KJ = -C5*F5/R5^-2 I = 79, 68
                s(68) = s(68) + f(85)*(2.0*q(68)*V5 - C5*q(79))/R5^2;
                
                % from line 86: V6      = C6*F6/R6;
                %//F6 J = 80 KJ = C5/R6 I = 69
                s(80) = s(80) - f(86)*q(69)*C6/R6^2;
                %//R6 J = 69 KJ = -C6*F6/R6^-2 I = 80, 69
                s(69) = s(69) + f(86)*(2.0*q(69)*V6 - C6*q(80))/R6^2;
                
                % All Second Derivative Vanish
                % from line 87: V       = mu*V1;
                % from line 88: D       = V2 + V3 + V4 + V5 + V6;
                
                % %from line 89: R       = mu*D;
                % %q(89) = q(89) + q(88)*mu;
                % % from line 90: Phi     = V + R;
                % %q(90) = q(90) + q(87) + q(89);
                
                % from line 90: Phi     = V*(CK + D);
                %//V J = 87  KJ = (CK + D)  I = 88
                s(87) = s(87) + f(90)*q(88);
                %//D J = 88  KJ = V  I = 87
                s(88) = s(88) + f(90)*q(87);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if LineCheck > 90
                    fac1 = this.fac1;
                    fac2 = this.fac2;
                    fac3 = this.fac3;
                    fac4 = this.fac4;
                    fac5 = this.fac5;
                    Mwe  = this.Mwe;
                    MeM  = this.MeM;
                    MaM  = this.MaM;
                    MWi  = this.MWi;
                    Miw  = this.Miw;
                    MLagrange = this.MLagrange;
                    % from line 91: fac1 = n*a
                    % K = 91
                    %//n J = 16, KJ = a I = 3
                    s(16) = s(16) + f(91)*q( 3);
                    %//a J = 3, KJ = n I = 16
                    s( 3) = s( 3) + f(91)*q(16);
                    
                    % from line 92: fac2 = fac1*a % line 92
                    % K = 92
                    %//a J = 3, KJ = fac1, I = 91)
                    s( 3) = s( 3) + f(92)*q(91);
                    %// fac1 J = 91 KJ = a I = 3
                    s(91) = s(91) + f(92)*q( 3);
                    
                    % from line 93: fac3 = fac2*e % line 93
                    % K = 93
                    %//e J = 2, KJ = fac2 I = 92
                    s( 2) = s( 2) + f(93)*q(92);
                    %//fac2 J = 92, KJ = e I = 2
                    s(92) = s(92) + f(93)*q( 2);
                    
                    % from line 94: fac4 = fac2*rootfac  % line 94
                    % K = 94
                    %//rootfac J = 14, KJ = fac2 I = 92
                    s(14) = s(14) + f(94)*q(92);
                    %//fac2 J = 92, KJ = rootfac I = 14
                    s(92) = s(92) + f(94)*q(14);
                    
                    % from line 95: fac5 = fac4*sinI        % line 95
                    % K = 95
                    %//I J = 4, KJ = fac4*cosI, I = 4, 94
                    s( 4) = s( 4) + f(95)*(q(94)*cosI - q( 4)*fac4*sinI);
                    %// fac4 J = 94, KJ = sinI, I = 4
                    s(94) = s(94) + f(95)*q( 4)*cosI;
                    
                    %from line 96: Mwe  = rootfac/fac3      % line 96
                    %K = 96
                    %//rootfac J = 14, KJ = 1/fac3 I = 93
                    s(14) = s(14) - f(96)*q(93)/fac3^2;
                    %//fac3 J = 93, KJ = -rootfac/fac3^2, I = 14, 93
                    s(93) = s(93)  + f(96)*(2.0*q(93)*Mwe - q(14))/fac3^2;
                    
                    % from line 97: MeM  = fac/fac3    % line 97
                    %// K = 97
                    %//fac J = 13, KJ = 1/fac3
                    s(13) = s(13) - f(97)*q(93)/fac3^2;
                    %//fac3 J = 93, KJ = -fac/fac3^2 I = 13, 93
                    s(93) = s(93) + f(97)*(2.0*q(93)*MeM - q(13))/fac3^2;
                    
                    % from line 98: MaM  = 2.0/fac1         % line 98
                    %//fac1 J = 91, KJ = -2.0/fac1^2 I = 91
                    s(91) = s(91) + f(98)*2.0*q(91)*MaM/fac1^2;
                    
                    % from line 99: MWi  = 1/fac5       % line 99
                    %//fac5 J = 95, KJ = -1/fac5^2 I = 95
                    s(95) = s(95) + 2.0*f(99)*q(95)/fac5^3;
                    
                    % from line 100: Miw  = cosI*MWi    % line 100
                    %//I J = 4, KJ = -sinI*MWi I = 4, 99
                    s( 4) = s( 4) - f(100)*(cosI*MWi*q( 4) + sinI*q(99));
                    %//MWi J = 99, KJ = cosI, I = 4
                    s(99) = s(99) - f(100)*sinI*q( 4);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % All Second Derivative Vanish
                    if LineCheck > 100
                        %los       = Rextrap - state_vector_sat(1:3);    % lines 101-103
                        %lx = los(1);                                    % line 101
                        %ly = los(2);                                    % line 102
                        %lz = los(3);                                    % line 103
                        % from line K = 85: V5      = C5*F5/R5;
                        %//F5 J = 79 KJ = C5/R5 I = 68
                        %s(J) = s(J) - f(K)*q(I)*C5/R5^2;
                        
                        % from line 104: lperpSq  = lx^2 + ly^2;         % line 104
                        % K = 104:
                        %//lx J = 101, KJ = 2.0*lx, I = 101
                        s(101)  = s(101) + 2.0*f(104)*q(101);
                        %//ly J = 102, KJ = 2.0*ly, I = 102
                        s(102)  = s(102) + 2.0*f(104)*q(102);
                        
                        % from line: 105: rangeSq  = lperpSq + lz^2;     % line 105
                        % K = 105
                        %//lz J = 103, KJ = 2.0*lz, I = 103
                        s(103)  = s(103) + 2.0*f(105)*q(103);
                        %//lperpSq J = 104, KJ = 1
                        % Second Derivative vanishes
                        
                        % from line 106: range    = sqrt(rangeSq);       % line 106
                        % K = 106
                        %// rangeSq J = 105, KJ = 0.5/range, I = 106
                        s(105)    = s(105) - 0.25*f(106)*q(105)/range^3;
                        
                        % from line 107: zLos     = lz/range;
                        % K = 107
                        %//lz J = 103, KJ = 1/range I = 106
                        s(103) = s(103) - f(107)*q(106)/rangeSq;
                        %//range J = 106, KJ = -lz/range^2 I = 103, 106
                        s(106) = s(106) + f(107)*(2.0*q(106)*zLos - q(103))/rangeSq;
                        
                        % from line 108: thetaLos = acos(zLos);          % line 108
                        % K = 108
                        %//zLos J = 107, KJ = -1/sqrt(1-zLos^2) I = 107
                        denom = sqrt(1.0-zLos^2);
                        s(107) = s(107) - f(108)*q(107)*zLos/denom^3;
                        
                        %phiLos   = atan2(ly,lx);                        % line 109
                        % K = 109
                        %//lx J = 101, KJ = -ly/(lx^2 + ly^2), I = 101, 102
                        s(101) = s(101) + f(109)*(2.0*q(101)*lx*ly + q(102)*(ly^2 - lx^2))/lperpSq^2;
                        %//ly J = 102, KJ = lx/(lx^2 + ly^2), I = 101, 102
                        s(102) = s(102) + f(109)*(q(101)*(ly^2 - lx^2) - 2.0*q(102)*lx*ly)/lperpSq^2;
                        
                        % All Second Derivative Vanish
                        %vlos     = Vextrap - state_vector_sat(4:6);       % lines 110-112
                        %vx       = vlos(1);                               % line 110
                        %vy       = vlos(2);                               % line 111
                        %vz       = vlos(3);                               % line 112
                        
                        % from line 113: ldotv    = lx*vx + ly*vy + lz*vz; % line 113
                        % K = 113
                        %//lx J = 101, KJ = vx I = 110
                        s(101)    = s(101) + f(113)*q(110);
                        %//vx J = 110, KJ = lx, I = 101
                        s(110)    = s(110) + f(113)*q(101);
                        
                        %//ly J = 102, KJ = vy I = 111
                        s(102)    = s(102) + f(113)*q(111);
                        %//vy J = 111, KJ = ly, I = 102
                        s(111)    = s(111) + f(113)*q(102);
                        
                        %//lz J = 103, KJ = vz I = 112
                        s(103)    = s(103) + f(113)*q(112);
                        %//vz J = 112, KJ = lz, I = 102
                        s(112)    = s(112) + f(113)*q(103);
                        
                        % from line 114: ldot     = ldotv/range;         % line 114
                        % K = 114
                        %//ldotv J = 113, KJ = 1/range I = 106
                        s(113) = s(113) - f(114)*q(106)/rangeSq;
                        %//range J = 106, KJ = -ldotv/rangeSq I = 106, 113
                        s(106) = s(106) + f(114)*(2.0*q(106)*ldot - q(113))/rangeSq;
                        
                        % from line 115: lperp    = sqrt(lperpSq);       % line 115
                        % K = 115
                        %//lperpSq J = 104, KJ = 0.5/sqrt(lperpSq) I = 104
                        s(104) = s(104) - 0.25*f(115)*q(104)/(lperpSq*lperp);
                        
                        % from line 116: thdotfac1 = range*vz - lz*ldot; % line 116
                        % K = 116
                        %//range J = 106, KJ = vz, I = 103
                        s(106) = s(106) + f(116)*q(112);
                        %//vz J = 112, KJ = range, I = 106
                        s(112) = s(112) + f(116)*q(106);
                        
                        %//lz J = 103, KJ = -ldot, I = 114
                        s(103) = s(103) - f(116)*q(114);
                        %//ldot J = 114, KJ = -vz, I = 103
                        s(114) = s(114) - f(116)*q(103);
                        
                        % from line 117: thdotfac2 = thdotfac1/range;    % line 117
                        %//thdotfac1 J = 116, KJ = 1/range I = 106
                        s(116) = s(116) - f(117)*q(106)/rangeSq;
                        %//range J = 106, KJ = -thdotfac1/rangeSq I = 106, 116
                        s(106) = s(106) + f(117)*(2.0*q(106)*thdotfac2 - q(116))/rangeSq;
                        
                        % from line 118: thdot     = -thdotfac2/lperp;   % line 118
                        %//thdotfac2 J = 117, KJ = -1/lperp I = 115
                        s(117) = s(117) + f(118)*q(115)/lperpSq;
                        %//lperp J = 115, KJ = thdotfac2/lperpSq I = 115, 117
                        s(115) = s(115) + f(118)*(q(117) + 2.0*q(115)*thdot)/lperpSq;
                        
                        % from line 119: phiDotNum = lx*vy - ly*vx;      % line 119
                        %K = 119
                        %//lx J = 101, KJ = vy I = 111
                        s(101) = s(101) + f(119)*q(111);
                        %//vy J = 111, KJ = lx I = 101
                        s(111) = s(111) + f(119)*q(101);
                        
                        %//ly J = 102, KJ = -vx I = 110
                        s(102) = s(102) - f(119)*q(110);
                        %//vx J = 110, KJ = -ly I = 102
                        s(110) = s(110) - f(119)*q(102);
                        
                        %phidot    = phiDotNum/lperpSq;                  % line 120
                        %K = 120
                        %//phiDotNum J = 119, KJ = 1/lperpSq I = 104
                        s(119) = s(119) - f(120)*q(104)/lperpSq^2;
                        %//lperpSq J = 104, KJ = -phiDotNum/lperpSq^2 I = 104, 119
                        s(104) = s(104) + f(120)*(2.0*q(104)*phidot - q(119))/lperpSq^2;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %cosThetaLos = cos(thetaLos);                        % line 121
                        % K = 121
                        %//thetaLos J = 108, KJ = -sin(thetaLos), I = 108
                        s(108) = s(108) - f(121)*q(108)*cosThetaLos;
                        %sinThetaLos = sin(thetaLos);                        % line 122
                        % K = 122
                        %//thetaLos J = 108, KJ = cos(thetaLos), I = 108
                        s(108) = s(108) - f(122)*q(108)*sinThetaLos;
                        %phiDiff     = phiLos - phiM;                        % line 123
                        % K = 123
                        %//phiLos J = 109, KJ = 1 -- All Second Derivatives vanish
                        %cosPhiDiff  = cos(phiDiff);                         % line 124
                        % K = 124
                        %//phiDiff J = 123, KJ = -sin(phiDiff) I = 123
                        s(123) = s(123) - f(124)*q(123)*cosPhiDiff;
                        %sinPhiDiff  = sin(phiDiff);                         % line 125
                        % K = 125
                        %//phiDiff J = 123, KJ = cos(phiDiff), I = 123
                        s(123) = s(123) - f(125)*q(123)*sinPhiDiff;
                        %delUV       = cosThetaLos*sinThetaM*cosPhiDiff;    % line 126
                        % K = 126
                        %cosThetaLos J = 121 KJ = sinThetaM*cosPhiDiff I = 124
                        s(121) = s(121) + f(126)*q(124)*sinThetaM;
                        %//cosPhiDiff J = 124 KJ = cosThetaLos*sinThetaM I = 121
                        s(124) = s(124) + f(126)*q(121)*sinThetaM;
                        %Uview       = sinThetaLos*cosThetaM - delUV;        % line 127
                        % K = 127
                        %//delUV K = 126, KJ = 1 -- All Second Derivatives vanish
                        %//sinThetaLos J = 122 KJ = cosThetaM -- All second derivatives vanish
                        %Vview       = sinThetaM*sinPhiDiff;                 % line 128
                        % K = 128
                        %//sinPhiDiff J = 123, KJ = sinThetaM -- All second derivatives vanish
                        
                        % from line 129: xLos     = lx/range;
                        % K = 129
                        %//lx J = 101, KJ = 1/range I = 106
                        s(101) = s(101) - f(129)*q(106)/rangeSq;
                        %//range J = 106, KJ = -lz/range^2 I = 101, 106
                        s(106) = s(106) + f(129)*(2.0*q(106)*xLos - q(101))/rangeSq;
                        
                        % from line 130: yLos     = ly/range;
                        % K = 130
                        %//ly J = 102, KJ = 1/range I = 106
                        s(102) = s(102) - f(130)*q(106)/rangeSq;
                        %//range J = 106, KJ = -ly/range^2 I = 102, 106
                        s(106) = s(106) + f(130)*(2.0*q(106)*yLos - q(102))/rangeSq;
                        
                        % from line 131: Collinearity = 1 - xLos*xLosM - yLos*yLosM - zLos*zLosM;
                        %  all second derivatives vanish
                        % from line 133: ChiSqLos = losDiff'*InvCovM*losDiff;
                        % K = 133
                        % xLos  J = 129 KJ = xUnit'*InvCovM*losDiff
                        % I = xLos 129, yLos 130, zLos 107
                        s(129)  = s(129) + f(133)*(  q(129)*(xUnit'*InvCovM*xUnit) ...
                            + q(130)*(xUnit'*InvCovM*yUnit) ...
                            + q(107)*(xUnit'*InvCovM*zUnit));
                        
                        % yLos  J = 130 KJ = yUnit'*InvCovM*losDiff
                        % I = xLos 129, yLos 130, zLos 107
                        s(130)  = s(130) + f(133)*(  q(129)*(yUnit'*InvCovM*xUnit) ...
                            + q(130)*(yUnit'*InvCovM*yUnit) ...
                            + q(107)*(yUnit'*InvCovM*zUnit));
                        
                        % zLos  J = 107 KJ = xUnit'*InvCovM*losDiff
                        % I = xLos 129, yLos 130, zLos 107
                        s(107)  = s(107) + f(133)*(  q(129)*(zUnit'*InvCovM*xUnit) ...
                            + q(130)*(zUnit'*InvCovM*yUnit) ...
                            + q(107)*(zUnit'*InvCovM*zUnit));
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if LineCheck > 133
                            % from line 142: Aeq      = omega + Omega; % line 142
                            % All second derivatives vanish
                            % from line 143: cosAeq   = cos(Aeq);
                            %//Aeq J = 142: KJ = -sin(Aeq)
                            s(142) = s(142) - f(143)*q(82)*cosAeq;
                            % from line 144: sinAeq   = sin(Aeq);  % line 144
                            %//Aeq J = 142: KJ = cos(Aeq)
                            s(142) = s(142) - f(144)*q(142)*sinAeq;
                            % from line 145: tanHalf  = tan(Inclination/2.0); % line 145
                            %//Inclination J = 4: KJ = (1 + tanHalf^2)/2
                            s( 4) = s( 4) + 0.5*f(145)*q( 4)*tanHalf*(1.0 + tanHalf^2);
                            % from line 146: Feq      = e*cosAeq;             % line 146
                            %// e J =   2: KJ = cosAeq I = 143
                            s( 2) = s( 2) + f(146)*q(143);
                            %//cosAeq J = 143: KJ = e
                            s(143) = s(143) + f(146)*q( 2);
                            % from line 147: Geq      = e*sinAeq;             % line 147
                            %//e J =  2: KJ = sinAeq
                            s( 2) = s( 2) + f(147)*q(144);
                            %//sinAeq J = 144: KJ = e
                            s(144) = s(144) + f(147)*q( 2);
                            
                            % from line 148: Heq      = tanHalf*cos(Omega);    % line 148
                            %q(148) = q(148) + q(145)*cosO - q(47)*tanHalf*sinO;
                            % from line 149: Keq      = tanHalf*sin(Omega);    % line 149
                            %q(149) = q(149) + q(145)*sinO + q(47)*tanHalf*cosO;
                            
                            % from line 148: Heq      = tanHalf*cosO;     % line 148
                            %//Omega J = 6: KJ = -tanHalf*sinO
                            s( 6) = s( 6) - f(148)*(q( 6)*Heq + q(145)*sinO);
                            %//tanHalf J = 145: KJ = cosO
                            s(145) = s(145) - f(148)*q( 6)*sinO;
                            % from line 149: Keq      = tanHalf*sinO;     % line 149
                            %//Omega J =  6: KJ = tanHalf*cosO
                            s( 6) = s( 6) + f(149)*(-q( 6)*Keq + q(145)*cosO);
                            %//tanHalf J = 145: KJ = sinO
                            s(145) = s(145) + f(149)*q( 6)*cosO;
                            % from line 150: Leq    = Aeq + nu;               % line 150
                            % All second derivatives vanish
                            % from line 151: CosL     = cos(Leq);             % line 151
                            %//Leq J = 150: KJ = -sin(Leg) I = 150
                            s(150) = s(150) - f(151)*q(150)*CosL;
                            % from line 152: SinL     = sin(Leq);             % line 152
                            %//LEQ J = 150: KJ = cos(Leq), I = 150
                            s(150) = s(150) - f(152)*q(150)*SinL;
                            % from line 153: alphaSq = Heq^2 - Keq^2;         % line 153
                            %//Heq J = 148: KJ = 2.0*Heq, I = 148
                            s(148) = s(148) + 2.0*f(153)*q(148);
                            %//Keq J = 89: KJ = -2.0*Keq, I = 149
                            s(149) = s(149) - 2.0*f(153)*q(149);
                            % from line 154: Seq     = 1.0 + Heq^2 + Keq^2;   % line 154
                            %//Heq J = 148: KJ = 2.0*Heq, I = 148
                            s(148) = s(148) + 2.0*f(154)*q(148);
                            %//Keq J = 149: KJ = 2.0*Keq, I = 149
                            s(149) = s(149) + 2.0*f(154)*q(149);
                            % from line 155: Weq = 1.0 + Feq*CosL + Geq*SinL; % line 155
                            %//Feq J = 146: KJ = CosL, I = 151
                            s(146) = s(146) + f(155)*q(151);
                            %//Geq J = 147: KJ = SinL, I = 152
                            s(147) = s(147) + f(155)*q(152);
                            %//CosL J = 151: KJ = Feq, I = 146
                            s(151) = s(151) + f(155)*q(146);
                            %//SinL J = 152: KJ = Geq, I = 147
                            s(152) = s(152) + f(155)*q(147);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % from line 156: Req     = p/Weq;              % line 156
                            %//p J = 15: KJ = 1/Weq
                            s(15) = s(15) -f(156)*q(155)/Weq^2;
                            %//Weq J = 155: KJ = -p/Weq^2, I = 15, 155
                            s(155) = s(155) + f(156)*(2.0*q(155)*Req - q(15))/Weq^2;
                            
                            % from line 157: RovS    = Req/Seq;            % line 157
                            %//Req J = 156: KJ = 1/Seq
                            s(156) = s(156) -f(157)*q(154)/Seq^2;
                            %//Seq J = 154: KJ = -Req/Seq^2, I = 156, 154
                            s(154) = s(154) + f(157)*(2.0*q(154)*RovS - q(156))/Seq^2;
                            %    q(157) = q(157) + q(156)/Seq - q(154)*RovS/Seq;
                            
                            % from line 158: srtpinv = rtpinv/Seq;         % line 158
                            %//rtpinv J = 53: KJ = 1/Seq, I = 154
                            s(53) = s(53) - f(158)*q(154)/Seq^2;
                            %//Seq J = 154: KJ = -rtpinv/Seq^2, I = 53, 154
                            s(154) = s(154) + f(158)*(2.0*q(154)*srtpinv - q(53))/Seq^2;
                            % q(158) = q(158) + q(53)/Seq - q(154)*srtpinv/Seq;
                            
                            % from line 159: HK  = Heq*Keq;            % line 159
                            %//Heq J = 148: KJ = Keq
                            s(148) = s(148) + f(159)*q(149);
                            s(149) = s(149) + f(159)*q(148);
                            %  q(159) = q(159) + q(148)*Keq + q(149)*Heq;
                            
                            % from line 160: OnePalphaSqCosL = (1.+alphaSq)*CosL;% line 160
                            %//alphaSq J = 153: KJ = CosL, I = 151
                            s(153) = s(153) + f(160)*q(151);
                            %//CosL J = 151: KJ = (1 + alphaSq), I = 153
                            s(151) = s(151) + f(160)*q(153);
                            % q(160) = q(160) + q(151)*(1.0 + alphaSq) + q(153)*CosL;
                            
                            % from line 161: OneMalphaSqCosL = (1.-alphaSq)*CosL;% line 161
                            %//alphaSq J = 153: KJ =-CosL, I = 151
                            s(153) = s(153) - f(161)*q(151);
                            %//CosL J = 151: KJ = (1 - alphaSq), I = 153
                            s(151) = s(151) - f(161)*q(153);
                            % q(161) = q(161) + q(151)*(1.0 - alphaSq) - q(153)*CosL;
                            
                            % from line 162: OnePalphaSqSinL = (1.+alphaSq)*SinL;% line 162
                            %//alphaSq J = 153: KJ = SinL, I = 152
                            s(153) = s(153) + f(162)*q(152);
                            %//SinL J = 152: KJ = (1 + alphaSq), I = 153
                            s(152) = s(152) + f(162)*q(153);
                            % q(102) = q(102) + q(92)*(1.0 + alphaSq) + q(93)*SinL;
                            
                            % from line 163: OneMalphaSqSinL = (1.-alphaSq)*SinL;% line 163
                            %//alphaSq J = 153: KJ =-SinL, I = 152
                            s(153) = s(153) - f(163)*q(152);
                            %//SinL J = 152: KJ = (1 - alphaSq), I = 153
                            s(152) = s(152) - f(163)*q(153);
                            % q(163) = q(163) + q(152)*(1.0 - alphaSq) - q(153)*SinL;
                            
                            % from line 164: Xfac    = OnePalphaSqCosL + 2.0*HK*SinL; % line 164
                            %//HK J = 159: KJ = 2.0*SinL
                            s(159) = s(159) + f(164)*2.0*q(152);
                            %//SinL J = 152: KJ = 2.0*HK
                            s(152) = s(152) + f(164)*2.0*q(159);
                            % OnePalphaSqCosL J = 100 -- All second derivatives vanish
                            % q(164) = q(164) + 2.0*q(152)*HK + 2.0*q(159)*SinL + q(100);
                            
                            % from line 165: Yfac    = OneMalphaSqSinL + 2.0*HK*CosL; % line 165
                            %//HK J = 159: KJ = 2.0*CosL
                            s(159) = s(159) + f(165)*2.0*q(151);
                            %//CosL J = 151: KJ = 2.0*HK
                            s(151) = s(151) + f(165)*2.0*q(159);
                            % OnePalphaSqCosL J = 100 -- All second derivatives vanish
                            % q(165) = q(165) + 2.0*q(151)*HK + 2.0*q(159)*CosL + q(163);
                            
                            % from line 166: Zfac    = Heq*SinL - Keq*CosL;          % line 166
                            s(148) = s(148) + f(166)*q(152);
                            s(149) = s(149) - f(166)*q(151);
                            s(151) = s(151) - f(166)*q(149);
                            s(152) = s(152) + f(166)*q(148);
                            % q(166) = q(148) + q(148)*SinL - q(149)*CosL  - q(151)*Keq  + q(152)*Heq;
                            
                            % from line 167: VXfac =  OnePalphaSqSinL - 2.0*HK*(Feq+CosL) + Geq*(1+alphaSq);% line 167
                            %//Feq J = 146: KJ = -2.0*HK, I = 159
                            s(146) = s(146) - 2.0*f(167)*q(159);
                            %//Geq J = 147: KJ = (1+alphaSq), I = 153
                            s(147) = s(147) + f(167)*q(153);
                            %//CosL J = 151: KJ = -2.0*HK, I = 159
                            s(151) = s(151) - 2.0*f(167)*q(159);
                            %//alphaSq J = 153: KJ = Geq, I =147
                            s(153) = s(153) + f(167)*q(147);
                            %//HK J = 159: KJ = -2.0*(Feq+CosL)
                            s(159) = s(159) -2.0*f(167)*(q(146) + q(151));
                            %//OnePalphaSqSinL J=102 -- All second derivatives vanish
                            % q(167) = q(167) -2.0*q(146)*HK +q(147)*(1.0 + alphaSq) -2.0*q(151)*HK +q(153)*Geq -2.0*q(159)*(Feq+CosL) +q(102);
                            
                            % from line 168: VYfac = -OneMalphaSqCosL + 2.0*HK*(Geq+SinL) + Feq*(alphaSq-1) ;% line 168
                            %//Feq J = 146: KJ = (alphaSq-1), I = 153
                            s(146) = s(146) + f(168)*q(153);
                            %//Geq J = 147: KJ = 2.0*HK, I = 159
                            s(147) = s(147) + 2.0*f(168)*q(159);
                            %//SinL J = 152: KJ = 2.0*HK, I = 159
                            s(152) = s(152) + 2.0*f(168)*q(159);
                            %//alphaSq J = 153: KJ = Feq, I =146
                            s(153) = s(153) + f(168)*q(146);
                            %//HK J = 159: KJ = 2.0*(Geq+SinL)
                            s(159) = s(159) +2.0*f(168)*(q(147) + q(152));
                            %//OneMalphaSqCosL J=101 -- All second derivatives vanish
                            %q(168) = q(168) +2.0*q(147)*HK +q(146)*(alphaSq - 1.0) +2.0*q(152)*HK +q(153)*Feq +2.0*q(159)*(Geq+SinL) -q(101);
                            
                            % from line 169: VZfac   =  Heq*(Feq + CosL) + Keq*(Geq + SinL) % line 169
                            %//Feq J = 146: KJ = Heq, I = 148
                            s(146) = s(146) + f(169)*q(148);
                            %//Geq J = 147: KJ = Keq, I = 149
                            s(147) = s(147) + f(169)*q(149);
                            %//Heq J = 148: KJ = (Feq + CosL)
                            s(148) = s(148) + f(169)*(q(146) + q(151));
                            %//Keq J = 149: KJ = (Geq + SinL)
                            s(149) = s(149) + f(169)*(q(147) + q(152));
                            %//CosL J = 151: KJ = Heq, I = 148
                            s(151) = s(151) + f(169)*q(148);
                            %//SinL J = 152: KJ = Keq, I = 149
                            s(152) = s(152) + f(169)*q(149);
                            %q(169) = q(169) +q(146)*Heq +q(147)*Keq +q(148)*(Feq+CosL) +q(149)*(Geq+SinL) +q(151)*Heq +q(152)*Keq;
                            
                            % from line 170: Xeq     =     RovS*Xfac;     % line 170
                            %//RovS J = 157: KJ = Xfac, I = 164
                            s( 157) = s( 157) + f(170)*q(164);
                            %//Xfac J = 164: KJ = RovS, I = 157
                            s(164) = s(164) + f(170)*q(157);
                            %q(170) = q(170) + q(157)*Xfac + q(164)*RovS;
                            
                            % from line 171: Yeq     =     RovS*Yfac;     % line 171
                            %//RovS J = 157: KJ = Yfac, I = 165
                            s( 157) = s( 157) + f(171)*q(165);
                            %//Yfac J = 165: KJ = RovS, I = 157
                            s(165) = s(165) + f(171)*q(157);
                            %q(171) = q(171) + q(157)*Yfac + q(165)*RovS;
                            
                            %from line 172: Zeq      = 2.0*RovS*Zfac;     % line 172
                            %//RovS J = 157: KJ = 2.0*Zfac, I = 166
                            s( 157) = s( 157) + 2.0*f(172)*q(166);
                            %//Zfac J = 166: KJ = 2.0*RovS, I = 157
                            s(166) = s(166) + 2.0*f(172)*q(157);
                            %q(172) = q(172) + 2.0*(q(157)*Zfac + q(166)*RovS);
                            
                            % from line 173: VXeq    =    -srtpinv*VXfac; % line 173
                            %//srtpinv J = 158: KJ = -VXfac, I = 167
                            s( 158) = s( 158) - f(173)*q(167);
                            %//VXfac J = 167: KJ = -srtpinv, I = 158
                            s(167) = s(167) - f(173)*q( 158);
                            %q(173) = q(173) - q(158)*VXfac - q(167)*srtpinv;
                            
                            % from line 174: VYeq    =    -srtpinv*VYfac; % line 174
                            %//srtpinv J = 158: KJ = -VYfac, I = 168
                            s( 158) = s( 158) - f(174)*q(168);
                            %//VYfac J = 168: KJ = -srtpinv, I = 158
                            s(168) = s(168) - f(174)*q( 158);
                            %q(174) = q(174) - q(158)*VYfac - q(168)*srtpinv;
                            
                            % from line 175: VZeq    = 2.0*srtpinv*VZfac; % line 175
                            %//srtpinv J = 158: KJ = -VZfac, I = 169
                            s( 158) = s( 158) - f(175)*q(169);
                            %//VZfac J = 169: KJ = -srtpinv, I = 158
                            s(169) = s(169) - f(175)*q( 158);
                            %q(175) = q(175) + 2.0*q(158)*VZfac + 2.0*q(169)*srtpinv;
                        end
                        
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            S = s;
        end
                
        function dJdt = exactwJ2_odefun(~,Inhomogeneous, Perturbation, y)
            
            % % mu = 1.0; %gravitational constant
            % % J2 = 1082.6269e-6;   %J2 perturbation constant
            % % Re = 1.0; %Radius of Earth (m)
            % %
            % % R = sqrt(y(1)^2 + y(2)^2 + y(3)^2);
            % % c = (3*J2*mu*Re^2)/(2*R^5);
            % Kvec = [y(1) y(2) y(3)];
            % [J6 Jperturb Jhom] = J6Gravity(Kvec);
            % %J6
            %
            % % J2K=     [(-mu*y(1)/(R^(3)) + c*(5*(y(3)^2/R^2) - 1)*y(1)),...
            % %        (-mu*y(2)/(R^(3)) + c*(5*(y(3)^2/R^2) - 1)*y(2)),...
            % %        (-mu*y(3)/(R^(3)) + c*(5*(y(3)^2/R^2) - 3)*y(3))];
            % % J2K
            % %
            %
            % % dJdt = [y(4);
            % %         y(5);
            % %         y(6);
            % %         -mu*y(1)/(R^(3)) + c*(5*(y(3)^2/R^2) - 1)*y(1);
            % %         -mu*y(2)/(R^(3)) + c*(5*(y(3)^2/R^2) - 1)*y(2);
            % %         -mu*y(3)/(R^(3)) + c*(5*(y(3)^2/R^2) - 3)*y(3)];
            % dJdt = [ y(4);
            %     y(5);
            %     y(6);
            %     J6(1);
            %     J6(2);
            %     J6(3)];
            
            
            dJdt = Inhomogenous + Perturbation*y;
            
            
        end
        
        function [H_I, Jacob, Hess] = WorkOrder(this, I)
            
            H_I          = this.DataList(I);
            
            Nvar         = 7;
            LineCheck    = I;
            F            = zeros(175,1);
            F(LineCheck) = 1;
            F = backdiff(this, F, LineCheck);
            Jacob        = [F(1) F(2) F(3) F(4) F(5) F(6) F(7)];
            % e,a,...         t    e    a    I    w    W    Mp
            Hess   = [];
            for j = 1: Nvar
                Q = zeros(175,1);
                S = zeros(175,1);
                %Q(ipoint(j)) = 1.0;
                Q(j) = 1.0;
                Q = fordiff(this, Q, LineCheck);
                S = secdiff(this, F, Q, S, LineCheck);
                S = backdiff(this, S, LineCheck);
                row = [];
                for k = 1 : Nvar
                    row = [row, S(k)];
                end
                Hess  = [Hess; row];
            end
            
        end
        
        function [H_I, Jacob, Hess] = GravityDerivs(this, I)
            
            H_I          = this.DataList(I);
            
            Nvar         = 3;
            LineCheck    = I;
            F            = zeros(175,1);
            F(LineCheck) = 1;
            F = this.backdiff(F, LineCheck);
            Jacob        = [F(61) F(62) F(63)];
            % e,a,...         X      Y     Z
            Hess   = [];
            for j = 1: Nvar
                Q = zeros(175,1);
                S = zeros(175,1);
                %Q(ipoint(j)) = 1.0;
                Q(j) = 1.0;
                Q = this.fordiff(Q, LineCheck);
                S = this.secdiff(F, Q, S, LineCheck);
                S = this.backdiff(S, LineCheck);
                row = [];
                for k = 61 : 63
                    row = [row, S(k)];
                end
                Hess  = [Hess; row];
            end
            
        end
        
        function [Jacob, JacobFinite, Hess, HessFinite] = DigitalJustice(this, I, Sensor, Observed)
            
            InvCovUV      = eye(2)*10^(8);
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
            [~,~,~,~,~,~] = Kdelta.OrbitDerivatives(time+epsilon, Sensor, Observed, InvCovUV);
            [HT, JacobT, HessT] = Kdelta.WorkOrder(I);
            
            Kepler1     = Kepler + KeplerDelta(1,:);
            Kdelta1     = LagrangePlanetary(Kepler1);
            [~,~,~,~,~,~] = Kdelta1.OrbitDerivatives(time, Sensor, Observed, InvCovUV);
            [H1, Jacob1, Hess1] = Kdelta1.WorkOrder(I);
            
            Kepler2     = Kepler + KeplerDelta(2,:);
            Kdelta2     = LagrangePlanetary(Kepler2);
            [~,~,~,~,~,~] = Kdelta2.OrbitDerivatives(time, Sensor, Observed, InvCovUV);
            [H2, Jacob2, Hess2] = Kdelta2.WorkOrder(I);
            
            Kepler3     = Kepler + KeplerDelta(3,:);
            Kdelta3     = LagrangePlanetary(Kepler3);
            [~,~,~,~,~,~] = Kdelta3.OrbitDerivatives(time, Sensor, Observed, InvCovUV);
            [H3, Jacob3, Hess3] = Kdelta3.WorkOrder(I);
            
            Kepler4     = Kepler + KeplerDelta(4,:);
            Kdelta4     = LagrangePlanetary(Kepler4);
            [~,~,~,~,~,~] = Kdelta4.OrbitDerivatives(time, Sensor, Observed, InvCovUV);
            [H4, Jacob4, Hess4] = Kdelta4.WorkOrder(I);
            
            Kepler5     = Kepler + KeplerDelta(5,:);
            Kdelta5     = LagrangePlanetary(Kepler5);
            [~,~,~,~,~,~] = Kdelta5.OrbitDerivatives(time, Sensor, Observed, InvCovUV);
            [H5, Jacob5, Hess5] = Kdelta5.WorkOrder(I);
            
            Kepler6     = Kepler + KeplerDelta(6,:);
            Kdelta6     = LagrangePlanetary(Kepler6);
            [~,~,~,~,~,~] = Kdelta6.OrbitDerivatives(time, Sensor, Observed, InvCovUV);
            [H6, Jacob6, Hess6] = Kdelta6.WorkOrder(I);
            
            JacobFinite = [HT-H0,H1-H0, H2-H0, H3-H0, H4-H0, H5-H0, H6-H0]/epsilon;
            
            HessFinite = [(JacobT-Jacob);...
                (Jacob1-Jacob);...
                (Jacob2-Jacob);...
                (Jacob3-Jacob);...
                (Jacob4-Jacob);...
                (Jacob5-Jacob);...
                (Jacob6-Jacob)]/epsilon;
        end
        
        function this = InitialValues()
        
            % this.time        = 0.0;
            % this.t           = 0.0;                             % line   1
            % this.e           = 0.0;                   % line   2
            % %this.p          = 0.0;                  % line   3
            % this.a           = 0.0;                   % line   3
            % this.Inclination = 0.0;                   % line   4
            % this.omega       = 0.0;                   % line   5
            % this.Omega       = 0.0;                   % line   6
            % this.Mp          = 0.0;                   % line   7
            % this.ClassicalElements = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0];
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Computations based on the above input variables begins ---
            % % start with line 11 to leave room for more input parameters in
            % % lines 1-10:
            %
            % this.onePe   = 0.0;                           % line  11
            % this.oneMe   = 0.0;                           % line  12
            % this.fac     = 0.0;                  % line  13
            % this.rootfac = 0.0;                   % line  14
            % %this.a          = 0.0;                 % line  15
            % this.p           = 0.0;                  % line  15
            % this.meanMotion  = 0.0;      % line  16
            %
            % this.Px   =  0.0;           % line  23
            % this.Py   =  0.0;           % line  24
            % this.Pz   =  0.0;                             % line  25
            %
            % this.Qx   = 0.0;           % line  26
            % this.Qy   = 0.0;           % line  27
            % this.Qz   = 0.0;                             % line  28
            %
            % % For Pure "Kepler" orbit we don't need Wvec, but it might play
            % % a role later when we include Perturbations, so leave it in
            % this.Wx   = 0.0;                              % line  29
            % this.Wy   = 0.0;                              % line  30
            % this.Wz   = 0.0;                                   % line  31
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Class Variable Storage Assignments -- not "computations"
            % this.cosI   = 0.0;
            % this.sinI   = 0.0;
            % this.cosom  = 0.0;
            % this.sinom  = 0.0;
            % this.cosO   = 0.0;
            % this.sinO   = 0.0;
            % this.Pvec   = [0.0; 0.0; 0.0];
            % this.Qvec   = [0.0; 0.0; 0.0];
            % this.Wvec   = [0.0; 0.0; 0.0];
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % this.M        = 0.;                                  % line  40
            % this.E        = 0.;                                  % line  41
            % this.cosK     = 0.;                                  % line  42
            % this.sinK     = 0.;                                  % line  43
            % %%%%%%%%%%%%%%%%%%% Implicit Function Theorem Derivatives
            % this.tanX     = 0.;                                  % line  44
            % this.tanY     = 0.;                                  % line  45
            % this.nu       = 0.;                                  % line  46
            % this.coss     = 0.;                                  % line  47
            % this.sins     = 0.;                                  % line  48
            % % OLD: rorbit   = p/(1.0 + e*coss);                  % line  49
            % this.rorbit   = 0.;                                  % line  49
            % this.rorbitx  = 0.;                                  % line  50
            % this.rorbity  = 0.;                                  % line  51
            % this.rorbitz  = 0.;                                  % line  52
            % this.rtpinv   = 0.;                                  % line  53
            % this.vorbitx  = 0.;                                  % line  54
            % this.vorbity  = 0.;                                  % line  55
            % this.vorbitz  = 0.;                                  % line  56
            % %  rmag
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Escobal Gravity Terms
            % this.X    = 0.;                                      % line  61
            % this.Y    = 0.;                                      % line  62
            % this.Z    = 0.;                                      % line  63
            % this.R2   = 0.;                                      % line  64
            % this.R1   = 0.;                                      % line  65
            % this.R3   = 0.;                                      % line  66
            % this.R4   = 0.;                                      % line  67
            % this.R5   = 0.;                                      % line  68
            % this.R6   = 0.;                                      % line  69
            % this.sd1  = 0.;                                      % line  70
            % this.sd2  = 0.;                                      % line  71
            % this.sd3  = 0.;                                      % line  72
            % this.sd4  = 0.;                                      % line  73
            % this.sd5  = 0.;                                      % line  74
            % this.sd6  = 0.;                                      % line  75
            % this.F2   = 0.;                                      % line  76
            % this.F3   = 0.;                                      % line  77
            % this.F4   = 0.;                                      % line  78
            % this.F5   = 0.;                                      % line  79
            % this.F6   = 0.;                                      % line  80
            % this.V1   = 0.;                                      % line  81
            % this.V2   = 0.;                                      % line  82
            % this.V3   = 0.;                                      % line  83
            % this.V4   = 0.;                                      % line  84
            % this.V5   = 0.;                                      % line  85
            % this.V6   = 0.;                                      % line  86
            % this.V    = 0.;           % Kepler Potential         % line  87
            % this.D    = 0.;                                      % line  88
            % this.Phi  = 0.;                                      % line  90
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % this.fac1 = 10^(-10);                                % line  91
            % this.fac2 = 10^(-10);                                % line  92
            % this.fac3 = 10^(-10);                                % line  93
            % this.fac4 = 10^(-10);                                % line  94
            % this.fac5 = 10^(-10);                                % line  95
            % this.Mwe  = 0.;                                      % line  96
            % this.MeM  = 0.;                                      % line  97
            % this.MaM  = 0.;                                      % line  98
            % this.MWi  = 0.;                                      % line  99
            % this.Miw  = 0.;                                      % line 100
            % this.MLagrange = zeros(6);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            this.los      = [0.;0.;0.];
            this.lx       = 0.;                                  % line 101
            this.ly       = 0.;                                  % line 102
            this.lz       = 0.;                                  % line 103
            this.lperpSq  = 10^(-20);                            % line 104
            this.rangeSq  = 0.;                                  % line 105
            this.range    = 10^(-10);                            % line 106
            this.zLos     = 0.;                                  % line 107
            this.thetaLos = 0.;                                  % line 108
            this.phiLos   = 0.;                                  % line 109
            this.vlos     = [0;0;0];
            this.vx       = 0.;                                  % line 110
            this.vy       = 0.;                                  % line 111
            this.vz       = 0.;                                  % line 112
            this.ldotv    = 0.;                                  % line 113
            this.ldot     = 0.;                                  % line 114
            this.lperp    = 10^(-10);                            % line 115
            this.thdotfac1= 0.;                                  % line 116
            this.thdotfac2= 0.;                                  % line 117
            this.thdot    = 0.;                                  % line 118
            this.phiDotNum= 0.;                                  % line 119
            this.phidot   = 0.;                                  % line 120
            this.thetaM      = 0.;
            this.cosThetaM   = 0.;
            this.sinThetaM   = 0.;
            this.phiM        = 0.;
            this.cosPhiM     = 0.;
            this.sinPhiM     = 0.;
            
            this.cosThetaLos = 0.;                               % line 121
            this.sinThetaLos = 0.;                               % line 122
            this.phiDiff     = 0.;                               % line 123
            this.cosPhiDiff  = 0.;                               % line 124
            this.sinPhiDiff  = 0.;                               % line 125
            this.delUV       = 0.;                               % line 126
            this.Uview       = 0.;                               % line 127
            this.Vview       = 0.;                               % line 128
            this.xLos        = 0.;                               % line 129
            this.yLos        = 0.;                               % line 130
            this.xLosM       = 0.;
            this.yLosM       = 0.;
            this.zLosM       = 0.;
            this.Collinearity = 0.;                              % line 131
            this.Mdelta      = zeros(2,3);
            this.InvCovM     = zeros(3,3);
            this.losDiff     = zeros(3,1);                       % line 132
            this.xUnit       = [1.0; 0.0; 0.0];
            this.yUnit       = [0.0; 1.0; 0.0];
            this.zUnit       = [0.0; 0.0; 1.0];
            this.ChiSqLos    = 0.0;                              % line 133
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            this.Aeq             = 0;                  % line 142 Aeq
            this.cosAeq          = 0;                  % line 143 cosAeq
            this.sinAeq          = 0;                  % line 144 sinAeq
            this.tanHalf         = 0;                  % line 145 tanHalf
            this.Feq             = 0;                  % line 146 Feq
            this.Geq             = 0;                  % line 147 Geq
            this.Heq             = 0;                  % line 148 Heq
            this.Keq             = 0;                  % line 149 Keq
            this.Leq             = 0;                  % line 150 Leq
            this.CosL            = 0;                  % line 151 CosL
            this.SinL            = 0;                  % line 152 SinL
            this.alphaSq         = 0;                  % line 153 alphaSq
            this.Seq             = 0;                  % line 154
            this.Weq             = 0;                  % line 155
            this.Req             = 0;                  % line 156
            this.RovS            = 0;                  % line 157
            this.srtpinv         = 0;                  % line 158
            this.HK              = 0;                  % line 159
            this.OnePalphaSqCosL = 0;                  % line 160
            this.OneMalphaSqCosL = 0;                  % line 161
            this.OnePalphaSqSinL = 0;                  % line 162
            this.OneMalphaSqSinL = 0;                  % line 163
            this.Xfac            = 0;                  % line 164
            this.Yfac            = 0;                  % line 165
            this.Zfac            = 0;                  % line 166
            this.VXfac           = 0;                  % line 167
            this.VYfac           = 0;                  % line 168
            this.VZfac           = 0;                  % line 169
            this.Xeq             = 0;                  % line 170
            this.Yeq             = 0;                  % line 171
            this.Zeq             = 0;                  % line 172
            this.VXeq            = 0;                  % line 173
            this.VYeq            = 0;                  % line 174
            this.VZeq            = 0;                  % line 175
        end
        
        function this = SetDataSet(this)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            this.DataList(  1) = this.time;                      % line   1
            this.DataList(  2) = this.e;                         % line   2
            this.DataList(  3) = this.a;                         % line   3
            this.DataList(  4) = this.Inclination;               % line   4
            this.DataList(  5) = this.omega;                     % line   5
            this.DataList(  6) = this.Omega;                     % line   6
            this.DataList(  7) = this.Mp;                        % line   7
            this.DataList( 11) = this.onePe;                     % line  11
            this.DataList( 12) = this.oneMe;                     % line  12
            this.DataList( 13) = this.fac;                       % line  13
            if this.fac > 0.0
                this.DataList( 14) = this.rootfac;               % line  14
            else
                this.DataList( 14) = 0;                          % line  14
            end
            this.DataList( 15) = this.p;                         % line  15
            this.DataList( 16) = this.meanMotion;                % line  16
            this.DataList( 17) = this.cosI;                      % line  17
            this.DataList( 18) = this.sinI;                      % line  18
            this.DataList( 19) = this.cosom;                     % line  19
            this.DataList( 20) = this.sinom;                     % line  20
            this.DataList( 21) = this.cosO;                      % line  21
            this.DataList( 22) = this.sinO;                      % line  22
            
            this.DataList( 23) = this.Px;                        % line  23
            this.DataList( 24) = this.Py;                        % line  24
            this.DataList( 25) = this.Pz;                        % line  25
            
            this.DataList( 26) = this.Qx;                        % line  26
            this.DataList( 27) = this.Qy;                        % line  27
            this.DataList( 28) = this.Qz;                        % line  28
            
            this.DataList( 29) = this.Wx;                        % line  29
            this.DataList( 30) = this.Wy;                        % line  30
            this.DataList( 31) = this.Wz;                        % line  31
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Class Variable Storage Assignments -- not "computations"
            %this.DataList() = this.Pvec;
            %this.DataList() = this.Qvec;
            %this.DataList() = this.Wvec;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % variables for use in Derivatives based on Implicit Function Theorem
            % on M = E - e*sin(E)
            %this.DataList() = this.eDenom;
            %this.DataList() = this.Eprime;
            %this.DataList() = this.dE_dM;
            %this.DataList() = this.dE_de;
            %this.DataList() = this.d2E_dMdM;
            %this.DataList() = this.d2E_dMde;
            %this.DataList() = this.d2E_dedM;
            %this.DataList() = this.d2E_dede;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            this.DataList( 40) = this.M;                         % line  40
            this.DataList( 41) = this.E;                         % line  41
            this.DataList( 42) = this.cosK;                      % line  42
            this.DataList( 43) = this.sinK;                      % line  43
            %%%%%%%%%%%%%%%%%%% Implicit Function Theorem Derivatives
            this.DataList( 44) = this.tanX;                      % line  44
            this.DataList( 45) = this.tanY;                      % line  45
            this.DataList( 46) = this.nu;                        % line  46
            this.DataList( 47) = this.coss;                      % line  47
            this.DataList( 49) = this.sins;                      % line  48
            % OLD: rorbit   = p/(1.0 + e*coss);                  % line  49
            this.DataList( 49) = this.rorbit;                    % line  49
            this.DataList( 50) = this.rorbitx;                   % line  50
            this.DataList( 51) = this.rorbity;                   % line  51
            this.DataList( 52) = this.rorbitz;                   % line  52
            this.DataList( 53) = this.rtpinv;                    % line  53
            this.DataList( 54) = this.vorbitx;                   % line  54
            this.DataList( 55) = this.vorbity;                   % line  55
            this.DataList( 56) = this.vorbitz;                   % line  56
            %  rmag
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Escobal Gravity Terms
            this.DataList( 61) = this.X;                         % line  61
            this.DataList( 62) = this.Y;                         % line  62
            this.DataList( 63) = this.Z;                         % line  63
            this.DataList( 64) = this.R2;                        % line  64
            this.DataList( 65) = this.R1;                        % line  65
            this.DataList( 66) = this.R3;                        % line  66
            this.DataList( 67) = this.R4;                        % line  67
            this.DataList( 68) = this.R5;                        % line  68
            this.DataList( 69) = this.R6;                        % line  69
            this.DataList( 70) = this.sd1;                       % line  70
            this.DataList( 71) = this.sd2;                       % line  71
            this.DataList( 72) = this.sd3;                       % line  72
            this.DataList( 73) = this.sd4;                       % line  73
            this.DataList( 74) = this.sd5;                       % line  74
            this.DataList( 75) = this.sd6;                       % line  75
            this.DataList( 76) = this.F2;                        % line  76
            this.DataList( 77) = this.F3;                        % line  77
            this.DataList( 78) = this.F4;                        % line  78
            this.DataList( 79) = this.F5;                        % line  79
            this.DataList( 80) = this.F6;                        % line  80
            this.DataList( 81) = this.V1;                        % line  81
            this.DataList( 82) = this.V2;                        % line  82
            this.DataList( 83) = this.V3;                        % line  83
            this.DataList( 84) = this.V4;                        % line  84
            this.DataList( 85) = this.V5;                        % line  85
            this.DataList( 86) = this.V6;                        % line  86
            this.DataList( 87) = this.V;   % Kepler Potential    % line  87
            this.DataList( 88) = this.D;                         % line  88
            
            this.DataList( 90) = this.Phi;                       % line  90
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            this.DataList( 91) = this.fac1;                      % line  91
            this.DataList( 92) = this.fac2;                      % line  92
            this.DataList( 93) = this.fac3;                      % line  93
            this.DataList( 94) = this.fac4;                      % line  94
            this.DataList( 95) = this.fac5;                      % line  95
            this.DataList( 96) = this.Mwe;                       % line  96
            this.DataList( 97) = this.MeM;                       % line  97
            this.DataList( 98) = this.MaM;                       % line  98
            this.DataList( 99) = this.MWi;                       % line  99
            this.DataList(100) = this.Miw;                       % line 100
            %this.MLagrange = zeros(6);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % this.DataList(101) = this.lx;                        % line 101
            % this.DataList(102) = this.ly;                        % line 102
            % this.DataList(103) = this.lz;                        % line 103
            % this.DataList(104) = this.lperpSq;                   % line 104
            % this.DataList(105) = this.rangeSq;                   % line 105
            % this.DataList(106) = this.range;                     % line 106
            % this.DataList(107) = this.zLos;                      % line 107
            % this.DataList(108) = this.thetaLos;                  % line 108
            % this.DataList(109) = this.phiLos;                    % line 109
            % this.DataList(110) = this.vx;                        % line 110
            % this.DataList(111) = this.vy;                        % line 111
            % this.DataList(112) = this.vz;                        % line 112
            % this.DataList(113) = this.ldotv;                     % line 113
            % this.DataList(114) = this.ldot;                      % line 114
            % this.DataList(115) = this.lperp;                     % line 115
            % this.DataList(116) = this.thdotfac1;                 % line 116
            % this.DataList(117) = this.thdotfac2;                 % line 117
            % this.DataList(118) = this.thdot;                     % line 118
            % this.DataList(119) = this.phiDotNum;                 % line 119
            % this.DataList(120) = this.phidot;                    % line 120
            % this.DataList(121) = this.cosThetaLos;               % line 121
            % this.DataList(122) = this.sinThetaLos;               % line 122
            % this.DataList(123) = this.phiDiff;                   % line 123
            % this.DataList(124) = this.cosPhiDiff;                % line 124
            % this.DataList(125) = this.sinPhiDiff;                % line 125
            % this.DataList(126) = this.delUV;                     % line 126
            % this.DataList(127) = this.Uview;                     % line 127
            % this.DataList(128) = this.Vview;                     % line 128
            % this.DataList(129) = this.xLos;                      % line 129
            % this.DataList(130) = this.yLos;                      % line 130
            % this.DataList(131) = this.Collinearity;              % line 131
            % %
            % % Mdelta      =[cosThetaM*cosPhiM, cosThetaM*sinPhiM, -sinThetaM;...
            % %                        -sinPhiM,           cosPhiM,       0.0 ];
            % % InvCovM  = Mdelta'*Mdelta;
            % % xUnit    = [1.0; 0.0; 0.0];
            % % yUnit    = [0.0; 1.0; 0.0];
            % % zUnit    = [0.0; 0.0; 1.0];
            % %
            % % losDiff  = [xLos - xLosM; yLos - yLosM; zLos - zLosM];% line 132
            % this.DataList(133) = this.ChiSqLos;                  % line 133
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % this.DataList(142)   = this.Aeq;                      % line 142 Aeq
            % this.DataList(143)   = this.cosAeq;                   % line 143 cosAeq
            % this.DataList(144)   = this.sinAeq;                   % line 144 sinAeq
            % this.DataList(145)   = this.tanHalf;                  % line 145 tanHalf
            % this.DataList(146)   = this.Feq;                      % line 146 Feq
            % this.DataList(147)   = this.Geq;                      % line 147 Geq
            % this.DataList(148)   = this.Heq;                      % line 148 Heq
            % this.DataList(149)   = this.Keq;                      % line 149 Keq
            % this.DataList(150)   = this.Leq;                      % line 150 Leq
            % this.DataList(151)   = this.CosL;                     % line 151 CosL
            % this.DataList(152)   = this.SinL;                     % line 152 SinL
            % this.DataList(153)   = this.alphaSq;                  % line 153 alphaSq
            % this.DataList(154)   = this.Seq;                      % line 154
            % this.DataList(155)   = this.Weq;                      % line 155
            % this.DataList(156)   = this.Req;                      % line 156
            % this.DataList(157)   = this.RovS;                     % line 157
            % this.DataList(158)   = this.srtpinv;                  % line 158
            % this.DataList(159)   = this.HK;                       % line 159
            % this.DataList(160)   = this.OnePalphaSqCosL;          % line 160
            % this.DataList(161)   = this.OneMalphaSqCosL;          % line 161
            % this.DataList(161)   = this.OnePalphaSqSinL;          % line 162
            % this.DataList(163)   = this.OneMalphaSqSinL;          % line 163
            % this.DataList(163)   = this.Xfac;                     % line 164
            % this.DataList(165)   = this.Yfac;                     % line 165
            % this.DataList(166)   = this.Zfac;                     % line 166
            % this.DataList(167)   = this.VXfac;                    % line 167
            % this.DataList(168)   = this.VYfac;                    % line 168
            % this.DataList(169)   = this.VZfac;                    % line 169
            % this.DataList(170)   = this.Xeq;                      % line 170
            % this.DataList(171)   = this.Yeq;                      % line 171
            % this.DataList(172)   = this.Zeq;                      % line 172
            % this.DataList(173)   = this.VXeq;                     % line 173
            % this.DataList(174)   = this.VYeq;                     % line 174
            % this.DataList(175)   = this.VZeq;                     % line 175
        end