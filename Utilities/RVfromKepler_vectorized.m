function [SV7Extrapolated] = RVfromKepler_vectorized(time, Kepler, units)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Canonical Coordinates assumed throughout -- length in DU, time in TU
    % Modifications can be made to the choice of units in wgs84Constants!
    % For Choice of units
    TU = units.TU;
    DU = units.DU;
    VU = units.VU;
    AU = units.AU;
    mu = units.mu;
    sqrtmu = sqrt(mu);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    e           = Kepler(:,1);                   % line   2
    a           = Kepler(:,2);                   % line   3
    Inclination = Kepler(:,3);                   % line   4
    omega       = Kepler(:,4);                   % line   5
    Omega       = Kepler(:,5);                   % line   6
    Mp          = Kepler(:,6);                   % line   7
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Computations based on the above input variables begins ---
    % start with line 11 to leave room for more input parameters in
    % lines 1-10:
    
    % Semi-Major Axis:
    onePe = 1.0 + e;                    % 39) 
    oneMe = 1.0 - e;                    % 40)
    fac   = onePe.*oneMe;               % 41)
    rootfac = zeros(size(fac,1),1);
    Inds = fac>0.0;
    rootfac(Inds) = sqrt(fac(Inds));            % 42)
    p           = a.*fac;                  % line  15
    meanMotion  = sqrtmu*a.^(-1.5);      % line  16
    
    cosI             = cos(Inclination);            % line  17
    sinI             = sin(Inclination);            % line  18
    cosom            = cos(omega);                  % line  19
    sinom            = sin(omega);                  % line  20
    cosO             = cos(Omega);                  % line  21
    sinO             = sin(Omega);                  % line  22
    
    Px   =  cosO.*cosom - sinO.*sinom.*cosI;           % line  23
    Py   =  sinO.*cosom + cosO.*sinom.*cosI;           % line  24
    Pz   =  sinom.*sinI;                             % line  25
    
    Qx   = -cosO.*sinom - sinO.*cosom.*cosI;           % line  26
    Qy   = -sinO.*sinom + cosO.*cosom.*cosI;           % line  27
    Qz   = cosom.*sinI;                             % line  28
    
    % For Pure "Kepler" orbit we don't need Wvec, but it might play
    % a role later when we include Perturbations, so leave it in
    Wx   =  sinO.*sinI;                              % line  29
    Wy   = -cosO.*sinI;                              % line  30
    Wz   =  cosI;                                   % line  31
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % variables for use in Derivatives based on Implicit Function Theorem
    % on M = E - e*sin(E)
    % invert Kepler's Equation
    % Brute Force iteration (good way to seed Newton's method which follows)
    % M = MeanAnomalyEpoch;
    %  meanAnomaly M to eccentric Anomaly E
    M = meanMotion.*time + Mp;
    %  meanAnomaly M to eccentric Anomaly EM
    EM = M + e;                                        
    %if (M > pi)
    %    EM = M - e;
    %elseif (-pi < M && M < 0)
    %    EM = M - e;
    %end
    Inds = M>pi;
    EM(Inds) = M(Inds) - e(Inds);
    Inds = (-pi<M) & (M<0);
    EM(Inds) = M(Inds) - e(Inds);
    %         //cout << " Input meanAnomaly = " << meanAnomaly << endl;
    for i=1:1:10
        EM  = M + e.*sin(EM);
        %//cout << " Mean Anomaly Solution " << E << endl;
    end
    % //      10 rounds of Newton's root finding method based on the above "seed". 
    for i=1:1:10
        Eprime      = 1.0 - e.*cos(EM);
        EM          = EM + (M - EM + e.*sin(EM))./Eprime;
    end
    % Solve Kepler's Equation for E:  M = EM - e*sin(EM);   % 66) EM = E(e,M)
    % Need to differentiate Implicitly here!
    eDenom   = Eprime.*Eprime.*Eprime;
    
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
    dE_dM    = 1.0./Eprime;
    dE_de    = sinK./Eprime;
    d2E_dMdM = -e.*sinK./eDenom;
    %// Mike Cain Corrections! 
    d2E_dMde = (cosK - e)./eDenom;
    d2E_dedM = (cosK - e)./eDenom;
    d2E_dede = ((2.0 - e.*cosK).*cosK.*sinK - e.*sinK)./eDenom;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % kfac  = 1.0 - e*cosK;                                 % 69)
    % coss  = (cosK - e)/kfac;                              % 70)                               
    % sins  = rootfac*sinK/kfac;                            % 71)
    % 	s     = atan2(sins, coss);                          % 72 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   % for convenience -- you can eliminate kfac; 
    tanX  = cosK - e;                                       % 70)                               
    tanY  = rootfac.*sinK;                                   % 71)
    nu     = atan2(tanY, tanX);                             % 72)
    coss  = cos(nu);                                        % 73)
    sins  = sin(nu);                                        % 74) 
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
    rorbit   = p./(1.0 + e.*coss);                        % line 75
    rorbitx  = rorbit.*(coss.*Px + sins.*Qx);              % line 76
    rorbity  = rorbit.*(coss.*Py + sins.*Qy);              % line 77
    rorbitz  = rorbit.*(coss.*Pz + sins.*Qz);              % line 78
    rtpinv   = sqrt(mu./p);
    vorbitx  = rtpinv.*(-sins.*Px + (e + coss).*Qx);  % line 79
    vorbity  = rtpinv.*(-sins.*Py + (e + coss).*Qy);  % line 80
    vorbitz  = rtpinv.*(-sins.*Pz + (e + coss).*Qz);  % line 81
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SV7Extrapolated = [rorbitx, rorbity, rorbitz, vorbitx, vorbity, vorbitz]';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end