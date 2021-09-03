function [ XEquin, XECI, JacKepler_2_Equinoctial] =  Kepler_2_All(XKepler, units)
    mu     = units.mu;
    sqrtmu = units.sqrtmu

    e           = XKepler(1);
    a           = XKepler(2);
    Inclination = XKepler(3);
    omega       = XKepler(4);
    Omega       = XKepler(5);
    %nu          = XKepler(6);
    Mp          = XKepler(6);
    tEpoch      = XKepler(7);
    
    onePe = 1.0 + e;                           % line  11
    oneMe = 1.0 - e;                           % line  12
    fac   = onePe*oneMe;                  % line  13
    if fac > 0.0
        rootfac = sqrt(fac);                   % line  14
    else
        rootfac = 0;                                % line  14
    end
    %a           = p/fac;                 % line  15
    p           = a*fac;                  % line  15
    n  = sqrtmu*(a)^(-1.5);      % line  16

    cosI             = cos(Inclination);            % line  17
    sinI             = sin(Inclination);            % line  18
    cosom            = cos(omega);                  % line  19
    sinom            = sin(omega);                  % line  20
    cosO             = cos(Omega);                  % line  21
    sinO             = sin(Omega);                  % line  22

    Px   =  cosO*cosom - sinO*sinom*cosI;           % line  23
    Py   =  sinO*cosom + cosO*sinom*cosI;           % line  24
    Pz   =  sinom*sinI;                             % line  25

    Qx   = -cosO*sinom - sinO*cosom*cosI;           % line  26
    Qy   = -sinO*sinom + cosO*cosom*cosI;           % line  27
    Qz   =  cosom*sinI;                             % line  28

    % For Pure "Kepler" orbit we don't need Wvec, but it might play
    % a role later when we include Perturbations, so leave it in
    Wx   =  sinO*sinI;                              % line  29
    Wy   = -cosO*sinI;                              % line  30
    Wz   =  cosI;                                   % line  31
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % invert Kepler's Equation
    % Brute Force iteration (good way to seed Newton's method which follows)
    % M = MeanAnomalyEpoch;
    %  meanAnomaly M to eccentric Anomaly E
    M = n*tEpoch + Mp;
    E = M + e;
    if (M > pi)
        E = M - e;
    elseif (-pi < M)
        if (M < 0)
            E = M - e;
        end
    end
    %         //cout << " Input meanAnomaly = " << meanAnomaly << endl;
    for k=1:1:10
        E  = M + e*sin(E);
        %//cout << " Mean Anomaly Solution " << E << endl;
    end
    % //      10 rounds of Newton's root finding method based on the above "seed".
    for k=1:1:10
        Eprime      = 1.0 - e*cos(E);
        E           = E + (M - E + e*sin(E))./Eprime;
    end
    KeplerInv = E;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cosK  = cos(E);                                   % line 42
    sinK  = sin(E);                                   % line 43
    %%%%%%%%%%%%%%%%%%% Implicit Function Theorem Derivatives
    eDenom   = Eprime*Eprime*Eprime;
    dE_dM    = 1.0/Eprime;
    dE_de    = sinK/Eprime;
    d2E_dMdM = -e*sinK/eDenom;
    %// Mike Cain Corrections!
    d2E_dMde = (cosK - e)/eDenom;
    d2E_dedM = (cosK - e)/eDenom;
    d2E_dede = ((2.0 - e*cosK)*cosK*sinK - e*sinK)/eDenom;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   % for convenience -- you can eliminate kfac;
    tanX  = cosK - e;                                 % line 44
    tanY  = rootfac*sinK;                             % line 45
    nu     = atan2(tanY, tanX);                       % line 46
    coss  = cos(nu);                                  % line 47
    sins  = sin(nu);                                  % line 48
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Let us propagate this Orbit to time = t
    %rorbit   = p/(1.0 + e*coss);                        % line  49
    rorbit   = a*(1.0 - e*cosK);                    % line  49
    rorbitx  = rorbit*(coss*Px + sins*Qx);               % line  50
    rorbity  = rorbit*(coss*Py + sins*Qy);               % line  51
    rorbitz  = rorbit*(coss*Pz + sins*Qz);               % line  52
    rtpinv   = sqrt(mu/p);                               % line  53
    vorbitx  = rtpinv*(-sins*Px + (e + coss)*Qx);        % line  54
    vorbity  = rtpinv*(-sins*Py + (e + coss)*Qy);        % line  55
    vorbitz  = rtpinv*(-sins*Pz + (e + coss)*Qz);        % line  56
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    XECI     = [rorbitx; rorbity; rorbitz; vorbitx; vorbity; vorbitz; tEpoch];

    tanHalf     = tan(Inclination/2.0);
    p           = a*(1-e^2);
    f           = e*cos(omega + Omega);
    g           = e*sin(omega + Omega);
    h           = tanHalf*cos(Omega);
    k           = tanHalf*sin(Omega);
    L           = omega + Omega + nu;

    XEquin = [p; f; g; h; k; L; tEpoch];
    
    JacKepler_2_Equinoctial(1,1 ) =  -(2.0*a*e);
    JacKepler_2_Equinoctial(1,2 ) =  -(e^2-1.0);
    JacKepler_2_Equinoctial(1,3 ) =  0.0;
    JacKepler_2_Equinoctial(1,4 ) =  0.0;
    JacKepler_2_Equinoctial(1,5 ) =  0.0;
    JacKepler_2_Equinoctial(1,6 ) =  0.0;

    JacKepler_2_Equinoctial(2,1 ) =  cos(Omega+omega);
    JacKepler_2_Equinoctial(2,2 ) =  0.0;
    JacKepler_2_Equinoctial(2,3 ) =  0.0;
    JacKepler_2_Equinoctial(2,4 ) =  -(sin(Omega+omega)*e);
    JacKepler_2_Equinoctial(2,5 ) =  -(sin(Omega+omega)*e);
    JacKepler_2_Equinoctial(2,6 ) =  0.0;

    JacKepler_2_Equinoctial(3,1 ) =  sin(Omega+omega);
    JacKepler_2_Equinoctial(3,2 ) =  0.0;
    JacKepler_2_Equinoctial(3,3 ) =  0.0;
    JacKepler_2_Equinoctial(3,4 ) =  cos(Omega+omega)*e;
    JacKepler_2_Equinoctial(3,5 ) =  cos(Omega+omega)*e;
    JacKepler_2_Equinoctial(3,6 ) =  0.0;

    JacKepler_2_Equinoctial(4,1 ) =  0.0;
    JacKepler_2_Equinoctial(4,2 ) =  0.0;
    JacKepler_2_Equinoctial(4,3 ) =  (tanHalf*tanHalf+1.0)*cos(Omega)/ 2.0;
    JacKepler_2_Equinoctial(4,4 ) =  0.0;
    JacKepler_2_Equinoctial(4,5 ) =  -(sin(Omega)*tanHalf);
    JacKepler_2_Equinoctial(4,6 ) =  0.0;

    JacKepler_2_Equinoctial(5,1 ) =  0.0;
    JacKepler_2_Equinoctial(5,2 ) =  0.0;
    JacKepler_2_Equinoctial(5,3 ) =  (tanHalf*tanHalf+1.0)*sin(Omega)/ 2.0;
    JacKepler_2_Equinoctial(5,4 ) =  0.0;
    JacKepler_2_Equinoctial(5,5 ) =  cos(Omega)*tanHalf;
    JacKepler_2_Equinoctial(5,6 ) =  0.0;

    JacKepler_2_Equinoctial(6,1 ) =  0.0;
    JacKepler_2_Equinoctial(6,2 ) =  0.0;
    JacKepler_2_Equinoctial(6,3 ) =  0.0;
    JacKepler_2_Equinoctial(6,4 ) =  1.0;
    JacKepler_2_Equinoctial(6,5 ) =  1.0;
    JacKepler_2_Equinoctial(6,6 ) =  1.0;
end