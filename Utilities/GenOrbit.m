function [Orbit, Kepler] = GenOrbit(tEpoch, units)

    %units = wgs84Constants;
    TU = units.TU;
    DU = units.DU;
    VU = units.VU;
    AU = units.AU;    
    mu = units.mu; 
    twopi = units.twopi;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     //  T0 = Epoch Time of the Chief = 0 = delt1 = 0.0;   // = (day1 - day1)*XMNPDA;
    GoldenEye1VelocityMin = 6500; % // m/s
    %     double GoldenEye1PropagationTimeAfterEpoch = 1.5; // min
    %     //double DaughterEpoch = 1.66666667;     // 100 seconds = 1.66667 minutes
    %     double DaughterEpoch = 0.5;     // 30 seconds = 1.66667 minutes
    %     double vsmear = 0.002;     // units = km/s 
    %     //double vsmear = 0.0;     // units = km/s 
    %     double rsmear = 0.1;       // units = km 
    %     //double vmag = 0.1;       // units = km/s
    %     //double vmag = 0.1;       // 100m/s
    %     double vmag = 0.01;        // 10 m/s

    %  // generate a radius vector
    %   // we will work with r between 100 to 400 km above the surface of the earth
    %   //double rmag = 1.0 + (100.0 + 300.0*r->Rndm())/XKMPER; // r in DU
    rmag =  1.0 + (100000.0 + 300000.0*rand)/DU; % // r in DU
    rcos = -1.0 + 2.0*rand;
    rsin = sqrt(1.0 - rcos*rcos);
    rphi = twopi*rand;
    rGoldenEye1 = rmag*[rsin*cos(rphi), rsin*sin(rphi), rcos];
    rVec  = rGoldenEye1;

    % rGoldenEye1(1)
    % rGoldenEye1(2)
    % rGoldenEye1(3)
    % 
    % // find vectors perpendicular to rGoldenEye1:
    %   double rx, ry, rz;
    if(abs(rGoldenEye1(1))>abs(rGoldenEye1(2)) && abs(rGoldenEye1(1))>abs(rGoldenEye1(3)))
        rx =  0.0;
        ry =  rGoldenEye1(3);
        rz = -rGoldenEye1(2);
    end
    if(abs(rGoldenEye1(2))>abs(rGoldenEye1(1)) && abs(rGoldenEye1(2))>abs(rGoldenEye1(3)))
        rx =  rGoldenEye1(3);
        ry =  0.0;
        rz = -rGoldenEye1(1);
    end
    if(abs(rGoldenEye1(3))>abs(rGoldenEye1(1)) && abs(rGoldenEye1(3))>abs(rGoldenEye1(2)))
        rx =  rGoldenEye1(2);
        ry = -rGoldenEye1(1);
        rz =  0.0;
    end

    rUnit   = rGoldenEye1;
    rPerp1  = [rx, ry, rz];
    rPerp2  = cross(rGoldenEye1, rPerp1);
    rUnit  = (1.0/ norm(rUnit))*rUnit;
    rPerp1 = (1.0/norm(rPerp1))*rPerp1;
    rPerp2 = (1.0/norm(rPerp2))*rPerp2;

    % check1 = dot(rUnit, rPerp1)
    % check2 = dot(rUnit, rPerp2)
    % check3 = dot(rPerp1,rPerp2)
    % check1 = dot(rUnit,  rUnit)
    % check2 = dot(rPerp1,rPerp1)
    % check3 = dot(rPerp2,rPerp2)

    % // generate the GoldenEye1's velocity vector in km/second [4 to 8 km/s] and convert to canonical units DU/TU
    vmagGoldenEye1 = GoldenEye1VelocityMin + 1000.0*rand; %  // vmagGoldenEye1 = 5.5 to 6.5 km/s
    % //          double vmagGoldenEye1 = 6.5 + 1.0*r->Rndm();   // vmagGoldenEye1 = 5.5 to 6.5 km/s
    % //        1 ER = XKMPER km, 1 time_unit = 806.81 seconds
    % //        conversion factor vfac: 1 = (1DU/XKMPER-km)*(time_unit-seconds/ 1 TU)
    vmagGoldenEye1 = vmagGoldenEye1/VU;   %// v in DU/TU
    %           //double vcos = -1.0 + 2.0*r->Rndm();
    %           //double vcos = 0.2 + 0.8*r->Rndm();  // produce v with components along the initial r-direction
    vcos = 0.0 + 1.0*rand; % // produce v with components along the initial r-direction
    vsin = sqrt(1.0 - vcos*vcos);
    vphi = twopi*rand;
    vGoldenEye1 = vmagGoldenEye1*(vsin*cos(vphi)*rPerp1 + vsin*sin(vphi)*rPerp2 + vcos*rUnit);
    vVec = vGoldenEye1;

    % Work on Kepler Components

    rmag = norm(rVec);
    vsq  = dot(vVec,vVec);
    er   = vsq - 1.0/rmag;
    ev   = dot(rVec, vVec);
    eVec = er*rVec - ev*vVec;
    e    = norm(eVec);
    PVec = eVec/e;
    hVec = cross(rVec, vVec);
    p    = dot(hVec,hVec);
    WVec = hVec/norm(hVec);
    QVec = cross(WVec, PVec);

    % Dot1 = dot(PVec, QVec)
    % Dot2 = dot(PVec, WVec)
    % Dot3 = dot(QVec, WVec)
    % Dot4 = dot(PVec, PVec)
    % Dot5 = dot(QVec, QVec)
    % Dot6 = dot(WVec, WVec)

    onePe = 1.0 + e;
    oneMe = 1.0 - e;
    fac   = onePe*oneMe;
    if(fac > 0)
        rootfac = sqrt(fac);
    else
        rootfac = 0.0;
    end
    a = p/fac;
    PeriodTU = twopi*a^(1.5);
    Period   = TU*PeriodTU;     % Period in Seconds

    meanMotion = twopi/Period;     % Radians per Second

    cosnu  = dot(rUnit, PVec);
    sinnu  = dot(WVec',cross(PVec, rUnit));
    %trig = cosnu*cosnu + sinnu*sinnu
    nuEpoch = atan2(sinnu, cosnu);
    while(nuEpoch < 0.0)
        nuEpoch = nuEpoch + twopi;
    end
    while(nuEpoch > twopi)
        nuEpoch = nuEpoch - twopi;
    end
    denom1 = 1.0 + e*cosnu;
    cosE = (e +   cosnu)/denom1;
    sinE = rootfac*sinnu/denom1;
    %trig = cosE*cosE + sinE*sinE
    EccentricAnomalyEpoch = atan2(sinE, cosE);

    % Find the Time at Periapsis
    % Solve Kepler's Equation for the MeanAnomaly at Epoch -- this side is easy

    MeanAnomalyEpoch              = EccentricAnomalyEpoch - e*sinE;
    TimeSincePeriapsis            = MeanAnomalyEpoch/meanMotion;
    FracionOfPeriodSincePeriapsis = TimeSincePeriapsis/Period;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Sanity Check on Eccentric Anomaly
    % % // 1)   Find Eccentric Anomaly from Mean Anomaly
    % % //      Solve Kepler's Equation:  M = E - e*sin(E);
    % % //      write as E = M + e*sin(E) and iterate 10 times.
    % %
    % % //  Brute Force iteration (good way to seed Newton's method which follows)
    % M  = MeanAnomalyEpoch;
    % % meanAnomaly M to eccentric Anomaly E
    % E = M + e;
    % if(M > pi)
    %     E = M - e;
    % elseif (-pi < M && M < 0)
    %     E = M - e;
    % end
    % %         //cout << " Input meanAnomaly = " << meanAnomaly << endl;
    % for i=1:10
    %     E  = M + e*sin(E);
    %     %//cout << " Mean Anomaly Solution " << E << endl;
    % end
    % % //      10 rounds of Newton's root finding method based on the above "seed".
    % for i=1:10;
    %     Eprime      = 1.0 - e*cos(E);
    %     E           = E + (M - E + e*sin(E))/Eprime;
    % end
    % KeplerInv = E;
    % %disp([" Newton-Raphson Kepler's Equation Eccentric Anomaly = ",num2str(EccentricAnomalyEpoch),"  ",num2str(KeplerInv)]);
    % %eDenom   = Eprime*Eprime*Eprime;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    QReco = rmag*cosnu*vVec + (1.0/sqrt(p))*sinnu*rVec;
    QReco = QReco/norm(QReco);
    [e a Inclination omega Omega Mp Period Jacobian] = kepler(tEpoch, rVec, vVec, mu);
    PeriodGoldenEye1 = Period*TU;
    %disp(' Kepler iSat 1 ')
    Kepler = [e a Inclination omega Omega Mp];
    Orbit  = LagrangePlanetary(Kepler, units);
end