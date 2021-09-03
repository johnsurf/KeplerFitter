function wgs84 = wgs84Constants()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%These Constants Need to be Collected and Stored in a Common Global Area
%           Constants
% Constants used to convert to Canonical Units
% Bulk parameters
% Mass (1024 kg)                  5.9723
% Volume (1010 km^3)               108.321
% Equatorial radius (km)	      6378.137
% Polar radius (km)               6356.752
% Volumetric mean radius (km)     6371.000
% Core radius (km)                3485
% Ellipticity (Flattening)        0.003353
% Mean density (kg/m^3)            5514
% Surface gravity (m/s^2)          9.798
% Surface acceleration (m/s^2)     9.780
% Escape velocity (km/s)          11.186
% GM (x 10^6 km^3/s^2)               0.3986012
% Bond albedo                     0.306
% Geometric albedo                0.434
% V-band magnitude V(1,0)         -3.99
% Solar irradiance (W/m^2)         1361.0
% Black-body temperature (K)       254.0
% Topographic range (km)            20.4
% Moment of inertia (I/MR^2)        0.3308
% J2 [x 10^(-6)]                     1082.63
% Number of natural satellites       1
% Planetary ring system             No
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Units
wgs84.twopi  = 2.0*pi;
%Re = 6378.1363*10^3;          %Radius of Earth (m)
wgs84.Grav   = 3.986004418e14; %mu Gravitational constant (m^3/s^2)
%wgs84.GM_m3ps2        = 3986004.418e8;
%wgs84.J2     = 1082.6269e-6;   %J2 perturbation constant
%wgs84.J3     =     -2.30e-6;
%wgs84.J4     =     -2.12e-6;
%wgs84.J5     =     -0.20e-6;
%wgs84.J6     =      1.00e-6;

%Gravity Models
wgs84.J2     =  1.08262668355e-3;
wgs84.J3     = -2.53265648533e-6;
wgs84.J4     = -1.61962159137e-6;
wgs84.J5     = -2.27296082869e-7;
wgs84.J6     =  5.40681239107e-7;

% Constants from gravityJ4
% %Earth Gravitational Constant
% mu = 3.986005e14; % m^3/s^2	
% %Earth Zonal Harmonics
% J = [0.0 0.108263e-2 -0.2532e-5 -0.1611e-5 -0.23578564879393e-6];	
% %Earth Equatorial Radius
% Re = 6378137.0; % meters
    
wgs84.XMNPDA = 1440.;          % 1440 minutes per day
%wgs84.XKE    = 0.0743669161;   % Hoot's ke Gravitational constant for Earth in RE^(1.5)/minutes
%wgs84.XKE   = 7.43668477395274e-2;   % Hoot's ke Gravitational constant for Earth in RE^(1.5)/minutes
wgs84.XKE    = 7.43668599821643e-2;   % Hoot's ke Gravitational constant for Earth in RE^(1.5)/minutes
%    mu = 3.986004415e14; %gravitational constant
%    J2 = 1082.6269e-6;   %J2 perturbation constant
%    Re = 6378.1363*10^3; %Radius of Earth (m)
%1 ER = DU km, 1 time unit (1 TU) = 806.81 seconds
wgs84.DU     = 6378136.3;      % Mean Equatorial Radius of the Earth in meters
wgs84.TU     = 60.0/wgs84.XKE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adam's numbers: 
% %Loads WGS84 Properties
% wgs84.a_m             = 6378137.0; %m
% wgs84.b_m             = 6356752.3; %m
% wgs84.f               = (wgs84.a_m-wgs84.b_m)/wgs84.a_m;
% wgs84.e               = sqrt(wgs84.a_m^2 - wgs84.b_m^2)/wgs84.a_m;
% wgs84.omega_earth_rps = 7.2921150e-5;
% wgs84.GM_m3ps2        = 3986004.418e8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modification of Units
%wgs84.DU     = wgs84.DU/8.0;
%wgs84.TU     = wgs84.TU/4.0; 
%cout << " Earth Time Unit (seconds) = " << TU << endl;
% COMPUTATION OF UNITS CONVERSION: 
% CAVIAT -- IF YOU MODIFY Length and Time from the above Canonical Settings
% Then "mu" will no longer be equal to one and you will have to make
% adjustments all over the code!
% conversion factor VU: 1 = (DU-m/1DU)*(1TU/TU-seconds)
% Application vConv = Vel(m/s)/VU
wgs84.VU      = wgs84.DU/wgs84.TU;
wgs84.AU      = wgs84.VU/wgs84.TU;
wgs84.GPS1980 = 1.0874984267e+09;  % Time offset in seconds from 1980 
%wgs84.mu      = 1.0;               % Canonical Gravitation Parameter = GM
TUsquared    = wgs84.TU^2; 
DUcubed      = wgs84.DU^3;
wgs84.mu     = wgs84.Grav*TUsquared/DUcubed;
wgs84.sqrtmu = sqrt(wgs84.mu);
wgs84.Rearth = 6378136.3/wgs84.DU;      % Mean Equatorial Radius of the Earth in meters/DU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%