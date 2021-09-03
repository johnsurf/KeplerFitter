function g = gravityJ4(p)	
	%4th Order Gravity Model (J4)
	%Based on "Fundamentals of Astrodynamics" Bate, Mueller, and White
	%  NOTE: valid for points above or near the earth's surface
	%  NOTE: fixed divide by zero problem above equator
	%  p = input position [X Y Z] in meters
	%
	%  g = acceleration due to gravity at the position p
	%      [Gx Gy Gz] in m/s^2
	
	%Earth Gravitational Constant
	mu = 3.986005e14; % m^3/s^2
	
	%Earth Zonal Harmonics
	J = [0.0 0.108263e-2 -0.2532e-5 -0.1611e-5 -0.23578564879393e-6];
    %------------------------------------------------------------------
    % wgs84 = wgs84Constants;
    % J2     = wgs84.J2;        % J2 perturbation constant
    % J3     = wgs84.J3;
    % J4     = wgs84.J4;
    % J5     = wgs84.J5;
    % J6     = wgs84.J6;
    % J = [0.0 J2 J3 J4 J5];
   %-------------------------------------------------------------------
	%Earth Equatorial Radius
	Re = 6378137.0; % meters
	
	x = p(1);
	y = p(2);
	z = p(3);
	r = sqrt(x*x + y*y + z*z);
	
	r_inv = 1.0/r;
	
	Redr = Re * r_inv;
	Redr_pow2 = Redr * Redr;
	Redr_pow3 = Redr_pow2 * Redr;
	Redr_pow4 = Redr_pow3 * Redr;
	
	
	zdr = z * r_inv;
	zdr_pow2 = zdr * zdr;
	zdr_pow3 = zdr_pow2 * zdr;
	zdr_pow4 = zdr_pow3 * zdr;
	
	%Terms for x and y
	xJ2 = J(2) *   1.5*Redr_pow2 * (5*zdr_pow2 - 1);
	xJ3 = J(3) *   2.5*Redr_pow3 * (3*zdr - 7*zdr_pow3);
	xJ4 = J(4) * 0.625*Redr_pow4 * (3 - 42*zdr_pow2 + 63*zdr_pow4);
	
	%Term for z
	zzJ2 = z *  J(2) * 1.5*Redr_pow2 * (3 - 5*zdr_pow2);
	zzJ3 =      J(3) * 1.5*Redr_pow3 * (z * (10*zdr - 11.66667*zdr_pow3) - r);
	zzJ4 = z *  J(4) * 0.625*Redr_pow4 * (15 - 70*zdr_pow2 + 63*zdr_pow4);
	
	g = zeros(3,1);
	
	%g(1) = -mu * x/r^3 * (1 - xJ2 + xJ3 - xJ4);
	%g(2) = -mu * y/r^3 * (1 - xJ2 + xJ3 - xJ4);
	%g(3) = -mu * z/r^3 * (1 + zJ2 + zJ3 - zJ4);
	
	mmudr_pow3 = -mu * r_inv * r_inv * r_inv;
	CXY = mmudr_pow3 * (1 - xJ2 + xJ3 - xJ4);
	
	g(1) = CXY * x;
	g(2) = CXY * y;
	g(3) = mmudr_pow3 * (z + zzJ2 + zzJ3 - zzJ4);
	
end