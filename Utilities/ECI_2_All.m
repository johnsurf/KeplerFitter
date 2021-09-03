function [XEquin, XKepler, XECI] =  ECI_2_All(XECI, units)

twopi  = units.twopi;
mu     = units.mu
sqrtmu = units.sqrtmu;

rx = XECI(1);
ry = XECI(2);
rz = XECI(3);
vx = XECI(4);
vy = XECI(5);
vz = XECI(6);
t  = XECI(7);

SignOrbit = 1;

rmag       = sqrt(rx*rx + ry*ry + rz*rz);
vsq        = vx*vx + vy*vy + vz*vz;
er         = vsq - mu/rmag;
ev         = rx*vx + ry*vy + rz*vz;
ex         = (er*rx - ev*vx)/mu;
ey         = (er*ry - ev*vy)/mu;
ez         = (er*rz - ev*vz)/mu;
ecc        = sqrt(ex*ex + ey*ey + ez*ez);
%deriv      = df(ecc,rx);
Px         = ex/ecc;
Py         = ey/ecc;
Pz         = ez/ecc;
%Pvec       = avec(Px,  Py,  Pz);
hx         = ry*vz - rz*vy;
hy         = rz*vx - rx*vz;
hz         = rx*vy - ry*vx;
%hVec       = avec(hx hy hz);
%p          = dot(hVec,hVec)/mu;
p          = (hx*hx + hy*hy + hz*hz)/mu;
hmag       = sqrt(mu*p);
%deriv      = df(hmag,rx);

Wx         = hx/hmag;
Wy         = hy/hmag;
Wz         = hz/hmag;
%Wvec       = avec(Wx, Wy, Wz);
Qx         = Wy*Pz - Wz*Py;
Qy         = Wz*Px - Wx*Pz;
Qz         = Wx*Py - Wy*Px;
%Qvec       = avec(Qx, Qy, Qz);
rUnitx     = rx/rmag;
rUnity     = ry/rmag;
rUnitz     = rz/rmag;
%rUnit      = [rUnitx rUnity rUnitz]
% Ax         = -Wy;
% Ay         =  Wx;
% normN      = sqrt(Ax*Ax + Ay*Ay);
% Nx         = Ax/normN;
% Ny         = Ay/normN;
Nx         = -Wy
Ny         =  Wx;
%Nvec       = [Nx,  Ny, 0.0]
onePe      = 1.0 + ecc;
oneMe      = 1.0 - ecc;
fac        = onePe*oneMe;
if fac < 0
    SignOrbit  = -1;
end
rootfac    = sqrt(SignOrbit*fac);
a          = p/fac;
Period     = twopi*(SignOrbit*a)*sqrt(SignOrbit*a)/sqrtmu;
meanMotion = twopi/Period;

Omega      = atan2(Ny, Nx);
while Omega<0.0
    Omega = Omega + twopi;
end
while Omega>twopi
    Omega = Omega - twopi;
end
cosO       = cos(Omega);
sinO       = sin(Omega);

cosP       = Nx*Px + Ny*Py;
sinP       = (Wx*Ny - Wy*Nx)*Pz + (Nx*Py - Ny*Px)*Wz;
omega      = atan2(sinP,cosP);
while omega < 0
    omega = omega + twopi;
end
while omega > twopi
    omega = omega - twopi;
end
omega      
cosom      = cos(omega);
sinom      = sin(omega);
cosnu      = rUnitx*Px + rUnity*Py + rUnitz*Pz;
sinnu      = Wx*(Py*rUnitz - Pz*rUnity) + Wy*(Pz*rUnitx - Px*rUnitz) +  Wz*(Px*rUnity - Py*rUnitx);
nu         = atan2(sinnu,cosnu);
cosT       = ecc + cosnu;
sinT       = sinnu*rootfac;
EccentricAnomalyEpoch = atan2(sinT,cosT);
cosE       = cos(EccentricAnomalyEpoch);
sinE       = sin(EccentricAnomalyEpoch);
MeanAnomalyEpoch = EccentricAnomalyEpoch - ecc*sinE;
Mp = MeanAnomalyEpoch - meanMotion*t;
while Mp < 0
    Mp = Mp + twopi;
end
while Mp > twopi
    Mp = Mp - twopi;
end
Inclination= acos(Wz);
tanHalf = tan(Inclination/2.0);
sumAngles = omega + Omega;
while sumAngles < 0
    sumAngles = sumAngles + twopi;
end
while sumAngles > twopi
    sumAngles = sumAngles - twopi;
end

cosSum    = cos(sumAngles);
sinSum    = sin(sumAngles);

p = a*(1. + ecc)*(1. - ecc);
f = ecc*cosSum;
g = ecc*sinSum;
h = tanHalf*cosO;
k = tanHalf*sinO;
L = sumAngles + nu;
while L < 0
    L = L + twopi;
end
while L > twopi
    L = L - twopi;
end

XEquin  = [p; f; g; h; k; L; t]

XKepler = [ecc; a; Inclination; omega; Omega; Mp; t]

XECI    = [rx; ry; rz; vx; vy; vz; t]

end
