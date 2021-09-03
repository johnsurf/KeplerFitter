function [XECI, JacEquinoctial_2_ECI] =  Equinoctial_2_ECI(XEquin, units)
    %orbitparams = wgs84Constants;
    %mu = orbitparams.mu;
    mu = units.mu;

    p = XEquin(1);
    f = XEquin(2);
    g = XEquin(3);
    h = XEquin(4);
    k = XEquin(5);
    l = XEquin(6);
    t = XEquin(7);
    
    cosL = cos(l);
    sinL = sin(l);
    rtpinv  = sqrt(mu/p);
    
    
    alf2    =  h^2 - k^2;
    s2      =  1 + h^2 + k^2;
    w       =  1 + f*cosL + g*sinL;
    r       =  p/w;
    vcoeff  =  (1./s2)*rtpinv;
    

    rx   =    (r/s2)*(cosL*(1. + alf2) + 2.*h*k*sinL);
    ry   =    (r/s2)*(sinL*(1. - alf2) + 2.*h*k*cosL);
    rz   =  2.*(r/s2)*(h*sinL - k*cosL);

    vx   =    -vcoeff*( sinL*(1. + alf2)  - 2.*h*k*(f + cosL) + (1. + alf2)*g);
    vy   =    -vcoeff*( cosL*(alf2 - 1.)  + 2.*h*k*(g + sinL) + (alf2 - 1.)*f);
    vz   =  2.*vcoeff*( h*(f+cosL) + k*(g+sinL) );

    XECI = [rx; ry; rz; vx; vy; vz; t];
    
    JacEquinoctial_2_ECI(1,1) =  -((k^2-1.0-h^2)*cosL-(2.0*h*k*sinL))/((g*sinL+1.0+cosL*f)*(k^2+1.0+h^2));

    JacEquinoctial_2_ECI(1,2) =  (((k^2-1.0-h^2)*cosL-(2.0*h*k*sinL))*cosL*p)/(((g*sinL+1.0+cosL*f)*(g* sinL+1.0+cosL*f))*(k^2+1.0+h^2));

    JacEquinoctial_2_ECI(1,3) =  (((k^2-1.0-h^2)*cosL-(2.0*h*k*sinL))*p*sinL)/(((g*sinL+1.0+cosL*f)*(g*  sinL+1.0+cosL*f))*(k^2+1.0+h^2));

    JacEquinoctial_2_ECI(1,4) =  (2.0*((k^2+1.0-h^2)*sinL+2.0*cosL*h*k)*k*p)/((g*sinL+1.0+cosL*f)*((k^2+1.0+h^2)*(k^2+1.0+h^2)));

    JacEquinoctial_2_ECI(1,5) =  -(2.0*((k^2-1.0-h^2)*h*sinL+2.0*(h^2+1.0)*cosL*k)*p)/((g*sinL+1.0+cosL*f)*( (k^2+1.0+h^2)*(k^2+1.0+h^2)));

    JacEquinoctial_2_ECI(1,6) =  ((((k^2-1.0-h^2)*cosL-(2.0*h*k*sinL))*(cosL*g-(f*sinL))+((k^2-1.0-h^2)* sinL+2.0*cosL*h*k)*(g*sinL+1.0+cosL*f))*p)/(((g*sinL+1.0+cosL*f)*(g*sinL+1.0+ cosL*f))*(k^2+1.0+h^2));

    JacEquinoctial_2_ECI(2,1) =  ((k^2+1.0-h^2)*sinL+2.0*cosL*h*k)/((g*sinL+1.0+cosL*f)*(k^2+1.0+h^2));

    JacEquinoctial_2_ECI(2,2) =  -(((k^2+1.0-h^2)*sinL+2.0*cosL*h*k)*cosL*p)/(((g*sinL+1.0+cosL*f)*(g*sinL+ 1.0+cosL*f))*(k^2+1.0+h^2));

    JacEquinoctial_2_ECI(2,3) =  -(((k^2+1.0-h^2)*sinL+2.0*cosL*h*k)*p*sinL)/(((g*sinL+1.0+cosL*f)*(g*sinL+ 1.0+cosL*f))*(k^2+1.0+h^2));

    JacEquinoctial_2_ECI(2,4) =  (2.0*((k^2+1.0-h^2)*cosL*k-(2.0*(k^2+1.0)*h*sinL))*p)/((g*sinL+1.0+cosL*f )*((k^2+1.0+h^2)*(k^2+1.0+h^2)));

    JacEquinoctial_2_ECI(2,5) =  -(2.0*((k^2-1.0-h^2)*cosL-(2.0*h*k*sinL))*h*p)/((g*sinL+1.0+cosL*f)*((k^2+ 1.0+h^2)*(k^2+1.0+h^2)));

    JacEquinoctial_2_ECI(2,6) =  ((((k^2+1.0-h^2)*cosL-(2.0*h*k*sinL))*(g*sinL+1.0+cosL*f)-(((k^2+1.0-(h^2 ))*sinL+2.0*cosL*h*k)*(cosL*g-(f*sinL))))*p)/(((g*sinL+1.0+cosL*f)*(g*sinL+1.0+ cosL*f))*(k^2+1.0+h^2));

    JacEquinoctial_2_ECI(3,1) =  -(2.0*(cosL*k-(h*sinL)))/((g*sinL+1.0+cosL*f)*(k^2+1.0+h^2));

    JacEquinoctial_2_ECI(3,2) =  (2.0*(cosL*k-(h*sinL))*cosL*p)/(((g*sinL+1.0+cosL*f)*(g*sinL+1.0+cosL*f))*(k *k+1.0+h^2));

    JacEquinoctial_2_ECI(3,3) =  (2.0*(cosL*k-(h*sinL))*p*sinL)/(((g*sinL+1.0+cosL*f)*(g*sinL+1.0+cosL*f))*(k *k+1.0+h^2));

    JacEquinoctial_2_ECI(3,4) =  (2.0*((k^2+1.0-h^2)*sinL+2.0*cosL*h*k)*p)/((g*sinL+1.0+cosL*f)*((k^2+1.0+h *h)*(k^2+1.0+h^2)));

    JacEquinoctial_2_ECI(3,5) =  (2.0*((k^2-1.0-h^2)*cosL-(2.0*h*k*sinL))*p)/((g*sinL+1.0+cosL*f)*((k^2+1.0 +h^2)*(k^2+1.0+h^2)));

    JacEquinoctial_2_ECI(3,6) =  (2.0*(((f*h+g*k)*cosL+h)*cosL+((g*sinL+1.0)*k+f*h*sinL)*sinL)*p)/(((g*sinL+ 1.0+cosL*f)*(g*sinL+1.0+cosL*f))*(k^2+1.0+h^2));

    JacEquinoctial_2_ECI(4,1) =  -(((k^2-1.0-h^2)*(g+sinL)+2.0*(cosL+f)*h*k)*mu)/(2.0*rtpinv* (k^2+1.0+h^2)*(p*p));

    JacEquinoctial_2_ECI(4,2) =  (2.0*rtpinv*h*k)/(k^2+1.0+h^2);

    JacEquinoctial_2_ECI(4,3) =  (rtpinv*(k^2-1.0-h^2))/(k^2+1.0+h^2);

    JacEquinoctial_2_ECI(4,4) =  (2.0*rtpinv*((k^2+1.0-h^2)*f-(2.0*(g+sinL)*h*k)+(k^2+1.0-(h* h))*cosL)*k)/((k^2+1.0+h^2)*(k^2+1.0+h^2));

    JacEquinoctial_2_ECI(4,5) =  -(2.0*rtpinv*(((k^2-1.0-h^2)*(g+sinL)+2.0*(cosL+f)*h*k)*k-(( (cosL+f)*h+(g+sinL)*k)*(k^2+1.0+h^2))))/((k^2+1.0+h^2)*(k^2+1.0+h^2));

    JacEquinoctial_2_ECI(4,6) =  (rtpinv*((k^2-1.0-h^2)*cosL-(2.0*h*k*sinL)))/(k^2+1.0+h^2);

    JacEquinoctial_2_ECI(5,1) =  -(((k^2+1.0-h^2)*(cosL+f)-(2.0*(g+sinL)*h*k))*mu)/(2.0*sqrt((mu/p) )*(k^2+1.0+h^2)*(p*p));

    JacEquinoctial_2_ECI(5,2) =  (rtpinv*(k^2+1.0-h^2))/(k^2+1.0+h^2);

    JacEquinoctial_2_ECI(5,3) =  -(2.0*rtpinv*h*k)/(k^2+1.0+h^2);

    JacEquinoctial_2_ECI(5,4) =  -(2.0*rtpinv*(((k^2+1.0-h^2)*(cosL+f)-(2.0*(g+sinL)*h*k))*h+ ((cosL+f)*h+(g+sinL)*k)*(k^2+1.0+h^2)))/((k^2+1.0+h^2)*(k^2+1.0+h^2));

    JacEquinoctial_2_ECI(5,5) =  (2.0*rtpinv*((k^2-1.0-h^2)*(g+sinL)+2.0*f*h*k+2.0*cosL*h*k)* h)/((k^2+1.0+h^2)*(k^2+1.0+h^2));

    JacEquinoctial_2_ECI(5,6) =  -(rtpinv*((k^2+1.0-h^2)*sinL+2.0*cosL*h*k))/(k^2+1.0+h^2);

    JacEquinoctial_2_ECI(6,1) =  -(((cosL+f)*h+(g+sinL)*k)*mu)/(rtpinv*(k^2+1.0+h^2)*(p*p));

    JacEquinoctial_2_ECI(6,2) =  (2.0*rtpinv*h)/(k^2+1.0+h^2);

    JacEquinoctial_2_ECI(6,3) =  (2.0*rtpinv*k)/(k^2+1.0+h^2);

    JacEquinoctial_2_ECI(6,4) =  -(2.0*rtpinv*(2.0*((cosL+f)*h+(g+sinL)*k)*h-((k^2+1.0+h^2)*( cosL+f))))/((k^2+1.0+h^2)*(k^2+1.0+h^2));

    JacEquinoctial_2_ECI(6,5) =  -(2.0*rtpinv*(2.0*((cosL+f)*h+(g+sinL)*k)*k-((k^2+1.0+h^2)*(g+ sinL))))/((k^2+1.0+h^2)*(k^2+1.0+h^2));

    JacEquinoctial_2_ECI(6,6) =  (2.0*rtpinv*(cosL*k-(h*sinL)))/(k^2+1.0+h^2);
    
end