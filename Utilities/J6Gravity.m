function [J6 Jperturb Jhom] = J6Gravity(Kvec)
    % mu = 3.986004415e14;   %Graviational constant (m^3/s^2);
    mu     =  1.0;         % Canonical Graviational Parameter = GM
    Rearth =  1.0;
    % J2     =  1082.6269e-6;   %J2 pertubative constant
    % J3     =      -2.30e-6;
    % J4     =      -2.12e-6;
    % J5     =      -0.20e-6;
    % J6     =       1.00e-6;
    %------------------------------------------------------------------
    wgs84 = wgs84Constants;
    J2     = wgs84.J2;        % J2 perturbation constant
    J3     = wgs84.J3;
    J4     = wgs84.J4;
    J5     = wgs84.J5;
    J6     = wgs84.J6;
    %------------------------------------------------------------------

    X    = Kvec(1);
    Y    = Kvec(2);
    Z    = Kvec(3);
    
    K2   = dot(Kvec,Kvec);
    K    = sqrt(K2);
    K3   = K*K2;
    K4   = K2*K2;
    K5   = K*K4;
    K6   = K2*K4;
    K7   = K2*K5;
    K9   = K2*K7;
    K11  = K2*K9;
    
    %Kx   = [Kx X];
    %Ky   = [Ky Y];
    %Kz   = [Kz Z];
    
    XX   = X*X;
    XY   = X*Y;
    XZ   = X*Z;
    YY   = Y*Y;
    YZ   = Y*Z;
    ZZ   = Z*Z;
    
    XYZ  = X*Y*Z;
    
    X4   = XX*XX;
    Y4   = YY*YY;
    Z4   = ZZ*ZZ;
    XXYY = XX*YY;
    XXZZ = XX*ZZ;
    YYZZ = YY*ZZ;

    F1Fact =     XX +     YY - 6.0*ZZ;
    F2Fact = 3.0*XX + 3.0*YY - 4.0*ZZ;
    F3Fact =     XX +     YY - 4.0*ZZ;
    F4Fact = 3.0*XX +     YY - 4.0*ZZ;
    F5Fact =     XX + 3.0*YY - 4.0*ZZ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    JKepler = mu/K3;
    sdel    = Z/K;
    J0sec   = [X, Y, Z];
    J2sec   = 1.5*(J2/K2)*[X*(1-5*sdel^2), Y*(1-5*sdel^2), Z*(3-5*sdel^2)]; 
    J3sec   = 2.5*(J3/K3)*sdel*[X*(3-7*sdel^2), Y*(3-7*sdel^2), Z*(6-7*sdel^2)];
    J3sec   = J3sec - [0, 0, 1.5*J3/K2];                           
    J4sec   =-0.625*(J4/K4)*[X*(3  - 42*sdel^2+63*sdel^4),... 
                             Y*(3  - 42*sdel^2+63*sdel^4),... 
                             Z*(15 - 70*sdel^2+63*sdel^4)]; 
    J5sec   =-0.375*(J5/K5)*sdel*[X*(35  - 210*sdel^2 + 231*sdel^4),... 
                                  Y*(35  - 210*sdel^2 + 231*sdel^4),... 
                                  Z*(105 - 315*sdel^2 + 231*sdel^4)];
    J5sec   = J5sec + [0,0,1.875*J5/K4];
    J6sec   = 0.0625*(J6/K6)*[X*(35  -  945*sdel^2 + 3465*sdel^4 - 3003*sdel^6),... 
                              Y*(35  -  945*sdel^2 + 3465*sdel^4 - 3003*sdel^6),... 
                              Z*(245 - 2205*sdel^2 + 4851*sdel^4 - 3003*sdel^6)];
    Jsec    = -JKepler*J0sec; 
    Jadd    = mu/K3*[0,0, (1.5*J3/K2 - 1.875*J5/K4)];                           
    %Jperturb= -JKepler*(J2sec+J3sec+J4sec+J5sec+J6sec) + Jadd;
    %Jperturb= -JKepler*(J2sec+J3sec+J4sec+J5sec+J6sec);
    % Kepler Hessian:
    % J0hom1   = (1/K2)*[3.0*XX-K2, 3.0*XY,    3.0*XZ;...
    %                3.0*XY,    3.0*YY-K2, 3.0*YZ;...
    %                3.0*XZ,    3.0*YZ,    3.0*ZZ-K2];               
    J0hom11 = 3.0*XX - K2;
    J0hom12 = 3.0*XY;
    J0hom13 = 3.0*XZ;
    J0hom21 = 3.0*XY;
    J0hom22 = 3.0*YY - K2;
    J0hom23 = 3.0*YZ;
    J0hom31 = 3.0*XZ;    
    J0hom32 = 3.0*YZ;
    J0hom33 = 3.0*ZZ - K2;
    J0hom   = (1/K2)*[J0hom11, J0hom12, J0hom13;...
                      J0hom21, J0hom22, J0hom23;...
                      J0hom31, J0hom32, J0hom33];
    % The determinant of this matrix is almost = 2 a constant. 
    % format long;
    % determinant = det(J0hom)
    % eigenvalues = eig(J0hom)
    % format short;
    J0hom = (-mu/K3)*J0hom;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perturbations due to Oblateness: 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % J2-Gravity Model
    J2pert  = 3.0*J2*mu*Rearth/2.0;
    % Non-Homogeneous Part due to J2 (evaluated at R = Kvec)
    % factor out 3 / (2.0*K7)
    %J2sec  = (J2pert/K7)*[(5*ZZ - K2)*X, (5*ZZ - K2)*Y, (5*ZZ - 3.0*K2)*Z];
    J2sec1 = -(XX+YY-4.0*ZZ)*X;
    J2sec2 = -(XX+YY-4.0*ZZ)*Y;
    J2sec3 = -(3.0*XX+3.0*YY-2.0*ZZ)*Z;    
    J2secR = (J2pert/K7)*[J2sec1, J2sec2, J2sec3];    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Old Stuff
    % %Jhom = (J2pert/K9)*[-(K2*(YY - 4.0*ZZ + 3.0*XX) - 7.0*F3Fact*XX), 5.0*F1Fact*XY, 5.0*F1Fact*XZ;...
    % %                     5.0*F1Fact*XY, -((3.0*YY - 4.0*ZZ + XX)*K2 - 7.0*F3Fact*YY), 5.0*F1Fact*Y*Z;...
    % %                     5.0*F2Fact*XZ, 5.0*F2Fact*Y*Z, -(3.0*X4 + 6.0*XX*YY - 14.0*XX*ZZ + 3.0*Y4 - 14.0*YY*ZZ + 18.0*Z4)];   
    % %  
    % %Jhom = [7.0*F3Fact*XX-K2*F4Fact, 5.0*F1Fact*XY,                 5.0*F2Fact*XZ;...
    % %        5.0*F1Fact*XY,           7.0*F3Fact*YY - 3.0*F3Fact*K2, 5.0*F2Fact*YZ;...
    % %        5.0*F2Fact*XZ,           5.0*F2Fact*YZ,                 -3.0*X4-6.0*XX*YY+24.0*XX*ZZ-3.0*Y4+24.0*YY*ZZ-8.0*Z4];
    % %m22 = -( 3.0*X4 + 6.0*XX*YY - (24.0*XX*ZZ) + 3.0*Y4 - (24.0*YY*ZZ) + 8.0*Z4);

    % J2hom = [7.0*F3Fact*XX-K2*F4Fact, 5.0*F1Fact*XY,                 5.0*F2Fact*XZ;...
    %          5.0*F1Fact*XY,           7.0*F3Fact*YY - K2*F5Fact,     5.0*F2Fact*YZ;...
    %          5.0*F2Fact*XZ,           5.0*F2Fact*YZ,                -3.0*(XX+YY)^2+24.0*(XX+YY)*ZZ-8.0*Z4];
    % J2hom1 = (J2pert/K9)*J2hom;
    % 
    % factor out 3 / (2.0*K9)
    J2hom11 =   7.0*(XX+YY-4.0*ZZ)*XX - K2*(3.0*XX+YY-4.0*ZZ);
    J2hom12 =   5.0*(XX+YY-6.0*ZZ)*XY;
    J2hom13 =   5.0*(3.0*XX+3.0*YY - 4.0*ZZ)*XZ;

    J2hom21 =   5.0*(XX+YY-6.0*ZZ)*XY;
    J2hom22 =   7.0*(XX+YY-4.0*ZZ)*YY - K2*(XX + 3.0*YY-4.0*ZZ);
    J2hom23 =   5.0*(3.0*XX+3.0*YY - 4.0*ZZ)*YZ;

    J2hom31 =   5.0*(3.0*XX+3.0*YY - 4.0*ZZ)*XZ;
    J2hom32 =   5.0*(3.0*XX+3.0*YY - 4.0*ZZ)*YZ;
    J2hom33 =  -(3.0*X4+6.0*XXYY - (24.0*XXZZ)+3.0*Y4-(24.0*YYZZ)+8.0*Z4);

    J2hom    = [ J2hom11, J2hom12, J2hom13;...
                 J2hom21, J2hom22, J2hom23;...
                 J2hom31, J2hom32, J2hom33];

    J2hom = (J2pert/K9)*J2hom;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % J3-Gravity Model
    J3pert  = J3*mu*Rearth/2.0;
    % factor out 1 / (2.0*K9)    
    J3sec1 = -(5.0*(3.0*XX+3.0*YY-4.0*ZZ)*XZ);
    J3sec2 = -(5.0*(3.0*XX+3.0*YY-4.0*ZZ)*YZ);
    J3sec3 = (3.0*X4+6.0*XX*YY-(24.0*XX*ZZ)+3.0*Y4-(24.0*YY*ZZ)+8.0*Z4);
    J3secR = (J3pert/K9)*[J3sec1, J3sec2, J3sec3];

    % factor out 1 / (2.0*K11)
    J3hom11 =  (5.0*(18.0*X4+15.0*XX*YY-(41.0*XX*ZZ)-(3.0*Y4)+YY*ZZ+4.0*Z4)*Z);
    J3hom12 =  (105.0*(XX+YY-(2.0*ZZ))*XYZ);
    J3hom13 = -(15.0*(X4+2.0*XX*YY-(12.0*XX*ZZ)+Y4-(12.0*YY*ZZ)+8.0*Z4)*X);
    
    J3hom21 =  (105.0*(XX+YY-(2.0*ZZ))*XYZ);
    J3hom22 = -(5.0*(3.0*X4-(15.0*XX*YY)-(XX*ZZ)-(18.0*Y4)+41.0*YY*ZZ-(4.0*Z4))*Z);
    J3hom23 = -(15.0*(X4+2.0*XX*YY-(12.0*XX*ZZ)+Y4-(12.0*YY*ZZ)+8.0*Z4)*Y);
    
    J3hom31 = -(15.0*(X4+2.0*XX*YY-(12.0*XX*ZZ)+Y4-(12.0*YY*ZZ)+8.0*Z4)*X);
    J3hom32 = -(15.0*(X4+2.0*XX*YY-(12.0*XX*ZZ)+Y4-(12.0*YY*ZZ)+8.0*Z4)*Y);
    J3hom33 = -(5.0*(15.0*X4+30.0*XX*YY-(40.0*XX*ZZ)+15.0*Y4-(40.0*YY*ZZ)+8.0*Z4)*Z);
    
    J3hom    = [ J3hom11, J3hom12, J3hom13;...
                 J3hom21, J3hom22, J3hom23;...
                 J3hom31, J3hom32, J3hom33];

    J3hom = (J3pert/K11)*J3hom;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % J4-Gravity Model
    X6   = X4*XX;
    Y6   = Y4*YY;
    Z6   = Z4*ZZ;
    % Factor out 5 / (8.0 * K11)
    J4pert = J4*mu*Rearth/8.0;
    J4sec1 = 3.0*(X4+2.0*XX*YY-(12.0*XX*ZZ)+Y4-(12.0*YY*ZZ)+8.0*Z4)*X; 
    J4sec2 = 3.0*(X4+2.0*XX*YY-(12.0*XX*ZZ)+Y4-(12.0*YY*ZZ)+8.0*Z4)*Y;
    J4sec3 = (15.0*X4+30.0*XX*YY-(40.0*XX*ZZ)+15.0*Y4-(40.0*YY*ZZ)+8.0*Z4)*Z;
    J4secR = (5.0*J4pert/K11)*[J4sec1, J4sec2, J4sec3];

    % factor out 1 /  (8.0 * K13)
    J4hom11 = -(15.0*(6.0*X6+11.0*X4*YY-(101.0*X4*ZZ)+4.0*XX*Y4-(90.0*XX*YY*ZZ)+116.0*XX*Z4-Y6+11.0*Y4*ZZ+4.0*YY*Z4-(8.0*Z6)));
    J4hom12 = -(105.0*(X4+2.0*XX*YY-(16.0*XX*ZZ)+Y4-(16.0*YY*ZZ)+16.0*Z4)*XY);
    J4hom13 = -(105.0*(5.0*X4+10.0*XX*YY-(20.0*XX*ZZ)+5.0*Y4-(20.0*YY*ZZ)+8.0*Z4)*XZ);

    J4hom21 = -(105.0*(X4+2.0*XX*YY-(16.0*XX*ZZ)+Y4-(16.0*YY*ZZ)+16.0*Z4)*XY);
    J4hom22 = (15.0*(X6-(4.0*X4*YY)-(11.0*X4*ZZ)-(11.0*XX*Y4)+90.0*XX*YY*ZZ-(4.0*XX*Z4)-(6.0*Y6)+101.0*Y4*ZZ-(116.0*YY*Z4)+8.0*Z6));
    J4hom23 = -(105.0*(5.0*X4+10.0*XX*YY-(20.0*XX*ZZ)+5.0*Y4-(20.0*YY*ZZ)+8.0*Z4)*YZ);

    J4hom31 = -(105.0*(5.0*X4+10.0*XX*YY-(20.0*XX*ZZ)+5.0*Y4-(20.0*YY*ZZ)+8.0*Z4)*XZ);
    J4hom32 = -(105.0*(5.0*X4+10.0*XX*YY-(20.0*XX*ZZ)+5.0*Y4-(20.0*YY*ZZ)+8.0*Z4)*YZ);
    J4hom33 = (15.0*(5.0*X6+15.0*X4*YY-(90.0*X4*ZZ)+15.0*XX*Y4-(180.0*XX*YY*ZZ)+120.0*XX*Z4+5.0*Y6-(90.0*Y4*ZZ)+120.0*YY*Z4-(16.0*Z6)));
    
    J4hom    = [ J4hom11, J4hom12, J4hom13;...
                 J4hom21, J4hom22, J4hom23;...
                 J4hom31, J4hom32, J4hom33];
    K13  = K2*K11;
    J4hom = (J4pert/K13)*J4hom;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % J5-Gravity Model
    % Factor out 3 / (8.0 * K13)
    J5pert = 3.0*J5*mu*Rearth/8.0;
    
    J5sec1 = 7.0*(5.0*X4+10.0*XX*YY-(20.0*XX*ZZ)+5.0*Y4-(20.0*YY*ZZ)+8.0*Z4)*XZ;
    J5sec2 = 7.0*(5.0*X4+10.0*XX*YY-(20.0*XX*ZZ)+5.0*Y4-(20.0*YY*ZZ)+8.0*Z4)*YZ;
    J5sec3 = -(5.0*X6+15.0*X4*YY-(90.0*X4*ZZ)+15.0*XX*Y4-(180.0*XX*YY*ZZ)+120.0*XX*Z4+5.0*Y6-(90.0*Y4*ZZ)+120.0*YY*Z4-(16.0*Z6));
    J5secR = (J5pert/K13)*[J5sec1, J5sec2, J5sec3];

    % Factor out 21/ (8.0 * K15)
    J5hom11 = -(40.0* X6+75.0*X4*YY-(225.0*X4*ZZ)+30.0*XX*Y4-(210.0*XX*YY*ZZ)+156.0*XX*Z4-(5.0* Y6)+15.0*Y4*ZZ+12.0*YY*Z4-(8.0*Z6))*Z;
    J5hom12 = -3*(15.0*X4+30.0*XX*YY-(80.0*XX*ZZ)+15.0*Y4-(80.0*YY*ZZ)+48.0*Z4)*XYZ;
    J5hom13 = (5.0* X6+15.0*X4*YY-(120.0*X4*ZZ)+15.0*XX*Y4-(240.0*XX*YY*ZZ)+240.0*XX*Z4+5.0* Y6-(120.0*Y4*ZZ)+240.0*YY*Z4-(64.0*Z6))*X;

    J5hom21 = -3*(15.0*X4+30.0*XX*YY-(80.0*XX*ZZ)+15.0*Y4-(80.0*YY*ZZ)+48.0*Z4)*XYZ;
    J5hom22 = (5.0* X6-(30.0*X4*YY)-(15.0*X4*ZZ)-(75.0*XX*Y4)+210.0*XX*YY*ZZ-(12.0*XX*Z4)-(40.0* Y6)+225.0*Y4*ZZ-(156.0*YY*Z4)+8.0*Z6)*Z;
    J5hom23 = (5.0* X6+15.0*X4*YY-(120.0*X4*ZZ)+15.0*XX*Y4-(240.0*XX*YY*ZZ)+240.0*XX*Z4+5.0* Y6-(120.0*Y4*ZZ)+240.0*YY*Z4-(64.0*Z6))*Y;

    J5hom31 = (5.0* X6+15.0*X4*YY-(120.0*X4*ZZ)+15.0*XX*Y4-(240.0*XX*YY*ZZ)+240.0*XX*Z4+5.0* Y6-(120.0*Y4*ZZ)+240.0*YY*Z4-(64.0*Z6))*X;
    J5hom32 = (5.0* X6+15.0*X4*YY-(120.0*X4*ZZ)+15.0*XX*Y4-(240.0*XX*YY*ZZ)+240.0*XX*Z4+5.0* Z6-(120.0*Y4*ZZ)+240.0*YY*Z4-(64.0*Z6))*Y ;
    J5hom33 = (35.0* Y6+105.0*X4*YY-(210.0*X4*ZZ)+105.0*XX*Y4-(420.0*XX*YY*ZZ)+168.0*X4*Z4+35.0* Y6-(210.0*Y4*ZZ)+168.0*YY*Z4-(16.0*Z6))*Z;
    
    J5hom    = [ J5hom11, J5hom12, J5hom13;...
                 J5hom21, J5hom22, J5hom23;...
                 J5hom31, J5hom32, J5hom33];
    K15  = K2*K13;
    J5hom = (7.0*J5pert/K15)*J5hom;
    Jperturb1 = J2secR + J3secR + J4secR + J5secR;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % J6-Gravity Model
    % factored out 7.0 / (16.0*K15)
    X8   = XX*X6;
    Y8   = YY*Y6;
    Z8   = ZZ*Z6;
    J6pert = 7.0*J6*mu*Rearth/16.0;
    J6sec1 = -(5.0*X6+15.0*X4*YY-(120.0*X4*ZZ)+15.0*XX*Y4-(240.0*XX*YY*ZZ)+240.0*XX*Z4+5.0*Y6-(120.0*Y4*ZZ)+240.0*YY*Z4-(64.0*Z6))*X;
    J6sec2 = -(5.0*X6+15.0*X4*YY-(120.0*X4*ZZ)+15.0*XX*Y4-(240.0*XX*YY*ZZ)+240.0*XX*Z4+5.0*Y6-(120.0*Y4*ZZ)+240.0*YY*Z4-(64.0*Z6))*Y;
    J6sec3 = -(35.0*X6+105.0*X4*YY-(210.0*X4*ZZ)+105.0*XX*Y4-(420.0*XX*YY*ZZ)+168.0*XX*Z4+35.0*Y6-(210.0*Y4*ZZ)+168.0*YY*Z4-(16.0*Z6))*Z;
    J6secR =  (J6pert/K15)*[J6sec1, J6sec2, J6sec3];
    
    % factored out 7.0 / (16.0*K17)
    J6hom11 = (40.0*X8+115.0*X6*YY-(1235.0*X6*ZZ)+105.0*X4*Y4-(2355.0*X4*YY*ZZ)+3480.0*X4*Z4+25.0*XX*Y6-(1005.0*XX*Y4*ZZ)+3360.0*XX*YY*Z4 -(1616.0*XX*Z6)-(5.0*Y8)+115.0*Y6*ZZ-(120.0*Y4*Z4)-(176.0*YY*Z6)+64.0*Z8);
    J6hom12 = 45.0*(X6+3.0*X4*YY-(30.0*X4*ZZ)+3.0*XX*Y4-(60.0*XX*YY*ZZ)+80.0*XX*Z4+Y6-(30.0*Y4*ZZ)+80.0*YY*Z4-(32.0*Z6))*XY;
    J6hom13 =  9.0*(35.0*X6+105.0*X4*YY-(280.0*X4*ZZ)+105.0*XX*Y4-(560.0*XX*YY*ZZ)+336.0*XX*Z4+35.0*Y6-(280.0*Y4*ZZ)+336.0*YY*Z4-(64.0*Z6))*XZ;
    
    J6hom21 = 45.0*(X6+3.0*X4*YY-(30.0*X4*ZZ)+3.0*XX*Y4-(60.0*XX*YY*ZZ)+80.0*XX*Z4+Y6-(30.0*Y4*ZZ)+80.0*YY*Z4-(32.0*Z6))*XY;
    J6hom22 = -(5.0*X8-(25.0*X6*YY)-(115.0*X6*ZZ)-(105.0*X4*Y4)+1005.0*X4*YY*ZZ+120.0*X4*Z4-(115.0*XX*Y6)+2355.0*XX*Y4*ZZ-(3360.0*XX*YY*Z4)+176.0*XX*Z6-(40.0*Y8)+1235.0*Y6*ZZ-(3480.0*Y4*Z4)+1616.0*YY*Z6-(64.0*Z8));
    J6hom23 = 9.0*(35.0*X6+105.0*X4*YY-(280.0*X4*ZZ)+105.0*XX*Y4-(560.0*XX*YY*ZZ)+336.0*XX*Z4+35.0*Y6-(280.0*Y4*ZZ)+336.0*YY*Z4-(64.0*Z6))*YZ;
    
    J6hom31 = 9.0*(35.0*X6+105.0*X4*YY-(280.0*X4*ZZ)+105.0*XX*Y4-(560.0*XX*YY*ZZ)+336.0*XX*Z4+35.0*Y6-(280.0*Y4*ZZ)+336.0*YY*Z4-(64.0*Z6))*XZ;
    J6hom32 = 9.0*(35.0*X6+105.0*X4*YY-(280.0*X4*ZZ)+105.0*XX*Y4-(560.0*XX*YY*ZZ)+336.0*XX*Z4+35.0*Y6-(280.0*Y4*ZZ)+336.0*YY*Z4-(64.0*Z6))*YZ;
    J6hom33 = -(35.0*X8+140.0*X6*YY-(1120.0*X6*ZZ)+210.0*X4*Y4-(3360.0*X4*YY*ZZ)+3360.0*X4*Z4+140.0*XX*Y6-(3360.0*XX*Y4*ZZ)+6720.0*XX*YY*Z4-(1792.0*XX*Z6)+35.0*Y8-(1120.0*Y6*ZZ)+3360.0*Y4*Z4-(1792.0*YY*Z6)+128.0*Z8);
    
    J6hom    = [J6hom11, J6hom12, J6hom13;...
                J6hom21, J6hom22, J6hom23;...
                J6hom31, J6hom32, J6hom33];            
    K17  = K2*K15;
    J6hom = (J6pert/K17)*J6hom;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Outputs
    Jperturb = J2secR + J3secR + J4secR + J5secR + J6secR;  % Just the Perturbations
    J6       = -JKepler*J0sec + Jperturb;                   % Full J6 Sectoral
    Jhom     = J2hom + J3hom + J4hom + J5hom + J6hom;       % Full Hessian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end