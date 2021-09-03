function [Gravity Phi, InHomogenous, Homogenous] = Escobal(Rvec, mu)

    DataList = zeros(30,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Constants
    % twopi  = 2.0*pi;
    wgs84  = wgs84Constants(); 
    % DU     = wgs84.DU;
    % TU     = wgs84.TU;
    % mu     = wgs84.mu;        % Canonical Gravitation Parameter    = GM = k^2*m
    % Rearth = wgs84.Rearth;   % Radius of earth in Canonical Units = 1
    J2     = wgs84.J2;        % J2 perturbation constant
    J3     = wgs84.J3;
    J4     = wgs84.J4;
    J5     = wgs84.J5;
    J6     = wgs84.J6;
    % Truncate Gravity Model
    %J2     = 0.0; 
    %J3     = 0.0; 
    %J4     = 0.0; 
    %J5     = 0.0; 
    %J6     = 0.0; 
    
    % Include Kepler Term CK = 1; 
    %CK     = 1.0;
    CK     = 0.0;
    
    C2     = 0.5000*J2;
    C3     = 0.5000*J3;
    C4     =-0.1250*J4;
    C5     =-0.1250*J5;
    C6     = 0.0625*J6;
    % Just use J4 for now. 
    %C5     = 0.0;
    %C6     = 0.0; 
    %Kvec   = Rvec(1:3)/DU;
    Kvec   = Rvec(1:3);   % For Rvec in Canonical Units
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X    = Kvec(1);                                              % line  1
    Y    = Kvec(2);                                              % line  2
    Z    = Kvec(3);                                              % line  3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    R2   = X^2 + Y^2 + Z^2;                                      % line  4
    R1   = sqrt(R2);                                             % line  5 
    R3   = R1*R2;                                                % line  6
    R4   = R1*R3;                                                % line  7
    R5   = R1*R4;                                                % line  8
    R6   = R1*R5;                                                % line  9
    
    sd1  = Z/R1;                                                 % line 10
    sd2  = sd1^2;                                                % line 11
    sd3  = sd1*sd2;                                              % line 12
    sd4  = sd1*sd3;                                              % line 13
    sd5  = sd1*sd4;                                              % line 14
    sd6  = sd1*sd5;                                              % line 15
    
    F2   =  1.0      -  3.0*sd2;                                 % line 16
    F3   =  3.0*sd1  -  5.0*sd3;                                 % line 17
    F4   =  3.0      -  30.*sd2  + 35*sd4;                       % line 18
    F5   =  15.*sd1  -  70.*sd3  + 63*sd5;                       % line 19
    F6   =  5.       - 105.*sd2  + 315*sd4 - 231*sd6;            % line 20
    
    V1      = 1.0/R1;                                            % line 21
    V2      = C2*F2/R2;                                          % line 22
    V3      = C3*F3/R3;                                          % line 23
    V4      = C4*F4/R4;                                          % line 24
    V5      = C5*F5/R5;                                          % line 25
    V6      = C6*F6/R6;                                          % line 26
    
    V       = mu*V1;           % Kepler Potential                % line 27
    % J6:
    D       = V2 + V3 + V4 + V5 + V6;                            % line 28
    % J4: 
    %D       = V2 + V3 + V4;                                      % line 28
    % J2: 
    %D       = V3;                                                % line 28
    %R       = mu*D;                                              % line 29
    %Phi     = V + R;                                             % line 30
    Phi     = V*(CK + D);                                         % line 30
    
    DataList( 1) = X;  % line  1
    DataList( 2) = Y;  % line  2
    DataList( 3) = Z;  % line  3
    DataList( 4) = R2; % line  4
    DataList( 5) = R1; % line  5 
    DataList( 6) = R3; % line  6
    DataList( 7) = R4; % line  7
    DataList( 8) = R5; % line  8
    DataList( 9) = R6; % line  9
    DataList(10) = sd1;% line 10
    DataList(11) = sd2;% line 11
    DataList(12) = sd3;% line 12
    DataList(13) = sd4;% line 13
    DataList(14) = sd5;% line 14
    DataList(15) = sd6;% line 15
    DataList(16) = F2; % line 16
    DataList(17) = F3; % line 17
    DataList(18) = F4; % line 18
    DataList(19) = F5; % line 19
    DataList(20) = F6; % line 20
    DataList(21) = V1; % line 21
    DataList(22) = V2; % line 22
    DataList(23) = V3; % line 23
    DataList(24) = V4; % line 24
    DataList(25) = V5; % line 25
    DataList(26) = V6; % line 26
    DataList(27) = V;  % Kepler Potential % line 27
    DataList(28) = D;  % line 28
    %R       = mu*D;   % line 29
    %Phi     = V + R;  % line 30
    DataList(30) = Phi;% V*(CK + D); % line 30
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %Jperturb = J2secR + J3secR + J4secR + J5secR + J6secR;  % Just the Perturbations
    %J6       = -JKepler*J0sec + Jperturb;                   % Full J6 Sectoral
    %Jhom     = J2hom + J3hom + J4hom + J5hom + J6hom;       % Full Hessian
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Debugging First Derivatives: 
    % format short
    % for LineCheck = 1:30
    % ForDerivs  = [];
    % F            = zeros(30,1);
    % F(LineCheck) = 1;
    % F            = backdiff(F);
    % forward      = [];
    % for j = 1:30 % check range of derivatives wrt forwards differentiation
    %     Q        = zeros(30,1);
    %     Q(j)     = 1;
    %     Q        = fordiff(Q);
    %     % Q(LineCheck) = dK_LineCheck/dLj
    %     forward  = [forward; j Q(LineCheck)];
    % end
    % ForDerivs    = [ForDerivs; forward];
    % Derivs       = [ForDerivs F (F-ForDerivs(:,2))]
    % %Derivs(1:LineCheck,:)
    % 'Debug Point'
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      Begin Finite Differences Check
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for LineCheck = 1:30
    %     [H_I, Jacob, Hess] = WorkOrder(LineCheck)
    %     Hess - Hess'
    %     % Sum = sum(sum(Hess - Hess'));
    %     % disp([' Line = ', num2str(LineCheck),'  sum = ', num2str(Sum)])
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % %%% TEST SECOND DERIVATIVES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %Debug Second Derivatives by checking symmetry:
    % %Begin Checking Out Second Derivatives:
    % Nvar     = 7;
    % % LineCheck    = 56;  % Kepler Part
    % % LineCheck    = 90;  % Escobal Gravitation to J6
    % % F            = zeros(90,1);
    % % F(LineCheck) = 1;
    % % F = this.backdiff(F);
    % Nlines = 133;
    % format short
    % for LineCheck = 1:Nlines
    %     ListMax   = LineCheck + 1;
    %     %Nvar = LineCheck;
    %     Hessian      = [];
    %     F            = zeros(Nlines,1);
    %     F(LineCheck) = 1;
    %     F = this.backdiff(F);
    %     %F(ListMax:Nlines) = 0;
    %     % for j = 1:Nvar
    %     for j = 1:LineCheck
    %         Q = zeros(Nlines,1);
    %         S = zeros(Nlines,1);
    %         Q(j) = 1.0;
    %         Q = this.fordiff(Q);
    %         %Q(ListMax:Nlines) = 0;
    %         S = this.secdiff(F,Q,S);
    %         %S(ListMax:Nlines) = 0;
    %         S = this.backdiff(S);
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
    %     max(max(Symmetry))
    %     disp([' Line Check = ', num2str(LineCheck)]);
    %     %disp('Debug Point')
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F     = zeros(30,1);
    % Let's computed the total J6 Gravity
    LineCheck = 30; 
    F(LineCheck) = 1.0;
    % Let's get the gradient of the disturbing potential
    % F(29) = 1.0;
    F           = backdiff(F);     
    Gravity     = [F(1); F(2); F(3)];  
    %GravityMKS  = Gravity*DU/TU/TU;
    [H_I, Jacob, Hess] = WorkOrder(LineCheck);
    InHomogenous = [0.0; 0.0; 0.0; Jacob'];
     Homogenous  = [zeros(3), eye(3)   ;...
                    Hess,     zeros(3)];

    
    function F = backdiff(F)
        % from line 30: Phi     = V + R;
        %F(27) = F(27) + F(30);
        %F(29) = F(29) + F(30);
        % from line 30: Phi     = V*(CK + D);
        F(27) = F(27) + F(30)*(CK + D);
        F(28) = F(28) + F(30)*V;
        % R is the disturbing potential
        % from line 29: R       = mu*D;
        % F(28) = F(28) + F(29)*mu;
        % Up to J6: 
        % from line 28: D       = V2 + V3 + V4 + V5 + V6;
        F(22) = F(22) + F(28);
        F(23) = F(23) + F(28);
        F(24) = F(24) + F(28);
        F(25) = F(25) + F(28);
        F(26) = F(26) + F(28);
        % from line 27: V       = mu*V1;
        F(21) = F(21) + F(27)*mu;
        % from line 26: V6      = C6*F6/R6;
        F(20) = F(20) + F(26)*C6/R6;
        F( 9) = F( 9) - F(26)*V6/R6;
        % from line 25: V5      = C5*F5/R5;
        F(19) = F(19) + F(25)*C5/R5;
        F( 8) = F( 8) - F(25)*V5/R5;
        % from line 24: V4      = C4*F4/R4;
        F(18) = F(18) + F(24)*C4/R4;
        F( 7) = F( 7) - F(24)*V4/R4;
        % from line 23: V3      = C3*F3/R3;
        F(17) = F(17) + F(23)*C3/R3;
        F( 6) = F( 6) - F(23)*V3/R3;
        % from line 22: V2      = C2*F2/R2;
        F(16) = F(16) + F(22)*C2/R2;
        F( 4) = F( 4) - F(22)*V2/R2;
        % from line 21: V1      = 1/R1;
        F( 5) = F( 5) - F(21)*V1/R1;
        % from line 20: F6   =  5.       - 105.*sd2  + 315*sd4 - 231*sd6;
        F(11) = F(11) - F(20)*105.;
        F(13) = F(13) + F(20)*315.;
        F(15) = F(15) - F(20)*231.;
        % from line 19: F5   =  15.*sd1  -  70.*sd3  + 63*sd5;
        F(10) = F(10) + F(19)*15.;
        F(12) = F(12) - F(19)*70.;
        F(14) = F(14) + F(19)*63.;
        % from line 18: F4   =  3.0      -  30.*sd2  + 35*sd4;
        F(11) = F(11) - F(18)*30.;
        F(13) = F(13) + F(18)*35.;
        % from line 17: F3   =  3.0*sd1  -  5.0*sd3;
        F(10) = F(10) + F(17)*3.;
        F(12) = F(12) - F(17)*5.;
        % from line 16: F2   =  1.0      -  3.0*sd2;
        F(11) = F(11) - F(16)*3.;
        % from line 15: sd6  = sd1*sd5;
        F(10) = F(10) + F(15)*sd5;
        F(14) = F(14) + F(15)*sd1;
        % from line 14: sd5  = sd1*sd4;
        F(10) = F(10) + F(14)*sd4;
        F(13) = F(13) + F(14)*sd1;
        % from line 13: sd4  = sd1*sd3;
        F(10) = F(10) + F(13)*sd3;
        F(12) = F(12) + F(13)*sd1;
        % from line 12: sd3  = sd1*sd2;
        F(10) = F(10) + F(12)*sd2;
        F(11) = F(11) + F(12)*sd1;
        % from line 11: sd2  = sd1^2;
        F(10) = F(10) + F(11)*2.0*sd1;
        % from line 10: sd1  = Z/R1;
        F( 3) = F( 3) + F(10)/R1;
        F( 5) = F( 5) - F(10)*sd1/R1;
        % from line  9: R6   = R1*R5;
        F( 5) = F( 5) + F( 9)*R5;
        F( 8) = F( 8) + F( 9)*R1;
        % from line  8: R5   = R1*R4;
        F( 5) = F( 5) + F( 8)*R4;
        F( 7) = F( 7) + F( 8)*R1;
        % from line  7: R4   = R1*R3;
        F( 5) = F( 5) + F( 7)*R3;
        F( 6) = F( 6) + F( 7)*R1;
        % from line  6: R3   = R1*R2;
        F( 5) = F( 5) + F( 6)*R2;
        F( 4) = F( 4) + F( 6)*R1;
        % from line  5: R1   = sqrt(R2);
        F( 4) = F( 4) + F( 5)*0.5/R1;
        % from line  4: R2   = X^2 + Y^2 + Z^2;
        F( 1) = F( 1) + F( 4)*2.0*X;
        F( 2) = F( 2) + F( 4)*2.0*Y;
        F( 3) = F( 3) + F( 4)*2.0*Z;
    end
    function Q = fordiff(Q)
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % X    = Kvec(1);                                              % line  1
        % Y    = Kvec(2);                                              % line  2
        % Z    = Kvec(3);                                              % line  3
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % from line 4: R2   = X^2 + Y^2 + Z^2;                         % line  4
        Q( 4) = Q( 4) + 2.0*(Q( 1)*X + Q( 2)*Y + Q( 3)*Z);
        % from line 5: R1   = sqrt(R2);                                % line  5
        Q( 5) = Q( 5) + Q( 4)*0.5/R1;
        % from line 6: R3   = R1*R2;
        Q( 6) = Q( 6) + Q( 5)*R2 + Q( 4)*R1;
        % from line 7: R4   = R1*R3;                                   % line  7
        Q( 7) = Q( 7) + Q( 5)*R3 + Q( 6)*R1;
        % from line 8: R5   = R1*R4;                                   % line  8
        Q( 8) = Q( 8) + Q( 5)*R4 + Q( 7)*R1;
        % from line 9: R6   = R1*R5;                                   % line  9
        Q( 9) = Q( 9) + Q( 5)*R5 + Q( 8)*R1;
        % 
        % from line 10: sd1  = Z/R1;
        Q(10) = Q(10) + Q( 3)/R1 - Q( 5)*sd1/R1;
        % from line 11: sd2  = sd1^2;
        Q(11) = Q(11) + Q(10)*2.0*sd1;
        % from line 12: sd3  = sd1*sd2;
        Q(12) = Q(12) + Q(10)*sd2 + Q(11)*sd1;
        % from line 13: sd4  = sd1*sd3;
        Q(13) = Q(13) + Q(10)*sd3 + Q(12)*sd1;
        % from line 14: sd5  = sd1*sd4;
        Q(14) = Q(14) + Q(10)*sd4 + Q(13)*sd1;
        % from line 15: sd6  = sd1*sd5;
        Q(15) = Q(15) + Q(10)*sd5 + Q(14)*sd1;
        % 
        % from line 16: F2   =  1.0      -  3.0*sd2;
        Q(16) = Q(16) - Q(11)*3.;
        % from line 17: F3   =  3.0*sd1  -  5.0*sd3;
        Q(17) = Q(17) + Q(10)*3.0 - Q(12)*5.0;
        % from line 18: F4   =  3.0      -  30.*sd2  + 35*sd4;
        Q(18) = Q(18) - Q(11)*30.0 + Q(13)*35.0;
        % from line 19: F5   =  15.*sd1  -  70.*sd3  + 63*sd5;
        Q(19) = Q(19) + Q(10)*15. - Q(12)*70.  + Q(14)*63.;
        % from line 20: F6   =  5.       - 105.*sd2  + 315*sd4 - 231*sd6;
        Q(20) = Q(20) - Q(11)*105. + Q(13)*315. - Q(15)*231.;
        % 
        % from line 21: V1    = 1/R1; Kepler Potential                 % line 21
        Q(21) = Q(21) - Q( 5)*V1/R1;
        % from line 22: V2      = C2*F2/R2;
        Q(22) = Q(22) + Q(16)*C2/R2 - Q( 4)*V2/R2;
        % % from line 23: V3      = C3*F3/R3;
        Q(23) = Q(23) + Q(17)*C3/R3 - Q(16)*V3/R3;
        % from line 24: V4      = C4*F4/R4;
        Q(24) = Q(24) + Q(18)*C4/R4 - Q( 7)*V4/R4;
        % from line 25: V5      = C5*F5/R5;
        Q(25) = Q(25) + Q(19)*C5/R5 - Q( 8)*V5/R5;
        % from line 26: V6      = C6*F6/R6;
        Q(26) = Q(26) + Q(20)*C6/R6 - Q( 9)*V6/R6;
        % from line 27: V       = mu*V1;
        Q(27) = Q(27) + Q(21)*mu;
        % Up to J6
        % from line 28: D       = V2 + V3 + V4 + V5 + V6;
        Q(28) = Q(28) + Q(22) + Q(23) + Q(24) + Q(25) + Q(26);
        % from line 29: R       = mu*D;
        %Q(29) = Q(29) + Q(28)*mu;
        % % from line 30: Phi     = V + R;
        %Q(30) = Q(30) + Q(27) + Q(29);
        % from line 30: Phi     = V*(CK + D);
        Q(30) = Q(30) + Q(27)*(CK + D) + Q(28)*V;

    end

        function S = secdiff(f, q, s)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            % Escobal Potential starts here
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % X    = Kvec(1);                                     % line  1
            % Y    = Kvec(2);                                     % line  2
            % Z    = Kvec(3);                                     % line  3
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % from line  1: X    = rorbitx = Kvec(1);             % line  1
            % All 2nd Derivatives vanish
            % from line  2: Y    = rorbity = Kvec(2);             % line  2
            % All 2nd Derivatives vanish
            % from line  3: Z    = rorbitz = Kvec(3);             % line  3
            % All 2nd Derivatives vanish
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % from line  4: R2   = X^2 + Y^2 + Z^2;               % line  4
            % K =  4
            %// X: J =  1  KJ = 2.0*X
            s( 1) = s( 1) + 2.0*f( 4)*q( 1);
            %// Y: J =  2  KJ = 2.0*Y
            s( 2) = s( 2) + 2.0*f( 4)*q( 2);
            %// Z: J =  3  KJ = 2.0*Z
            s( 3) = s( 3) + 2.0*f( 4)*q( 3);
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % from line 6 : R1   = rorbit;                        % line  4
            % All 2nd Derivatives vanish
            % from line  5: R1   = sqrt(R2);                      % line  5
            %// R2: J =  4  KJ = 0.5/sqrt(R2)
            s( 4) = s( 4) - f( 5)*q( 4)*0.25/R3;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % from line  4: R2   = R1^2;                         % line   5
            %// R1: J =  5  KJ = 2.0*R1   I =  5
            %s( 5) = s( 5) + 2.0*f( 4)*q( 5);
            
            % from line  6: R3   = R1*R2;                        % line   6
            %// R1: J =  5 KJ = R2  I =  4
            s( 5) = s( 5) + f( 6)*q( 4); 
            %// R2: J =  4 KJ = R1 I =  5
            s( 4) = s( 4) + f( 6)*q( 5); 
            
            % from line  7: R4   = R1*R3;                        % line   7
            %// R1: J =  5 KJ = R3  I =  6
            s( 5) = s( 5) + f( 7)*q( 6);
            %// R3: J =  6 KJ = R1
            s( 6) = s( 6) + f( 7)*q( 5);
            
            % from line 8: R5   = R1*R4;                          % line  8
            % K = 8
            %// R1: J = 5  KJ = R4 I = 7
            s( 5) = s( 5) + f( 8)*q( 7);
            %// R4: J = 7  KJ = RI I = 5
            s( 7) = s( 7) + f( 8)*q( 5);
            
            % from line 9: R6   = R1*R5;                          % line  9
            %// R1: J = 5  KJ = R5 I = 8
            s( 5) = s( 5) + f( 9)*q( 8);
            %// R5: J = 8  KJ = R1 I = 5
            s( 8) = s( 8) + f( 9)*q( 5);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % from line 10: sd1  = rorbitz/R1;
            %//rorbitz: J = 3 KJ = 1/R1 I = 5
            s( 3) = s( 3) - f(10)*q( 5)/R2;
            %// R1: J =  5 KJ = -Z*R1^(-2) I = 3,  5
            s( 5) = s( 5) + f(10)*(2.0*q( 5)*sd1 - q( 3))/R2; 
            
            % from line 10: sd1  = Z/R1;
            %// Z: J =  3 KJ = 1/R1 I =  5
            %s(63) = s(63) - f(70)*q(65)/R2;
            %// R1: J =  5 KJ = -Z*R1^(-2) I =  3,  5
            %s( 5) = s( 5) + f(10)*(2.0*q( 5)*sd1 - q( 3))/R2; 
            % from line 10: sd1  = sin(I)*sin(omega + nu);
            %q(10)  = q(10) + q(18)*sinomnu + sinI*cosomnu*(q(5) + q(46));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % from line 11: sd2  = sd1^2;
            %//sd1: J = 10 KJ = 2.0*sd1 I = 10
            s(10) = s(10) + f(11)*q(10)*2.0;
            
            % from line 12: sd3  = sd1*sd2;
            %//sd1: J = 10 KJ = sd2 I = 11
            s(10) = s(10) + f(12)*q(11);
            %//sd2: J = 11 KJ = sd1 I = 10
            s(11) = s(11) + f(12)*q(10);
            
            % from line 13: sd4  = sd1*sd3;
            %//sd1: J = 10 KJ = sd3 I = 12
            s(10) = s(10) + f(13)*q(12);
            %//sd3: J = 12 KJ = sd1 I = 10
            s(12) = s(12) + f(13)*q(10);
            
            % from line 14: sd5  = sd1*sd4;
            %//sd1: J = 10 KJ = sd4 I = 13
            s(10) = s(10) + f(14)*q(13);
            %//sd4: J = 13 KJ = sd1 I = 10
            s(13) = s(13) + f(14)*q(10);
            
            % from line 15: sd6  = sd1*sd5;
            %//sd1: J = 10 KJ = sd5 I = 14
            s(10) = s(10) + f(15)*q(14);
            %//sd5: J = 13 KJ = sd1 I = 10
            s(14) = s(14) + f(15)*q(10);
            
            % All Second Derivatives Vanish here
            % from line 16: F2   =  1.0      -  3.0*sd2;
            % from line 17: F3   =  3.0*sd1  -  5.0*sd3;
            % from line 18: F4   =  3.0      -  30.*sd2  + 35*sd4;
            % from line 19: F5   =  15.*sd1  -  70.*sd3  + 63*sd5;
            % from line 20: F6   =  5.       - 105.*sd2  + 315*sd4 - 231*sd6;
            
            % from line 21: V1    = 1/R1;          Kepler Potential                 % line 21
            %//R1: J = 15 KJ = -1/R1^2 I = 5
            s( 5) = s( 5) + 2.0*f(21)*q( 5)/R3;
            
            % from line 22: V2      = C2*F2/R2;
            %//F2 J = 16 KJ = C2/R2 I =  4
            s(16) = s(16) - f(22)*q( 4)*C2/R4; 
            %//R2 J =  4 KJ = -C2*F2/R2^2 I =  16,  7
            s( 4) = s( 4) + f(22)*(2.0*q( 4)*V2 - C2*q(16))/R4;
            
            % from line 23: V3      = C3*F3/R3;
            %//F3 J = 17 KJ = C3/R3 I =  6
            s(17) = s(17) - f(23)*q( 6)*C3/R6; 
            %//R3 J =  6 KJ = -C3*F3/R3^2 I = 17,  6
            s( 6) = s( 6) + f(23)*(2.0*q( 6)*V3 - C3*q(17))/R6;
            
            % from line 24: V4      = C4*F4/R4;
            %//F4 J = 18 KJ = C4/R4 I =  7
            s(18) = s(18) - f(24)*q( 7)*C4/R4^2; 
            %//R4 J =  7 KJ = -C4*F4/R4^-2 I = 18,  7
            s( 7) = s( 7) + f(24)*(2.0*q( 7)*V4 - C4*q(18))/R4^2;
            
            % from line 25: V5      = C5*F5/R5;
            %//F5 J = 19 KJ = C5/R5 I =  8
            s(19) = s(19) - f(25)*q( 8)*C5/R5^2; 
            %//R5 J =  8 KJ = -C5*F5/R5^-2 I = 19,  8
            s( 8) = s( 8) + f(25)*(2.0*q( 8)*V5 - C5*q(19))/R5^2;
            
            % from line 26: V6      = C6*F6/R6;
            %//F6 J = 20 KJ = C5/R6 I =  9
            s(20) = s(20) - f(26)*q( 9)*C6/R6^2; 
            %//R6 J =  9 KJ = -C6*F6/R6^-2 I = 20,  9
            s( 9) = s( 9) + f(26)*(2.0*q( 9)*V6 - C6*q(20))/R6^2;

            % All Second Derivative Vanish
            % from line 87: V       = mu*V1;
            % from line 88: D       = V2 + V3 + V4 + V5 + V6;
            
            % %from line 89: R       = mu*D;
            % %q(89) = q(89) + q(88)*mu;
            % % from line 90: Phi     = V + R;
            % %q(90) = q(90) + q(87) + q(89);
            
            % from line 30: Phi     = V*(CK + D);
            %//V J = 27  KJ = (CK + D)  I = 28
            s(27) = s(27) + f(30)*q(28);
            %//D J = 28  KJ = V  I = 87
            s(28) = s(28) + f(30)*q(27); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
            S = s;
        end
    
        function [H_I, Jacob, Hess] = WorkOrder(I)
            
            H_I          = DataList(I);
            Nvar         = 3;
            LineCheck    = I;
            F            = zeros(30,1);
            F(LineCheck) = 1;
            F = backdiff(F);
            Jacob        = [F(1) F(2) F(3)];
            % e,a,...         X    Y    Z
            Hess   = [];
            for j = 1: Nvar
                Q = zeros(30,1);
                S = zeros(30,1);
                %Q(ipoint(j)) = 1.0;
                Q(j) = 1.0;
                Q = fordiff(Q);
                S = secdiff(F, Q, S);
                S = backdiff(S);
                row = [];
                for k = 1 : Nvar
                    row = [row, S(k)];
                end
                Hess  = [Hess; row];
            end
        end
    % function [Jacob, JacobFinite, Hess, HessFinite] = DigitalJustice(I)
    %
    %     Kdelta        = [X, Y, Z];
    %
    %     [H0, Jacob, Hess] = WorkOrder(I);
    %     %Hess
    %
    %     % Do time Increment Separately
    %
    %     epsilon     = 0.000000001;
    %     KeplerDelta = epsilon*eye(3);
    %     [~,~,~,~,~,~] = Kdelta.OrbitDerivatives(time+epsilon, Sensor, Observed);
    %     [HT, JacobT, HessT] = Kdelta.WorkOrder(I);
    %
    %     Kepler1     = Kepler + KeplerDelta(1,:);
    %     Kdelta1     = LagrangePlanetary(Kepler1);
    %     [~,~,~,~,~,~] = Kdelta1.OrbitDerivatives(time, Sensor, Observed);
    %     [H1, Jacob1, Hess1] = Kdelta1.WorkOrder(I);
    %
    %     Kepler2     = Kepler + KeplerDelta(2,:);
    %     Kdelta2     = LagrangePlanetary(Kepler2);
    %     [~,~,~,~,~,~] = Kdelta2.OrbitDerivatives(time, Sensor, Observed);
    %     [H2, Jacob2, Hess2] = Kdelta2.WorkOrder(I);
    %
    %     Kepler3     = Kepler + KeplerDelta(3,:);
    %     Kdelta3     = LagrangePlanetary(Kepler3);
    %     [~,~,~,~,~,~] = Kdelta3.OrbitDerivatives(time, Sensor, Observed);
    %     [H3, Jacob3, Hess3] = Kdelta3.WorkOrder(I);
    %
    %     Kepler4     = Kepler + KeplerDelta(4,:);
    %     Kdelta4     = LagrangePlanetary(Kepler4);
    %     [~,~,~,~,~,~] = Kdelta4.OrbitDerivatives(time, Sensor, Observed);
    %     [H4, Jacob4, Hess4] = Kdelta4.WorkOrder(I);
    %
    %     Kepler5     = Kepler + KeplerDelta(5,:);
    %     Kdelta5     = LagrangePlanetary(Kepler5);
    %     [~,~,~,~,~,~] = Kdelta5.OrbitDerivatives(time, Sensor, Observed);
    %     [H5, Jacob5, Hess5] = Kdelta5.WorkOrder(I);
    %
    %     Kepler6     = Kepler + KeplerDelta(6,:);
    %     Kdelta6     = LagrangePlanetary(Kepler6);
    %     [~,~,~,~,~,~] = Kdelta6.OrbitDerivatives(time, Sensor, Observed);
    %     [H6, Jacob6, Hess6] = Kdelta6.WorkOrder(I);
    %
    %     JacobFinite = [HT-H0,H1-H0, H2-H0, H3-H0, H4-H0, H5-H0, H6-H0]/epsilon;
    %
    %     HessFinite = [(JacobT-Jacob);...
    %                   (Jacob1-Jacob);...
    %                   (Jacob2-Jacob);...
    %                   (Jacob3-Jacob);...
    %                   (Jacob4-Jacob);...
    %                   (Jacob5-Jacob);...
    %                   (Jacob6-Jacob)]/epsilon;
    % end

end