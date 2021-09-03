    function [F] = losDerivs(los,f)
        % On the fly Backwards Differentiation: 
        % Prediction: 
        % lx = los(1)    Line 1
        % ly = los(2)    Line 2
        % lz = los(3)    Line 3
        %range   = norm(los)
        %theta   = acos(los(3)/range)
        %phi     = atan2(los(2),los(1))
        rangeSq  = los(1)^2 + los(2)^2 + los(3)^2;   % line 4  
        range    = sqrt(rangeSq);                    % line 5
        u        = los(3)/range;                     % line 6
        theta   = acos(u);                           % line 7
        phi     = atan2(los(2),los(1));              % line 8
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        denom   = los(1)^2 + los(2)^2 ;
        % from line 8: phi     = atan2(los(2),los(1))
        f(1)   = f(1) - f(8)*los(2)/denom; 
        f(2)   = f(2) + f(8)*los(1)/denom ;
        % from line 7: theta   = acos(u)
        f(6)   = f(6) - f(7)/sqrt(1 - u^2); 
        %from line 6: u        = los(3)/range
        f(3)   = f(3) + f(6)/range; 
        f(5)   = f(5) - f(6)*u/range;
        % from line 5: range    = sqrt(rangeSq)
        f(4)   = f(4) + 0.5*f(5)/range;
        % from line 4: rangeSq  = los(1)^2 + los(2)^2 + los(3)^2
        f(1)  = f(1) + 2.0*f(4)*los(1); 
        f(2)  = f(2) + 2.0*f(4)*los(2); 
        f(3)  = f(3) + 2.0*f(4)*los(3);
        F = f; 
    end