function [Iorder, fit, CovFit, ChiSq] = PolyCorrFit(Object, track_data, tFit, Iorder, iType, units, InvCovUV)
CovFit  = [];
ChiSq   = 0.0;
[Matrix, Minc]    =  BuildPolyCorrMatrix(Object, track_data, tFit, Iorder, iType, units, InvCovUV);
MatrixRank  = rank(Minc);
while MatrixRank < size(Matrix,1)
    Iorder         = Iorder - 1;
    [Matrix, Minc] = BuildPolyCorrMatrix(Object, track_data, tFit, Iorder, iType, units, InvCovUV);
    MatrixRank  = rank(Minc);
end
if MatrixRank < size(Matrix,1)
    disp(' Polynomial Fit Fails: Matrix rank is too low')
else
    A      = Matrix(:,end);
    %CovFit = inv(Minc)
    fit    = Minc\A;
    CovFit = inv(Minc);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rows   = size(track_data,1);
    for ii = 1:rows
        trackID   = track_data(ii,1);
        tinc      = (track_data(ii,2) - tFit)/units.TU;
        Sii  = track_data(ii,3:5)'/units.DU;   % Satellite
        Uii  = track_data(ii,6:8)';   % Line of Site
        Uii  = Uii/norm(Uii);
        Uprojection = eye(3) - Uii*Uii'; % N.B. This matrix is symmetric
        thetaM      = acos(Uii(3));
        phiM        = atan2(Uii(2), Uii(1));
        cosThetaM   = cos(thetaM);
        sinThetaM   = sin(thetaM);
        cosPhiM     = cos(phiM);
        sinPhiM     = sin(phiM);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Mdelta      =[cosThetaM*cosPhiM, cosThetaM*sinPhiM, -sinThetaM;...
            -sinPhiM,           cosPhiM,       0.0 ];
        CovInv      = Mdelta'*InvCovUV*Mdelta;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Polynomial Interpolation for the Sensor
        X           = [];   % Incidence/Design Matrix
        Y           = [];
        Column      = eye(3);
        InvCov      = Uprojection*eye(3)*CovInv*Uprojection;
        CL          = InvCov*Sii;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Polynomials in time
        if iType == 1
            for j      = 1:Iorder
                X      = [X, Column];
                %Y      = [Y; CL];
                Column = Column*tinc;
                %CL     = CL*tinc;
            end
        end
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Chebyshev Polynomials of the First Kind
        % if iType == 2
        %     tscaled = tinc/tMag;
        %     T       = [1, tscaled];
        %     for j   = 1:Iorder
        %         if j >= 3
        %             T  = [T, (2*tscaled*T(j-1) - T(j-2))];
        %         end
        %         Column = Column*T(j);
        %         CL     = CL*T(j);
        %         X      = [X, Column];
        %         Y      = [Y; CL];
        %     end
        % end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %IncMat    = X'*InvCov*X;
        %Yerror    = X'*InvCov*X;
        %IncMat    = X'*InvCov*X;
        YValue    = X*fit;
        Delta     = YValue - Sii;
        ChiSq = ChiSq + 0.5*Delta'*InvCov*Delta;
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Reshape into matrix form
    end
    fit     = reshape(fit,3,Iorder);
    
    %[state_vector, ChiSq] = PolyExtrapArray(1,track_data,fit, tFit);
    %state_vector(4:6) = state_vector(4:6)/units.TU;
    %state_vector(7:9) = state_vector(7:9)/units.TU^2;
    
end