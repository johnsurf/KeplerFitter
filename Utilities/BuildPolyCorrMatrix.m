    function [Matrix, Minc] = BuildPolyCorrMatrix(Object, track_data, tFit, Iorder, iType, units, InvCovUV)
    % Called by PolyCorrFit.m 
        rows = size(track_data,1);

        % Modfied Polynomial Fit Tools
        MatRow = Iorder*3;
        MatCol = MatRow + 1;
        IncMat = zeros(MatRow, MatRow);
        Minc   = zeros(MatRow, MatRow);
        Matrix = zeros(MatRow, MatCol);
        CovInv = eye(3);

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
                    Y      = [Y; CL];
                    Column = Column*tinc;
                    CL     = CL*tinc;
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
            IncMat    = X'*InvCov*X;
            %Yerror    = X'*InvCov*X;
            Minc      = Minc   + IncMat;
            Matrix    = Matrix + [IncMat, Y];
        end

        MatrixRank  = rank(Minc);

    end