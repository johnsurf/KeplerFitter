    function [Matrix, Minc] = BuildMatrix(Object, track_data, tFit, Iorder, iType)
        rows       = size(track_data,1);

        % Modfied Polynomial Fit Tools
        MatRow = Iorder*3;
        MatCol = MatRow + 1;
        IncMat = zeros(MatRow, MatRow);
        Minc   = zeros(MatRow, MatRow);
        Matrix = zeros(MatRow, MatCol);

        for ii = 1:rows
            trackID   = track_data(ii,1);
            tinc      = track_data(ii,2) - tFit;
            if Object == 1
                states_used  = track_data(ii,3:5)';   % Satellite
            else
                states_used  = track_data(ii,6:8)';   % Line of Site
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Polynomial Interpolation for the Sensor
            X           = [];   % Incidence/Design Matrix
            Y           = [];
            Column      = eye(3);
            InvCov      = eye(3);
            CL          = InvCov*states_used;
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
            Yerror    = X'*InvCov*X;
            Minc      = Minc   + IncMat;
            Matrix    = Matrix + [IncMat, Y];
        end

        MatrixRank  = rank(Minc);

    end