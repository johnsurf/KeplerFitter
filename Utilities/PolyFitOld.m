    function [state_vector, fit] = PolyFitOld(Object, track_data, tFit, Iorder, iType)
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
            Covariance  = eye(3);
            CL          = states_used;
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
            IncMat    = X'*Covariance*X;
            Minc      = Minc   + IncMat;
            Matrix    = Matrix + [IncMat, Y];
        end

        MatrixRank  = rank(Minc);
        if MatrixRank < size(Matrix,1)
            Iorder      = Iorder - 1;
            MatrixRank  = rank(Minc);
            disp(' Matrix rank is too low')
        else
            A      = Matrix(:,end);
            fit    = Minc\A;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Reshape into matrix form
            fit     = reshape(fit,3,Iorder);
        end
        [state_vector, ChiSq] = PolyExtrapArray(Object,track_data,fit, tFit);
    end