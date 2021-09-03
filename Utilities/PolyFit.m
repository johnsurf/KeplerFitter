    function [state_vector, fit, Iorder, ChiSq] = PolyFit(Object, track_data, tFit, Iorder, iType)
        ChiSq  = [];
        [Matrix, Minc]   =  BuildMatrix(Object, track_data, tFit, Iorder, iType);
        MatrixRank  = rank(Minc);
        while MatrixRank < size(Matrix,1)
            Iorder      = Iorder - 1;
            [Matrix, Minc]   =  BuildMatrix(Object, track_data, tFit, Iorder, iType);
            MatrixRank  = rank(Minc);
        end
        if MatrixRank < size(Matrix,1)
            disp(' Polynomial Fit Fails: Matrix rank is too low')
        else
            A      = Matrix(:,end);
            %CovFit = inv(Minc)
            fit    = Minc\A;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Reshape into matrix form
            fit     = reshape(fit,3,Iorder);
        end

        [state_vector, ChiSq] = PolyExtrapArray(Object,track_data,fit, tFit);
        
    end