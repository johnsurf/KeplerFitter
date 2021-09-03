function [SVextrapolated, CovSV, fit, CovFit, ChiSq] = BuildPolyKeplerMatrix(SVInit, tExtrap, TPolyFit, track_data, Iorder, iType, units, InvCovUV)
    rows = size(track_data,1);
    mu = units.mu;

    % Modfied Polynomial Fit Tools
    MatRow = Iorder*3;
    MatCol = MatRow + 3 + 1;
    Column = eye(3);
    X1     = [  Column, zeros(3),   zeros(3)];
    X2     = [zeros(3),   Column,   zeros(3)];
    X3     = [zeros(3), zeros(3), 2.0*Column];
    X      = [X1; X2; X3];
    
    SVextrapolated   = SVInit;    
    Tinit            = SVInit(10);
    SVInitCanonical  = [SVInit(1:3)/units.DU; SVInit(4:6)/units.VU; SVInit(7:9)/units.AU; Tinit/units.TU];
    % Initialize to the incoming state: 
    betaInit         = inv(X)*SVInitCanonical(1:9);
    beta             = [betaInit; zeros(MatRow-9,1)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %lambda = zeros(3,1);
    lambda = [1; 1; 1];

    for NR = 1:10

        Minc   = zeros(MatRow+3, MatRow+3);
        Matrix = zeros(MatRow+3, MatCol);
        Gradient = zeros(MatRow+3,1);
        CovInv = eye(3);
        InvFitCov = zeros(MatRow,MatRow);
        ChiSq  = 0.0;

        for ii = 1:rows
            trackID   = track_data(ii,1);
            tinc      = (track_data(ii,2) - TPolyFit)/units.TU;

            %tinc = 1.0

            Sii  = track_data(ii,3:5)'/units.DU;   % Satellite
            Uii  = track_data(ii,6:8)';   % Line of Site
            Uii  = Uii/norm(Uii);
            Uprojection = eye(3) - Uii*Uii'; % N.B. This matrix is symmetric
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %InvCov      = (10^8)*Uprojection*eye(3)*Uprojection;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            thetaM      = acos(Uii(3));
            phiM        = atan2(Uii(2), Uii(1));
            cosThetaM   = cos(thetaM);
            sinThetaM   = sin(thetaM);
            cosPhiM     = cos(phiM);
            sinPhiM     = sin(phiM);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Mdelta      =[cosThetaM*cosPhiM, cosThetaM*sinPhiM, -sinThetaM;...
                                   -sinPhiM,           cosPhiM,       0.0 ];
            InvCov      = Mdelta'*InvCovUV*Mdelta;
            % %InvCov   = Uprojection*eye(3)*CovInv*Uprojection;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Polynomial Interpolation for the Sensor
            X1          = [];   % Incidence/Design Matrix
            X2          = [zeros(3)];
            X3          = [zeros(3), zeros(3)];
            Column1     = eye(3);
            Column2     = eye(3);
            Column3     = eye(3);
            CL          = InvCov*Sii;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Polynomials in time
            for j      = 1:Iorder
                deriv1st    = j-1;
                deriv2nd    = (j-2)*deriv1st;
                X1      = [X1, Column1];
                Column1 = Column1*tinc;
                if j > 1
                    X2      = [X2, deriv1st*Column2];
                    Column2 = Column2*tinc;
                end
                if j > 2
                    X3     = [X3, deriv2nd*Column3];
                    Column3      = Column3.*tinc;
                end
            end
            %IncMat = zeros(MatRow+3, MatRow+3);
            %IncMat    = X1'*InvCov*X1;

            Y = X1*beta;
            Displacement = Uprojection*Sii;
            DistSq = Displacement'*Displacement;
            InvFitCov = InvFitCov + (X1'*InvCov*X1)/DistSq;

            dChiSq  = 2.0*X1'*InvCov*(Y - Sii);
            ChiSq   = ChiSq + 0.5*(Y - Sii)'*InvCov*(Y - Sii);

            Fbeta     = Y'*Y;
            dFbeta    = 2.0*X1'*Y;
            d2Fbeta   = 2.0*X1'*X1;

            Lagrange    = X3*beta + mu*Fbeta^(-1.5)*X1*beta;
            dLagrange   = X3' + mu*Fbeta^(-1.5)*X1' - 1.5*mu*Fbeta^(-2.5)*dFbeta*beta'*X1';
            Constraint  = Lagrange'*lambda;
            dConstraint  = dLagrange*lambda;

            Gradient  = Gradient + [dChiSq + dConstraint; Lagrange];
            d2ChiSq   = 2.0*X1'*InvCov*X1;
            d2Constraint = mu*( 3.75*(lambda'*Y)*Fbeta^(-3.5)*dFbeta*dFbeta'...
                -3*Fbeta^(-2.5)*dFbeta*(lambda'*X1)...
                -1.5*(lambda'*Y)*Fbeta^(-2.5)*d2Fbeta);

            IncMat = [d2ChiSq + d2Constraint, dLagrange;...
                dLagrange',                zeros(3)];

            %Yerror    = X'*InvCov*X;
            Minc      = Minc   + IncMat;
            Matrix    = Matrix + [IncMat, Gradient];
        end
        delta = Minc\Gradient;
        beta = beta - delta(1:Iorder*3);
        lambda = lambda - delta(Iorder*3+1:end);
        %Gradient
        %Constraint
    end
    MatrixRank  = rank(Minc);
    % CovFit = inv(Minc);
    CovFit = inv(InvFitCov);
    %CovFit  = PosDefInverse(InvFitCov); 
    fit    = reshape(beta,3,Iorder);
    
    % Compute at TPolyFit time: tinc = 0.0
    %tinc        = 0.0;
    tinc      = (tExtrap - TPolyFit)/units.TU;
    X1          = [];   % Incidence/Design Matrix
    X2          = [zeros(3)];
    X3          = [zeros(3), zeros(3)];
    Column1     = eye(3);
    Column2     = eye(3);
    Column3     = eye(3);
    
    % Polynomials in time
    for j      = 1:Iorder
        deriv1st    = j-1;
        deriv2nd    = (j-2)*deriv1st;
        X1      = [X1, Column1];
        Column1 = Column1*tinc;
        if j > 1
            X2      = [X2, deriv1st*Column2];
            Column2 = Column2*tinc;
        end
        if j > 2
            X3     = [X3, deriv2nd*Column3];
            Column3      = Column3.*tinc;
        end
    end
    X      = [X1; X2; X3];
    CovSV = X*CovFit(1:MatRow,1:MatRow)*X';
    CovSV = CovSV./units.COV9Ref;
    SVextrapolated = X*beta;
    SVextrapolated = [units.DU*SVextrapolated(1:3); units.VU*SVextrapolated(4:6); units.AU*SVextrapolated(7:9); tExtrap];
    
end