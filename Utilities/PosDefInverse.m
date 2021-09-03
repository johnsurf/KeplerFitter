function [ PDInv, SqrtPDInv, SqrtPD ] = PosDefInverse( PD )
%% Inverse of Positive Definite Matrix (i.e., a covariance)

    DiagMinVal = 1.0e-10;
    [L, D, Perm] = ldl(PD);

    SqrtD    = D;
    SqrtDInv = D;

    N = size(D,1);

    for k = 1:N
        if (D(k,k) < DiagMinVal) % Repair any negaitve (or too small)diag values
            D(k,k) = DiagMinVal;
        end

        SqrtD(k,k) = sqrt(D(k,k));
        SqrtDinv(k,k) = 1/SqrtD(k,k);
    end

    LInv = Inverse_Lower_Triangle(L);
    
    %L * LInv

    SqrtPDInv = Perm * LInv' * SqrtDInv;

    PDInv     = SqrtPDInv * SqrtPDInv';

    SqrtPD = Perm * L * SqrtD; % Replaces RCHol

end