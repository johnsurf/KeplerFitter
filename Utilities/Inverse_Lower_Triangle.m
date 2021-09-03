function [ LInv ] = Inverse_Lower_Triangle ( L );
    % Compute inverse of unit lower triangular matrix
    N    = size(L,1);
    LInv = eye(N);

    for jj = 1:N-1
        for ii = jj+1:N
            s = L(ii,jj);
            for k = (jj+1):(ii-1)
                s = s + L(ii,k) * LInv(k,jj);
            end
            LInv(ii,jj) = -s;
        end
    end
    
end