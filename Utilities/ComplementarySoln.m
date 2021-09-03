    function [expAt] = ComplementarySoln(A)
             dim = size(A,2);
             Inc = eye(dim);
             expAt = zeros(6);
             for i = 1:20
                 expAt = expAt + Inc; 
                 Inc   = Inc*A/i;
             end
    end