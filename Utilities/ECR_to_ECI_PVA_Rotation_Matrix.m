function [ M ] = ECR_to_ECI_PVA_Rotation_Matrix( time )
%Coordinate transformation for covariance matrix from ECR to ECI 

    [ W T ] = Earth_Rotation_Matrix( time );

    T = T';
    Z = zeros(3,3);

    M = [ T      Z      Z;
         -W*T    T      Z;
          W*W*T -2*W*T  T ];
  
end