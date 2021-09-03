function [ W T ] = Earth_Rotation_Matrix( time )
% rotation matrices for ECI - ECR conversions

    %Figure the delta longitude due to earth rotation
    We = 7.292115E-5; %earth rotation rate (rad/s)

    del_lon = time * We;

    %Euler Rotation Matrix (CCW about z-axis)
    cosl = cos(del_lon);
    sinl = sin(del_lon);

    T = [ cosl  sinl  0;
         -sinl  cosl  0;
          0     0     1  ];


    W =  [ 0   We 0;
          -We  0  0;
           0   0  0  ];
     
end