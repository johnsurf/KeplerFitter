function MyUnits = Units(I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MyUnits.twopi  = 2.0*pi;         % 2*pi
    orbitparams = wgs84Constants;
    if I == 0
        % For Canonical Units
        MyUnits.TU = orbitparams.TU;
        MyUnits.DU = orbitparams.DU;
        MyUnits.VU = orbitparams.VU;
        MyUnits.AU = orbitparams.AU;
        MyUnits.mu = orbitparams.mu;
        MyUnits.sqrtmu = sqrt(MyUnits.mu);
    else
        % For MKS units:
        MyUnits.TU = 1.0;
        MyUnits.DU = 1.0;
        MyUnits.VU = 1.0;
        MyUnits.AU = 1.0;
        MyUnits.mu = orbitparams.Grav;
        MyUnits.sqrtmu = sqrt(MyUnits.mu);
    end
    TRef_s    = MyUnits.TU;
    LRef_m    = MyUnits.DU;
    VRef_mps  = MyUnits.VU;
    ARef_mpss = MyUnits.AU;
    % Units Conversion for vectors and Covariances Bill's method
    % Converson from Internal Scaled Units (Canonical) back to MKS
    MyUnits.JacRef6  = diag([LRef_m,LRef_m,LRef_m,VRef_mps,VRef_mps,VRef_mps,]);
    MyUnits.JacRef9  = diag([LRef_m,LRef_m,LRef_m,VRef_mps,VRef_mps,VRef_mps,...
                            ARef_mpss,ARef_mpss,ARef_mpss]);
    MyUnits.X_PVRef  = [1/MyUnits.DU,1/MyUnits.DU,1/MyUnits.DU,1/MyUnits.VU,1/MyUnits.VU,1/MyUnits.VU];
    MyUnits.X_PVARef = [1/MyUnits.DU,1/MyUnits.DU,1/MyUnits.DU,1/MyUnits.VU,1/MyUnits.VU,...
                        1/MyUnits.VU,1/MyUnits.AU,1/MyUnits.AU,1/MyUnits.AU];
    MyUnits.COV6Ref  = MyUnits.X_PVRef' * MyUnits.X_PVRef;
    MyUnits.COV9Ref  = MyUnits.X_PVARef'* MyUnits.X_PVARef;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%