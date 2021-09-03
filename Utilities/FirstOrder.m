    function [KeplerDot] = FirstOrder(t,K)             
        Kepler     = LagrangePlanetary(K);
        %J6Check    = Kepler.J6Gravity(time, EscobalGravity*TU^2/DU);
        J6Check    = Kepler.J6Gravity(t);
        KeplerDot  = Kepler.FirstOrder(t);
    end