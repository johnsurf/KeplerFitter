function  InitialState = InitialStateFit(Iorder, TPolyFit, fit, CovFit, SVInitial)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    InitialState.Iorder      = Iorder;    
    InitialState.TPolyFit    = TPolyFit;
    InitialState.fit         = fit; 
    InitialState.CovFit      = CovFit;
    
    %InitialState.TExtrap     = TExtrap;    
    %InitialState.TEnd        = TEnd;
    InitialState.SV10Initial = SVInitial;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end