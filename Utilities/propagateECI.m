function SV7Prop = propagateECI(SV7, toTime)
    %{
    Ballistic Object Propagation
    ---------------------------------------------------------------------------
    SV7     = beginning ECI state vector column [x y z vx vy vz validity_time]
    toTime  = validity time of requested propagated state
    SV7Prop = propagated ECI state vector column [x y z vx vy vz validity_time]
    ---------------------------------------------------------------------------
    %}
    
	% Define variables
	f1 = double(zeros(6, 1));
	f2 = double(zeros(6, 1));
	f3 = double(zeros(6, 1));
	f4 = double(zeros(6, 1));
	
    integTimeStepMax = 5.0;
    
    propTime = toTime - SV7(7,1);
    SV7Prop = SV7;
    
    nIterations = ceil( abs(propTime) / integTimeStepMax );
    if nIterations < 10 
        nIterations = 10;
    end
    
    if( nIterations == 0 ); return; end;
    
    delt = propTime / nIterations;
    
    for i = 1 : nIterations
        
        %4th Order Runge Kutta Integration
        SV7Old = SV7Prop;
        
        f1(1:6,1) = [SV7Old(4:6); gravityJ4(SV7Old)];
        %f1(1:6,1) = [SV7Old(4:6); EscobalMKS(SV7Old)];
        
        temp = SV7Old(1:6) + 0.5*delt*f1;
        f2(1:6,1) = [temp(4:6); gravityJ4(temp)];
        %f2(1:6,1) = [temp(4:6); EscobalMKS(temp)];
        
        temp = SV7Old(1:6) + 0.5*delt*f2;
        f3(1:6,1) = [temp(4:6); gravityJ4(temp)];
        %f3(1:6,1) = [temp(4:6); EscobalMKS(temp)];
        
        temp = SV7Old(1:6) + delt*f3;
        f4(1:6,1) = [temp(4:6); gravityJ4(temp)];
        %f4(1:6,1) = [temp(4:6); EscobalMKS(temp)];
        
        SV7Prop(1:6,1) = SV7Old(1:6) + (delt/6)*(f1 + 2*f2 + 2*f3 + f4);
        SV7Prop(7,1)   = SV7Old(7)   + delt;
        
    end
    
end