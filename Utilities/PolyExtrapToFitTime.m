    function [state_vector] = PolyExtrapToFitTime(fit, tFit)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %elements         = zeros(3,rows);
        %derivElements    = zeros(3,rows);
        %secDerivElements = zeros(3,rows);
        elements      = zeros(3,1);
        derivElements = zeros(3,1);
        secDerivElements = zeros(3,1);
        state_vector_sat = [];
        %tIndex = find(track_data(:,2) == tFit);
        %for ii = 1:rows
            ii = 1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Tinc     = 0.0;
            %Tinc     = track_data(ii,2) - tFit;
            %timeFac  = ones(1,rows);
            %derivFac = ones(1,rows);
            timeFac     = 1;
            derivFac    = 1;
            secDerivFac = 1;
            for i = 1:size(fit,2)
                %elements = elements + fit(:,i)*Tinc.^(i-1);
                deriv1st    = i-1;
                deriv2nd    = (i-2)*deriv1st;
                elements(:,ii) = elements(:,ii) + fit(:,i)*timeFac;
                timeFac        = timeFac.*Tinc;
                if i > 1
                    derivElements(:,ii) = derivElements(:,ii) + fit(:,i)*deriv1st*derivFac;
                    derivFac      = derivFac.*Tinc;
                end
                if i > 2
                    secDerivElements(:,ii) = secDerivElements(:,ii) + fit(:,i)*deriv2nd*secDerivFac;
                    secDerivFac      = secDerivFac.*Tinc;
                end
            end
        %end
        %%%%  Satellite 9-State Vector at the median time
        %state_vector_sat = [elements/DU; derivElements/VU; secDerivElements/AU]
        state_vector = [elements; derivElements; secDerivElements; tFit];
    end