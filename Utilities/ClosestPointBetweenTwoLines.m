    function [SVFit] = ClosestPointBetweenTwoLines(track_data,tFit)
        % Method 1: 
        % First do a simple position at closest approach of two line of
        % sight vectors
        % make initial guess using first two sightings -- assume same
        % point is described by both sightings
        r1 = track_data(1,3:5)';
        u1 = track_data(1,6:8)';
        t1 = track_data(1,2);

        r2 = track_data(2,3:5)';
        u2 = track_data(2,6:8)';
        t2 = track_data(2,2);

        r3 = track_data(3,3:5)';
        u3 = track_data(3,6:8)';
        t3 = track_data(3,2);

        CoefMatrix = [u1'*u1,  -u1'*u2;...
            -u2'*u1,  u2'*u2];
        bPoint     = [(r2 - r1)'*u1;...
            (r1 - r2)'*u2];
        tValues    = CoefMatrix\bPoint;
        l1         = r1 + tValues(1)*u1;
        l2         = r2 + tValues(2)*u2;
        v1         = (l2 - l1)/(t2 - t1);

        CoefMatrix = [u2'*u2,  -u2'*u3;...
            -u3'*u2,  u3'*u3];
        bPoint     = [(r3 - r2)'*u2;...
            (r2 - r3)'*u3];
        tValues    = CoefMatrix\bPoint;
        l2a        = r2 + tValues(1)*u2;
        l3         = r3 + tValues(2)*u3;
        v2         = (l3 - l2a)/(t3 - t2);

        SV1        = [ (l2 + l2)/2.0; v1; (t1+t2)/2.0];
        SV2        = [ (l3 + l2)/2.0; v2; (t3+t2)/2.0];
        toTime     = (t1 + t2)/2.0;
        SV2Prop    = propagateECI(SV2, toTime);            
        SVStart    = (SV1 + SV2Prop)/2.0;
        SVFit    = propagateECI(SVStart, tFit);
    end