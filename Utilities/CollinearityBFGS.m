function [KeplerFit, ChiSq, CovFit] = CollinearityBFGS(track_data, KeplerInit, InvCovUV, units)

rows       = size(track_data,1);

%InvCovUV

% For Choice of Units
twopi = units.twopi;
TU = units.TU;
DU = units.DU;
%VU = units.VU;
%AU = units.AU;
mu = units.mu;
%sqrtmu = units.sqrtmu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TU = units.TU;
% DU = units.DU;
% VU = units.VU;
% AU = units.AU;

delta      = zeros(6,1);
delp       = zeros(6,1);

TangentArray   = [];
ZhatArray      = [];
DelZhatArray   = [];

ChiSqArray     = [];
DelChiSqArray  = [];

DChiSqArray    = [];
DiffChiSqArray = [];

PhiFundamental = eye(6);
InvCovMeas   = eye(2);
InvCovM      = eye(2);
deltaPar     = zeros(6,1);
epsilon      = 0.000000001*ones(6,1);
Kepler       = KeplerInit

angle         = Kepler(4)
while(angle <   0.0)
    angle = angle + twopi;
end
while(angle > twopi)
    angle = angle - twopi;
end
Kepler(4) = angle

angle         = Kepler(5)
while(angle <   0.0)
    angle = angle + twopi;
end
while(angle > twopi)
    angle = angle - twopi;
end
Kepler(5) = angle


Mind  = 6
oparm = zeros(Mind,1);
Jmatrix   = zeros(6);

% Stephen's Paper:
%       W  = HU
%       U  = JW   where J is the update for H^(-1)
%  [J_k + VV^T]W = U
%   (VV^T)W  = U - J_k*W
%    V       = [U - J_k*W]*(V^T*W)^(-1)
%    W^T VV^T W    = W^T[U - J_k*W]*[W^T VV^T W]^(-1)*[U - J_k*W]^TW
%    A = W^T VV^T W
%    B = W^T[U - J_k]W

Umatrix   = zeros(6);
delp      = zeros(6,1);

for NR = 1:30
    
    % nfree    = 0;
    % oparm(1) = true;
    % oparm(2) = true;
    % oparm(3) = true;
    % oparm(4) = true;
    % oparm(5) = true;
    % oparm(6) = true;
    
    % deltaPar = zeros(6,1);
    % delta     = zeros(6,1);
    
    IndexRows = find(oparm == 1);
    IndexCols = [IndexRows; Mind+1];
    
    % for ipar=1:Mind
    %     if oparm(ipar) == true
    %         nfree = nfree+1;
    %     end
    % end
    
    KeplerObj = LagrangePlanetary(Kepler, units);
    ChiSq     = 0.0
    %delp      = ones(6,1);
    %delp      = delp/norm(delp);
    gradChiSq = zeros(6,1);
    %Umatrix   = zeros(6,6);
    
    for ii =1:rows
        tRecord = track_data(ii,2)/TU;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %losM    = state_vector_losArray(1:3,ii);
        losM     = track_data(ii,6:8)';
        losM     = losM/norm(losM);   % los
        % rangeM   = norm(losM);
        % thetaM   = acos(losM(3)/rangeM);
        % phiM     = atan2(losM(2), losM(1));
        % Observed = [thetaM; phiM];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %Sat Position at measured position
        Rsensor  = [track_data(ii,3:5)']/DU;
        Sensor   = [Rsensor; zeros(6,1)];
        
        % %Sat Position from Polynomial at tFit
        % Sensor = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU]
        
        % %Sat Position from Orbit Fit to Kepler elements:
        % [Rsensor, Vsensor, ParamList, Jacobian, Hessian, GravityCan] = KeplerSat.OrbitDerivatives(tRecord, Sensor, Observed);
        % Sensor   = [Rsensor; Vsensor; zeros(3,1)];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Rinit,Vinit,ParamList] = OrbitAtTime(KeplerObj, tRecord, Sensor, losM, InvCovUV);
        % LineCheck    = 131 Collinearity
        ChiSq        = ChiSq + abs(KeplerObj.DataList(131));
        LineCheck    = 131;
        F            = zeros(133,1);
        F(LineCheck) = 1;
        F = backdiff(KeplerObj, F, LineCheck);
        Jacob        = [F(2) F(3) F(4) F(5) F(6) F(7)];
        gradChiSq    = gradChiSq + Jacob';
        %Umatrix      = Umatrix  + eye(6);
        %Umatrix(:,1)      = Umatrix(:,1) + Jacob';
        % if NR == 1
        %     Umatrix      = Umatrix + [Jacob', 0.001*ones(6,1)];
        % else
        %     Umatrix      = Umatrix + [Jacob', delp];
        % end
        % Q = zeros(133,1);
        % S = zeros(133,1);
        % %Q(ipoint(j)) = 1.0;
        % Q(2:7) = Jacob';
        % Q = fordiff(KeplerObj, Q, LineCheck);
        % S = secdiff(KeplerObj, F, Q, S, LineCheck);
        % S = backdiff(KeplerObj, S, LineCheck);
        % Wmatrix = Wmatrix + S(2:7);
    end
    
    %Umatrix(:,1) = Umatrix(:,1)/norm(Umatrix(:,1));
    Umatrix = eye(6);
    Umatrix(:,1) = gradChiSq/norm(gradChiSq);
    if NR == 1
        Umatrix(:,2) = [0;1;0;0;0;0];
    else
        Umatrix(:,2) =-delp/norm(delp);
    end
    
    %Wmatrix = zeros(6,2);
    Wmatrix = zeros(6);
    for ii =1:rows
        tRecord = track_data(ii,2)/TU;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %losM    = state_vector_losArray(1:3,ii);
        losM     = track_data(ii,6:8)';
        losM     = losM/norm(losM);   % los
        % rangeM   = norm(losM);
        % thetaM   = acos(losM(3)/rangeM);
        % phiM     = atan2(losM(2), losM(1));
        % Observed = [thetaM; phiM];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %Sat Position at measured position
        Rsensor  = [track_data(ii,3:5)']/DU;
        Sensor   = [Rsensor; zeros(6,1)];
        
        % %Sat Position from Polynomial at tFit
        % Sensor = [state_vector_sat(1:3)/DU; state_vector_sat(4:6)/VU; state_vector_sat(7:9)/AU]
        
        % %Sat Position from Orbit Fit to Kepler elements:
        % [Rsensor, Vsensor, ParamList, Jacobian, Hessian, GravityCan] = KeplerSat.OrbitDerivatives(tRecord, Sensor, Observed);
        % Sensor   = [Rsensor; Vsensor; zeros(3,1)];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Rinit,Vinit,ParamList] = OrbitAtTime(KeplerObj, tRecord, Sensor, losM, InvCovUV);
        % LineCheck    = 131 Collinearity
        %H_I          = KeplerObj.DataList(131);
        LineCheck    = 131;
        F            = zeros(133,1);
        F(LineCheck) = 1;
        F = backdiff(KeplerObj, F, LineCheck);
        
        Q = zeros(133,1);
        S = zeros(133,1);
        %Q(ipoint(j)) = 1.0;
        Q(2:7) = Umatrix(:,1);
        Q = fordiff(KeplerObj, Q, LineCheck);
        S = secdiff(KeplerObj, F, Q, S, LineCheck);
        S = backdiff(KeplerObj, S, LineCheck);
        W1 = S(2:7);
        
        Q = zeros(133,1);
        S = zeros(133,1);
        %Q(ipoint(j)) = 1.0;
        Q(2:7) = Umatrix(:,2);
        Q = fordiff(KeplerObj, Q, LineCheck);
        S = secdiff(KeplerObj, F, Q, S, LineCheck);
        S = backdiff(KeplerObj, S, LineCheck);
        W2 = S(2:7);
        
        Q = zeros(133,1);
        S = zeros(133,1);
        %Q(ipoint(j)) = 1.0;
        Q(2:7) = Umatrix(:,3);
        Q = fordiff(KeplerObj, Q, LineCheck);
        S = secdiff(KeplerObj, F, Q, S, LineCheck);
        S = backdiff(KeplerObj, S, LineCheck);
        W3 = S(2:7);
        
        
        Q = zeros(133,1);
        S = zeros(133,1);
        %Q(ipoint(j)) = 1.0;
        Q(2:7) = Umatrix(:,4);
        Q = fordiff(KeplerObj, Q, LineCheck);
        S = secdiff(KeplerObj, F, Q, S, LineCheck);
        S = backdiff(KeplerObj, S, LineCheck);
        W4 = S(2:7);
        
        
        Q = zeros(133,1);
        S = zeros(133,1);
        %Q(ipoint(j)) = 1.0;
        Q(2:7) = Umatrix(:,5);
        Q = fordiff(KeplerObj, Q, LineCheck);
        S = secdiff(KeplerObj, F, Q, S, LineCheck);
        S = backdiff(KeplerObj, S, LineCheck);
        W5 = S(2:7);
        
        Q = zeros(133,1);
        S = zeros(133,1);
        %Q(ipoint(j)) = 1.0;
        Q(2:7) = Umatrix(:,6);
        Q = fordiff(KeplerObj, Q, LineCheck);
        S = secdiff(KeplerObj, F, Q, S, LineCheck);
        S = backdiff(KeplerObj, S, LineCheck);
        W6 = S(2:7);
        
        Wmatrix = Wmatrix + [W1, W2, W3, W4, W5, W6];
        %Wmatrix = Wmatrix + [W1'; W2'; W3'; W4'; W5'; W6'];
    end
    
    if NR < 3
        Wmatrix(:,3:6) = [];
        Umatrix(:,3:6) = [];
    elseif NR < 6
        Wmatrix(:,4:6) = [];
        Umatrix(:,4:6) = [];
    elseif NR < 9
        Wmatrix(:,5:6) = [];
        Umatrix(:,5:6) = [];
    end
    
    WU    = Wmatrix'*Umatrix;
    % WJW   = Wmatrix'*Jmatrix*Wmatrix;
    % Amatrix = WU - WJW;
    % Nmatrix = Umatrix - Jmatrix*Wmatrix
    % Jmatrix = Jmatrix + Nmatrix*inv(Amatrix)*Nmatrix'
    % delp = -Jmatrix*Umatrix(:,1)
    
    % Try Direction Set  HU = W
    UHU  = Umatrix'*Wmatrix;
    delp = -Umatrix*inv(UHU)*Umatrix'*gradChiSq
    
    
    % if det(WU) > 0 && det(Amatrix) > 0
    %     Jmatrix = Jmatrix + Nmatrix*inv(Amatrix)*Nmatrix'
    % end
    %
    % if det(WU) > 0 && det(Amatrix) < 0
    %     Jmatrix = zeros(6);
    % end
    
    %Wmatrix = Wmatrix/norm(Wmatrix);
    
    %[deltaM, JacFac, Hess] = WorkOrder(KeplerObj, 131);
    %alphaSq    = Umatrix'*(Wmatrix - Hess(2:7,2:7)*Umatrix)
    %alphaSq    = Umatrix'*(Wmatrix - Hessian*Umatrix)
    %alpha      = sqrt(alphaSq)
    %Update1    = (Wmatrix - Hessian*Umatrix)/alpha
    %denom      = 1 + Update1'*Jmatrix*Update1
    %InvUpdate  = Jmatrix*Update1*Update1'*Jmatrix/denom
    %Jmatrix = Jmatrix - InvUpdate
    %Hessian    = inv(Jmatrix)
    
    ChiSqArray   = [ChiSqArray, log(ChiSq)];
    %delta        = -(HessianSyst)\Umatrix
    %deltaPar     = -hessNR\betaNR
    %TangentArray = [TangentArray, log(norm(betaNR))]
    %Kepler       =  Kepler + deltaPar'
    %Kepler       =  Kepler + delta';
    
    % AugmentedMatrix = [Hessian,-Umatrix];
    % AugmentedMatrix = AugmentedMatrix(IndexRows, IndexCols);
    % hitzero   = 0;
    % opzero = 0.0000000000001;
    % disp("Entering NRMin ");
    % pinc = NRMin(AugmentedMatrix, nfree, opzero, hitzero);
    % if(hitzero == 1)
    %     disp("  hit operational zero it NRMin  ")
    % end
    % ii = 1;
    % stepSize = 0.0;
    % for i=1:Mind
    %     if oparm(i)
    %         delp(i)   = pinc(ii);
    %         ii = ii + 1;
    %         stepSize = stepSize + delp(i)*delp(i);
    %     end
    % end
    delp
    % if norm(delp)<10^(-12)
    %     break
    % end
    % if NR < 3
    %     delp = -0.15*Umatrix/rows;
    % end
    Kepler = Kepler + delp'
    
    TangentArray = [TangentArray, log(norm(gradChiSq))]
    angle         = Kepler(4);
    while(angle <   0.0)
        angle = angle + twopi;
    end
    while(angle > twopi)
        angle = angle - twopi;
    end
    Kepler(4) = angle;
    
    angle         = Kepler(5);
    while(angle <   0.0)
        angle = angle + twopi;
    end
    while(angle > twopi)
        angle = angle - twopi;
    end
    Kepler(5) = angle;
    
end

KeplerFit = Kepler;
CovFit    = Umatrix*inv(UHU)*Umatrix';
% Adjust for Spec Values for cov_UV = 10^(-8)*eye*2(
CovFit    = 10^(-8)*CovFit;

figure
plot([1:numel(ChiSqArray)], ChiSqArray)
title('log(\chi^2) per iteration')
ylabel('log(\chi^2)')
xlabel('Newton-Raphson Iteration Number')

figure
plot([1:numel(TangentArray)], TangentArray)
title('Log norm of the Derivatives of \chi^2 per iteration')
ylabel('Log norm \chi^2 Derivatives')
xlabel('Newton-Raphson Iteration Number')
end
