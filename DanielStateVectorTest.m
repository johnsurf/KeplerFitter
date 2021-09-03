clear all;
close all;
clc;

restoredefaultpath;
CurrentPath = pwd;
Utilities = fullfile(CurrentPath,'Utilities')
%genpath(Utilities)
addpath(Utilities)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iUnits = 0
units = Units(iUnits);
twopi = units.twopi;
RadianToDegree = 360/twopi;
% For Canonical Units
TU = units.TU;
DU = units.DU;
VU = units.VU;
AU = units.AU;
mu = units.mu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SVInitial = 10^6 * [-4.050066411885000; ...
                     1.251143732427000; ...
                     6.659743546495000; ...
                     0.000046211881000; ...
                    -0.004678656648000; ...
                     0.001815955236000; ...
                     0.0];

TBeg = SVInitial(7);
tBeg = TBeg/TU;
time = tBeg;
%TFit = 2000.0;
TFit = 3.795452450038239e+03
TEnd = TFit;

tFit = TFit/TU
tEnd = TEnd/TU

timeArray = [TBeg:TFit]/TU;

Rpos = SVInitial(1:3)/DU;
Rdot = SVInitial(4:6)/VU;
[e a Inclination omega Omega Mp Period Jacobian] = kepler(time, Rpos, Rdot, mu)
PeriodSeconds = Period*TU
KeplerBest = [e a Inclination omega Omega Mp]
KeplerObject = KeplerFromECI(time, Rpos, Rdot, units);
InMotion = KeplerObject.Extrapolator(time);
KeplerObject.JacECI_2_Kepler*KeplerObject.JacKepler_2_ECI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
losM     = ones(3,1);
losM     = losM/norm(losM);
InvCovUV = eye(2)*10^(8);
epsilon = 0.001*ones(6,1);
Sensor = [Rpos + 1000.0*epsilon(1:3); zeros(6,1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

KLagrange = LagrangePlanetary(KeplerBest, units);
[Rinit,   Vinit,   ParamList] = OrbitAtTime(KLagrange, tFit, Sensor, losM, InvCovUV);
%options   = odeset('RelTol', 1e-10, 'AbsTol', 1e-13);
options   = odeset('RelTol', 1e-11, 'AbsTol', 1e-14);
dJdt      = perturbations_odefun(KLagrange,tFit,zeros(6,1),units)
dJdtLP    = LagrangePlanetary_odefun(KeplerBest, tFit, zeros(6,1), units)
% Y         = zeros(6,1);
% % Propagate from fit Reference time tFit back to tBeg
% [Tout, Yout] = ode45(@(t,y) perturbations_odefun(KLagrange, t, y, units), [tFit, tBeg], Y, options);
% %Yout(end,:)
% Perturbed = [Tout, Yout];
% Perturbed = sortrows(Perturbed,1);
% Y = Perturbed(end,2:7);
% Propagate from fit Reference time tFit back to tEnd
if tBeg < tEnd
    % tic
    %   [Tout, Yout] = ode45(@(t,y) perturbations_odefun(KLagrange, t, y, units), timeArray, zeros(6,1), options);
    % toc
   tic
     [Tout, Yout] = ode45(@(t,y) LagrangePlanetary_odefun(KeplerBest, t, y, units), timeArray, zeros(6,1), options);
   toc
    %[Tout, Yout] = ode45(@(t,y) perturbations_odefun(KLagrange, t, y, units), [tBeg, tEnd], Y, options);
    % Yout(end,:)
    % Perturbed(end,:) = [];
    % Perturbed = [Perturbed; Tout, Yout];
    PerturbedKepler = [Tout, Yout];
    PerturbedKepler = sortrows(PerturbedKepler,1);
end
KeplerInterpolated  = interp1(PerturbedKepler(:,1),PerturbedKepler(:,1:7),timeArray);
KeplerInstantaneous = KeplerBest + KeplerInterpolated(:,2:7);
SVExtrapolated      = RVfromKepler_vectorized(Tout, KeplerInstantaneous, units);

% % Convert Kepler States to Position and Velocity States
% SVExtrapolated = [];
% for Interval = 1:numel(timeArray)
%     time = timeArray(Interval);
%     % Instantaneouls Kepler State at time = time;
%     %KeplerInstantaneous = KeplerBets + PerturbedKepler(Interval:2:7);
%     KeplerInstantaneous = KeplerBest + KeplerInterpolated(Interval,2:7);
%     %KLagrange = LagrangePlanetary(KeplerInstantaneous, units);
%     %[Rinit,  Vinit, ParamList] = OrbitAtTime(KLagrange, time, Sensor, losM, InvCovUV);
%     %SVExtrapolated = [SVExtrapolated, [Rinit; Vinit]];
%     [dJdt, this] = LagrangePlanetary_odefun(KeplerInstantaneous, time, zeros(6,1), units);
%     Rextrap = [this.rorbitx; this.rorbity; this.rorbitz]; 
%     Vextrap = [this.vorbitx; this.vorbity; this.vorbitz]; 
%     SVExtrapolated = [SVExtrapolated, [Rextrap; Vextrap]];
% end

% %Gravity Comparioson
% format long
% j4MKS = gravityJ4(Rpos*DU)
% [Gravity Phi, InHomogenous, Homogenous] = Escobal(Rpos, mu);
% EcobalMKS = Gravity*DU/TU/TU
% [Rextrap,Vextrap,ParamList,Jacobian,Hessian,GravityCan] = OrbitDerivatives(KLagrange, time, Sensor, losM, InvCovUV);
% LP_GravityMKS = GravityCan*DU/TU/TU
% [J6 Jperturb Jhom] = J6Gravity(Rpos);
% J6GravityMKS = J6'*DU/TU/TU

% Caution --> the following vector is supposed to have 6 rows and 1 column
SV6Vector = [Rpos; Rdot]
dJdt = SVIntegration_odefun(SV6Vector);
Y         = SV6Vector;
tic
    % [Tout, Yout] = ode45(@(t,y) SVIntegration_odefun(y), [tBeg, tEnd], Y, options);
    [Tout, Yout] = ode45(@(t,y) SVIntegration_odefun(y), timeArray, Y, options);
    PerturbedSV = [Tout, Yout];
    PerturbedSV = sortrows(PerturbedSV,1);
    %SV6Interpolated = PerturbedSV(:,2:7);
    SV6Interpolated  = interp1(PerturbedSV(:,1),PerturbedSV(:,2:7),timeArray);
    RposVector       = SV6Interpolated(:,1:3);
    RdotVector       = SV6Interpolated(:,4:6);
    Distance         = sum(RposVector.*RposVector,2);
    Istop            = find(Distance<1.0);
    OrbitIndex       = 1:Istop(1)-1;
    [e, a, Inclination, omega, Omega, Mp, hVec, Latitude, Jacobian] = kepler_vectorized(timeArray, RposVector, RdotVector, mu);
    KeplerElements   = [e a Inclination omega Omega Mp] - KeplerBest;
    %DifferenceKepler = [DifferenceKepler; KeplerElements - PerturbedKepler(Interval,2:7)];
    DifferenceKepler = KeplerElements - KeplerInterpolated(:,2:7);
    KeplerReconstructed = [timeArray', KeplerElements];
    IndexGood        = OrbitIndex;
    TimeActive       = timeArray(IndexGood); 
    hVec             = hVec';
toc

% KeplerReconstructed = [];
% DifferenceKepler    = [];
% Latitude            = [];
% Distance            = [];
% OrbitIndex          = [];
% hVec                = [];
% EarthIntersection   = true;
% TimeActive          = 0.0;
% for Interval = 1:numel(timeArray)
%     time = timeArray(Interval);
%     Rpos = SV6Interpolated(Interval,1:3)';
%     RTransverseSQ = Rpos(1)*Rpos(1) + Rpos(2)*Rpos(2);
%     RadialDistance = sqrt(RTransverseSQ + Rpos(3)*Rpos(3));
%     % Upon dist encounter with the surface of the Earth -- turn it off
%     if RadialDistance < 1.0
%         EarthIntersection = false;
%     end
%     % Only include up to the first encounter with the surface of the Earth
%     if EarthIntersection
%         OrbitIndex = [OrbitIndex, Interval];
%         TimeActive = time*TU;
%     end
%     Rtransverse = sqrt(RTransverseSQ);
%     Distance = [Distance, RadialDistance];
%     Theta = RadianToDegree*atan2(Rpos(3), Rtransverse);
%     Latitude = [Latitude, Theta];
%     Rdot = SV6Interpolated(Interval,4:6)';
%     hVec = [hVec, cross(Rpos,Rdot)];
%     [e a Inclination omega Omega Mp Period Jacobian] = kepler(time, Rpos, Rdot, mu);
%     KeplerElements  = [e a Inclination omega Omega Mp] - KeplerBest;
%     %DifferenceKepler = [DifferenceKepler; KeplerElements - PerturbedKepler(Interval,2:7)];
%     DifferenceKepler = [DifferenceKepler; KeplerElements - KeplerInterpolated(Interval,2:7)];
%     KeplerReconstructed = [KeplerReconstructed; time, KeplerElements];
% end
% 
% IndexGood = OrbitIndex;
% TimeActive

hVecMean = mean(hVec(:,IndexGood),2);

[U,S,V] = svd((hVec(:,IndexGood) - hVecMean)*(hVec(:,IndexGood) - hVecMean)');
EigenHVec = sort(eig(S),'descend');

% Begin Figures: 
figure 
hx = hVec(1,IndexGood);
hy = hVec(2,IndexGood);
hz = hVec(3,IndexGood);
h1 = hVec(:,IndexGood)'*U(:,1);
h2 = hVec(:,IndexGood)'*U(:,2);
h3 = hVec(:,IndexGood)'*U(:,3);
colormat = rand(1,3);
scatter3(hx, hy, hx,'MarkerEdgeColor','k','MarkerFaceColor',colormat);
%scatter3(h1, h2, h3,'MarkerEdgeColor','k','MarkerFaceColor',colormat);
xlabel('\it h1')
ylabel('\it h2')
zlabel('\it h3')
%daspect([1 1 1])

figure
%plot(PerturbedKepler(:,1), PerturbedKepler(:,2:7),'linewidth',2)
%plot(timeArray, KeplerInterpolated(:,1:6),'linewidth',2)
%plot(timeArray, Distance,'linewidth',2)
plot(Latitude, Distance, 'linewidth', 2)
hold on
thresh = 1.0;
%plot([timeArray(1), timeArray(end)], [thresh, thresh],'-', 'linewidth', 2)
plot([min(Latitude), max(Latitude)], [thresh, thresh],'-', 'linewidth', 2)
hold on
scatter(Latitude(IndexGood(1)), Distance(IndexGood(1)), 'Marker', 'x', 'MarkerEdgeColor', 'r', 'MarkerFaceColor','r','linewidth', 2);
xlabel('Latitude [degrees]')
ylabel('Distance from the Center of the Earth [DU]')
grid on; 

figure
%plot(PerturbedKepler(:,1), PerturbedKepler(:,2:7),'linewidth',2)
%plot(timeArray, KeplerInterpolated(:,1:6),'linewidth',2)
%plot(timeArray(indexGood), Distance(IndexGood),'linewidth',2)
plot(Latitude(IndexGood), Distance(IndexGood), 'linewidth', 2)
hold on
scatter(Latitude(IndexGood(1)), Distance(IndexGood(1)), 'Marker', 'x', 'MarkerEdgeColor', 'r', 'MarkerFaceColor','r','linewidth', 2);
xlabel('Latitude [degrees]')
ylabel('Altitude [DU]')
grid on; 

figure
%plot(PerturbedKepler(:,1), PerturbedKepler(:,2:7),'linewidth',2)
%plot(timeArray, KeplerInterpolated(:,1:6),'linewidth',2)
plot(Latitude(IndexGood), PerturbedKepler(IndexGood,2:7), 'linewidth', 2)
legend('\delta e','\delta a','\delta I','\delta\omega','\delta\Omega','\delta Mp')
xlabel('Kepler Elements [time] - Initial Kepler Elements')
ylabel('Latitude [degrees]')
grid on; 

figure
plot(Latitude(IndexGood), KeplerReconstructed(IndexGood,2:7),'linewidth',2)
%plot(KeplerReconstructed(:,1), KeplerReconstructed(:,2:7),'linewidth',2)
legend('\delta e','\delta a','\delta I','\delta\omega','\delta\Omega','\delta Mp')
xlabel('Kepler Elements [time] - Initial Kepler Elements')
ylabel('Latitude [degrees]')
grid on; 

figure
plot(Latitude(IndexGood), DifferenceKepler(IndexGood,1:6),'linewidth',2)
%plot(KeplerReconstructed(:,1), KeplerReconstructed(:,2:7),'linewidth',2)
legend('\delta e','\delta a','\delta I','\delta\omega','\delta\Omega','\delta Mp')
%legend('e','a','I','\omega','\Omega','Mp')
ylabel('Difference in Kepler Elements Calculations')
xlabel('Latitude [degrees]')
grid on; 

% Plot Positions
figure
plot(timeArray(IndexGood), SVExtrapolated(1, IndexGood),'linewidth',2)
hold on
plot(timeArray(IndexGood), PerturbedSV(IndexGood,2),'linewidth',2)
ylabel(' X component [DU]')
xlabel(' Time [TU]')

figure
plot(timeArray(IndexGood), SVExtrapolated(2, IndexGood),'linewidth',2)
hold on
plot(timeArray(IndexGood), PerturbedSV(IndexGood,3),'linewidth',2)
ylabel(' Y component [DU]')
xlabel(' Time [TU]')

figure
plot(timeArray(IndexGood), SVExtrapolated(3, IndexGood),'linewidth',2)
hold on
plot(timeArray(IndexGood), PerturbedSV(IndexGood,4),'linewidth',2)
ylabel(' Z component [DU]')
xlabel(' Time [TU]')

figure
plot(timeArray(IndexGood), SVExtrapolated(4, IndexGood),'linewidth',2)
hold on
plot(timeArray(IndexGood), PerturbedSV(IndexGood,5),'linewidth',2)
ylabel(' Vx component [DU]')
xlabel(' Time [TU]')

figure
plot(timeArray(IndexGood), SVExtrapolated(5, IndexGood),'linewidth',2)
hold on
plot(timeArray(IndexGood), PerturbedSV(IndexGood,6),'linewidth',2)
ylabel(' Vy component [DU]')
xlabel(' Time [TU]')

figure
plot(timeArray(IndexGood), SVExtrapolated(6, IndexGood),'linewidth',2)
hold on
plot(timeArray(IndexGood), PerturbedSV(IndexGood,7),'linewidth',2)
ylabel(' Vz component [DU]')
xlabel(' Time [TU]')

% Plot Residuals Kepler - Position and Velocity
%SV6  = PerturbedSV(:,2:7)';
SV6  = SV6Interpolated';
Differences = SVExtrapolated - SV6;
diffPos = sqrt(sum(Differences(1:3, IndexGood).*Differences(1:3,IndexGood)));
figure
plot(timeArray(IndexGood), diffPos*DU,'linewidth',2)
ylabel(' Residuals of Position Differences [m]')
xlabel(' Time [TU]')

diffVel = sqrt(sum(Differences(4:6, IndexGood).*Differences(4:6,IndexGood)));
figure
plot(timeArray(IndexGood), diffVel*VU,'linewidth',2)
ylabel(' Residuals of Velocity Differences [m/s]')
xlabel(' Time [TU]')

figure
plot(timeArray(IndexGood), Differences(1,IndexGood)*DU,'linewidth',2)
ylabel(' X component Differences [m]')
xlabel(' Time [TU]')

figure
plot(timeArray(IndexGood), Differences(2,IndexGood)*DU,'linewidth',2)
ylabel(' Y component Differences [m]')
xlabel(' Time [TU]')

figure
plot(timeArray(IndexGood), Differences(3,IndexGood)*DU,'linewidth',2)
ylabel(' Z component Differences [m]')
xlabel(' Time [TU]')

figure
plot(timeArray(IndexGood), Differences(4,IndexGood)*VU,'linewidth',2)
ylabel(' Vx component Differences [m/s]')
xlabel(' Time [TU]')

figure
plot(timeArray(IndexGood), Differences(5,IndexGood)*VU,'linewidth',2)
ylabel(' Vy component Differences [m/s]')
xlabel(' Time [TU]')

figure
plot(timeArray(IndexGood), Differences(6,IndexGood)*VU,'linewidth',2)
ylabel(' Vz component Differences [m/s]')
xlabel(' Time [TU]')


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dJdt = SVIntegration_odefun(SV6Vector)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Test out the Lagrange Planetary Equations:
    [J6 Jperturb Jhom] = J6Gravity(SV6Vector(1:3));
    dJdt = [SV6Vector(4:6); J6'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%
% % All these State Vectors have been Extrapolated to the latest data point
% %StateVectorECI = state_vector
% %StateVectorECI
% %StateVectorECI = SVestimated1
% SVBest
% ChiSqMin
% %[KeplerFit1, CovFit1, TFit, SV7Fit1] = LeastSquares(track_data, tFirst, 1, 0, zeros(7,1) )
% [KeplerFit1, CovFit1, TFit, SV7Fit1] = LeastSquaresAllSatellites(track_data, TFit, 1, 0, InvCovUV, units, InitialState)