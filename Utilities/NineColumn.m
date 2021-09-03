function [row] = NineColumn(Orbit, Thief, timeTU, ISensor) 

    orbitparams = wgs84Constants;
    TU = orbitparams.TU;
    DU = orbitparams.DU;
    VU = orbitparams.VU;
    AU = orbitparams.AU;

    [rSensor, vSensor] = extrapolate(Orbit, timeTU);
    [rThief, vThief]   = extrapolate(Thief, timeTU);
    los   = rThief - rSensor;
    range = norm(los);
    los   = los/range;
    row   = [ISensor, timeTU*TU, rSensor'*DU, los', range*DU, rThief'*DU, vThief*VU];
end