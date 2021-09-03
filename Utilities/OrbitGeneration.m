twopi = 2.0*pi;

orbitparams = wgs84Constants;
TU = orbitparams.TU;
DU = orbitparams.DU;
VU = orbitparams.VU;
AU = orbitparams.AU;

%mu     =  1.0;         % Canonical Graviational Parameter = GM
%Rearth =  1.0;
%J2     =  1082.6269e-6;   %J2 pertubative constant
%J3     =      -2.30e-6;
%J4     =      -2.12e-6;
%J5     =      -0.20e-6;
%J6     =       1.00e-6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input / Output
%Input_filepath = 'InputDirectory';
%filename       = 'ThreeGoldenEyes'
%filename       = 'FourGoldenEyes'
filename       = 'FiveGoldenEyes'
in_dir         = pwd
out_dir        = fullfile(in_dir,'los_files')
%fid            = fopen(fullfile(Input_filepath, filename),'r');
if ~isfolder(out_dir)
    mkdir(out_dir);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Master Clock Start Time: 
%tEpoch          = 200/TU;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time   = [];
iFirst = 0;
for tSeconds = 600:900
    time = [time; tSeconds];    
end
TimeSize  = size(time);
timeSorted = sortrows(time);
if iFirst == 0
    tFirst   = timeSorted(1);
end
iMid           = floor(numel(time)/2);
TFit           = time(iMid)
%TFit            = median(time); 
tEpoch          = TFit/TU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iKepler = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a Sensor Orbit called Golden Eye1
[GoldenEye1, Kepler] = GenOrbit(tEpoch);
iKepler = [iKepler; Kepler];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a Sensor Orbit called Golden Eye 2
[GoldenEye2, Kepler] = GenOrbit(tEpoch);
iKepler = [iKepler; Kepler];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a Sensor Orbit called Golden Eye 3
[GoldenEye3, Kepler] = GenOrbit(tEpoch);
iKepler    = [iKepler; Kepler];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a Sensor Orbit called Golden Eye 4
[GoldenEye4, Kepler] = GenOrbit(tEpoch);
iKepler    = [iKepler; Kepler];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a Sensor Orbit called Golden Eye 4
[GoldenEye5, Kepler] = GenOrbit(tEpoch);
iKepler    = [iKepler; Kepler];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a Target called "The Thief"
[Thief, KeplerThief] = GenOrbit(tEpoch);
iKepler = [iKepler; KeplerThief];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Let us Digitize the Theif's Orbit's LOS vector with respect to the GoldenEye Sensors  

icnt = 0;
imax = 0;
imin = 0;

time = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ISensor    = 1:10;   % Allows for Sensor Satellite Numbering. 
track_data = [];
iFirst     = 0;

%for tSeconds = 600:PeriodGoldenEye1/4000.0:900
%for tSeconds = 0:PeriodGoldenEye1/4000.0:400
for tSeconds = 600:900
    
    %timeArray = tEpoch*TU + 300*rand(3,1);
    %timeArray = 300*rand(3,1);
    %time      = timeArray(1); 
    
    time = [time, tSeconds];
    icnt = icnt + 1;
    ISat = 0;
    
    ISat      = ISat + 1;
    timeMeas = tSeconds;
    %timeMeas = timeArray(1);
    %timeMeas = tSeconds + 2.0*(0.5 - rand(1))*50.;
    
    [row] = NineColumn(GoldenEye1, Thief, timeMeas/TU, ISensor(ISat)); 
    % %[rSensor, vSensor] = extrapolate(GoldenEye1, timeMeas/TU);
    % [rSensor, vSensor] = extrapolate(GoldenEye1, timeMeas/TU);
    % [rThief, vThief]   = extrapolate(Thief,   timeMeas/TU);
    % los   = rThief - rSensor;
    % range = norm(los);
    % los   = los/range;
    % row   = [ISensor(ISat), timeMeas, rSensor*DU, los, range*DU];
    track_data = [track_data; row];
    
    ISat     = ISat     + 1;
    timeMeas = tSeconds + 5;
    %timeMeas = timeArray(2);
    %timeMeas = tSeconds + 2.0*(0.5 - rand(1))*50.;
    
    [row] = NineColumn(GoldenEye2, Thief, timeMeas/TU, ISensor(ISat)); 
    track_data = [track_data; row];
    
    ISat     = ISat     + 1;
    timeMeas = tSeconds + 10;
    %timeMeas = timeArray(3);
    %timeMeas = tSeconds + 2.0*(0.5 - rand(1))*50.;    
    [row] = NineColumn(GoldenEye3, Thief, timeMeas/TU, ISensor(ISat)); 
    track_data = [track_data; row];

    ISat     = ISat     + 1;
    timeMeas = tSeconds + 10;
    %timeMeas = timeArray(3);
    %timeMeas = tSeconds + 2.0*(0.5 - rand(1))*50.;    
    [row] = NineColumn(GoldenEye4, Thief, timeMeas/TU, ISensor(ISat));
    track_data = [track_data; row];
    
    ISat     = ISat     + 1;
    timeMeas = tSeconds - 10;
    %timeMeas = timeArray(3);
    %timeMeas = tSeconds + 2.0*(0.5 - rand(1))*50.;    
    [row] = NineColumn(GoldenEye5, Thief, timeMeas/TU, ISensor(ISat));
    track_data = [track_data; row];
    
end

TimeArray = size(time);
timeSorted = sortrows(track_data,2);
if iFirst == 0
    tFirst   = timeSorted(1,2);
    iFirst = -1;
end
track_data      = timeSorted;
%track_data(:,2) = track_data(:,2) - tFirst;
%iMid           = floor(numel(track_data(:,2))/2);
%TFit           = tEpoch*TU
%TFit           = track_data(iMid,2)
%TFit            = median(track_data(:,2))

saveName = sprintf('%s_los_data.txt',strrep(filename,'.txt',''));
save(fullfile(out_dir,saveName),'track_data','-ascii','-double');
%save(fullfile(out_dir,'track_data.txt'),'track_data','-ascii','-double');
