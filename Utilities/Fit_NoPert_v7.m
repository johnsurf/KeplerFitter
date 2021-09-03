clc, clear all
restoredefaultpath
file_path  = mfilename('fullpath')
in_dir  = fileparts(file_path)
out_dir = strcat(in_dir,'\Converted')
orbitparams = wgs84Constants;
TU = orbitparams.TU;
DU = orbitparams.DU;
VU = orbitparams.VU;
AU = orbitparams.AU;
%mu = orbitparams.mu;
mu = orbitparams.Grav;
twopi  = 2.0*pi;         % 2*pi
Set = 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Input_filepath = 'InputDirectory';
filename       = 'InputFile'
out_dir        = fullfile(in_dir,'Converted')
fid            = fopen(fullfile(Input_filepath, filename),'r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[message] = null_los_message_v1();
trackID   = 1

if ~isfolder(out_dir)
    mkdir(out_dir);
end
%%
count     = 0;
SVInitial = zeros(9,1);
iFirst    = 0; 
while ~feof(fid)
    tline = fgetl(fid);
    line = textscan(tline,'%s','Delimiter',' ','MultipleDelimsAsOne',1);
    newline = line{1,1}';
    % If header row, start new block read
    if contains(tline,'Track ID') | contains(tline,'Sequence')
        % if isempty(message(1).trackID)
        %     count = 1;
        % end
        % % % Header Row Data Fields
        % message(count).trackID   = str2double(newline{1,4});
        % %message(count).blockID   = str2double(newline{1,6});
        % %message(count).availTime = str2double(newline{1,9});
        if count ~= 0
            disp(' analyze block of LOS data ')
            % Convert messages into Our format
            rows       = size(message,2);
            %track_data = zeros(rows,121);
            track_data = zeros(rows,8);
            for ii = 1:rows
                time                 = message(ii).timeUTC;
                if Set == 1
                    state                = message(ii).satVec
                elseif Set == 2
                    % Adding in conversion from km to m if necessary!!!
                    state                = message(ii).satVec';
                    state                = 1e3*state;
                    los                  = message(ii).losVec';
                    % no point scaling the Cov9x9 matrix since it is "truth" and we're
                    % using unit matrices for the covariances.
                end
                % Convert ECEF to ECI at some arbitrary time:
                [W T] = Earth_Rotation_Matrix(time);
                T     = T';
                %new_state         = state;
                %new_los           = los;
                new_state           = T*state;
                new_los             = T*los;
                new_los             = new_los/norm(new_los);
                %C.data(ii,2:4)       = new_9state(1:3);
                %C.data(ii,5:7)       = new_9state(4:6);
                %C.data(ii,8:10)      = new_9state(7:9);
                %Cov9x9_ECI           = ECR_to_ECI_PVA_cov(Cov9x9,time);
                
                %track_data(ii,1)     = trackID;
                track_data(ii,1)     = message(ii).satNum;
                track_data(ii,2)     = time;
                % Adding in conversion from km to m
                track_data(ii,3:5)   = new_state(1:3);
                track_data(ii,6:8)   = new_los(1:3);
                %track_data(ii,9:11)  = new_9state(7:9);
                %track_data(ii,12:92) = reshape(Cov9x9_ECI,[1 81]);
                       
            end
            timeSorted = sortrows(track_data,2)
            if iFirst == 0
                tFirst   = timeSorted(1,2)
                iFirst = -1;
            end
            track_data(:,2) = track_data(:,2) - tFirst
            
            %[Zfit, tFit] = LeastSquares(track_data, tFirst, 402, 0 )
            
            disp(' reset for next block ')
            count = 0
        end
        
    continue
    else
        count = count + 1;
        % Message Block, row 2
        message(count).timeUTC     = str2double(newline{1,1});
        message(count).pad1        = str2double(newline{1,2});
        message(count).satNum      = str2double(newline{1,3});
        message(count).losVec(1,1) = str2double(newline{1,4});
        message(count).losVec(1,2) = str2double(newline{1,5});
        message(count).losVec(1,3) = str2double(newline{1,6});
        message(count).satVec(1,1) = str2double(newline{1,7});
        message(count).satVec(1,2) = str2double(newline{1,8});
        message(count).satVec(1,3) = str2double(newline{1,9});
        message(count).radioInt  = str2double(newline{1,10});
        message(count).band      = str2double(newline{1,11});
        message(count).pad2      = str2double(newline{1,12});
    end

end
Cov9x9 = eye(9);

%%
% Convert messages into Our format
rows       = size(message,2);
M          = unique(track_data,'rows');
track_data = M;

saveName = sprintf('%s_los_data.txt',strrep(filename,'.txt',''));
save(fullfile(out_dir,saveName),'track_data','-ascii','-double');
save(fullfile(out_dir,'track_data.txt'),'track_data','-ascii','-double');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check Out Polynomial Fits   Sensor Orbits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%plot3(track_data(:,3), track_data(:,4),track_data(:,4))
Satellites = unique(track_data(:,1))
ChiSqSat      = [];
for iSat=1:numel(Satellites)
    Index      = find(track_data(:,1) == Satellites(iSat));
    TrackOrbit = track_data(Index,:);
    plot3(TrackOrbit(:,3),TrackOrbit(:,4),TrackOrbit(:,5),'linewidth',2);
    hold on
    rows       = size(TrackOrbit,1);
    Iorder     = 5;
    %Iorder     = 1;
    Object     = 1; % for Sensor
    iType      = 1;  % iType is hard-wired at the moment
    iMid       = floor(numel(track_data(:,2))/2);
    tFit       = median(TrackOrbit(:,2))
    %tFit       = 0.;
    [state_vector_fit, fit, Iorder, ChiSq] = PolyFit(Object, TrackOrbit,tFit,Iorder,iType);
    plot3(state_vector_fit(1,:),state_vector_fit(2,:),state_vector_fit(3,:),'linewidth',2);
    hold on
    Diff = state_vector_fit(1:3,:) - TrackOrbit(:,3:5)';
    ChiSqSat = [ChiSqSat, sum(sum(Diff.*Diff))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check Out Polynomial Fits   Line of Sight Polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%plot3(track_data(:,3), track_data(:,4),track_data(:,4))
Satellites = unique(track_data(:,1))
ChiSqSat      = [];
for iSat=1:numel(Satellites)
    Index      = find(track_data(:,1) == Satellites(iSat));
    TrackOrbit = track_data(Index,:);
    plot3(TrackOrbit(:,6),TrackOrbit(:,7),TrackOrbit(:,8),'linewidth',2);
    hold on
    rows       = size(TrackOrbit,1);
    Iorder     = 4;
    %Iorder     = 1;
    Object     = 2; % for los
    iType      = 1;  % iType is hard-wired at the moment
    iMid       = floor(numel(track_data(:,2))/2);
    tFit       = median(TrackOrbit(:,2))
    %tFit       = 0.;
    [state_vector_fit, fitSat, Iorder, ChiSq] = PolyFit(Object, TrackOrbit,tFit,Iorder,iType);
    plot3(state_vector_fit(1,:),state_vector_fit(2,:),state_vector_fit(3,:),'linewidth',2);
    hold on
    Diff = state_vector_fit(1:3,:) - TrackOrbit(:,6:8)';
    ChiSqSat = [ChiSqSat, sum(sum(Diff.*Diff))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[KeplerFit, CovFit, tFit] = LeastSquaresAll(track_data, tFirst, 0, 0, zeros(7,1))
[KeplerFit402, CovFit402, tFit, SV7Fit] = LeastSquares(track_data, tFirst, 402, 0, zeros(7,1) )
%[KeplerFit401, CovFit401, tFit, SV7Fit] = LeastSquares(track_data, tFirst, 401, 0, zeros(7,1) )
%[KeplerFit401, CovFit401, tFit] = LeastSquares(track_data, tFirst, 401, 1, SV7Fit )
InvComb     = inv(CovFit401 + CovFit402)
KeplerF     = CovFit402*InvComb*KeplerFit401 + CovFit401*InvComb*KeplerFit402
CovFit1      = CovFit401*InvComb*CovFit402
CovFit2      = CovFit402*InvComb*CovFit401
