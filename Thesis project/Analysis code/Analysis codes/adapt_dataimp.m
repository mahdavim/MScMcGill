function [edfStruct1, cat_var, lbl_names, datestr, gx, gy, gstx, gsty, genx, geny, st_times, en_times] = ...
    adapt_dataimp(path, cue_pos, eye, cue_positions)

% Outputs:

% gx and gy are in dva
%gstx, gsty, genx, gny gaze start and end points of saccades

% Cat_var is returned in the order that you provide the cue_positions
addpath(genpath(path{1})) % Main codes


edfStruct1 = edfmex(path{2}); % EDF file
psych_dat = readtable(path{3}); %psychopy csv file

encoded_positions = {'1','2', '3', '4', '5', '6', '7', '8'}; % Code id for que positions. CHANGE TO DEGREES


%---- import cue positions and convert them to 1:8 categorical values
tmp_q = psych_dat.q_position(~cellfun(@isempty,psych_dat.q_position)); % Que positions.only choose non-empty cells
cat_var = double(categorical(tmp_q,cue_positions, encoded_positions )); 
lbl_names = unique(tmp_q);



gx = double(edfStruct1.FSAMPLE.gx(eye,1:end)); % For the right eye(2) or left (1)
gx(gx == max(gx)) = nan;
gy = double(edfStruct1.FSAMPLE.gy(eye,1:end));
gy(gy  == max(gy)) = nan;

% Align gx and gy to center
gx = gx - 1920/2;
gy = 1080/2 - gy;

rx = double(edfStruct1.FSAMPLE.rx); % For the right eye(2) or left (1)
ry = double(edfStruct1.FSAMPLE.ry);

gx = gx ./ rx;
gy = gy ./ ry;

% Extract the date from CSV file name (The third entry in pth)
datePattern = '\d{4}-\d{2}-\d{2}';

datestr = regexp(path{3}, datePattern, 'match', 'once');

% -----------Extract data for saccades
 event_type = double([edfStruct1.FEVENT.type]);% Event type of code strings STARTBLINK=3, ENDBLINK=4, STARTSACC=5, ENDSACC=6,
%STARTFIX=7, ENDFIX=8, FIXUPDATE=9)type event type STARTFIX=7, ENDFIX=8, FIXUPDATE=9
 
gstx = (double([edfStruct1.FEVENT.gstx])- 1920/2) ./ (double([edfStruct1.FEVENT.supd_x])); % In DVA
gstx = gstx([edfStruct1.FEVENT.eye] == 1 & event_type == 6); % Filter for right eye events and saccade ends only
gsty = (1080/2 - double([edfStruct1.FEVENT.gsty])) ./ (double([edfStruct1.FEVENT.supd_y]));
gsty = gsty([edfStruct1.FEVENT.eye] == 1 & event_type == 6);


genx = (double([edfStruct1.FEVENT.genx]) - 1920/2) ./ (double([edfStruct1.FEVENT.eupd_x]));
genx = genx([edfStruct1.FEVENT.eye] == 1 & event_type == 6);
geny = (1080/2 - double([edfStruct1.FEVENT.geny])) ./ (double([edfStruct1.FEVENT.eupd_y]));
geny = geny([edfStruct1.FEVENT.eye] == 1 & event_type == 6);


st_times = double([edfStruct1.FEVENT.sttime]);
st_times = st_times([edfStruct1.FEVENT.eye] == 1 & event_type == 6);
en_times = double([edfStruct1.FEVENT.entime]);
en_times = en_times([edfStruct1.FEVENT.eye] == 1 & event_type == 6);


end