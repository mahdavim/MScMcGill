function [edfStruct1, datestr, gx, gy, gstx, gsty, genx,...
    geny,pupil_a, st_times, en_times, tmp_dvnt, tmp_adapt, tmp_resp, tmp_side, tmp_q, tmp_rt] = ...
    adapt_dataimp_main(path, eye, typ)
% To do: Add adaptor train, deviant, deviant side, reaction time
%-----Inputs
% typ: experiment type (e.g., nwm, nwm2)


%-------- Outputs:

% tmp_dvnt: Chosen deviant slope
% tmp_adapt: Chosen  adaptro train
% tmp_resp: Key response (Ascedning is 1 descending is 0 no response is -1)
% tmp_side: Deviant side for experiment 2 (For experiment 1 its nan) (0 is
% left and 1 is right same as tmp_q)
% tmp_q: Cue position for working memory experiments (Nan for non working
% memory) (34, 2.5)' => 1 (Right on screen), (-34, 2.5)' => 0 (Left on screen)

% gx and gy are in dva
%gstx, gsty, genx, gny gaze start and end points of saccades
% st_times, en_times gaze start and end times

% typ = input("What is the epxeriment type? ");

% Cat_var is returned in the order that you provide the cue_positions
addpath(genpath(path{1})) % Main codes


edfStruct1 = edfmex(path{2}); % EDF file
psych_dat = readtable(path{3}); %psychopy csv file

%CHANGE
% encoded_positions = {'1','2'}; % Code id for que positions. CHANGE TO DEGREES


%---- import cue positions and convert them to 1:8 categorical values
% if cue_pos > 0
%     tmp_q = psych_dat.q_position(~cellfun(@isempty,psych_dat.q_position)); % Que positions.only choose non-empty cells
%     cat_var = double(categorical(tmp_q,cue_positions, encoded_positions ));
%     lbl_names = unique(tmp_q);
% else
%     cat_var = NaN;
%     lbl_names=NaN;
% end




gx = double(edfStruct1.FSAMPLE.gx(eye,1:end)); % For the right eye(2) or left (1)
gx(gx == max(gx)) = nan;
gy = double(edfStruct1.FSAMPLE.gy(eye,1:end));
gy(gy  == max(gy)) = nan;

% Align gx and gy to center
gx = gx - 1920/2;
gy = 1080/2 - gy;

rx = double(edfStruct1.FSAMPLE.rx); % For the right eye(2) or left (1)
ry = double(edfStruct1.FSAMPLE.ry);

gx = gx ./ rx; %Convert to DVA
gy = gy ./ ry;
pupil_a = double(edfStruct1.FSAMPLE.pa(eye, :)); % @Hadi

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


%----------  adaptor train, deviant, deviant side, reaction time, cue position ----

% if strcmpi(typ, 'nwm') | strcmpi(typ, 'nwm2')

if strcmpi(typ, 'nwm') | strcmpi(typ, 'wm')
    dvnts = [-8, -6, -4, -2,0, 2, 4, 6, 8];


elseif strcmpi(typ, 'nwm2') | strcmpi(typ, 'wm2')
    dvnts = [-12, -9, -6, -3,0, 3, 6, 9, 12];

end

%------- Deviant type-------
tmp_dvnt = psych_dat.selected_dvnt;% Selected deviant
tmp_dvnt = cell2mat(cellfun(@str2num, tmp_dvnt, 'UniformOutput', false));
dvs = unique(tmp_dvnt);
% Enter true deviant values
for kk=1:length(dvs)
    tmp_dvnt(tmp_dvnt==dvs(kk)) = dvnts(kk);
end
%---------- Adaptor train ------
tmp_adapt = psych_dat.adptr_train_list; % Adaptor train
tmp_adapt = cell2mat(cellfun(@str2num, tmp_adapt, 'UniformOutput', false));

% ------------- Key response ------------
tmp_resp = psych_dat.key_index; % Key response
tmp_resp = tmp_resp(~cellfun( @isempty, tmp_resp)); % Remove the empty character cell
tmp_resp = cellfun(@(x) ...
    (strcmp(x, '[[''right'']]') * 1) + ...   % 'right' -> 1
    (strcmp(x, '[[''left'']]') * 0) + ...
    (strcmp(x, '[]') * -1), tmp_resp); %code key response

%----------- Deviant side (For experiment 2)
if strcmpi(typ, 'nwm2') | strcmpi(typ, 'wm2') % Store deviant side
    tmp_side = psych_dat.deviant_side; % Key response
    tmp_side = tmp_side(~cellfun( @isempty, tmp_side)); % Remove the empty character cell
    tmp_side = cell2mat(cellfun(@str2num, tmp_side, 'UniformOutput', false));

    tmp_side = tmp_side - 2; % HADI REMOVE
else
    tmp_side = nan(size(tmp_resp));

end

%   ----------------- q position

if strcmpi(typ, 'wm') | strcmpi(typ, 'wm2')
    tmp_q = psych_dat.q_position;
    tmp_q = tmp_q(~cellfun( @isempty, tmp_q)); % Remove the empty character cell
    tmp_q = cellfun(@(x) ...
        (strcmp(x, '(-34, 2.5)') * 0) + ...   % 'right' -> 1
        (strcmp(x, '(34, 2.5)') * 1) ...
        , tmp_q);
else
    tmp_q = nan((size(tmp_resp)));
end

%----------- Reaction time-----------------
tmp_rt = psych_dat.reaction_time; % Adaptor train
tmp_rt = tmp_rt(2:end-1);

numVector = NaN(size(tmp_rt));
for ii = 1:length(tmp_rt)
    % Check if the cell contains '[]' and assign NaN
    if strcmp(tmp_rt{ii}, '[]')
        numVector(ii) = NaN;
    else
        % Convert the number to double and assign
        numVector(ii) = str2double(tmp_rt{ii}(2:end-1));
    end
end

tmp_rt = numVector;


end