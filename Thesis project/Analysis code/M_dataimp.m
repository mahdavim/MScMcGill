function [edfStruct1, cat_var, lbl_names, tmp_deviant, tmp_key, tmp_adapt,tmp_side,tmp_time, gx, gy] = ...
M_dataimp(path, cue_pos, eye, cue_positions, keys)

addpath(genpath(path{1})) % Main codes


edfStruct1 = edfmex(path{2}); % EDF file
psych_dat = readtable(path{3}); %psychopy csv file

%---- import cue positions and convert them to 1:8 categorical values
tmp_q = psych_dat.q_position(~cellfun(@isempty,psych_dat.q_position)); % Que positions.only choose non-empty cells
cat_var = double(categorical(tmp_q,cue_positions, {'1','2'} )); 
lbl_names = unique(tmp_q);
% cats = split(num2str(1:cue_pos), ' '); % Number of different cue positions in the task
% cats = cats(strlength(cats) > 0); %categories
% 
% cat_var = renamecats(cat_var, cats); % Rename trials to dummy category names
% cat_var = cat_var(2:end); % Discard the first trial
%------ Extract deviant, adaptor train number, deviant side, reaction time for key, key data

tmp_deviant = psych_dat.selected_dvnt(~cellfun(@isempty,psych_dat.selected_dvnt));
tmp_deviant = cellfun(@(x) regexp(x, '\d+', 'match'), tmp_deviant); % Remove[] from numbers
tmp_deviant = cellfun(@str2double, tmp_deviant); % Convert strings to numbers and the cell to a matrix

tmp_adapt = psych_dat.adptr_train_list(~cellfun(@isempty,psych_dat.adptr_train_list));
tmp_adapt = cellfun(@(x) regexp(x, '\d+', 'match'), tmp_adapt); % Remove[] from numbers
tmp_adapt = cellfun(@str2double, tmp_adapt); 

try
tmp_side = psych_dat.deviant_side(~cellfun(@isempty,psych_dat.deviant_side));
tmp_side = cellfun(@(x) regexp(x, '\d+', 'match'), tmp_side); % Remove[] from numbers
tmp_side = cellfun(@str2double, tmp_side); 
catch
tmp_side = NaN;
end

tmp_time = psych_dat.reaction_time(~cellfun(@isempty,psych_dat.reaction_time));
tmp_time = cellfun(@(x) regexp(x, '[\d.]+', 'match'), tmp_time); % Remove[] from numbers
tmp_time = cellfun(@str2double, tmp_time); 

tmp_key = psych_dat.key_index(~cellfun(@isempty,psych_dat.key_index)); % Que positions.only choose non-empty cells
tmp_key = double(categorical(tmp_key,keys, {'1','2'} )); 

%-----------


gx = double(edfStruct1.FSAMPLE.gx(eye,1:end)); % For the right eye(2) or left (1)
gx(gx == max(gx)) = nan;
gy = double(edfStruct1.FSAMPLE.gy(eye,1:end));
gy(gy  == max(gy)) = nan;


end