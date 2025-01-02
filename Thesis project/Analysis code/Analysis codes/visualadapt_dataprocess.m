clear all
close all
clc

% pth = 'C:\Users\mahdi\Desktop\Current Tray\Eyedata\Eyetracker analysis codes\edfImport-master\edfImport-master';
% cd(pth)

cnd = input("What is the input condition? ");

addpath ('C:\Users\mahdi\Desktop\Current Tray\Eyedata\Eyetracker analysis codes\edfImport-master\edfmex')
pth = genpath('C:\Users\mahdi\Desktop\Current Tray\Eyedata\Eyetracker analysis codes\edfImport-master\edfImport-master');

addpath(pth)

pth = {'C:\Users\mahdi\Desktop\Current Tray\Eyedata\Eyetracker analysis codes\preCompiled_edfmex\edfmex',...
    'C:\Users\mahdi\Desktop\Current Tray\Eyedata\wm1_nk_2_2024_08_21_15_56.EDF', ...
    'C:\Users\mahdi\Desktop\Current Tray\Eyedata\wm1_nk_2_sample 1_2024-08-21_15h56.57.491.csv'};

%---------- Select wd folder based on experiment type

if strcmpi(cnd, 'nwm')
    wd = 'C:\Users\mahdi\Desktop\Current Tray\Eyedata\Hadi data\Nwm';
elseif strcmpi(cnd, 'nwm2')
    wd = 'C:\Users\mahdi\Desktop\Current Tray\Eyedata\Hadi data\Nwm2';

elseif strcmpi(cnd, 'wm')
    wd = 'C:\Users\mahdi\Desktop\Current Tray\Eyedata\Hadi data\Wm';

else
    wd = 'C:\Users\mahdi\Desktop\Current Tray\Eyedata\Hadi data\Wm2';
end

%-------------


items = dir(wd);
items = items(~ismember({items.name}, {'.', '..'}));

% Number of cue positions !!!!!! THIS SHOULD BE CHANGED FOR ADAPTATION
eye=2;  %1 left 2 right

dist = 116.8;

if strcmpi(cnd, 'nwm') | strcmpi(cnd, 'nwm2') % Main events will be coded based on their
    
    msg_lst = [];% CHECK THIS WHENEVER THE RESPONSE IS EMPTY !!! CHANGE FOR ADAPTAION

else
    msg_lst = ["Cue Cue_ons", "Response Resp_ons"];
end


%% Read and process the data for each participant and save it in a mat file

subject_id = 0;
for ii=1: length(items) % Loop over all participant folders



    matchingFiles = {};
    if items(ii).isdir && ~strcmp(items(ii).name, '.') && ~strcmp(items(ii).name, '..')
        subject_id = subject_id + 1;
        disp(["Subject: "+ subject_id])
        folderPath = fullfile(wd, items(ii).name); % Create oath for the item

        filesInFolder = dir(folderPath);

        for j = 1:length(filesInFolder) %Loop over files for a participants

            [~, ~, ext] = fileparts(filesInFolder(j).name); % Get file extension

            % Check if the file has either '.xlsx' or '.EDF' extension
            if strcmpi(ext, '.csv')

                pth{3} = fullfile(folderPath, filesInFolder(j).name);
            end
            if strcmpi(ext, '.EDF')

                pth{2} = fullfile(folderPath, filesInFolder(j).name);

            end

        end
        % Read the data
        %====== Process the data for each participant ===============
        [edfStruct1, datestr, gx, gy, gstx, gsty, genx,...
            geny,pupil_a, st_times, en_times, slc_dvnt, key_resp,...
            dvnt_side, cue_side, key_rt] = ...
            visualadapt_dataimp_main(pth, eye, cnd);
        all_time = double(edfStruct1.FSAMPLE.time); % Vector of all time points from all trials



        %---------Initial process of messages

        msgs = {edfStruct1.FEVENT.message}; % the difference between the first trialid time and start time of eyetracker is due to
        % the initial lag and drift correction. This is ok.


        timestmps = double([edfStruct1.FEVENT.sttime]);


        idx = cellfun(@isempty,msgs); % Find the indexes of empty cell
        msgs(idx) = {''};


        %---------End initial process of messages------

        %---- sometimes there is a trial result -1 which needs to be excluded
        if any(contains(msgs,'TRIAL_RESULT -1'))

            idc_end = find(contains(msgs,'TRIAL_RESULT') & ~contains(msgs,'TRIAL_RESULT -1')); % Find cells that have TRIALID in them

        else
            idc_end = find(contains(msgs,'TRIAL_RESULT') ); % Find cells that have TRIALID in them
        end

        %---------------
        %-------------Find start and end of trials

        trial_end = timestmps(idc_end);% Find trial end time (ms)

        idc = find(contains(msgs,'TRIALID')); % Find cells that have TRIALID in them
        trial_sttime = timestmps(idc); % Find the starting time of trials (in ms)

        %----------------- Start of the main process of trials -------------

        idx = dsearchn(all_time', trial_sttime');
        idx_end = dsearchn(all_time', trial_end'); % Find end time points

        gaze_struc = struct('gx', [],'gy', [], 'pupilA', []);

        trial_msgs = {};

        trial_msgs_time = {};
        trial_msgs_time_trialref = {}; % Trial times with respect to the start time of each trial
        main_evnts = {};
        sacade_struc = struct('gstx',[],  'gsty',[], 'genx',[], 'geny',[], 'st_times', [],'en_times',[]);
        trial_code_time = {};

        sacade_events = struct('start_t', 'end_t', 'gstx', 'gsty', 'genx', 'geny');

        for kk=1:length(trial_end)
            tmp_timestmp = timestmps(idc(kk):idc_end(kk));
            trial_msgs{kk} = msgs(idc(kk):idc_end(kk)); %get all the messages between the start and end trial messages
            trial_msgs_time{kk} = timestmps(idc(kk):idc_end(kk)); % get trial message time stamps
            tmp_idx = dsearchn(all_time', trial_msgs_time{kk}');

            trial_msgs_time{kk} = all_time(tmp_idx);
            trial_msgs_time_trialref{kk} =  trial_msgs_time{kk} - all_time(idx(kk));

            % Extract and store saccade events
            scd_msk = tmp_timestmp(1) <= st_times & st_times <= tmp_timestmp(end);% Saccade mask to choose points whose start time is within the trial
            sacade_struc(kk).start_t = st_times(scd_msk) - all_time(idx(kk)); %Convert time to trial time
            sacade_struc(kk).end_t = en_times(scd_msk)- all_time(idx(kk));

            sacade_struc(kk).gstx = gstx(scd_msk);
            sacade_struc(kk).gsty = gsty(scd_msk);
            sacade_struc(kk).genx = genx(scd_msk);
            sacade_struc(kk).geny = geny(scd_msk);

            tmp_main_envnts = [];
            for jj=1:length(msg_lst)

                tmp_idx = find(contains(trial_msgs{kk},msg_lst(jj)) );
                tmp_t = trial_msgs_time_trialref{kk}(tmp_idx);
                tmp_main_envnts = [tmp_main_envnts; [jj, tmp_t]];


            end
            main_evnts{kk} = tmp_main_envnts;

            disp("Trial:" + kk)
        end




        for kk=1:length(idx)



            gaze_struc(kk).gx = gx(idx(kk): idx_end(kk));
            gaze_struc(kk).gy = gy(idx(kk): idx_end(kk));
            gaze_struc(kk).pupilA = pupil_a(idx(kk): idx_end(kk));


        end

        %Create time vectors using sampling rate
        trial_times = {};


        for kk=1:length(trial_end)

            trial_times{kk} = [all_time(idx(kk): idx_end(kk))] - all_time(idx(kk));% Reference with respect to the start of trial
        end
        %--------------- End main process of each trial-----

        save("C:\Users\mahdi\Desktop\Current Tray\Eyedata\Adapt_proc_hd\"+cnd+"\" +subject_id+".mat", ...
            "main_evnts", "trial_times", "gaze_struc", "datestr",...
            "sacade_struc", "slc_dvnt", "key_resp",...
            "dvnt_side", "cue_side", "cnd", "key_rt")

    end





end


%% Loop over the saved data and read them to store all the info in a table

clear
cnd = input("What is the task condition that you'd like to access? ");
dat_dir = "C:\Users\mahdi\Desktop\Current Tray\Eyedata\Adapt_proc_hd\"+cnd;

if strcmpi(cnd, 'wm2')
    tsk_name = "Visual_Working_memory_exp2";
elseif strcmpi(cnd, 'wm')
    tsk_name = "Visual_Working_memory_exp1";
elseif strcmpi(cnd, 'nwm2')
    tsk_name = "Visual_Noworking_memory_exp2";
elseif strcmpi(cnd, 'nwm')
    tsk_name = "Visual_Noworking_memory_exp1";
end



items = dir(dat_dir);
items = items(~ismember({items.name}, {'.', '..'}));
theta_list = ["180", "0"];%111 THIS SHOULD BE CHANGED FOR ROUTINE TASKS (0 is 180 and 1 is 0)
keys = ["No response" ,"Left_tilt", "Right_tilt"]; %0 descending 1 ascending
dviant_sd = ["Left", "Right"];


if strcmpi(cnd, 'nwm') | strcmpi(cnd, 'nwm2') % Main events will be coded based on their
    cue_num = 0;
    msg_lst = ["Cue adap_ons", "Cue blank_ons"];% CHECK THIS WHENEVER THE RESPONSE IS EMPTY !!! CHANGE FOR ADAPTAION

else
    cue_num = 2;
    msg_lst = ["Cue Cue_ons", "Response Resp_ons"];
end

distnc = 116.8;
radius = rad2deg( atan(9.75/distnc)  );


subject_id = 0;
drift_trial = 5; %!! Change this

all_dat = table([], [], [], [], [], [], [], [], [],[],[], [], [],[],[], [], [], [],[],...
    'VariableNames',["Id", "date", "taskName", "trialNumber", "fixationOn", "targetRadius",...
    "targetTheta", "targetOn", "targetOff", "fixationOff", "eyeData",...
    "saccadeReactionTime", "saccadeVector", "saccadeLandingPoint", "trialTime",...
    "deviantSlope","keyResponse", "deviantSide", "KeyReactiontime"]);


% In saved main event, the first entry is for cue and the second one is
% fore response onset
for ii=1: length(items) % Loop over all participant folders
    subject_id = subject_id + 1;
    disp(["Subject: "+ subject_id])
    folderPath = fullfile(dat_dir, items(ii).name); % Create oath for the item
    load(folderPath)

    for jj=1:length(main_evnts)% Loop over each trial
        %--------calculate saccade reaction time based on fixation off time

        if strcmpi(cnd, 'wm') |  strcmpi(cnd, 'wm2')
            sacade_rt =  sacade_struc(jj).start_t - main_evnts{jj}(2,2); %Calculate reaction times of saccades
            sacade_ent = sacade_struc(jj).end_t - main_evnts{jj}(2,2); % Saccade end time (!NEW)
            %Reaction times
        else
            sacade_rt = NaN;
            sacade_ent = NaN;
        end


        if rem(jj, drift_trial)==0 % To account for drift correction times
            ofset = trial_times{jj}(2);
            trial_times{jj} = [0, trial_times{jj}(2:end)-ofset+1];

            %------For wm experiments -------
            try
                main_evnts{jj}(1,2) = main_evnts{jj}(1,2) - ofset + 1;
                main_evnts{jj}(2,2) = main_evnts{jj}(2,2) - ofset + 1;
            catch
            end
        end



        %-------Calculate saccade radious-------
        sacade_radious = sqrt( ( sacade_struc(jj).gstx-sacade_struc(jj).genx ).^2 +...
            ( sacade_struc(jj).gsty -  sacade_struc(jj).geny).^2);

        % A threshold radious of 2dva was used to sieve the saccades
        tmp_rt = sacade_rt(sacade_rt>0 & sacade_radious>2);

        if isempty(tmp_rt)
            tmp_rt = NaN;
        end

        tmp_ent = sacade_ent(sacade_rt>0 & sacade_radious>2); % Saccade end (!NEW)

        try % Try to see if there are any valid saccades
            scd_vectx = gaze_struc(jj).gx(tmp_rt(1)+main_evnts{jj}(2,2): tmp_ent(1) + main_evnts{jj}(2,2));
            scd_vecty = gaze_struc(jj).gy(tmp_rt(1)+main_evnts{jj}(2,2): tmp_ent(1) + main_evnts{jj}(2,2));
        catch
            scd_vectx = NaN;
            scd_vecty = NaN;
        end


        % Some trials the participants might not pay attention or forget to saccade
        %Define a threshold of 2 for sacade radious
        tmp_scdfinx = sacade_struc(jj).genx(sacade_rt>0 & sacade_radious>2);% Get sacade final position sand initial position for those that started after fixation off
        if isempty(tmp_scdfinx)
            tmp_scdfinx = NaN;
            scd.engx = NaN;
        end

        tmp_scdstx = sacade_struc(jj).gstx(sacade_rt>0& sacade_radious>2);
        if isempty(tmp_scdstx)
            tmp_scdstx = NaN;
            scd.stgx = NaN;
        end

        tmp_scdfiny = sacade_struc(jj).geny(sacade_rt>0& sacade_radious>2);
        if isempty(tmp_scdfiny)
            tmp_scdfiny = NaN;
            scd.engy = NaN;
        end


        tmp_scdsty = sacade_struc(jj).gsty(sacade_rt>0& sacade_radious>2);
        if isempty(tmp_scdsty)
            tmp_scdsty = NaN;
            scd.stgy = NaN;
        end




        %Saccade vector

        scd.stgx = tmp_scdstx(1);
        scd.engx = tmp_scdfinx(1);
        scd.stgy = tmp_scdsty(1);
        scd.engy = tmp_scdfiny(1);
        scd.r = sacade_radious(sacade_rt>0 & sacade_radious>2);

        if isempty(scd.r)
            scd.r = NaN;
        else
            scd.r = scd.r(1); % Choose the first saccade raious
        end
        scd.x = scd_vectx; %(!NEW)
        scd.y = scd_vecty;


        % If there is deviant side data, convert it into left/right
        if ~isnan(dvnt_side(jj))
            dvnt_sd = dviant_sd(dvnt_side(jj)+1);
        else
            dvnt_sd = NaN;
        end


        %no working memory trials dont have target on/off and fixation off
        if strcmpi(cnd, 'wm') | strcmpi(cnd, 'wm2')
            tmp_taron = main_evnts{jj}(1,2);
            tmp_tarof = main_evnts{jj}(1,2)+200;
            tmp_fixof = main_evnts{jj}(2,2);
            cue_theta = theta_list(cue_side(jj)+1);
        else
            tmp_taron = NaN;
            tmp_tarof = NaN;
            tmp_fixof = NaN;
            cue_theta = NaN;
            radius = NaN;
            scd = NaN;
        end

        tmp_entry = {subject_id, string(datestr), tsk_name, jj,0, radius, cue_theta ,...
            tmp_taron, tmp_tarof,tmp_fixof, gaze_struc(jj),tmp_rt(1),...
            scd, [tmp_scdfinx(1),tmp_scdfiny(1)] , trial_times(jj), slc_dvnt(jj), ...
            keys(key_resp(jj)+2), dvnt_sd, key_rt(jj)};
        all_dat = [all_dat; tmp_entry];
        clear scd

    end % End of loop over trials
    all_dat(1,: ) = [];


end % End of loop over participant data

%%
save('V_table_nwm.mat', 'all_dat')

load('adap_table.mat')


kkk = all_dat.saccadeVector;

