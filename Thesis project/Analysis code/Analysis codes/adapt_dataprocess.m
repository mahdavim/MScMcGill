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
    wd = 'C:\Users\mahdi\Desktop\Current Tray\Eyedata\Adapt data\Nwm';
elseif strcmpi(cnd, 'nwm2')
    wd = 'C:\Users\mahdi\Desktop\Current Tray\Eyedata\Adapt data\Nwm2';

elseif strcmpi(cnd, 'wm')
    wd = 'C:\Users\mahdi\Desktop\Current Tray\Eyedata\Adapt data\Wm';

else
    wd = 'C:\Users\mahdi\Desktop\Current Tray\Eyedata\Adapt data\Wm2';
end

%-------------


items = dir(wd);
items = items(~ismember({items.name}, {'.', '..'}));

% Number of cue positions !!!!!! THIS SHOULD BE CHANGED FOR ADAPTATION
eye=2;  %1 left 2 right

dist = 116.8;

if strcmpi(cnd, 'nwm') | strcmpi(cnd, 'nwm2') % Main events will be coded based on their

    msg_lst = ["Cue adap_ons", "Cue blank_ons"];% CHECK THIS WHENEVER THE RESPONSE IS EMPTY !!! CHANGE FOR ADAPTAION

else
    msg_lst = ["Cue cue_ons","Cue adap_ons", "Cue blank_ons", "Response resp_ons"];
end


%% Read and process the data for each participant and save it in a mat file

subject_id = 0;
sm_r = [];
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
            geny,pupil_a, st_times, en_times, slc_dvnt, adptr_len, key_resp,...
            dvnt_side, cue_side, key_rt] = ...
            adapt_dataimp_main(pth, eye, cnd);
        all_time = double(edfStruct1.FSAMPLE.time); % Vector of all time points from all trials

        disp("Recording freq: " + edfStruct1.RECORDINGS(1).sample_rate)
        sm_r = [sm_r, double(edfStruct1.RECORDINGS(1).sample_rate)];

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


            if double(edfStruct1.RECORDINGS(1).sample_rate) == 500
                gaze_struc(kk).gx = repelem(gaze_struc(kk).gx,2);
                gaze_struc(kk).gx(2) = gaze_struc(kk).gx(3);
                gaze_struc(kk).gy = repelem(gaze_struc(kk).gy,2);
                gaze_struc(kk).gy(2) = gaze_struc(kk).gy(3);
                gaze_struc(kk).pupilA = repelem(gaze_struc(kk).pupilA,2);
                gaze_struc(kk).pupilA(2) = gaze_struc(kk).pupilA(3);
            end


        end

        %Create time vectors using sampling rate
        trial_times = {};


        for kk=1:length(trial_end)

            trial_times{kk} = [all_time(idx(kk): idx_end(kk))] - all_time(idx(kk));% Reference with respect to the start of trial
            if double(edfStruct1.RECORDINGS(1).sample_rate) == 500
                if rem(kk,10)==0
                    trial_times{kk} = [0, trial_times{kk}(2):max(trial_times{kk})];
                else
                    trial_times{kk} = 0:max(trial_times{kk});
                end
            end
        end
        %--------------- End main process of each trial-----

        save("C:\Users\mahdi\Desktop\Current Tray\Eyedata\Adapt_proc\"+cnd+"\" +subject_id+".mat", ...
            "main_evnts", "trial_times", "gaze_struc", "datestr",...
            "sacade_struc", "slc_dvnt", "adptr_len", "key_resp",...
            "dvnt_side", "cue_side", "cnd", "key_rt")

    end





end


%% Loop over the saved data and read them to store all the info in a table

clear
cnd = input("What is the task condition that you'd like to access? ");
dat_dir = "C:\Users\mahdi\Desktop\Current Tray\Eyedata\Adapt_proc\"+cnd;

if strcmpi(cnd, 'wm2')
    tsk_name = "Auditory_Workingmemory_deviantOneside";
elseif strcmpi(cnd, 'wm')
    tsk_name = "Auditory_Workingmemory_deviantBothside";
elseif strcmpi(cnd, 'nwm2')
    tsk_name = "Auditory_passiveFixation_deviantOneside";
elseif strcmpi(cnd, 'nwm')
    tsk_name = "Auditory_passiveFixation_deviantBothside";
end



items = dir(dat_dir);
items = items(~ismember({items.name}, {'.', '..'}));
theta_list = [180, 0];%111 THIS SHOULD BE CHANGED FOR ROUTINE TASKS (0 is 180 and 1 is 0)
keys = ["No response", "Descending", "Ascending"]; %0 descending 1 ascending
dviant_sd = ["180", "0", "[0 & 180]"];


if strcmpi(cnd, 'nwm') | strcmpi(cnd, 'nwm2') % Main events will be coded based on their
    cue_num = 0;
    msg_lst = ["Cue adap_ons", "Cue blank_ons"];% CHECK THIS WHENEVER THE RESPONSE IS EMPTY !!! CHANGE FOR ADAPTAION

else
    cue_num = 2;
    msg_lst = ["Cue cue_ons","Cue adap_ons", "Cue blank_ons", "Response resp_ons"];
end

distnc = 116.8;
radius = rad2deg( atan(36/distnc)  );
dvnt_radius = radius;


subject_id = 0;
drift_trial = 10; %!! Change this


if strcmpi(cnd, 'wm') |strcmpi(cnd, 'nwm')

    all_dat = table([], [], [], [], [], [], [], [], [],[],[], [],[], [],[],[],[], ...
        [], [], [],[],[], [], [], [],[],[],...
        'VariableNames',["Id", "date", "taskName", "trialNumber", "fixationOn", "targetRadius",...
        "targetTheta", "targetOn", "targetOff", "fixationOff", "eyeData",...
        "saccadeReactionTime", "saccadeVector", "saccadeLandingPoint", "trialTime",...
        "deviantSlope", "adaptorNumber","keyResponse", "keyReactiontime", ...
        "adaptorOn", "adaptorOff", "deviantOn", "deviantOff", "deviantRadius1", ...
        "deviantTheta1", "deviantRadius2", "deviantTheta2"]);

else


    all_dat = table([], [], [], [], [], [], [], [], [],[],[], [],[], [],[],[],[], ...
        [], [], [],[],[], [],[],[],...
        'VariableNames',["Id", "date", "taskName", "trialNumber", "fixationOn", "targetRadius",...
        "targetTheta", "targetOn", "targetOff", "fixationOff", "eyeData",...
        "saccadeReactionTime", "saccadeVector", "saccadeLandingPoint", "trialTime",...
        "deviantSlope", "adaptorNumber","keyResponse", "keyReactiontime", ...
        "adaptorOn", "adaptorOff", "deviantOn", "deviantOff", "deviantRadius", ...
        "deviantTheta"]);
end

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
            sacade_rt =  sacade_struc(jj).start_t - main_evnts{jj}(4,2); %Calculate reaction times of saccades
            sacade_ent = sacade_struc(jj).end_t - main_evnts{jj}(4,2); % Saccade end time (!NEW)


        else
            sacade_rt = NaN;
            sacade_ent = NaN;
        end
        %-------Calculate saccade radious-------
        sacade_radious = sqrt( ( sacade_struc(jj).gstx-sacade_struc(jj).genx ).^2 +...
            ( sacade_struc(jj).gsty -  sacade_struc(jj).geny).^2);
        %--------------------
        sac_en_t = sacade_struc(jj).end_t(sacade_rt>0 & sacade_radious>2);
        sac_st_t = sacade_struc(jj).start_t(sacade_rt>0 & sacade_radious>2);

        try
        sac_st_t = sac_st_t(1);
        sac_en_t = sac_en_t(1);
        catch
            sac_st_t = NaN;
            sac_en_t = NaN;
        end

        if rem(jj, drift_trial)==0 % To account for drift correction times
            ofset = trial_times{jj}(2);
            trial_times{jj} = [0, trial_times{jj}(2:end)-ofset+1];
            main_evnts{jj}(1,2) = main_evnts{jj}(1,2) - ofset + 1;
            main_evnts{jj}(2,2) = main_evnts{jj}(2,2) -  ofset + 1;
            sac_st_t = sac_st_t-ofset;
            sac_en_t = sac_en_t- ofset;
            %------For wm experiments -------
            try
                main_evnts{jj}(3,2) = main_evnts{jj}(3,2) - ofset + 1;
                main_evnts{jj}(4,2) = main_evnts{jj}(4,2) - ofset + 1;
            catch
            end
        end





        % A threshold radious of 2dva was used to sieve the saccades
        tmp_rt = sacade_rt(sacade_rt>0 & sacade_radious>2);

        if isempty(tmp_rt)
            tmp_rt = NaN;
        end

        tmp_ent = sacade_ent(sacade_rt>0 & sacade_radious>2); % Saccade end (!NEW)

        try % Try to see if there are any valid saccades
            % scd_vectx = gaze_struc(jj).gx(tmp_rt(1)+main_evnts{jj}(2,2): tmp_ent(1) + main_evnts{jj}(2,2));
            % scd_vecty = gaze_struc(jj).gy(tmp_rt(1)+main_evnts{jj}(2,2): tmp_ent(1) + main_evnts{jj}(2,2));
            scd_vectx = gaze_struc(jj).gx(sac_st_t: sac_en_t);
            scd_vecty = gaze_struc(jj).gy(sac_st_t: sac_en_t);
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
            if strcmpi(dvnt_sd, "180")
                dvnt1_th = 180;
                dvnt1_r = dvnt_radius;

            else
                dvnt1_th = 0;
                dvnt1_r = dvnt_radius;
            end

        else
            dvnt_sd = dviant_sd(3);
            dvnt1_th = 180;
            dvnt2_th = 0;
            dvnt1_r = dvnt_radius;
            dvnt2_r = dvnt_radius;
        end


        %no working memory trials dont have target on/off and fixation off
        if strcmpi(cnd, 'wm') | strcmpi(cnd, 'wm2')
            tmp_taron = main_evnts{jj}(1,2);
            tmp_tarof = main_evnts{jj}(1,2)+200;
            tmp_fixof = main_evnts{jj}(4,2);
            cue_theta = theta_list(cue_side(jj)+1);

            tmp_adapon = main_evnts{jj}(2,2); %Adaptation onset wm
            tmp_adapoff = tmp_adapon + adptr_len(jj) * 810 - 10; %800 is tone length 10 ms between last adaoptor and test
            tmp_tston = tmp_adapoff + 10;
            tmp_tstoff = tmp_tston + 800;
        else
            tmp_taron = NaN;
            tmp_tarof = NaN;
            tmp_fixof = NaN;
            cue_theta = NaN;
            radius = NaN;
            scd = NaN;

            tmp_adapon = main_evnts{jj}(1,2); % Adaptation onset nwm
            tmp_adapoff = tmp_adapon + adptr_len(jj) * 810 - 10; %800 is tone length 10 ms between last adaoptor and test
            tmp_tston = tmp_adapoff + 10;
            tmp_tstoff = tmp_tston + 800;
        end



        if strcmpi(cnd, "wm")| strcmpi(cnd, "nwm")


            tmp_entry = {subject_id, string(datestr), tsk_name, jj,0, radius, cue_theta ,...
                tmp_taron, tmp_tarof,tmp_fixof, gaze_struc(jj),tmp_rt(1),...
                scd, [tmp_scdfinx(1),tmp_scdfiny(1)] , trial_times(jj), slc_dvnt(jj), adptr_len(jj), ...
                keys(key_resp(jj)+2), key_rt(jj), tmp_adapon, tmp_adapoff, ...
                tmp_tston, tmp_tstoff, dvnt1_r,  dvnt1_th, dvnt2_r, dvnt2_th};
            all_dat = [all_dat; tmp_entry];
            clear scd
        else
            tmp_entry = {subject_id, string(datestr), tsk_name, jj,0, radius, cue_theta ,...
                tmp_taron, tmp_tarof,tmp_fixof, gaze_struc(jj),tmp_rt(1),...
                scd, [tmp_scdfinx(1),tmp_scdfiny(1)] , trial_times(jj), slc_dvnt(jj), adptr_len(jj), ...
                keys(key_resp(jj)+2), key_rt(jj), tmp_adapon, tmp_adapoff, ...
                tmp_tston, tmp_tstoff, dvnt1_r,  dvnt1_th};
            all_dat = [all_dat; tmp_entry];
            clear scd
        end

    end % End of loop over trials



end % End of loop over participant data


all_dat(all_dat.trialNumber==1,: ) = [];

if strcmpi (cnd, 'wm') | strcmpi (cnd, 'nwm')
    all_dat = movevars(all_dat, ["adaptorOn", "adaptorOff", "deviantOn", "deviantOff",...
        "deviantSlope", "adaptorNumber", "deviantRadius1", "deviantTheta1", ...
        "deviantRadius2","deviantTheta2","keyResponse", "keyReactiontime"], "After","targetOff");

else
    all_dat = movevars(all_dat, ["adaptorOn", "adaptorOff", "deviantOn", "deviantOff",...
        "deviantSlope", "adaptorNumber", "deviantRadius", "deviantTheta", ...
        "keyResponse", "keyReactiontime"], "After","targetOff");
end


% save("C:\Users\mahdi\Desktop\Current Tray\Eyedata\Final tables\"+tsk_name+".mat", "all_dat")
%%

sv_dir = "C:\Users\mahdi\Desktop\Current Tray\Eyedata\Final tables\"+cnd+"\";

for ii=1:length(unique(all_dat.Id))

    tmp_dat = all_dat(all_dat.Id==ii, :);
    save(sv_dir+"subject"+ii+".mat", 'all_dat')

end

% save(sv_dir, 'all_dat')
%% Save all data fields in one structure



itms = dir("C:\Users\mahdi\Desktop\Current Tray\Eyedata\Final tables");
itms = itms(3:end);

all_audio_struc = struct();

for ii=1:length(itms)
    clear all_dat
    tmp_name = itms(ii).name;
    load("C:\Users\mahdi\Desktop\Current Tray\Eyedata\Final tables\"+ tmp_name)

    for jj=1:length(unique(all_dat.Id)) %Loop over IDs
        % all_audio_struc(jj).tmp_name(1:end-4) = all_dat(all_dat.Id==ii, :);
        eval("all_audio_struc(jj)."+tmp_name(1:end-4)+ "=all_dat(all_dat.Id==jj, :)");
    end
end
save("C:\Users\mahdi\Desktop\Current Tray\Eyedata\Final tables\Audio.mat", "all_audio_struc")
%% Separately save each data for participants

for jj=1:length(all_audio_struc)

    tmp_dat = all_audio_struc(jj);
    save("C:\Users\mahdi\Desktop\Current Tray\Eyedata\Final tables\strucs\subject"+...
        jj+".mat", "tmp_dat")
end




%%
load("C:\Users\mahdi\Desktop\Current Tray\Eyedata\Final tables\Data5.mat")

all_nms = fieldnames(sub01);

for ii= 1:length(all_nms)
    tmp_name = all_nms{ii};
    eval("all_audio_struc.Visual"+tmp_name+ "= sub01."+tmp_name);
end

save("C:\Users\mahdi\Desktop\Current Tray\Eyedata\Final tables\all_data.mat", "all_audio_struc")

%%
load("C:\Users\mahdi\Desktop\Current Tray\Eyedata\Final tables\Auditory_passiveFixation_deviantBothside.mat")


dd1 = all_dat;
dd1.keyResponse = strcmpi(dd1.keyResponse, "Ascending");
load("C:\Users\mahdi\Desktop\Current Tray\Eyedata\Final tables\Auditory_Workingmemory_deviantBothside.mat")

dd2 = all_dat;
dd2.keyResponse = strcmpi(dd2.keyResponse, "Ascending");


nwm_dd= [dd1.keyResponse, dd1.Id, dd1.deviantSlope];
nwm_dd = [nwm_dd, zeros(size(nwm_dd,1), 1)];
nwm_dd = array2table(nwm_dd, "VariableNames",["key", "id", "slope", "wm"]);
wm_dd = [dd2.keyResponse, dd2.Id, dd2.deviantSlope];
wm_dd = [wm_dd, ones(size(wm_dd,1), 1)];
wm_dd = array2table(wm_dd, "VariableNames",["key", "id", "slope", "wm"]);

all_dd = [nwm_dd;wm_dd];


forml = 'key ~ slope * wm + (1 + slope | id)';

glme = fitglme(all_dd, forml, 'Distribution', 'Binomial', 'Link', 'logit');

%% Yaser section
clear
load("C:\Users\mahdi\Desktop\Current Tray\Eyedata\Final tables\Audio.mat")

hlf_cong_all = [];
hlf_incong_all = [];

pc_x_cong = [];
pc_x_incong = [];
pc_y_cong = [];
pc_y_incong = [];

for ii=1:length(all_audio_struc)

    tmp_dat = all_audio_struc(ii).Auditory_Workingmemory_deviantOnesides;

    if isempty(tmp_dat)
        % continue
    else
        cong_dat = tmp_dat(tmp_dat.targetTheta==tmp_dat.deviantTheta, :);
        incong_dat = tmp_dat(tmp_dat.targetTheta~=tmp_dat.deviantTheta, :);

        cong_dv = cong_dat.deviantSlope;
        cong_k = cong_dat.keyResponse;
        cong_k = arrayfun(@(x) strcmpi(x, "Ascending"), cong_k);
        cong_dd = [cong_dv, cong_k];
        incong_dv = incong_dat.deviantSlope;
        incong_k = incong_dat.keyResponse;
        incong_k = arrayfun(@(x) strcmpi(x, "Ascending"), incong_k);
        incong_dd = [incong_dv, incong_k];

        [~,~,hlf_cong,~,cong_x, cong_y,cng_xx,cng_yy] =    fitting_log(cong_dd, 50, 1);
        [~,~,hlf_incong,~,incong_x, incong_y,incng_xx,incng_yy] =  fitting_log(incong_dd, 50, 1);

        hlf_cong_all = [hlf_cong_all, hlf_cong];
        hlf_incong_all = [hlf_incong_all, hlf_incong];
        pc_x_cong = [pc_x_cong, cong_x];
        pc_x_incong = [pc_x_incong, incong_x];
        pc_y_cong = [pc_y_cong, cong_y];
        pc_y_incong = [pc_y_incong, incong_y];
    end
end


save('C:\Users\mahdi\Desktop\Current Tray\Eyedata\Final tables\audio_fit.mat', ...
    "hlf_cong_all", "hlf_incong_all", "pc_x_cong", "pc_x_incong", "pc_y_cong", "pc_y_incong")

%% Hadi data
clear
dr = 'C:\Users\mahdi\Desktop\Current Tray\Eyedata\Hadi';
itms = dir(dr);
itms = itms(3:end);


hlf_cong_all = [];
hlf_incong_all = [];

pc_x_cong = [];
pc_x_incong = [];
pc_y_cong = [];
pc_y_incong = [];

pd_x_cong = [];
pd_x_incong = [];
pd_y_cong = [];
pd_y_incong = [];

pd_y_cong = nan(6, length(itms));
pd_y_incong = nan(6, length(itms));

for ii=1:length(itms)

    dat = load(dr+"\"+ itms(ii).name);
    nms =fieldnames(dat);
    eval("dat = dat."+nms{1}+";");
    dat = dat.workingMemoryOneside;

    dat.deviantTheta(dat.deviantTheta=="-180") = "180";
    if ii>1
        dat.targetTheta = string(dat.targetTheta);
    end
    dat.targetTheta(dat.targetTheta=="-180") = "180";

    tmp_key = (dat.keyResponse==1 & dat.deviantTheta=="0") | (dat.keyResponse==-1 & dat.deviantTheta=="180");
    cong_dat = dat(dat.targetTheta==dat.deviantTheta, :);
    cong_k = tmp_key(dat.targetTheta==dat.deviantTheta);

    incong_dat = dat(dat.targetTheta~=dat.deviantTheta, :);
    incong_k = tmp_key(dat.targetTheta~=dat.deviantTheta);

    cong_dv = cong_dat.deviantOri;



    cong_dd = [cong_dv, cong_k];
    incong_dv = incong_dat.deviantOri;

    incong_dd = [incong_dv, incong_k];

    [~,~,hlf_cong,~,cong_x, cong_y,cng_xx,cng_yy]  =    fitting_log(cong_dd, 50, 1);

    title(itms(ii).name + " Cong")
    xlim([-2 12])
    [~,~,hlf_incong,~,incong_x, incong_y,incng_xx,incng_yy] =  fitting_log(incong_dd, 50, 1);
    title(itms(ii).name+ " Incong")
    xlim([-2 12])

    hlf_cong_all = [hlf_cong_all, hlf_cong];
    hlf_incong_all = [hlf_incong_all, hlf_incong];
    pc_x_cong = [pc_x_cong, cong_x];
    pc_x_incong = [pc_x_incong, incong_x];
    pc_y_cong = [pc_y_cong, cong_y];
    pc_y_incong = [pc_y_incong, incong_y];

    pd_x_cong = [pd_x_cong, [0,2,4,6,8,10]'];
    pd_x_incong = [pd_x_incong, [0,2,4,6,8,10]'];
    pd_y_cong(1:length(cng_yy), ii) = cng_yy;
    pd_y_incong(1:length(incng_yy),ii) = incng_yy;
    save("C:\Users\mahdi\Desktop\Current Tray\Eyedata\Hadi data\new data\sub" + ii+".mat", "dat")

end
% save('C:\Users\mahdi\Desktop\Current Tray\Eyedata\Hadi\visual_fit.mat', ...
%     "hlf_cong_all", "hlf_incong_all", "pc_x_cong", "pc_x_incong", "pc_y_cong", "pc_y_incong", ...
%     "pd_x_cong", "pd_x_incong", "pd_y_cong", "pd_y_incong")

% save('C:\Users\mahdi\Desktop\Current Tray\Eyedata\Hadi\visual_fit.mat', ...
%     "hlf_cong_all", "hlf_incong_all", "pc_x_cong", "pc_x_incong", "pc_y_cong", "pc_y_incong")



%%
load("C:\Users\mahdi\Desktop\Current Tray\Eyedata\Final tables\Auditory_Workingmemory_deviantOnesides.mat")

dd1 = all_dat(all_dat.deviantTheta==all_dat.targetTheta, :);
dd1.keyResponse = strcmpi(dd1.keyResponse, "Ascending");

dd2 = all_dat(all_dat.deviantTheta~=all_dat.targetTheta, :);
dd2.keyResponse = strcmpi(dd2.keyResponse, "Ascending");


cng_all = [];
incng_all = [];

% for ii=1:length(unique(all_dat.deviantSlope))
%
%     tmp_cng = mean()

%% Remove blink artifacts
clear
dr = "C:\Users\mahdi\Desktop\Current Tray\Eyedata\Final tables\strucs";
% dr = "C:\Users\mahdi\Desktop\Current Tray\Eyedata\Hadi";

thr_y = 3;
thr_x = 3;
slpthr_y = 0.05;
slpthr_x = 0.05;

itms = dir(dr);
itms = itms(3:end);
for ii=1:length(itms)
    disp("Iter "+ii)
    tmp_dat = load(dr + "\"+itms(ii).name);
    tmp_nms = fieldnames(tmp_dat);
    eval("tmp_dat = tmp_dat."+ tmp_nms{1}+";") 
    tmp_nms = fieldnames(tmp_dat); 
    for jj=1:length(tmp_nms); % Loop over data structure field        
        tmp_f = eval("tmp_dat."+tmp_nms(jj)+";");
        if ~isempty(tmp_f)

        tmp_f.eyeDataFiltered = tmp_f.eyeData;
        for kk=1:size(tmp_f,1) % loop over table entries
            tmp_gy = tmp_f(kk,:).eyeData.gy;
            tmp_gx = tmp_f(kk,:).eyeData.gx;
            
            resp_t = tmp_f(kk,:).fixationOff;

            try % For passive fixation
            xx = tmp_gx(1:resp_t);
            catch
                xx = tmp_gx;
            end
            yy = tmp_gy;

            blinkMask_y = isnan(yy) | (abs(yy) > thr_y);% Threshold for extreme values
blinkMask_x = isnan(xx) | (abs(xx) > thr_x);% Threshold for extreme values

            % blinkMask = isnan(tmp_gy) | (abs(tmp_gy) > thr_y);% Threshold for extreme values

            dy = diff(yy);
            dx = diff(xx);

steep_regions_y = abs(dy) > slpthr_y;
steep_regions_x = abs(dx) > slpthr_x;

steep_indices_y = [find(steep_regions_y), find(steep_regions_y) +1];
steep_indices_x = [find(steep_regions_x), find(steep_regions_x) +1];

yy(steep_indices_y) = NaN;
yy(blinkMask_y) = NaN;
bas_yy = tmp_gy(abs(tmp_gy)<(abs(nanmean(tmp_gy))) + 0.4);

yy(isnan(yy)) = randsample(bas_yy, sum(isnan(yy)), 1); %Impute Y

xx(steep_indices_x) = NaN;
xx(blinkMask_x) = NaN;
bas_xx = tmp_gx(abs(tmp_gx)<(abs(nanmean(tmp_gx))) + 0.4);
xx(isnan(xx)) = randsample(bas_xx, sum(isnan(xx)), 1); %Impute X

yy = movmean(yy, 50);
xx = movmean(xx, 50);

try 
xx = [xx, movmean(tmp_gx(resp_t+1:end), 50)];
catch 
    xx = xx;
end


%--------------------- Interpolate 


% yData = struct("yData", yData);


tmp_f(kk,:).eyeDataFiltered.gy = yy; %Replace
tmp_f(kk,:).eyeDataFiltered.gx = xx;
        end % End of loop on table fields
        eval("tmp_dat."+tmp_nms(jj)+"=tmp_f"); %Replace the original data structure
        end% If statement to check if the field is empty
    end % End of loop over data structure fields
    save("C:\Users\mahdi\Desktop\Current Tray\Eyedata\Final tables\Modified_struc\" +...
        itms(ii).name, "tmp_dat");

        % save("C:\Users\mahdi\Desktop\Current Tray\Eyedata\Hadi data\Filtered data" +...
        % itms(ii).name, "tmp_dat"); %Hadi
end % End of loop over participants


%% Cherknevis

load ("C:\Users\mahdi\Desktop\Current Tray\Eyedata\Final tables\strucs\subject4.mat")


dd = tmp_dat.Auditory_Workingmemory_deviantOnesides;

yy = dd(10,:).eyeData.gy;
resp_t = dd(10,:).fixationOff;
xx  = dd(10,:).eyeData.gx;
xx = xx(1:resp_t);

figure
plot(yy)
figure
plot(xx)
thr_y = 3;
thr_x = 3;
slpthr_y = 0.01;
slpthr_x = 0.01;


blinkMask_y = isnan(yy) | (abs(yy) > thr_y);% Threshold for extreme values
blinkMask_x = isnan(xx) | (abs(xx) > thr_x);% Threshold for extreme values

dy = diff(yy);
dx = diff(xx);

steep_regions_y = abs(dy) > slpthr_y;
steep_regions_x = abs(dx) > slpthr_x;

steep_indices_y = [find(steep_regions_y), find(steep_regions_y) +1];
steep_indices_x = [find(steep_regions_x), find(steep_regions_x) +1];



