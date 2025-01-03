function [all_thr, all_dat, all_log, log_struc,response_ratios, cong_log,...
    incong_log, cong_resp_ratio, incong_resp_ratio] = read_psychophys(pth, log_win, plt)

% typ: thr, nwm, wm
%log_win: window for logistic fit
%{
 all_thr: deviant code and response key to plot logistic 
 all_dat: 'dviant', 'response', 'adaptor_train', 'reaction_time'

%}
%---------- If block to choose data table format
typ = input('Experiment type?');
if strcmpi(typ, 'nwm')
    all_dat = table([],[],[],[],[], [], 'VariableNames',{'dviant', 'response',...
        'adaptor_train', 'wmGroup', 'ID','reaction_time' }); % all data extracted
elseif  strcmpi(typ, 'nwm2')
    all_dat = table([],[],[],[],[],[], [], 'VariableNames',{'dviant', 'response',...
        'adaptor_train', 'wmGroup', 'ID','dviant_side', 'reaction_time'}); % all data extracted
elseif strcmpi(typ, 'wm')
    all_dat = table([],[],[],[],[],[], [], 'VariableNames',{'dviant', 'response', ...
        'adaptor_train', 'q_position','wmGroup', 'ID', 'reaction_time'}); % all data extracted
elseif strcmpi (typ, 'wm2')
    all_dat = table([],[],[],[],[],[],[], [],[], 'VariableNames',{'dviant', 'response', ...
        'adaptor_train','q_position','wmGroup', 'ID', 'congGroup', 'dviant_side', 'reaction_time', }); % all data extracted

end
%--------------

all_log = table([],[],[],[],[],'VariableNames',{'alpha', 'beta', 'half_max_x', 'amplitude', 'n_trials'}); % all data extracted

cong_log=struct();
incong_log=struct();


log_struc = struct();% original data points and x and y for the fited function

thr_itms = dir(pth);
thr_itms = thr_itms(~ismember({thr_itms.name}, {'.', '..'}));

if strcmpi(typ, 'thr')
    all_thr = [];
    dvnts = [-8, -6, -4, -2, -0.5, 0.5, 2, 4, 6, 8];
    for ii = 1:length(thr_itms)
        tmp_dat = readtable(fullfile(pth, thr_itms(ii).name));
        tmp_dvnt = tmp_dat.selected_dvnt;
        tmp_dvnt = cell2mat(cellfun(@str2num, tmp_dvnt, 'UniformOutput', false));
        dvs = unique(tmp_dvnt);

        % ------------ Enter true deviant values in octave multiplied by
        %  96
        for kk=1:length(dvs)
            tmp_dvnt(tmp_dvnt==dvs(kk)) = dvnts(kk);
        end


        tmp_resp = tmp_dat.key_index;
        tmp_resp = tmp_resp(~cellfun( @isempty, tmp_resp)); % Remove the empty character cell
        tmp_resp = cellfun(@(x) ...
            (strcmp(x, '[[''right'']]') * 1) + ...   % 'right' -> 1
            (strcmp(x, '[[''left'']]') * 0) + ...
            (strcmp(x, '[]') * -1), tmp_resp);

        tmp_dat = [tmp_dvnt, tmp_resp];
        tmp_dat(tmp_dat(:,2)==-1,:) = [];
        [tmp_alpha, tmp_beta, tmp_half_max_x, tmp_amplitude] = fitting_log(tmp_dat, log_win, plt);
        tmp_ntrials = size(tmp_dat,1);
        tmp_log = table(tmp_alpha,tmp_beta,tmp_half_max_x,tmp_amplitude,tmp_ntrials, ...
            'VariableNames',{'alpha', 'beta', 'half_max_x', 'amplitude', 'n_trials'});
        all_log = [all_log; tmp_log]; % Concat logistic properties


        fg = gcf;
        fg.Name = thr_itms(ii).name;
        hold on
        xline(7.5)

        all_thr = [all_thr; tmp_dat];
        all_dat = all_thr;
        log_struc = NaN;
    end




elseif strcmpi(typ, 'nwm') | strcmpi(typ, 'nwm2')
    all_thr = []; % Only deviant code and key response
    % all_dat = table([],[],[],[],'VariableNames',{'dviant', 'response', 'adaptor_train', 'reaction_time'}); % all data extracted
    response_ratios = [];

    if strcmpi(typ, 'nwm') %!!CHANGE
        dvnts = [-8, -6, -4, -2,0, 2, 4, 6, 8];


    elseif strcmpi(typ, 'nwm2')
        dvnts = [-12, -9, -6, -3,0, 3, 6, 9, 12];

    end
par_id = 0; % Participant ID
    for ii = 1:length(thr_itms) % Loop over participants
        par_id = par_id +1;
        tmp_dat = readtable(fullfile(pth, thr_itms(ii).name));
        tmp_dvnt = tmp_dat.selected_dvnt;% Selected deviant
        tmp_dvnt = cell2mat(cellfun(@str2num, tmp_dvnt, 'UniformOutput', false));
        dvs = unique(tmp_dvnt);

        % Enter true deviant values
        for kk=1:length(dvs)
            tmp_dvnt(tmp_dvnt==dvs(kk)) = dvnts(kk);
        end
        %-----------
        tmp_adapt = tmp_dat.adptr_train_list; % Adaptor train
        tmp_adapt = cell2mat(cellfun(@str2num, tmp_adapt, 'UniformOutput', false));

        %--------------response time
        tmp_rt = tmp_dat.reaction_time; %
        tmp_rt = cell2mat(cellfun(@str2num, tmp_rt, 'UniformOutput', false));

        % ------------- Key response ------------
        tmp_resp = tmp_dat.key_index; % Key response
        tmp_resp = tmp_resp(~cellfun( @isempty, tmp_resp)); % Remove the empty character cell
        tmp_resp = cellfun(@(x) ...
            (strcmp(x, '[[''right'']]') * 1) + ...   % 'right' -> 1
            (strcmp(x, '[[''left'']]') * 0) + ...
            (strcmp(x, '[]') * -1), tmp_resp); %code key response

        % ------------- choose dviant side if nwm2
        if strcmpi(typ, 'nwm2') % Store deviant side
            tmp_side = tmp_dat.deviant_side; % Key response
            tmp_side = tmp_side(~cellfun( @isempty, tmp_side)); % Remove the empty character cell
            tmp_side = cell2mat(cellfun(@str2num, tmp_side, 'UniformOutput', false));
            tmp_side = tmp_side - 2;

        end
        %-------------

        %--------------- Add WM and ID columns -------------
        tmp_id = par_id*ones(size(tmp_dvnt));
        tmp_wm = zeros(size(tmp_dvnt));




        tmp_dat = [tmp_dvnt, tmp_resp];
        tmp_all = [tmp_dvnt, tmp_resp, tmp_adapt, tmp_wm,tmp_id ];

        if strcmpi(typ, 'nwm2')
            tmp_all = [tmp_all, tmp_side];
        end
        tmp_all(tmp_dat(:,2)==-1,:) = []; % tmp dat will be overwriten in next line
        tmp_dat(tmp_dat(:,2)==-1,:) = [];
        tmp_all = [tmp_all, tmp_rt];% empty rts will be excluded previously in the code

        %=========== Extract ratio of responses for each deviant type
        tmp_ratio = [];
        for ll = 1:length(dvnts)
            tmp_ratio = [tmp_ratio, mean( tmp_all(tmp_all(:,1)==dvnts(ll), 2) )];

        end
        response_ratios = [response_ratios;tmp_ratio];
        %--------- For second experiment, add deviant side------------
        if strcmpi(typ, 'nwm')
            tmp_all = array2table(tmp_all, "VariableNames",{'dviant', 'response', 'adaptor_train',...
                'wmGroup', 'ID','reaction_time'});
        elseif strcmpi(typ, 'nwm2')
            tmp_all = array2table(tmp_all, "VariableNames",{'dviant', 'response', 'adaptor_train',...
                'wmGroup', 'ID','dviant_side', 'reaction_time'});
        end
        %--------- ------------

        all_dat = [all_dat; tmp_all];
        tmp_ntrials = size(tmp_all,1);% Number of trials

        [tmp_alpha, tmp_beta, tmp_half_max_x, tmp_amplitude, pc_x, pc_y, pd_x, pd_y] = fitting_log(tmp_dat, log_win, plt);
        log_struc(ii).pc_x = pc_x;
        log_struc(ii).pc_y = pc_y;
        log_struc(ii).pd_x = pd_x;
        log_struc(ii).pd_y = pd_y;

        tmp_log = table(tmp_alpha,tmp_beta,tmp_half_max_x,tmp_amplitude,tmp_ntrials,...
            'VariableNames',{'alpha', 'beta', 'half_max_x', 'amplitude', 'n_trials'});
        all_log = [all_log; tmp_log]; % Concat logistic properties

        fg = gcf;
        fg.Name = thr_itms(ii).name;
        hold on
        xline(7, 'k--')

        all_thr = [all_thr; tmp_dat];
    end % End of loop over participants
    response_ratios = array2table(response_ratios, 'VariableNames',string(dvnts));

    % ----------------- WM --------------------
elseif strcmpi(typ, 'wm') | strcmpi(typ, 'wm2')
    all_thr = []; % Only deviant code and key response
    % all_dat = table([],[],[],[],'VariableNames',{'dviant', 'response', 'adaptor_train', 'reaction_time'}); % all data extracted
    response_ratios = []
    if strcmpi(typ, 'wm')
        dvnts = [-8, -6, -4, -2,0, 2, 4, 6, 8];
    elseif strcmpi(typ, 'wm2')
        dvnts = [-12, -9, -6, -3,0, 3, 6, 9, 12];
        cong_resp_ratio = [];
        incong_resp_ratio = [];
    end
par_id = 0;%Participant ID
    for ii = 1:length(thr_itms) %Loop over participants
        par_id = par_id+1;
        tmp_dat = readtable(fullfile(pth, thr_itms(ii).name));
        tmp_dvnt = tmp_dat.selected_dvnt;% Selected deviant
        tmp_dvnt = cell2mat(cellfun(@str2num, tmp_dvnt, 'UniformOutput', false));
        dvs = unique(tmp_dvnt);

        % Enter true deviant values
        for kk=1:length(dvs)
            tmp_dvnt(tmp_dvnt==dvs(kk)) = dvnts(kk);
        end
        %------------- Adaptor train
        tmp_adapt = tmp_dat.adptr_train_list; % Adaptor train
        tmp_adapt = cell2mat(cellfun(@str2num, tmp_adapt, 'UniformOutput', false));

        %----------- Reaction time-----------------
        tmp_rt = tmp_dat.reaction_time; % Adaptor train
        tmp_rt = cell2mat(cellfun(@str2num, tmp_rt, 'UniformOutput', false));

        %   ---------- Key response-----------
        tmp_resp = tmp_dat.key_index; % Key response
        tmp_resp = tmp_resp(~cellfun( @isempty, tmp_resp)); % Remove the empty character cell
        tmp_resp = cellfun(@(x) ...
            (strcmp(x, '[[''right'']]') * 1) + ...   % 'right' -> 1
            (strcmp(x, '[[''left'']]') * 0) + ...
            (strcmp(x, '[]') * -1), tmp_resp); %code key response

        %   ----------------- q position
        tmp_q = tmp_dat.q_position;
        tmp_q = tmp_q(~cellfun( @isempty, tmp_q)); % Remove the empty character cell
        tmp_q = cellfun(@(x) ...
            (strcmp(x, '(-34, 2.5)') * 0) + ...   % 'right' -> 1
            (strcmp(x, '(34, 2.5)') * 1) ...
            , tmp_q); %que position
        % --------------
        % ------------- choose dviant side if wm2
        if strcmpi(typ, 'wm2') % Store deviant side
            tmp_side = tmp_dat.deviant_side; % Key response
            tmp_side = tmp_side(~cellfun( @isempty, tmp_side)); % Remove the empty character cell
            tmp_side = cell2mat(cellfun(@str2num, tmp_side, 'UniformOutput', false));
            tmp_side = tmp_side - 2;

            tmp_cong = (tmp_side == tmp_q);

        end
        %-------------
        % -------- Add Id and wm group 
        tmp_id = par_id*ones(size(tmp_dvnt));
        tmp_wm = ones(size(tmp_dvnt));






        tmp_dat = [tmp_dvnt, tmp_resp];
        tmp_all = [tmp_dvnt, tmp_resp, tmp_adapt, tmp_q, tmp_wm,tmp_id ];
        %------ add dviant side if wm2
        if strcmpi(typ, 'wm2')
            tmp_all = [tmp_all, tmp_cong, tmp_side];
        end
        %----------------
        tmp_all(tmp_dat(:,2)==-1,:) = []; % tmp dat will be overwriten in next line
        tmp_dat(tmp_dat(:,2)==-1,:) = [];
        tmp_all = [tmp_all, tmp_rt];% empty rts will be excluded previously in the code


        %--------------- Extract ratio of responses for each deviant type
        tmp_ratio = [];
        for ll = 1:length(dvnts)
            tmp_ratio = [tmp_ratio, mean( tmp_all(tmp_all(:,1)==dvnts(ll), 2) )];

        end
        response_ratios = [response_ratios;tmp_ratio];
        %-----------------
        %--------- For second experiment, add deviant side------------
        if strcmpi(typ, 'wm')
            tmp_all = array2table(tmp_all, "VariableNames",{'dviant', 'response', 'adaptor_train','q_position',...
                'wmGroup', 'ID','reaction_time'});
        elseif strcmpi(typ, 'wm2')
            tmp_all = array2table(tmp_all, "VariableNames",{'dviant', 'response',...
                'adaptor_train','q_position', 'wmGroup', 'ID', 'congGroup','dviant_side', 'reaction_time'});
        end
        %--------- ------------

        all_dat = [all_dat; tmp_all];
        tmp_ntrials = size(tmp_all,1);

        [tmp_alpha, tmp_beta, tmp_half_max_x, tmp_amplitude, pc_x, pc_y, pd_x, pd_y] = fitting_log(tmp_dat, log_win, plt);

        log_struc(ii).pc_x = pc_x;
        log_struc(ii).pc_y = pc_y;
        log_struc(ii).pd_x = pd_x;
        log_struc(ii).pd_y = pd_y;
        tmp_log = table(tmp_alpha,tmp_beta,tmp_half_max_x,tmp_amplitude,tmp_ntrials, ...
            'VariableNames',{'alpha', 'beta', 'half_max_x', 'amplitude', 'n_trials'});
        all_log = [all_log; tmp_log]; % Concat logistic properties

        fg = gcf;
        fg.Name = thr_itms(ii).name;
        hold on
        xline(7, 'k--')

        all_thr = [all_thr; tmp_dat];

        if strcmpi(typ, 'wm2') %
            tmp_cong = tmp_all( tmp_all.dviant_side == tmp_all.q_position, :);
            tmp_cong_rt = tmp_cong.reaction_time;

            tmp_cong = table2array(tmp_cong(:,1:2)); % Overwrite to remove other columns

            tmp_incong  = tmp_all( tmp_all.dviant_side ~= tmp_all.q_position, :);
            tmp_incong_rt = tmp_incong.reaction_time; %Reaction time
            tmp_incong = table2array(tmp_incong(:,1:2)); % Overwrite to remove other columns

            [cong_alpha, cong_beta, cong_half_max_x, cong_amplitude, cong_pc_x,...
                cong_pc_y, cong_pd_x, cong_pd_y] = fitting_log(tmp_cong, log_win, plt);
            title("Congruent")

            [incong_alpha, incong_beta, incong_half_max_x, incong_amplitude, incong_pc_x,...
                incong_pc_y, incong_pd_x, incong_pd_y] = fitting_log(tmp_incong, log_win, plt);
            title("Incongruent")

            cong_log(ii).alpha = cong_alpha;
            cong_log(ii).beta = cong_beta;
            cong_log(ii).hlf_max = cong_half_max_x;
            cong_log(ii).amplitude = cong_amplitude;
            cong_log(ii).pc_x = cong_pc_x;
            cong_log(ii).pc_y = cong_pc_y;
            cong_log(ii).pd_x = cong_pd_x;
            cong_log(ii).pd_y = cong_pd_y;
            cong_log(ii).data = tmp_cong; % Psychometry data
            cong_log(ii).rt = tmp_cong_rt;
            cong_log(ii).rt2 = mean(tmp_cong_rt);


            incong_log(ii).alpha = incong_alpha;
            incong_log(ii).beta = incong_beta;
            incong_log(ii).hlf_max = incong_half_max_x;
            incong_log(ii).amplitude = incong_amplitude;
            incong_log(ii).pc_x = incong_pc_x;
            incong_log(ii).pc_y = incong_pc_y;
            incong_log(ii).pd_x = incong_pd_x;
            incong_log(ii).pd_y = incong_pd_y;
            incong_log(ii).data = tmp_incong; % Psychometry data
            incong_log(ii).rt = tmp_incong_rt;
            incong_log(ii).rt2 = mean(tmp_incong_rt);


            % ---------Extract response ratios

            tmp_ratio1 = [];
            tmp_ratio2 = [];

            for ll = 1:length(dvnts)
                tmp_ratio1 = [tmp_ratio1, mean( tmp_cong(tmp_cong(:,1)==dvnts(ll), 2) )];
                tmp_ratio2 = [tmp_ratio2, mean( tmp_incong(tmp_incong(:,1)==dvnts(ll), 2) )];

            end
            cong_resp_ratio = [cong_resp_ratio;tmp_ratio1];
            incong_resp_ratio = [incong_resp_ratio; tmp_ratio2];
        end
    end %End loop over participants
    response_ratios = array2table(response_ratios, 'VariableNames',string(dvnts));




end








