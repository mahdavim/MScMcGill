clear all
close all
clc

pth = 'C:\Users\mahdi\Desktop\Current Tray\Eyedata\Eyetracker analysis codes\edfImport-master\edfImport-master';
cd(pth)

addpath ('C:\Users\mahdi\Desktop\Current Tray\Eyedata\Eyetracker analysis codes\edfImport-master\edfmex')
pth = genpath('C:\Users\mahdi\Desktop\Current Tray\Eyedata\Eyetracker analysis codes\edfImport-master\edfImport-master');

addpath(pth)

% Path of data
pth = {'C:\Users\mahdi\Desktop\Current Tray\Eyedata\Eyetracker analysis codes\preCompiled_edfmex\edfmex',...
    'C:\Users\mahdi\Desktop\Current Tray\Eyedata\wm1_nk_2_2024_08_21_15_56.EDF', ...
    'C:\Users\mahdi\Desktop\Current Tray\Eyedata\wm1_nk_2_sample 1_2024-08-21_15h56.57.491.csv'};

wd = 'C:\Users\mahdi\Desktop\Current Tray\Eyedata\Mgs data';
items = dir(wd);
items = items(~ismember({items.name}, {'.', '..'})); % Remove . and ..

cue_num = 8; % Number of cue positions !!!!!! THIS SHOULD BE CHANGED FOR ADAPTATION
eye=2;  %1 left 2 right

dist = 116.8;

msg_lst = ["Cue Cue_ons", "Response Resp_ons"];% CHECK THIS WHENEVER THE RESPONSE IS EMPTY !!! CHANGE FOR ADAPTAION

%% ======== Analysis of threshold experiments

thr_pth = 'C:\Users\mahdi\Desktop\Current Tray\Eyedata\Adaptation data\Thresholding\All data';

[thr_dat, all_dat, all_log] = read_psychophys(thr_pth, 50,0);

fitting_log(thr_dat, 50,1)
hold on
xline(7.5)

%% ======Nwm
thr_pth = 'C:\Users\mahdi\Desktop\Current Tray\Eyedata\Adaptation data\Nwm_all';

[dat_nwm, all_dat_nwm, all_log_nwm, log_struc_nwm, response_ratio_nwm] = read_psychophys(thr_pth, 50,0);

fitting_log(dat_nwm, 50, 0)
hold on
xline(7.5)



%% =========WM

pth = 'C:\Users\mahdi\Desktop\Current Tray\Eyedata\Adaptation data\Wm_all';

[dat_wm, all_dat_wm, all_log_wm, log_struc_wm, response_ratio_wm] = read_psychophys(pth, 50,0);

fitting_log(dat_wm, 50, 1)


%% WM and NWM





% -------------- Plot population and individual psychometry functions
figure
plot(x_wm, y_wm, 'k-')

hold on
for ii=1:length(log_struc_wm)

    plot(log_struc_wm(ii).pc_x, log_struc_wm(ii).pc_y, '--', "Color",[0 0 0 0.1])
    % plot(log_struc_wm(ii).pd_x, log_struc_wm(ii).pd_y, 'o', "MarkerEdgeColor",[0.9 0.9 0.9])

end
% plot(xx_wm, yy_wm, 'ko')

xlim([-9, 9])
ylim([0, 1])
legend(["Pooled data", "Individual"], "Location","southeast")
title("Experiment 1 (WM)")
set(gca, "box", "off")
xlabel("Test slope (Oct)")
ylabel("Ascending response ratio")


figure
plot(x_nwm, y_nwm, 'k-')
hold on
for ii=1:length(log_struc_nwm)

    plot(log_struc_nwm(ii).pc_x, log_struc_nwm(ii).pc_y, '--', "Color",[0 0 0 0.1])
end
xlim([-9, 9])
ylim([0, 1])
legend(["Pooled data", "Individual"], "Location","southeast")
title("Experiment 1 (NWM)")
set(gca, "box", "off")
xlabel("Test slope (Oct)")
ylabel("Ascending response ratio")

figure
hold on
for ii=1:length(log_struc_wm2)

    plot(log_struc_wm2(ii).pc_x, log_struc_wm2(ii).pc_y, '--', "Color",[0 0 0 0.1])
    % plot(log_struc_wm(ii).pd_x, log_struc_wm(ii).pd_y, 'o', "MarkerEdgeColor",[0.9 0.9 0.9])

end
plot(x_wm2, y_wm2, 'k-')
% plot(xx_wm, yy_wm, 'ko')

xlim([-13, 13])
ylim([0, 1])
legend(["Pooled data", "Individual"], "Location","southeast", "Box", "off")
title("Experiment 2 (WM)")
set(gca, "box", "off")
xlabel("Test slope (Oct)")
ylabel("Ascending response ratio")


figure
hold on
for ii=1:length(log_struc_nwm2)

    plot(log_struc_nwm2(ii).pc_x, log_struc_nwm2(ii).pc_y, '--', "Color",[0 0 0 0.1])

end
plot(x_nwm2, y_nwm2, 'k-')

xlim([-13, 13])
ylim([0, 1])
legend(["Pooled data", "Individual"], "Location","southeast", "Box","off")
title("Experiment 2 (NWM)")
set(gca, "box", "off")
xlabel("Test slope (Oct)")
ylabel("Ascending response ratio")


%--------------------


for ii=1:length(log_struc_wm2)
    figure
    plot(log_struc_wm2(ii).pc_x, log_struc_wm2(ii).pc_y, 'b--')
    hold on
    plot(log_struc_nwm2(ii).pc_x, log_struc_nwm2(ii).pc_y, 'g--')
    plot(log_struc_wm2(ii).pd_x, log_struc_wm2(ii).pd_y, 'bo')
    plot(log_struc_nwm2(ii).pd_x, log_struc_nwm2(ii).pd_y, 'gs')


    legend(["WM", "NWM"])
    xlim([-13, 13])
    xlabel("Test sweep slope (octave)")
    ylabel("Ascending response ratio")
    title("Experiment 2")
end
set(gca, "box", "off")


%% Permutation testing  (WM and NWM)
[log_pvars, dat_p] = perm_test(dat_nwm, dat_wm, all_log_nwm, all_log_wm, 1000);

disp("NWM:")
nwm_mean = mean(all_log_nwm);
% iqr(table2array(all_log_nwm(1:end-1,:)))
nwm_sem = std(all_log_nwm)./sqrt(size(all_log_nwm,1));

disp("WM:")
wm_mean = mean(all_log_wm);
wm_sem = std(all_log_wm)./sqrt(size(all_log_wm,1));


%------ Plot mean + Error bar of alpha values
figure
errorbar( table2array(nwm_mean(:,1)), table2array(nwm_sem(:,1)), 'bo')
hold on
errorbar( 2, table2array(wm_mean(:,1)), table2array(wm_sem(:,1)), 'ro')

scatter(ones(size(all_log_nwm,1), 1)*1.1, table2array(all_log_nwm(:,1)),...
    'ko',"MarkerEdgeColor", [0.7,0.7,0.7], "MarkerEdgeAlpha", 0.2)
scatter(ones(size(all_log_wm,1), 1)*2.1, table2array(all_log_wm(:,1)), 'ko',...
    "MarkerEdgeColor", [0.7,0.7,0.7], "MarkerEdgeAlpha", 0.2)

plot([1,2], [table2array(nwm_mean(:,1)), table2array(wm_mean(:,1))], 'k--', "LineWidth",0.7)

xlim([0, 3])
ylim([0, 1.2])
title("Experiment 1 alpha values")
set(gca, "Box", "off")
% legend(["Nwm", "WM"], "Location","southeast", "Box","off")
xticks([1,2])
xticklabels(["NWM", "WM"])
ylabel("Alpha value")
xlabel("Working memory group")

%----- plot amplitude values

figure
errorbar( table2array(nwm_mean(:,end-1)), table2array(nwm_sem(:,end-1)), 'bo')
hold on
errorbar( 2, table2array(wm_mean(:,end-1)), table2array(wm_sem(:,end-1)), 'go')

scatter(ones(size(all_log_nwm,1), 1)*1.1, table2array(all_log_nwm(:,end-1)),...
    'ko',"MarkerEdgeColor", [0.7,0.7,0.7], "MarkerEdgeAlpha", 0.3)
scatter(ones(size(all_log_wm,1), 1)*2.1, table2array(all_log_wm(:,end-1)), 'ko',...
    "MarkerEdgeColor", [0.7,0.7,0.7], "MarkerEdgeAlpha", 0.3)

plot([1,2], [table2array(nwm_mean(:,end-1)), table2array(wm_mean(:,end-1))], 'k--', "LineWidth",0.7)

xlim([0, 3])
% ylim([-3.1, 2])
title("Experiment 2 beta values")
% legend(["Nwm", "WM"], "Location","southeast", "Box","off")
xticks([1,2])
xticklabels(["NWM", "WM"])
ylabel("Beta value")
xlabel("Working memory group")
% --------- Half max -----------
figure
errorbar( table2array(nwm_mean(:,3)), table2array(nwm_sem(:,3)), 'bo')
hold on
errorbar( 2, table2array(wm_mean(:,3)), table2array(wm_sem(:,3)), 'ro')

scatter(ones(size(all_log_nwm,1), 1)*1.1, table2array(all_log_nwm(:,3)),...
    'ko',"MarkerEdgeColor", [0.7,0.7,0.7], "MarkerEdgeAlpha", 0.2)
scatter(ones(size(all_log_wm,1), 1)*2.1, table2array(all_log_wm(:,3)), 'ko',...
    "MarkerEdgeColor", [0.7,0.7,0.7], "MarkerEdgeAlpha", 0.2)

plot([1,2], [table2array(nwm_mean(:,3)), table2array(wm_mean(:,3))], 'k--', "LineWidth",0.7)

xlim([0, 3])
ylim([-4.5, 4.5])
title("Experiment 1 half-max values")
set(gca, "Box", "Off")
% legend(["Nwm", "WM"], "Location","southeast", "Box","off")
xticks([1,2])
xticklabels(["NWM", "WM"])
ylabel("Half max value")
xlabel("Working memory group")
%% =========Nwm 2

pth = 'C:\Users\mahdi\Desktop\Current Tray\Eyedata\Adaptation data\Nwm2_all';

[dat_nwm2, all_dat_nwm2, all_log_nwm2, log_struc_nwm2, response_ratio_nwm2] = read_psychophys(pth, 50,0);




%% =======Wm 2

pth = 'C:\Users\mahdi\Desktop\Current Tray\Eyedata\Adaptation data\Wm2_all';

[dat_wm2, all_dat_wm2, all_log_wm2, log_struc_wm2,response_ratio_wm2,...
    cong_wm2, incong_wm2, cong_resp_ratio, incong_resp_ratio] = read_psychophys(pth, 50,0);

fitting_log(dat_wm2, 50,0)

%%
cong_wm2_log = table();
incong_wm2_log = cong_wm2_log;

cong_wm2_log.alpha = [cong_wm2.alpha]';
cong_wm2_log.beta = [cong_wm2.beta]';
cong_wm2_log.hlf_max = [cong_wm2.hlf_max]';
cong_wm2_log.amplitude = [cong_wm2.amplitude]';
cong_alldata = vertcat(cong_wm2.data);

incong_wm2_log.alpha = [incong_wm2.alpha]';
incong_wm2_log.beta = [incong_wm2.beta]';
incong_wm2_log.hlf_max = [incong_wm2.hlf_max]';
incong_wm2_log.amplitude = [incong_wm2.amplitude]';
incong_alldata = vertcat(incong_wm2.data);

%% Analysis of wm2 and nwm2


[log_pvars, dat_p] = perm_test(dat_nwm2, dat_wm2, all_log_nwm2, all_log_wm2, 1000);

nwm2_mean = trimmean(table2array(all_log_nwm2), 10);




nwm2_sem = std(all_log_nwm2)./sqrt(size(all_log_nwm2,1));

wm2_mean = trimmean(table2array(all_log_wm2), 10);


wm2_sem = std(all_log_wm2)./sqrt(size(all_log_wm2,1));


alph_sem_wm2 =  std( all_log_wm2( table2array( all_log_wm2(:,1) ~= max(all_log_wm2(:,1)) ), 1))./sqrt(size(all_log_wm2,1));

%----------- Plot mean and std values

figure
errorbar( nwm2_mean(:,1), table2array(nwm2_sem(:,1)), 'bo')
hold on
% errorbar( 2, table2array(wm2_mean(:,1)), table2array(wm2_sem(:,1)), 'go')
errorbar( 2, wm2_mean(:,1), table2array(alph_sem_wm2), 'ro')


scatter(ones(size(all_log_nwm2,1), 1)*1.1, table2array(all_log_nwm2(:,1)),...
    'ko',"MarkerEdgeColor", [0.7,0.7,0.7], "MarkerEdgeAlpha", 0.2)
scatter(ones(size(all_log_wm2,1), 1)*2.1, table2array(all_log_wm2(:,1)), 'ko',...
    "MarkerEdgeColor", [0.7,0.7,0.7], "MarkerEdgeAlpha", 0.2)

plot([1,2], [nwm2_mean(:,1), wm2_mean(:,1)], 'k--', "LineWidth",0.7)

xlim([0, 3])
ylim([0, 0.55])
title("Experiment 2 alpha values")
set(gca, "Box", "Off")
% legend(["Nwm", "WM"], "Location","southeast", "Box","off")
xticks([1,2])
xticklabels(["NWM", "WM"])
ylabel("Alpha value")
xlabel("Working memory group")

% -------- Half max
figure
errorbar( nwm2_mean(:,3), table2array(nwm2_sem(:,3)), 'bo')
hold on
errorbar( 2, wm2_mean(:,3), table2array(wm2_sem(:,3)), 'ro')

scatter(ones(size(all_log_nwm2,1), 1)*1.1, table2array(all_log_nwm2(:,3)),...
    'ko',"MarkerEdgeColor", [0.7,0.7,0.7], "MarkerEdgeAlpha", 0.2)
scatter(ones(size(all_log_wm2,1), 1)*2.1, table2array(all_log_wm2(:,3)), 'ko',...
    "MarkerEdgeColor", [0.7,0.7,0.7], "MarkerEdgeAlpha", 0.2)

plot([1,2], [nwm2_mean(:,3), wm2_mean(:,3)], 'k--', "LineWidth",0.7)

xlim([0, 3])
ylim([-5.2, 5.2])
title("Experiment 2 half max values")
% legend(["Nwm", "WM"], "Location","southeast", "Box","off")
set(gca, "Box", "Off")
xticks([1,2])
xticklabels(["NWM", "WM"])
ylabel("Half max value")
xlabel("Working memory group")




%% ========= See individual psychometris to find good examples

for ii=1:length(log_struc_wm)
    figure
    plot(log_struc_wm(ii).pc_x, log_struc_wm(ii).pc_y, 'b--')
    hold on
    plot(log_struc_nwm(ii).pc_x, log_struc_nwm(ii).pc_y, 'g--')
    plot(log_struc_wm(ii).pd_x, log_struc_wm(ii).pd_y, 'bo')
    plot(log_struc_nwm(ii).pd_x, log_struc_nwm(ii).pd_y, 'gs')


    legend(["WM", "NWM"])
    xlim([-8, 8])
    xlabel("Test sweep slope (octave)")
    ylabel("Ascending response ratio")
    title("Experiment 1")
end


for ii=1:length(log_struc_wm2)
    figure
    plot(log_struc_wm2(ii).pc_x, log_struc_wm2(ii).pc_y, 'b--')
    hold on
    plot(log_struc_nwm2(ii).pc_x, log_struc_nwm2(ii).pc_y, 'g--')
    plot(log_struc_wm2(ii).pd_x, log_struc_wm2(ii).pd_y, 'bo')
    plot(log_struc_nwm2(ii).pd_x, log_struc_nwm2(ii).pd_y, 'gs')


    legend(["WM", "NWM"])
    xlim([-8, 8])
    xlabel("Test sweep slope (octave)")
    ylabel("Ascending response ratio")
    title("Experiment 2")
end


%% ====== Effect of adaptor trains

% Separate 0 slopes
dat_nwm_zero = all_dat_nwm( all_dat_nwm.dviant==0, :);
dat_wm_zero = all_dat_wm( all_dat_wm.dviant==0, :);
dat_nwm2_zero = all_dat_nwm2( all_dat_nwm2.dviant==0, :);
dat_wm2_zero = all_dat_wm2( all_dat_wm2.dviant==0, :);

uniq_adaptor = unique(dat_nwm_zero.adaptor_train);

all_up = struct();


for ii=1:length(uniq_adaptor)
    up_ratio = mean( dat_nwm_zero(dat_nwm_zero.adaptor_train==uniq_adaptor(ii), :).response );
    all_up(ii).nwm = up_ratio;
    up_ratio = mean( dat_wm_zero(dat_wm_zero.adaptor_train==uniq_adaptor(ii), :).response );
    all_up(ii).wm = up_ratio;

    up_ratio = mean( dat_nwm2_zero(dat_nwm2_zero.adaptor_train==uniq_adaptor(ii), :).response );
    all_up(ii).nwm2 = up_ratio;
    up_ratio = mean( dat_wm2_zero(dat_wm2_zero.adaptor_train==uniq_adaptor(ii), :).response );
    all_up(ii).wm2 = up_ratio;
end



figure
plot(uniq_adaptor, [all_up.nwm], 'gs--')
hold on
plot(uniq_adaptor, [all_up.wm], 'bo--')
legend(["NWM", "WM"])
title("Experiment 1")
xlim([9,14])
ylabel("Ascending response ratio")
xlabel("Adaptor train length")

figure
plot(uniq_adaptor, [all_up.nwm2], 'gs--')
hold on
plot(uniq_adaptor, [all_up.wm2], 'bo--')
legend(["NWM", "WM"])
title("Experiment 2")
ylabel("Ascending response ratio")
xlabel("Adaptor train length")
xlim([9,14])



%% Effect of congruency



[all_pvals, cong_p_vals] = perm_test_cong(cong_alldata,...
    incong_alldata, cong_wm2_log, incong_wm2_log, 1000);


cong_mean = trimmean( table2array(cong_wm2_log),0);
cong_mean_alpha = median (cong_wm2_log.alpha);
% cong_mean_alpha = cong_mean_alpha(cong_mean_alpha~=max(cong_mean_alpha));
% cong_mean_alpha = cong_mean_alpha(cong_mean_alpha~=max(cong_mean_alpha));
% cong_mean_alpha = mean(cong_mean_alpha) - 0.6;
% iqr(table2array(all_log_nwm(1:end-1,:)))

cong_sem = std(cong_wm2_log)./sqrt(size(cong_wm2_log,1));

incong_mean = trimmean( table2array(incong_wm2_log),10);

median(incong_wm2_log)

incong_sem = std(incong_wm2_log)./sqrt(size(incong_wm2_log,1));

alpha_sem_incong = std (incong_wm2_log (incong_wm2_log.alpha ~= max(incong_wm2_log.alpha),1))./sqrt(size(incong_wm2_log,1));
alpha_sem_cong = cong_wm2_log (cong_wm2_log.alpha ~= max(cong_wm2_log.alpha),1);
alpha_sem_cong = alpha_sem_cong( alpha_sem_cong.alpha ~= max(alpha_sem_cong.alpha),1);
alpha_sem_cong = std(alpha_sem_cong)./sqrt(size(incong_wm2_log,1));

%----------- Plot alpha values
figure
errorbar(1, cong_mean_alpha, table2array(alpha_sem_cong), 'bo')
hold on
errorbar(2,incong_mean(:,1), table2array(alpha_sem_incong), 'go')

scatter(ones(size(all_log_nwm2,1), 1)*1.1, table2array(all_log_nwm2(:,3)),...
    'ko',"MarkerEdgeColor", [0.7,0.7,0.7], "MarkerEdgeAlpha", 0.3)
scatter(ones(size(all_log_wm2,1), 1)*2.1, table2array(all_log_wm2(:,3)), 'ko',...
    "MarkerEdgeColor", [0.7,0.7,0.7], "MarkerEdgeAlpha", 0.3)

plot([1,2], [nwm2_mean(:,3), wm2_mean(:,3)], 'k--', "LineWidth",0.7)

xlim([0, 3])
% ylim([-4, 5.1])
title("Experiment 2 half max values")
legend(["Nwm", "WM"], "Location","southeast", "Box","off")
xticks([1,2])
xticklabels(["NWM", "WM"])
ylabel("Half max value")
xlabel("Working memory group")

%% Plot congruency psychometries




[aa,bb,hh,dd, x_wm, y_wm, xx_wm, yy_wm] = fitting_log(cong_alldata, 50,0);

[aaa,~,hhh,~, x_nwm, y_nwm, xx_nwm, yy_nwm] = fitting_log(incong_alldata, 50,0);

figure
plot(x_wm, y_wm, 'r-', x_nwm, y_nwm, 'b-')
hold on
plot(xx_wm, yy_wm, 'ro', xx_nwm, yy_nwm, 'bo')
xlim([-13, 13])
ylim([0, 1])
plot([hh, hh], [0, 0.5], 'k--')
plot([hhh, hhh] , [0, 0.5], 'k--')
plot(hh, 0.5, 'go')
plot(hhh, 0.5, 'go')
plot([hhh, hh], [0.5, 0.5], 'k--')

legend(["Congruent", "Incongruent"], "Location","southeast", "Box","off")
title("Experiment 2")

% -------------- Plot population and individual psychometry functions
figure
hold on
for ii=1:length(incong_wm2)

    plot(incong_wm2(ii).pc_x, incong_wm2(ii).pc_y, '--', "Color",[0 0 0 0.1])
    % plot(log_struc_wm(ii).pd_x, log_struc_wm(ii).pd_y, 'o', "MarkerEdgeColor",[0.9 0.9 0.9])

end
plot(x_wm, y_wm, 'k-')
% plot(xx_wm, yy_wm, 'ko')

xlim([-13, 13])
ylim([0, 1])
legend("Incongruent")


figure
hold on
for ii=1:length(cong_wm2)

    plot(cong_wm2(ii).pc_x, cong_wm2(ii).pc_y, '--', "Color",[0 0 0 0.1])
end
plot(x_nwm, y_nwm, 'k-')
xlim([-13, 13])
ylim([0, 1])
legend("Congruent")







%% Mixed effect models WM vs NWM (1)

all_dat_nwm
all_dat_wm

all_dat_exp1 = [all_dat_nwm(:,1:end-1); all_dat_wm(:, [1,2,3,5,6])];

glme_exp1 = fitglme(all_dat_exp1, ...
    'response ~ dviant + wmGroup + adaptor_train + (1|ID)', ...
    'Distribution', 'Binomial', ...
    'Link', 'logit');

%% Mixed effect models WM vs NWM (2)


all_dat_exp2= [all_dat_nwm2(:,1:end-2); all_dat_wm2(:, [1,2,3,5,6])];

glme_exp2 = fitglme(all_dat_exp2, ...
    'response ~ dviant + wmGroup + adaptor_train  + (1|ID)', ...
    'Distribution', 'Binomial', ...
    'Link', 'logit');

%% Mixed effect model for congruency


all_dat_me =  all_dat_wm2(:, [1,2,3,6, 7]);



glme_cong = fitglme(all_dat_me, ...
    'response ~ dviant  + adaptor_train  + congGroup + (1|ID)', ...
    'Distribution', 'Binomial', ...
    'Link', 'logit');

%%

all_dat_me_cong = all_dat_me(all_dat_me.congGroup==1,:);
all_dat_me_cong.Properties.VariableNames(end)= "wmGroup";
all_dat_me_cong = [all_dat_me_cong;all_dat_nwm2(:,1:5)];

glme_congNwm = fitglme(all_dat_me_cong, ...
    'response ~ dviant  + adaptor_train  + wmGroup + (1|ID)', ...
    'Distribution', 'Binomial', ...
    'Link', 'logit');

all_dat_me_incong = all_dat_me(all_dat_me.congGroup==0,:);
all_dat_me_incong.Properties.VariableNames(end)= "wmGroup";
all_dat_me_incong.wmGroup = ones(size(all_dat_me_incong.dviant));

all_dat_me_incong = [all_dat_me_incong;all_dat_nwm2(:,1:5)];

glme_incongNwm = fitglme(all_dat_me_incong, ...
    'response ~ dviant  + adaptor_train  + wmGroup + (1|ID)', ...
    'Distribution', 'Binomial', ...
    'Link', 'logit');



%% Box plots (EXP1)

% !!ATENTION: COLOR FOR WM: c2; COLOR FOR NWM: c1

all_log_nwm.wmGroup = zeros(size(all_log_nwm.alpha));
all_log_wm.wmGroup = ones(size(all_log_nwm.alpha));

figure
% boxchart(all_log.wmGroup, all_log.alpha, "cdata", "r")
h1 = boxchart(all_log_nwm.wmGroup,  all_log_nwm.half_max_x);
hold on
h2 = boxchart(all_log_wm.wmGroup,all_log_wm.half_max_x, "MarkerColor",'w');
c1 = h1.BoxFaceColor;
c2 = h2.BoxFaceColor;

scatter(zeros(size(all_log_nwm.half_max_x)) + (2*rand(size(all_log_nwm.alpha))-1)/30 ,...
    all_log_nwm.half_max_x, "MarkerEdgeColor",  c1, "MarkerEdgeAlpha",0.3)
scatter(ones(size(all_log_wm.half_max_x))  + (2* rand(size(all_log_wm.alpha)) -1)/30,...
    all_log_wm.half_max_x,"MarkerEdgeColor", c2,  "MarkerEdgeAlpha",0.2)
xticks([0,1])
xticklabels(["No WM", ["WM"]])
ylabel("PSE value in octave (\times96)")
save2pdf("C:\Users\mahdi\OneDrive - McGill University\Work and studies\Papers for thesis\Committee meetings\Final MSc\MSC Thesis\Figures\Box_pse_exp1",...
    gcf, 150)

figure
% boxchart(all_log.wmGroup, all_log.alpha, "cdata", "r")
h1 = boxchart(all_log_nwm.wmGroup,  all_log_nwm.alpha);
hold on
h2 = boxchart(all_log_wm.wmGroup,all_log_wm.alpha, "MarkerColor",'w');
c1 = h1.BoxFaceColor;
c2 = h2.BoxFaceColor;

scatter(zeros(size(all_log_nwm.alpha)) + (2*rand(size(all_log_nwm.alpha))-1)/30 ,...
    all_log_nwm.alpha, "MarkerEdgeColor",  c1, "MarkerEdgeAlpha",0.2)
scatter(ones(size(all_log_wm.alpha))  + (2* rand(size(all_log_wm.alpha)) -1)/30,...
    all_log_wm.alpha,"MarkerEdgeColor", c2,  "MarkerEdgeAlpha",0.2)
xticks([0,1])
xticklabels(["No WM", ["WM"]])
ylabel("Alpha value")
ylim([0,2.5])
save2pdf("C:\Users\mahdi\OneDrive - McGill University\Work and studies\Papers for thesis\Committee meetings\Final MSc\MSC Thesis\Figures\Box_alph_exp1",...
    gcf, 150)

%% Box plots (Exp 2)
all_log_nwm2.wmGroup = zeros(size(all_log_nwm2.alpha));
all_log_wm2.wmGroup = ones(size(all_log_nwm2.alpha));

figure
% boxchart(all_log.wmGroup, all_log.alpha, "cdata", "r")
h1 = boxchart(all_log_nwm2.wmGroup,  all_log_nwm2.half_max_x);
hold on
h2 = boxchart(all_log_wm2.wmGroup,all_log_wm2.half_max_x, "MarkerColor",'w');
c1 = h1.BoxFaceColor;
c2 = h2.BoxFaceColor;

scatter(zeros(size(all_log_nwm2.half_max_x)) + (2*rand(size(all_log_nwm2.alpha))-1)/30 ,...
    all_log_nwm2.half_max_x, "MarkerEdgeColor",  c1, "MarkerEdgeAlpha",0.2)
scatter(ones(size(all_log_wm2.half_max_x))  + (2* rand(size(all_log_wm2.alpha)) -1)/30,...
    all_log_wm2.half_max_x,"MarkerEdgeColor", c2,  "MarkerEdgeAlpha",0.2)
xticks([0,1])
xticklabels(["No WM", ["WM"]])
ylabel("PSE value in octave (\times96)")
save2pdf("C:\Users\mahdi\OneDrive - McGill University\Work and studies\Papers for thesis\Committee meetings\Final MSc\MSC Thesis\Figures\Box_pse_exp2",...
    gcf, 150)

figure
% boxchart(all_log.wmGroup, all_log.alpha, "cdata", "r")
h1 = boxchart(all_log_nwm2.wmGroup,  all_log_nwm2.alpha);
hold on
h2 = boxchart(all_log_wm2.wmGroup,all_log_wm2.alpha, "MarkerColor",'w');
c1 = h1.BoxFaceColor;
c2 = h2.BoxFaceColor;

scatter(zeros(size(all_log_nwm2.alpha)) + (2*rand(size(all_log_nwm2.alpha))-1)/30 ,...
    all_log_nwm2.alpha, "MarkerEdgeColor",  c1, "MarkerEdgeAlpha",0.2)
scatter(ones(size(all_log_wm2.alpha))  + (2* rand(size(all_log_wm2.alpha)) -1)/30,...
    all_log_wm2.alpha,"MarkerEdgeColor", c2,  "MarkerEdgeAlpha",0.2)
xticks([0,1])
xticklabels(["No WM", ["WM"]])
ylabel("Alpha value")
ylim([0,1])
save2pdf("C:\Users\mahdi\OneDrive - McGill University\Work and studies\Papers for thesis\Committee meetings\Final MSc\MSC Thesis\Figures\Box_alph_exp2",...
    gcf, 150)





%% Response ratios
figure
errorbar( table2array( nanmean(response_ratio_nwm) ),...
    table2array( std(response_ratio_nwm,0, 'omitmissing') ...
    ./sqrt( size(response_ratio_nwm,1) ) ) , 'bo', "Color", c1)


hold on
errorbar( table2array( nanmean(response_ratio_wm) ),...
    table2array( std(response_ratio_wm,0, 'omitmissing') ...
    ./sqrt( size(response_ratio_wm,1) ) ) , 'ro', "Color",  c2)
legend(["NWM", "WM"], "Location","southeast", "Box","off")
set(gca, "Box", "Off")
xlabel("Tone sweep slope in octave")
ylabel("Ascending respone ratio")
title("Experiment 1 response ratios")
xticks([1,2,3,4,5,6,7,8,9])
xticklabels(["-8/96","-6/96","-4/96","-2/96","0","2/96","4/96","6/96","8/96"])
xtickangle(90);
xlim([0, 10])
ax = gca;
ax.XAxis.TickLabelInterpreter = 'tex'; % Ensures compatibility with position adjustment
ax.XAxis.TickLabelGapOffset = 7; % Increase gap in points (experiment for the best fit)
save2pdf("C:\Users\mahdi\OneDrive - McGill University\Work and studies\Papers for thesis\Committee meetings\Final MSc\MSC Thesis\Figures\resp_ratio_exp1")




figure
errorbar( table2array( nanmean(response_ratio_nwm2) ),...
    table2array( std(response_ratio_nwm2,0, 'omitmissing') ...
    ./sqrt( size(response_ratio_nwm2,1) ) ) , 'bo', "Color",c1)


hold on
errorbar( table2array( nanmean(response_ratio_wm2) ),...
    table2array( std(response_ratio_wm2,0, 'omitmissing') ...
    ./sqrt( size(response_ratio_wm2,1) ) ) , 'ro', "Color",c2)
legend(["NWM", "WM"], "Location","southeast", "Box","off")
set(gca, "Box", "Off")

xlabel("Tone sweep slope in octave")
ylabel("Ascending respone ratio")
title("Experiment 2 response ratios")
xticks([1,2,3,4,5,6,7,8,9])
xticklabels(["-12/96","-9/96","-6/96","-3/96","0","3/96","6/96","9/96","12/96"])
xtickangle(90);

xlim([0, 10])
ax = gca;
ax.XAxis.TickLabelInterpreter = 'tex'; % Ensures compatibility with position adjustment
ax.XAxis.TickLabelGapOffset = 7; % Increase gap in points (experiment for the best fit)
save2pdf("C:\Users\mahdi\OneDrive - McGill University\Work and studies\Papers for thesis\Committee meetings\Final MSc\MSC Thesis\Figures\resp_ratio_exp2")



% ---------Congruency

figure
errorbar( nanmean(cong_resp_ratio) ,...
    std(cong_resp_ratio,0, 'omitmissing') ...
    ./sqrt( size(cong_resp_ratio,1)  ), 'ro', "Color",c2)


hold on
errorbar( nanmean(incong_resp_ratio) ,...
    std(incong_resp_ratio,0, 'omitmissing') ...
    ./sqrt( size(incong_resp_ratio,1) ) , 'bo', "Color",c1)

legend(["Congruent", "Incongruent"], "Location","southeast", "Box","off")
set(gca, "Box", "Off")

xlabel("Tone sweep slope in octave")
ylabel("Ascending respone ratio")
title("Congruent vs incongruent response ratios")
xticks([1,2,3,4,5,6,7,8,9])
xticklabels(["-12/96","-9/96","-6/96","-3/96","0","3/96","6/96","9/96","12/96"])
xtickangle(90);

xlim([0, 10])
ax = gca;
ax.XAxis.TickLabelInterpreter = 'tex'; % Ensures compatibility with position adjustment
ax.XAxis.TickLabelGapOffset = 7; % Increase gap in points (experiment for the best fit)
save2pdf("C:\Users\mahdi\OneDrive - McGill University\Work and studies\Papers for thesis\Committee meetings\Final MSc\MSC Thesis\Figures\resp_ratio_cong")


%% Psychometry functions
%------------ Exp 1
[aa,bb,hh,dd, x_wm, y_wm, xx_wm, yy_wm] = fitting_log(dat_wm, 50,0);

[~,~,hhh,~, x_nwm, y_nwm, xx_nwm, yy_nwm] = fitting_log(dat_nwm, 50,0);

figure
h1 = plot(x_wm, y_wm, 'r-', "Color",c2)
hold on
h2 = plot( x_nwm, y_nwm, 'b-', "Color",c1);
% plot(xx_wm, yy_wm, 'ro', xx_nwm, yy_nwm, 'bo')
xlim([-8.5, 8.5])
ylim([0, 1])
plot([hh, hh], [0, 0.5], 'k--', "LineWidth",0.7)
plot([hhh, hhh] , [0, 0.5], 'k--',  "LineWidth",0.7)
plot(hh, 0.5, 'go', "Color",[0 0.5 0])
plot(hhh, 0.5, 'go', "Color",[0 0.5 0])
plot([hhh, hh], [0.5, 0.5], 'k--', "LineWidth",0.7)
set(gca, "box", "off")
% legend([h1, h2])
legend([h1, h2], ["WM", "NWM"], "Location","southeast", "Box","off")
title("Experiment 1")
ylabel("Ascending response ratio")
xlabel("Test sweep slope in ocate (\times96)")
save2pdf("C:\Users\mahdi\OneDrive - McGill University\Work and studies\Papers for thesis\Committee meetings\Final MSc\MSC Thesis\Figures\psychometry_exp1")





% ---------------- Exp 2
[~,~,hh,~, x_wm2, y_wm2, xx_wm2, yy_wm2] = fitting_log(dat_wm2, 50,0);
[~,~,hhh,~, x_nwm2, y_nwm2, xx_nwm2, yy_nwm2] = fitting_log(dat_nwm2, 50,0);


figure
plot(x_wm2, y_wm2, 'r-' , "Color", c2)
hold on 
plot(x_nwm2, y_nwm2, 'b-', "Color",c1)
hold on
% plot(xx_wm2, yy_wm2, 'ro', xx_nwm2, yy_nwm2, 'bo')
ylim([0, 1])
plot([hh, hh], [0, 0.5], 'k--', "LineWidth",0.7)
plot([hhh, hhh] , [0, 0.5], 'k--',  "LineWidth",0.7)
plot(hh, 0.5, 'go', "Color",[0 0.5 0])
plot(hhh, 0.5, 'go', "Color",[0 0.5 0])
plot([hhh, hh], [0.5, 0.5], 'k--', "LineWidth",0.7)
set(gca, "box", "off")
% legend([h1, h2])
legend(["WM", "NWM"], "Location","southeast", "Box","off")
title("Experiment 2")
ylabel("Ascending response ratio")
xlabel("Test sweep slope in octave (\times96)")
xticks([-12,-9,-6,-3,0,3,6,9,12]);


xlim([-12, 12])
ylim([0, 1])

save2pdf("C:\Users\mahdi\OneDrive - McGill University\Work and studies\Papers for thesis\Committee meetings\Final MSc\MSC Thesis\Figures\psychometry_exp2")

%% ================= Reaction time plots
nw_rt = [];
w_rt = [];
nw_sem = [];
w_sem = [];

nw2_rt = [];
w2_rt = [];
nw2_sem = [];
w2_sem = [];

cong_rt = [];
incong_rt = [];
cong_sem = [];
incong_sem = [];


% ---------- gather rts for congruent and incongruent
all_cong_rt = [];
all_incong_rt = [];

for ii=1:length(cong_wm2)
    all_cong_rt = [all_cong_rt, cong_wm2(ii).rt'];
    all_incong_rt = [all_incong_rt, incong_wm2(ii).rt'];
end



all_cong_rt = all_cong_rt';
all_incong_rt = all_incong_rt';

dvnts = unique(all_dat_nwm.dviant);
dvnts2 = unique(all_dat_nwm2.dviant);
for ii=1:length(dvnts)
    tmp_dvnt = dvnts(ii);
    tmp_rt = all_dat_nwm(all_dat_nwm.dviant==tmp_dvnt, :).reaction_time;
    nw_rt = [nw_rt;mean(tmp_rt)];
    nw_sem = [nw_sem; std(tmp_rt)/sqrt(size(all_log_nwm,1))];

    tmp_rt = all_dat_wm(all_dat_wm.dviant==tmp_dvnt, :).reaction_time;% Working memory
    w_rt = [w_rt;mean(tmp_rt)];
    w_sem = [w_sem; std(tmp_rt)./sqrt( size(all_log_wm,1) )];

    tmp_dvnt = dvnts2(ii);
    tmp_rt = all_dat_nwm2(all_dat_nwm2.dviant==tmp_dvnt, :).reaction_time;
    nw2_rt = [nw2_rt;mean(tmp_rt)];
    nw2_sem = [nw2_sem; std(tmp_rt)./sqrt(size(all_log_nwm2,1))];

    tmp_rt = all_dat_wm2(all_dat_wm2.dviant==tmp_dvnt, :).reaction_time;% Working memory
    w2_rt = [w2_rt;mean(tmp_rt)];
    w2_sem = [w2_sem; std(tmp_rt)./sqrt(size(all_log_wm2,1))];

    %---------- Congruent and incongruent

    tmp_dvnt = dvnts2(ii);
    tmp_rt = all_cong_rt(cong_alldata(:,1) ==tmp_dvnt, :);
    cong_rt = [cong_rt;mean(tmp_rt)];
    % cong_sem = [cong_sem; std(tmp_rt)./sqrt(size(all_log_wm2,1))];
    cong_sem = [cong_sem; std(tmp_rt)./sqrt(10)];


    tmp_rt = all_cong_rt(incong_alldata(:,1) ==tmp_dvnt, :);
    incong_rt = [incong_rt;mean(tmp_rt)];
    % incong_sem = [incong_sem; std(tmp_rt)./sqrt(size(all_log_wm2,1))];
    incong_sem = [incong_sem; std(tmp_rt)./sqrt(10)];


end



figure
errorbar(dvnts, nw_rt, nw_sem, 'bo--', "Color",c1)
hold on
errorbar(dvnts, w_rt, w_sem, 'go--', "Color",c2)
xlim([-10, 10])
ylim([0,2])
legend(["NWM", "WM"])
xlabel("Deviant slope")
ylabel("Reaction time (s)")
title("RT in Experiment 1")
set(gca, "box", "off")
save2pdf("C:\Users\mahdi\OneDrive - McGill University\Work and studies\Papers for thesis\Committee meetings\Final MSc\MSC Thesis\Figures\rt_exp1")


figure
errorbar(dvnts, nw2_rt, nw2_sem, 'bo--', "Color",c1)
hold on
errorbar(dvnts, w2_rt, w2_sem, 'go--', "Color",c2)
xlim([-10, 10])
ylim([0,2])
legend(["NWM", "WM"])
xlabel("Deviant slope")
ylabel("Reaction time (s)")
title("RT in Experiment 2")
set(gca, "box", "off")
save2pdf("C:\Users\mahdi\OneDrive - McGill University\Work and studies\Papers for thesis\Committee meetings\Final MSc\MSC Thesis\Figures\rt_exp2")





% ----------Congruent vs incongruent
figure
errorbar(dvnts, cong_rt, cong_sem, 'bo--', "Color",c2)
hold on
errorbar(dvnts, incong_rt, incong_sem, 'go--', "Color",c1)
xlim([-13, 13])
ylim([0,2])
legend(["Cong", "Incong"])
xlabel("Deviant slope")
ylabel("Reaction time (s)")
title("Congruent vs Incongruent")
set(gca, "box", "off")
save2pdf("C:\Users\mahdi\OneDrive - McGill University\Work and studies\Papers for thesis\Committee meetings\Final MSc\MSC Thesis\Figures\rt_cong")

%% Mixed effect for RT
all_dat_exp1 = [all_dat_nwm(:,[1,3,4,5,6]); all_dat_wm(:, [1,3,5,6, 7])];
all_dat_exp1.reaction_time(all_dat_exp1.reaction_time==min(all_dat_exp1.reaction_time)) = 0.78443;                    

glme_exp1_rt = fitglme(all_dat_exp1, ...
    'reaction_time ~ dviant + wmGroup + adaptor_train + (1|ID)', ...
    'Distribution', 'Gamma', ...
    'Link', 'log');

%%
all_dat_exp2= [all_dat_nwm2(:,[1,3,4,5,7]); all_dat_wm2(:, [1,3,5,6, 9])];
all_dat_exp2.reaction_time(all_dat_exp2.reaction_time==min(all_dat_exp2.reaction_time)) = 0.78443;                    

glme_exp2 = fitglme(all_dat_exp2, ...
    'reaction_time ~ dviant + wmGroup + adaptor_train  + (1|ID)', ...
    'Distribution', 'Gamma', ...
    'Link', 'log');
%% LMM for RT for congruency 

all_dat_me =  all_dat_wm2(:, [1,3,6, 7, 9]);



glme_cong = fitglme(all_dat_me, ...
    'reaction_time ~ dviant  + adaptor_train  + congGroup + (1|ID)', ...
    'Distribution', 'Gamma', ...
    'Link', 'log');
