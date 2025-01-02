function [all_pval, p_vals] = perm_test_cong(data1, data2, log_vars1, log_vars2, n_perm)

% Log_vars 1 and 2 are results of logistic fit for each participant. The
% last column of them is number of valid trials for each participant
% log_pvals: permutation results from logistic variables
% p_vals: Permutation results from data


all_stats = table();
%------------- Loop to get t values for all logistic fit variables
for ii = 1:size(log_vars1,2) % Exclude number of trials;
    tmp_name = log_vars1.Properties.VariableNames{ii};
    [~,p,~, log_tstat] = ttest2(table2array( log_vars1(:, ii) ), table2array( log_vars2(:, ii) ));
    eval("all_stats." + tmp_name + "=log_tstat.tstat");


end


%------------- Permutation testing



perm_stats = [];
n_condition = 4; %Number of trials per condition in congruent or incongruent (e.g., -8, -6 etc)
uniq_condition = unique(data1(:,1)); %Number of unique conditions
data_pool = [data1; data2];

for ii = 1:n_perm
    disp("Iteration: "+ii)



    for jj = 1:size(log_vars1,1) % Number of subjects

        tmp_logs1 = [];
        tmp_logs2 = [];
        tmp_dat1 = [];
        tmp_dat2 = [];
        for kk = 1:length( uniq_condition )
            tmp_condition = uniq_condition(kk);
            tmp_pool = data_pool(data_pool(:,1) == tmp_condition, :);
            tmp_ind = randsample(size( tmp_pool, 1), n_condition);

            tmp_dat1 = [tmp_dat1; tmp_pool(tmp_ind, :)];
            tmp_ind = randsample(size( tmp_pool, 1), n_condition);

            tmp_dat2 = [tmp_dat2; tmp_pool(tmp_ind, :)];

        end

        %-------- fit logistics on two permuted datasets and store values
        [tmp_alph, tmp_beta, tmp_hlfmx, tmp_amp, ~, ~] = fitting_log(tmp_dat1, 50,0);
        tmp_logs1(jj,:) = [tmp_alph, tmp_beta, tmp_hlfmx, tmp_amp];
        [tmp_alph, tmp_beta, tmp_hlfmx, tmp_amp, ~, ~] = fitting_log(tmp_dat2, 50,0);
        tmp_logs2(jj,:) = [tmp_alph, tmp_beta, tmp_hlfmx, tmp_amp];
    end

    [~,~,~, tmp_tstat] =  ttest2(tmp_logs1, tmp_logs2);
    perm_stats = [perm_stats;tmp_tstat.tstat];

end


p_vals = [];
for ii= 1:size(all_stats,2)
    p_vals(ii) = mean( (abs(perm_stats(:,ii)) > abs(table2array(all_stats(1, ii))) ) );
    disp("P value for "+ all_stats.Properties.VariableNames{ii} + " is: "+ p_vals(ii))
    figure
    histogram(perm_stats(:,ii), 10)
    hold on
    xline(table2array(all_stats(1, ii)), 'r-')
    % fg = gcf;
    % fg.Name = all_stats.Properties.VariableNames{ii};
end
p_vals = array2table(p_vals,"VariableNames",log_vars1.Properties.VariableNames(1:end) );




%-------- Permute logistic fit features on all data---------

[alpha1,beta1,hlfmx1,amp1, ~, ~, ~, ~] = fitting_log(data1, 50,0);     

[alpha2,beta2,hlfmx2,amp2, ~, ~, ~, ~] = fitting_log(data2, 50,0);

ind = [alpha2 - alpha1, beta2 - beta1, hlfmx2 - hlfmx1, amp2 - amp1];

perm_stats = [];
uniq_condition = unique(data1(:,1)); %Number of unique conditions
data_pool = [data1; data2];

n_condition = ceil(size(data1,1) / length(uniq_condition));%Number of trials per condition in congruent or incongruent (e.g., -8, -6 etc)


all_ind = [];
for ii = 1:n_perm
    disp("Iteration: "+ii)



 % Number of subjects

        tmp_logs1 = [];
        tmp_logs2 = [];
        tmp_dat1 = [];
        tmp_dat2 = [];
        for kk = 1:length( uniq_condition )
            tmp_condition = uniq_condition(kk);
            tmp_pool = data_pool(data_pool(:,1) == tmp_condition, :);
            tmp_ind = randsample(size( tmp_pool, 1), n_condition);

            tmp_dat1 = [tmp_dat1; tmp_pool(tmp_ind, :)];
            tmp_ind = randsample(size( tmp_pool, 1), n_condition);

            tmp_dat2 = [tmp_dat2; tmp_pool(tmp_ind, :)];

        end

        %-------- fit logistics on two permuted datasets and store values
        [tmp_alph1, tmp_beta1, tmp_hlfmx1, tmp_amp1, ~, ~] = fitting_log(tmp_dat1, 50,0);
        [tmp_alph2, tmp_beta2, tmp_hlfmx2, tmp_amp2, ~, ~] = fitting_log(tmp_dat2, 50,0);
        % all_ind = [all_ind, tmp_alph1 - tmp_alph2];
        tmp_ind = [tmp_alph2-tmp_alph1, tmp_beta2 - tmp_beta1, tmp_hlfmx2-tmp_hlfmx1, tmp_amp2-tmp_amp1];
        all_ind = [all_ind; tmp_ind];
    

end

all_pval = mean(abs(all_ind) > abs(ind));





% Second permutation test
% Permute logistic regression features
% Fit classifiers for wm/nwm


