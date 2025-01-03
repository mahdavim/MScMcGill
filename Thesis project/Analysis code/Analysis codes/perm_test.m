function [log_pvals, p_vals] = perm_test(data1, data2, log_vars1, log_vars2, n_perm)

% Log_vars 1 and 2 are results of logistic fit for each participant. The
% last column of them is number of valid trials for each participant
% log_pvals: permutation results from logistic variables
% p_vals: Permutation results from data


all_stats = table();
%------------- Loop to get t values for all logistic fit variables
for ii = 1:size(log_vars1,2)-1 % Exclude number of trials;
    tmp_name = log_vars1.Properties.VariableNames{ii};
    [~,p,~, log_tstat] = ttest2(table2array( log_vars1(:, ii) ), table2array( log_vars2(:, ii) ));
    disp(tmp_name + " " +p)
    eval("all_stats." + tmp_name + "=log_tstat.tstat");


end


%------------- Permutation testing



perm_stats = [];
n_condition = 7; %Number of trials per condition (e.g., -8, -6 etc)
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
p_vals = array2table(p_vals,"VariableNames",log_vars1.Properties.VariableNames(1:end-1) );


%-------- Permute logistic fit features---------
log_t = nan(n_perm, size(log_vars1,2)-1);
subj_n = size(log_vars1,1) + size(log_vars2,1);
for ii=1:n_perm
    disp("Iteration "+ ii)
    for jj=1:size(log_vars1,2)-1
        tmp_all = [ table2array(log_vars1(:,jj)); table2array(log_vars2(:,jj)) ];
        tmp_all = randsample(tmp_all , length(tmp_all));% Shuffle
        % tmp_ind1 = randsample(subj_n, size(log_vars1,1), 0);
        % tmp_log1 = tmp_all(tmp_ind1);
        % tmp_log2 = setdiff(tmp_logs1,tmp_all);
        tmp_log1 = randsample(tmp_all, size(log_vars1,1));
        tmp_log2 = randsample(tmp_all, size(log_vars2,1));

        [~,~,~, tmp_stats] = ttest2(tmp_log1, tmp_log2);
        log_t(ii, jj) =  tmp_stats.tstat;

    end
end


log_pvals = mean(abs(log_t) > abs(table2array(all_stats)));
log_pvals = array2table(log_pvals, "VariableNames",log_vars1.Properties.VariableNames(1:end-1));



% Second permutation test
% Permute logistic regression features
% Fit classifiers for wm/nwm


