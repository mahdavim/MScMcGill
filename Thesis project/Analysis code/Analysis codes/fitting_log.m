function [alpha, beta, half_max, amplitude,pc_x,pc_y, pd_x, pd_y] = fitting_log(data, x_window, plt)

% x_window is the window for drawing the ligistic function

% Upper and lower asymptotes are calculated from the data NOT from the
% fited logistic function
% pc_x and pc_y are for plotting the main logistic fit
% pd_x and pd_y are the data points


global designmx

designmx = [ones(length(data),1),data];



guess=[-4.5,0.6];
guess=guess';
p = fminsearch(@efcn,guess);
beta=p(1);
alpha=p(2);


pd_x = unique(data(:,1));
pc_x = linspace((pd_x(1)-x_window)*1.1,(pd_x(end)+x_window*1.1),1000)'; %changed
design_pc = [ones(length(pc_x),1),pc_x];
pc_y = (1./( 1 + exp(-(design_pc*p))));

for j = 1 : length(pd_x),

   ind=find(data(:,1)==pd_x(j));
   pd_y(j) = mean(data(ind,2));

end;

%==============================================
% Define a finer grid for x values
% fine_grid_x = linspace(pd_x(1), pd_x(end), 1000);

% Find the x-value at which the predicted probability is 0.5 (half-max point)
% half_max_x = fminsearch(@(x) abs(interp1(fine_grid_x, logistic_function(x, design_pc, p), x) - 0.5), pd_x(1));
% =======================================

% Calculate the upper and lower asymptotes
PF_values = 1 ./ (1 + exp(-(designmx(:,1:end-1)*p)));
upper_asymptote = max(PF_values);
lower_asymptote = min(PF_values);
amplitude = upper_asymptote - lower_asymptote;
 half_max = -beta / alpha;

if plt
 figure
 plot(pd_x,pd_y,'b.',pc_x,pc_y,'r-')
 
 hold on
 
 % y=[0.51,0.6,0.7,0.8,0.9];
 % x=(log((1./(y))-1)+beta)./(-alpha);
 % xlim([-1  1])
 ylim([0 1])







% half_max = -alpha / beta;

% 
% % Plot the half-max point on the graph
plot(half_max, 0.5, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'm');
end
% 
% plot(x,y,'g*');
% disp(['Beta: ', num2str(beta)]);
% disp(['Alpha: ', num2str(alpha)]);