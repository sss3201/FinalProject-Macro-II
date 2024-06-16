/*
 * Este codigo es una adaptación e inspiración del codigo realizado por: 
 * Dr. Johannes Pfeifer  
 * donde replica el paper de Hansen, Gary D. (1985): "Invisible labor and the business cycle", Journal 
 * of Monetary Economics 16, pp.309-327.
 * 
 */

title_string='Economy with indivisible labor'

var c $c$ (long_name='consumption')
    w $w$ (long_name='real wage')
    r $r$ (long_name='real interest rate')
    y $y$ (long_name='output')
    h $h$ (long_name='hours')
    k $k$ (long_name='capital stock')
    invest $i$ (long_name='investment')
    lambda $\lambda$ (long_name='TFP')
    productivity ${\frac{y}{h}}$ (long_name='Productivity');
    
varexo eps_a;

parameters beta $\beta$ (long_name='discount factor')
    delta $\delta$ (long_name='depreciation rate')
    theta $\theta$ (long_name='capital share')
    gamma $\gamma$ (long_name='AR coefficient TFP')
    A $A$ (long_name='labor disutility parameter')
    h_0 ${h_0}$ (long_name='full time workers in steady state')
    sigma_eps $\sigma_e$ (long_name='TFP shock volatility')
    B $B$ (long_name='composite labor disutility parameter')
    ;

//Calibration, p. 319
beta = 0.99;
delta = 0.025;
theta = 0.36;
gamma = 0.95;
A = 2;
sigma_eps=0.00712;
h_0=0.53;

model;
//1. Euler Equation
1/c = beta*((1/c(+1))*(r(+1) +(1-delta)));
//2. Labor FOC
(1-theta)*(y/h) = B*c;
//3. Resource constraint
c = y +(1-delta)*k(-1) - k;
//4. LOM capital
k= (1-delta)*k(-1) + invest;
//5. Production function
y = lambda*k(-1)^(theta)*h^(1-theta);
//6. Real wage
r = theta*(y/k(-1));
//7. Real interest rate
w = (1-theta)*(y/h);
//8. LOM TFP
log(lambda)=gamma*log(lambda(-1))+eps_a;
//9. Productivity
productivity= y/h;
end;

steady_state_model;
//Follows footnote 15
B=-A*(log(1-h_0))/h_0; %called psi in footnote
lambda = 1;
h = (1-theta)*(1/beta -(1-delta))/(B*(1/beta -(1-delta)-theta*delta));
k = h*((1/beta -(1-delta))/(theta*lambda))^(1/(theta-1));
invest = delta*k;
y = lambda*k^(theta)*h^(1-theta);
c = y-delta*k;
r =  1/beta - (1-delta);
w = (1-theta)*(y/h);
productivity = y/h;
end;

steady;

shocks;
var eps_a; stderr sigma_eps;
end;

check;
steady;
stoch_simul(order=1,irf=20,loglinear,hp_filter=1600) y c invest k h productivity w;

stoch_simul(order=1,irf=20,loglinear,hp_filter=1600,simul_replic=10000,periods=200) y c invest h productivity w;

%read out simulations
simulated_series_raw=get_simul_replications(M_,options_);

%filter series
simulated_series_filtered=NaN(size(simulated_series_raw));
for ii=1:options_.simul_replic
    [trend, cycle]=sample_hp_filter(simulated_series_raw(:,:,ii)',1600);
    simulated_series_filtered(:,:,ii)=cycle';
end

%get variable positions
y_pos=strmatch('y',M_.endo_names,'exact');
c_pos=strmatch('c',M_.endo_names,'exact');
i_pos=strmatch('invest',M_.endo_names,'exact');
h_pos=strmatch('h',M_.endo_names,'exact');
productivity_pos=strmatch('productivity',M_.endo_names,'exact');
w_pos=strmatch('w',M_.endo_names,'exact');

var_positions=[y_pos; c_pos; i_pos; h_pos; productivity_pos];
var_names=M_.endo_names_long(var_positions,:);

%Compute standard deviations
std_mat=std(simulated_series_filtered(var_positions,:,:),0,2)*100;
std_h = std(simulated_series_filtered(h_pos,:,:),0,2)*100;
std_w = std(simulated_series_filtered(w_pos,:,:),0,2)*100;

%Compute relative standard deviations
rel_std_mat = std_mat ./ std_mat(1,:,:); % std(x) / std(y)
rel_std_h_w = std_h ./ std_w; % std(h) / std(w)

%Compute correlations
corr_mat = NaN(5, options_.simul_replic); % Adjust the size of corr_mat
for ii=1:options_.simul_replic
    corr_mat(1,ii)=corr(simulated_series_filtered(y_pos,:,ii)',simulated_series_filtered(y_pos,:,ii)');
    corr_mat(2,ii)=corr(simulated_series_filtered(y_pos,:,ii)',simulated_series_filtered(c_pos,:,ii)');
    corr_mat(3,ii)=corr(simulated_series_filtered(y_pos,:,ii)',simulated_series_filtered(i_pos,:,ii)');
    corr_mat(4,ii)=corr(simulated_series_filtered(y_pos,:,ii)',simulated_series_filtered(h_pos,:,ii)');
    corr_mat(5,ii)=corr(simulated_series_filtered(y_pos,:,ii)',simulated_series_filtered(productivity_pos,:,ii)');
end

%Compute correlation between h and w
corr_h_w = NaN(1, options_.simul_replic);
for ii=1:options_.simul_replic
    corr_h_w(1,ii) = corr(simulated_series_filtered(h_pos,:,ii)', simulated_series_filtered(w_pos,:,ii)');
end

%Print table with standard deviations and correlations
fprintf('\n%-40s \n',title_string)
fprintf('%-20s \t %11s \t %11s \n','','std(x)','corr(y,x)')
for ii=1:size(corr_mat,1)
    fprintf('%-20s \t %3.2f (%3.2f) \t %3.2f (%3.2f) \n', var_names{ii,:}, mean(std_mat(ii,:,:),3), std(std_mat(ii,:,:),0,3), mean(corr_mat(ii,:),2), std(corr_mat(ii,:),0,2))
end

%Print table with relative standard deviations
fprintf('\n%-40s \n','Relative Standard Deviations')
fprintf('%-20s \t %11s \n','','std(x)/std(y)')
for ii=2:size(corr_mat,1) % Start from 2 to skip y itself
    fprintf('%-20s \t %3.2f (%3.2f) \n', var_names{ii,:}, mean(rel_std_mat(ii,:,:),3), std(rel_std_mat(ii,:,:),0,3))
end

%Print table with correlation of h and w and relative std of h and w
fprintf('\n%-40s \n','Correlation of h and w, Relative Std of h and w')
fprintf('%-20s \t %11s \t %11s \n','Variable','corr(h,w)','std(h)/std(w)')
fprintf('%-20s \t %3.2f (%3.2f) \t %3.2f (%3.2f) \n', 'h, w', mean(corr_h_w,2), std(corr_h_w,0,2), mean(rel_std_h_w,3), std(rel_std_h_w,0,3))

% Generate histograms with mean lines
figure;

% Relative Std of Consumption
subplot(3,2,1);
histogram(rel_std_mat(2,:,:));
hold on;
mean_rel_std_c = mean(rel_std_mat(2,:,:), 'all');
xline(mean_rel_std_c, 'r', 'LineWidth', 2);
title('Relative Std of Consumption');
xlabel('std(c)/std(y)');
ylabel('Frequency');
hold off;

% Relative Std of Investment
subplot(3,2,2);
histogram(rel_std_mat(3,:,:));
hold on;
mean_rel_std_i = mean(rel_std_mat(3,:,:), 'all');
xline(mean_rel_std_i, 'r', 'LineWidth', 2);
title('Relative Std of Investment');
xlabel('std(i)/std(y)');
ylabel('Frequency');
hold off;

% Relative Std of Hours
subplot(3,2,3);
histogram(rel_std_mat(4,:,:));
hold on;
mean_rel_std_h = mean(rel_std_mat(4,:,:), 'all');
xline(mean_rel_std_h, 'r', 'LineWidth', 2);
title('Relative Std of Hours');
xlabel('std(h)/std(y)');
ylabel('Frequency');
hold off;

% Relative Std of Productivity
subplot(3,2,4);
histogram(rel_std_mat(5,:,:));
hold on;
mean_rel_std_productivity = mean(rel_std_mat(5,:,:), 'all');
xline(mean_rel_std_productivity, 'r', 'LineWidth', 2);
title('Relative Std of Productivity');
xlabel('std(productivity)/std(y)');
ylabel('Frequency');
hold off;

% Relative Std of Hours and Wage
subplot(3,2,5);
histogram(rel_std_h_w);
hold on;
mean_rel_std_h_w = mean(rel_std_h_w, 'all');
xline(mean_rel_std_h_w, 'r', 'LineWidth', 2);
title('Relative Std of Hours and Wage');
xlabel('std(h)/std(w)');
ylabel('Frequency');
hold off;

% Correlation of Hours and Wage
subplot(3,2,6);
histogram(corr_h_w);
hold on;
mean_corr_h_w = mean(corr_h_w, 'all');
xline(mean_corr_h_w, 'r', 'LineWidth', 2);
title('Correlation of Hours and Wage');
xlabel('corr(h,w)');
ylabel('Frequency');
hold off;

% Standard Deviation of Consumption
figure;
histogram(std_mat(2,:,:));
hold on;
mean_std_c = mean(std_mat(2,:,:), 'all');
xline(mean_std_c, 'r', 'LineWidth', 2);
title('Standard Deviation of Consumption');
xlabel('std(c)');
ylabel('Frequency');
hold off;