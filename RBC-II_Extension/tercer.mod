%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  LAMBDA GROUP %%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%% TOPICOS DSGE - RBC %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Codigo que describe un modelo RBC basico loglinealizado (expansion de
// Taylor de 1er orden). Las CPO se presentan de forma no lineal.
// El Estado Estacionario es obtenido manualmente.
// La funcion de utilidad corresponde a la forma KPR. 
// (c) Carlos Rojas Quiroz 

var lab c w r y kap innv z g;
predetermined_variables kap;
varexo e_z e_g;
parameters alpha delta betta sigma theta rho_z rho_g 
z_ss lab_ss r_ss  kap_ss w_ss y_ss c_ss inv_ss g_ss C_Y I_Y G_Y;

alpha  = 1-0.33;
delta  = 0.023;
betta  = 0.99;
sigma  = 1.0;
rho_z  = 0.95;
rho_g  = 0.75;
z_ss   = 1; 
G_Y    = 0.155;
lab_ss = 0.383409634485720;
y_ss   = z_ss^(1/alpha)*((1-alpha)*betta/(1-betta+betta*delta))^((1-alpha)/alpha)*lab_ss;
c_ss   = ((1-betta+alpha*betta*delta)/(1-betta+betta*delta)-G_Y)*y_ss;
theta  = 1/((1/lab_ss-1)*alpha*y_ss/c_ss+1);
w_ss   = alpha*y_ss/lab_ss;
kap_ss = (1-alpha)*betta/(1-betta+betta*delta)*y_ss;
inv_ss = delta*kap_ss;
r_ss   = (1-alpha)*y_ss/kap_ss-delta;
g_ss   = G_Y*y_ss;
C_Y    = c_ss/y_ss;
I_Y    = inv_ss/y_ss;

model;
exp(lab)        = 1-(1-theta)/theta*(exp(c)/exp(w));
exp(c)^(theta-1)*(1-exp(lab))^(1-theta)*(exp(c)^theta*(1-exp(lab)^(1-theta))^(-sigma)) = 
betta*(exp(c(+1))^(theta-1)*(1-exp(lab(+1)))^(1-theta)*(exp(c(+1))^theta*(1-exp(lab(+1))^(1-theta))^(-sigma)))*(1+exp(r(+1)));
exp(w)          = alpha*exp(y)/exp(lab);
exp(r)+delta    = (1-alpha)*exp(y)/exp(kap);
exp(y)          = exp(c)+exp(innv)+exp(g);
exp(kap(+1))    = (1-delta)*exp(kap)+exp(innv);
exp(y)          = exp(z)*exp(kap)^(1-alpha)*exp(lab)^alpha;
z               = (1-rho_z)*log(z_ss) + rho_z*z(-1) + e_z;
g               = (1-rho_g)*log(g_ss) + rho_g*g(-1) + e_g;
end;

steady_state_model;
lab =log(lab_ss);
c   =log(c_ss); 
w   =log(w_ss); 
r   =log(r_ss); 
y   =log(y_ss); 
kap =log(kap_ss); 
innv=log(inv_ss); 
z   =log(z_ss);
g   =log(g_ss);
end;

shocks;

var e_z; stderr 0.01;
var e_g; stderr 0.01;
end;
 
resid;
steady;
check;

stoch_simul(order=1, irf=20, loglinear, hp_filter=1600, simul_replic=100, periods=179) y c innv lab z w;

% Read out simulations
simulated_series_raw = get_simul_replications(M_, options_);

% Filter series
simulated_series_filtered = NaN(size(simulated_series_raw));
for ii = 1:options_.simul_replic
    [trend, cycle] = sample_hp_filter(simulated_series_raw(:,:,ii)', 1600);
    simulated_series_filtered(:,:,ii) = cycle';
end

% Get variable positions
y_pos = strmatch('y', M_.endo_names, 'exact');
c_pos = strmatch('c', M_.endo_names, 'exact');
i_pos = strmatch('innv', M_.endo_names, 'exact');
h_pos = strmatch('lab', M_.endo_names, 'exact');
productivity_pos = strmatch('z', M_.endo_names, 'exact');
w_pos = strmatch('w', M_.endo_names, 'exact');

var_positions = [y_pos; c_pos; i_pos; h_pos; productivity_pos];
var_names = M_.endo_names_long(var_positions,:);

% Compute standard deviations
std_mat = std(simulated_series_filtered(var_positions,:,:), 0, 2) * 100;
std_h = std(simulated_series_filtered(h_pos,:,:), 0, 2) * 100;
std_w = std(simulated_series_filtered(w_pos,:,:), 0, 2) * 100;

% Compute relative standard deviations
rel_std_mat = std_mat ./ std_mat(1,:,:); % std(x) / std(y)
rel_std_h_w = std_h ./ std_w; % std(h) / std(w)

% Compute correlations
corr_mat = NaN(5, options_.simul_replic);
for ii = 1:options_.simul_replic
    corr_mat(1,ii) = corr(simulated_series_filtered(y_pos,:,ii)', simulated_series_filtered(y_pos,:,ii)');
    corr_mat(2,ii) = corr(simulated_series_filtered(y_pos,:,ii)', simulated_series_filtered(c_pos,:,ii)');
    corr_mat(3,ii) = corr(simulated_series_filtered(y_pos,:,ii)', simulated_series_filtered(i_pos,:,ii)');
    corr_mat(4,ii) = corr(simulated_series_filtered(y_pos,:,ii)', simulated_series_filtered(h_pos,:,ii)');
    corr_mat(5,ii) = corr(simulated_series_filtered(y_pos,:,ii)', simulated_series_filtered(productivity_pos,:,ii)');
end

% Compute correlation between h and w
corr_h_w = NaN(1, options_.simul_replic);
for ii = 1:options_.simul_replic
    corr_h_w(1,ii) = corr(simulated_series_filtered(h_pos,:,ii)', simulated_series_filtered(w_pos,:,ii)');
end

% Print table with standard deviations and correlations
fprintf('\n%-40s \n', 'Standard Deviations and Correlations')
fprintf('%-20s \t %11s \t %11s \n', '', 'std(x)', 'corr(y,x)')
for ii = 1:size(corr_mat,1)
    fprintf('%-20s \t %3.2f (%3.2f) \t %3.2f (%3.2f) \n', var_names{ii,:}, mean(std_mat(ii,:,:), 3), std(std_mat(ii,:,:), 0, 3), mean(corr_mat(ii,:), 2), std(corr_mat(ii,:), 0, 2))
end

% Print table with relative standard deviations
fprintf('\n%-40s \n', 'Relative Standard Deviations')
fprintf('%-20s \t %11s \n', '', 'std(x)/std(y)')
for ii = 2:size(corr_mat,1)
    fprintf('%-20s \t %3.2f (%3.2f) \n', var_names{ii,:}, mean(rel_std_mat(ii,:,:), 3), std(rel_std_mat(ii,:,:), 0, 3))
end

% Print table with correlation of h and w and relative std of h and w
fprintf('\n%-40s \n', 'Correlation of h and w, Relative Std of h and w')
fprintf('%-20s \t %11s \t %11s \n', 'Variable', 'corr(h,w)', 'std(h)/std(w)')
fprintf('%-20s \t %3.2f (%3.2f) \t %3.2f (%3.2f) \n', 'h, w', mean(corr_h_w, 2), std(corr_h_w, 0, 2), mean(rel_std_h_w, 3), std(rel_std_h_w, 0, 3))
