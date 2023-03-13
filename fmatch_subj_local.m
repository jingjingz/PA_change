%% supp material to manuscript 
% "A Riemann Manifold Model Framework for Longitudinal Changes in Physical Activity"
% example Matlab code for estimating the deformation between baseline and follow-up PA curves (parallel computing for subject #id)

function [] = fmatch_subj_local(id)

addpath(genpath('/fshapesTK'))

%% read baseline (cluster)
subj_base_all = csvread('MENU_baseline.csv', 0, 0);
subj_target_all = csvread('MENU_month12.csv', 0, 0);

%% get data for subject i
subj_base = subj_base_all(id, :);
subj_target = subj_target_all(id, :);

%% construct data: source curves
total_time = size(subj_base_all, 2)/2;

f = zeros(total_time,1);
g = ([(1:1:(total_time-1));(2:1:total_time)])';
x = [subj_base(1:total_time)', subj_base((total_time+1):total_time*2)'];

x_target = [subj_target(1:total_time)', subj_target((total_time+1):total_time*2)'];

source = struct('x', x, 'f', f, 'G', g);
target = struct('x', x_target, 'f', f, 'G', g);


%% ESTIMATE DEFORMATION
comp_method = 'matlab'; % possible values are 'cuda' or 'matlab'

% Parameters for the deformations
defo.kernel_size_mom = [.08, .02]; 
% the kernel used to generate the deformations is a sum of 2 kernels
% kernel size is appropriate if both x and y coordinates of the curves are approximately within the range [-1,1]
% otherwise the researcher can apply scaling so that the ranges are in [-1,1]
% it is suggested the researcher sets the first kernel size to be about 1/10 of the actual data range
% and the second about 1/4 - 1/2 of the first kernel
% one can also use only one kernel value (~1/10 of data range)

defo.nb_euler_steps = 10; % nbr of steps in the (for||back)ward integration
% usually no. of step = 10 is sufficient (11 if including the final step)

defo.method =comp_method; % possible values are 'cuda' or 'matlab'

% the following parameters can be kept as they are for the estimation of deformation of PA curves
% Parameters for the matchterm
objfun.distance = 'kernel';
objfun.kernel_distance.distance = 'var';
objfun.kernel_distance.kernel_size_geom=.1;
% size of the geometric kernel  in the data attachment term  (2 runs : 1st at large scale and 2nd at smaller scale)
objfun.kernel_distance.kernel_size_signal=1;
% size of the functional kernel in the data attachment term
%(2 runs with infinite scale size so that signal do not affect the matching)
objfun.kernel_distance.method=comp_method; % possible values are 'cuda' or 'matlab'

objfun.signal_type = 'vertex';
objfun.weight_coef_dist = 3000; % weighting coeff in front of the data attachment term
objfun.weight_coef_pen_fr =1;% weighting coeff in front of funres penalization
objfun.weight_coef_pen_f = 1; % weighting coeff in front of fun penalization

% Parameters for the optimization
optim.method = 'gradDesc';
optim.gradDesc.max_nb_iter = 50; % !! it is better to have n_iter >= 40
optim.gradDesc.step_size_fr = 1e-6;
optim.gradDesc.step_size_p = 1e-6;
optim.gradDesc.kernel_size_signal_reg=0;
optim.gradDesc.min_fun_decrease=1e-10;

objfun.pen_signal = 'l2';
objfun.fem_type	= 'lump';


% ESTIMATE DEFORMATION
[p,pf,summary] = fsmatch_tan(source,target,defo,objfun,optim);

% shoot from template to get estimated target
momentum_e = p{1,1};

[points_evol,obj_evol,mom_evol] = shoot_and_flow_tan(source.x, momentum_e, defo);

% save resulting manifold of shoot_and_flow
data_e.x = points_evol{1, defo.nb_euler_steps + 1};
data_e.G = source.G;
data_e.f = source.f;

%% save estimated target
csvwrite(['MENU_defor_est_' num2str(id)...
    '.csv'], momentum_e)
csvwrite(['MENU_est_' num2str(id) '.csv'],...
    data_e.x)


end