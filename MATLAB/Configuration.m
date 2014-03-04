function [ option , sys_parameter ] = Configuration()
%This function includes general tuning parameters.
%   option is model and algorithm tuning parameters.
%   sys_parameter is plots and system configuration.
% Authour : Mahdi Jadaliha
% e-mail  : jadaliha@gmail.com
% =========================================================================
sys_parameter.plot_enabled = 1;
sys_parameter.color = {'b','g','r','c','m','y'};
sys_parameter.usesaveddata = '';%'filename.mat';
sys_parameter.savedata = 'test.mat';%'filename.mat';
sys_parameter.torus = 1; %The precision matrix Q introduced above is chosen with the regular lattices wrapped on a torus.
sys_parameter.realdata = 0; % 1 -> use real data instead of generating fields

%--------Dimension of the field--------------------------------------------
if isempty(sys_parameter.usesaveddata)
    Xmax = 15; Ymax =15; Tmax =10;
else
    load(sys_parameter.usesaveddata,'true_field');
    Xmax = (size(true_field,2)-1)/2;
    Ymax = (size(true_field,1)-1)/2; 
    Tmax =  size(true_field,3);
end
option.finegridsize = 1;
option.margin = 10;
option.X_mesh = (-Xmax:option.finegridsize:Xmax);                           % X axes cordinate
option.Y_mesh = (-Ymax:option.finegridsize:Ymax);                           % Y axes cordinate
option.T_mesh = (1:1:Tmax); % Time   cordinate
[tmp_S1 tmp_S2] = meshgrid(option.X_mesh,option.Y_mesh);
option.grids = [reshape(tmp_S1,[],1) reshape(tmp_S2,[],1)];


%--------Hyperparameters---------------------------------------------------
option.kappa = 50; option.alpha = 0.1;
% the last column is the probability for the rowset.
% option.hyperparameters_possibilities = ...
%     [    option.kappa,               option.alpha,       0.2;
%     0.1* option.kappa,               option.alpha,       0.2;
%     10 * option.kappa,               option.alpha,       0.2;
%          option.kappa,          0.1* option.alpha,       0.2;
%          option.kappa,          10 * option.alpha,       0.2];
%      0.1*option.kappa,      option.alpha,       0.12 ];% option.hyperparameters_possibilities = ...
%     [option.kappa,          option.alpha,       0.12 ;
%      0.1*option.kappa,      option.alpha,       0.12 ;
%      option.kappa,          0.1*option.alpha,   0.12 ;
%      10*option.kappa,       option.alpha,       0.12 ;
%      option.kappa,          10*option.alpha,    0.12 ;
%      0.1*option.kappa,      0.1*option.alpha,   0.1 ;
%      0.1*option.kappa,      10*option.alpha,    0.1 ;   
%      10*option.kappa,       0.1*option.alpha,   0.1 ;
%      10*option.kappa,       10*option.alpha,    0.1 ];
option.hyperparameters_possibilities = ...
    [option.kappa,          option.alpha,       1 ];
option.c_randomwallk = 0.0;
%--------Vehicle Dynamic---------------------------------------------------
%   q(t+1) = A q(t) + B u(t) + w(t)
%   ~q(t)  = Lt q(t)+ e(t)
option.vehicle.StateMatrix = 1;         %A
option.vehicle.InputMatrix = 1;         %B
option.vehicle.InputVariance = 1;       %Sigma(ut)
option.vehicle.ModelUncertanity = 0.25;    %Sigma(wt)
option.vehicle.OutputMatrix = 1;        %Lt
option.vehicle.ObservationNoise = 9.1;    %Sigma(et)
%--------Sensor network----------------------------------------------------
option.agentnumbers = 1;                %N
option.agentmobilityconstrain = 0;      %In the case that there was some mobility constrain
option.s_e2 = 0.01;                      %Sigma(epsilon^2)

%-------uncertain localization parameter-----------------------------------
option.roughgrid = 1;                                                       % This parameter determine the size of rogh grids with respect to finer grids
nx = size(option.X_mesh,2); ny = size(option.Y_mesh,2);                     % $nx$ and $ny$ are the size of mesh grids which we are interested to do prediction
centergrid = ceil(size(option.grids,1)/2);                                  % compute the index of center poins of the field
tmp1 = (floor(nx/option.roughgrid)-1)/2;
tmp2 = (floor(ny/option.roughgrid)-1)/2;
tmp3 = ones(tmp2*2+1,1)*(-tmp1:tmp1) * option.roughgrid * ny + ...
    (-tmp2:tmp2)'*ones(1,tmp1*2+1) * option.roughgrid + centergrid;
option.roughgrids = reshape(tmp3,[],1);                                     % The rough grids are all possible positions which we can get measurments.

option.number_uncertain_location = min(0, option.agentnumbers);

%%%%%%%%%%%%%%%%%%%
option.Sigmau = eye(option.agentnumbers);
%-------warnings-----------------------------------------------------------
if(size(option.roughgrids,1)<option.agentnumbers)
    disp('Number of agents is more than possible positions.')
end
end

