% =========================================================================
%       ---------------------------------------------------------
%        Test the effect of Tetha and its convergance properties
%       ---------------------------------------------------------
% Authour : Mahdi Jadaliha
% e-mail  : jadaliha@gmail.com
% =========================================================================
tic  
close all
clear all
clc
global globaloption  sys_parameter
%----- initialize system configuration-------------------------------------
[ globaloption , sys_parameter ] = Configuration();
option = globaloption;

%%% Case 1: Generate data
%----- generate or load the field------------------------------------------
% [true_field, VehicleStates, AgumentedData, ~]  =  ...                       % [field, VehicleStates, AgumentedData, numberofpossibilities] = Generate_new_field(savetofile)
%     Generate_new_field();                                                   % the data will be saved in 'savetofile'


%%% Case 2: Load data
importxls                                                                   % Load the excel sheet
%load data %for loading data quickly
noiseMagnitude=30;
location = [locationX, locationY];    % Assign the Location
dataSize=size(locationX,1);
noiseAdded=noiseMagnitude*rand(dataSize,1);
noiseAdded=noiseAdded-mean(noiseAdded);
location_noise=[locationX+noiseAdded,...
                locationY+noiseAdded];
features = AFFT1;                                                           % select the one or more features
%---- we normalize features here
features = features - ones(size(features,1),1)*mean(features,1);            % remove mean average
cov_features = cov(features);
[U,S,~] =  svd(cov_features);                                               % make feature orthogonal using SVD 
features = features*U*sqrt(S^-1);                                           %transfer top new coordinates
%---- Predict hyper parameters
nf = size(features,2);
hyper_p = zeros(4,nf);
for index = 1:nf
    f = features(:,index);
    p0 = [var(f) 1 1 0.5]'; % sig_f^2, sig_x sig_y sigma2w                  % Initial guess of hyper-parameters
    [hyper_p(:,index), ~] = HyperParameter(f,location(:,1:2), p0);          % extract hyper parameters for each layer sepratly
    %hyper_p(:,index) = [1;0.5;0.5];                                           % You may overwrite the hyper parameter values manually
    str = sprintf('hyper-parameters # %d: \t (sig_f^2: %0.2f, \t sig_x: %0.2f, \t sig_y: %0.2f, \t sig_w^2: %0.2f)',index,hyper_p(1,index),hyper_p(2,index),hyper_p(3,index),hyper_p(4,index));
    disp(str)
end

% Normalize coordinate and compute Precision hyperparameters from GP
% hyperparameters
index = 1; % this part should be changed for multi dimensional cases
option.x_scale = hyper_p(2,index);                                          % scale x coordinate with the x direction bandwidth
option.y_scale = hyper_p(3,index);
location       = [locationX/option.x_scale, locationY/option.y_scale];      % rescale cartesian coordinates
location_noise = [location_noise(:,1)/option.x_scale,...
                  location_noise(:,2)/option.y_scale];

h = 1; % discritization factor
option.alpha = h^(-2) * 2; % $\ell = \frac{1}{h} \sqrt{\frac{\alpha}{2}}$
option.kappa = (hyper_p(1,index)* 4 * pi * option.alpha)^-1;  % $\sigma_f^2 = \frac{1}{4 \pi \alpha \kappa}$
% because of normalization bandwith = 1 here
option.hyperparameters_possibilities = ...
    [option.kappa,          option.alpha,       1 ];
option.X_mesh = (150:5:450)/option.x_scale;
option.Y_mesh = (200:5:400)/option.y_scale;
option.T_mesh = (1:size(f,1));
[tmp_S1, tmp_S2] = meshgrid(option.X_mesh,option.Y_mesh);
option.grids = [reshape(tmp_S1,[],1) reshape(tmp_S2,[],1)];
nx = size(option.X_mesh,2) ;                                                % number of X grid
ny = size(option.Y_mesh,2) ;                                                % number of Y grid
nt = size(option.T_mesh,2) ;                                                % number of time steps

%regression loop:
for t=1:nt
%     AgumentedData(t).y = f(t);
%     dist_grid_from_continous = ...
%         sqrt((option.grids(:,1) - location(t,1)).^2 + ...
%         (option.grids(:,2) - location(t,2)).^2);
%     [IC,IX] = sort(dist_grid_from_continous);
%     AgumentedData(t).possible_q.true_q = IX(1);
%     AgumentedData(t).possible_q.support_qt{1} = IX(1);
%     AgumentedData(t).possible_q.prior_qt = 1;
%     AgumentedData(t).possible_q.N_possible_qt = 1;
%     AgumentedData(t).possible_q.measuredposition = IX(1);
    
    AgumentedData(t).y = f(t);
    
    dist_grid_from_continous_true= ...
        sqrt((option.grids(:,1) - location(t,1)).^2 + ...
        (option.grids(:,2) - location(t,2)).^2);
    [IC_true,IX_true] = sort(dist_grid_from_continous_true);
    
    dist_grid_from_continous_noise = ...
        sqrt((option.grids(:,1) - location_noise(t,1)).^2 + ...
        (option.grids(:,2) - location_noise(t,2)).^2);
    [IC_noise,IX_noise] = sort(dist_grid_from_continous_noise);
    omegaSet=omegaMaker(option.grids,IX_true,IX_noise,nx,ny);
    sizeOmg=size(omegaSet,1);
    AgumentedData(t).possible_q.true_q = IX_true(1);
    for i=1:sizeOmg
       AgumentedData(t).possible_q.support_qt{i}=omegaSet(i);
       AgumentedData(t).possible_q.prior_qt{i} = 1/sizeOmg;
    end
    %AgumentedData(t).possible_q.support_qt{1} = IX_true(1);
    %AgumentedData(t).possible_q.prior_qt = 1;
    AgumentedData(t).possible_q.N_possible_qt = sizeOmg; %??
    AgumentedData(t).possible_q.measuredposition = IX_noise(1);
end
ntheta= size(option.hyperparameters_possibilities,1);                       % number of possibilities for $theta$
n = size(option.grids,1);                                                   % number of spatial sites
N = option.agentnumbers;
x_width = option.X_mesh(end) - option.X_mesh(1) + option.finegridsize;
y_width = option.Y_mesh(end) - option.Y_mesh(1) + option.finegridsize;
phi = 0 ; option.vehicle.ModelUncertanity = 4;
%----- Construct precision matrix Qz---------------------------------------
muz_theta = zeros(n,ntheta);
Sigmaz_theta = zeros(n^2,ntheta);
for indtheta=1:ntheta 
    tmp1 = prec_mat_5by5(option.hyperparameters_possibilities(indtheta,1)...
        ,option.hyperparameters_possibilities(indtheta,2),nx,ny);
    Sigmaz_theta(:,indtheta) = reshape(tmp1^(-1),[],1);
end 



f_theta = log(option.hyperparameters_possibilities(:,end));
itmp = (1:2*N);
jtmp = kron((1:N)',[1;1]);
mux = zeros(2*N,1);
Sigmax = eye(2*N)*1000;
mux_ = mux;
Sigmax_ = Sigmax;

for t=option.T_mesh 
    t %#ok<NOPTS>
    xtilda = reshape(option.grids(...
        AgumentedData(1,t).possible_q.measuredposition,:)',[],1);
    ztilda = AgumentedData(t).y;

    
    qt= AgumentedData(1,t).possible_q.support_qt{1,1}';
    for indj = 2:N
        qtN = AgumentedData(1,t).possible_q.support_qt{1,indj}';
        qt= [kron(ones(1,size(qtN,2)),qt);...
            kron(qtN,ones(1,size(qt,2)))];
    end
    Sigma_e = eye(2*size(qt,1))* option.vehicle.ObservationNoise;


    
    numberofpossibleqt = size(qt,2);
    Sigma_epsilon = eye(size(qt,1))* option.s_e2;
    f_qtandtheta = zeros(numberofpossibleqt,ntheta);
    muz_qandtheta =    zeros(n, 1);
    Sigmaz_qandtheta = zeros(n^2,1);
    for indq = 1:numberofpossibleqt
        q = qt(:,indq);
        H = sparse(1:N,q',ones(1,N),N,n);                                   % find the map $H_{qt}$ from $q{t}$ to spacial sites S.
        x = reshape(option.grids(q,:)',[],1);
        f_qt  = logmvnpdf2(x,mux_,Sigmax_,option);
        f_qtildaGqt  = logmvnpdf2(xtilda,x,Sigma_e,option);
        for indtheta = 1:ntheta
            muztheta = muz_theta(:,indtheta);
            Sigmaztheta      = reshape(Sigmaz_theta(:,indtheta),n,n);
            muztildatheta    = H * muztheta;                                % $mu_{tilde{zt}|theta,D{t-1},\q{t}}  =  H_{qt} mu_{z| theta,D{t-1}}$            
            Sigmaztildatheta = Sigma_epsilon + H * Sigmaztheta * H';
            f_ztilda = logmvnpdf(ztilda,muztildatheta,Sigmaztildatheta);
            % approximate the distribution of $\q{t},\theta |\D{0}{t}$           
            f_qtandtheta(indq,indtheta) = ...
                f_qtildaGqt  + f_theta(indtheta) + f_qt + f_ztilda;         % pi(q{t},theta|D{t})
            

        end
    end
    tmp_c = log(sum(sum(exp(f_qtandtheta))));                               % Here we find normalization factor $tmp_c$
%     if (abs(tmp_c)>100)
%         pause
%     end
    f_qtandtheta = f_qtandtheta - tmp_c;                                    % After this line $f_qtandtheta$ is normalized, 
    pi_qtandtheta = exp(f_qtandtheta);                                      % $sum(pi_qtandtheta) = 1$.
    
    pi_theta = sum(pi_qtandtheta,1);                
    pi_theta = pi_theta/sum(pi_theta);                                      % Compensate computational error
    f_theta = log(pi_theta);
    pi_q = sum(pi_qtandtheta,2);
    
    
    
    
    % For the toures correction on the sampling positions measured  close
    % to the borders
    mux = zeros(2*N,1); %#ok<NASGU>
    xxx = zeros(2*N,numberofpossibleqt);
    for indj = 1:N
        xx = option.grids(qt(indj,:),:);
%         % uncomment for tourse
%         median_xx = median(xx);
%         tmp1 = find((xx(:,1) - median_xx(1))> x_width/2);
%         xx(tmp1,1) = xx(tmp1,1)- x_width;
%         tmp2 = find((xx(:,1) - median_xx(1))<-x_width/2);
%         xx(tmp2,1) = xx(tmp2,1)+ x_width;
%         tmp1 = find((xx(:,2) - median_xx(2))> y_width/2);
%         xx(tmp1,2) = xx(tmp1,2)- y_width;
%         tmp2 = find((xx(:,2) - median_xx(2))<-y_width/2);
%         xx(tmp2,2) = xx(tmp2,2)+ y_width;
        
        xxx((indj*2)-1:indj*2,:) = xx';
    end
    mux = (xxx * pi_q);
    Sigmax = reshape(...
        VectorSquare(xxx - mux * ones(1,numberofpossibleqt)) * pi_q,...
        2*N,2*N); 
    

    for indtheta = 1:ntheta
        pi_qtGtheta = (pi_qtandtheta(:,indtheta)/pi_theta(indtheta));
        muz_theta_temp = zeros(n,1);
        Sigmaz_theta_temp = zeros(n^2,1);
        for indq = 1:numberofpossibleqt
            q = qt(:,indq);
            H = sparse(1:N,q',ones(1,N),N,n);                               % find the map $H_{qt}$ from $q{t}$ to spacial sites S.
            muztheta = muz_theta(:,indtheta);
            Sigmaztheta      = reshape(Sigmaz_theta(:,indtheta),n,n);
            Sigmaztheta = 0.5 * (Sigmaztheta+Sigmaztheta');
            muztildatheta    = H * muztheta;
            tmp1 = Sigmaztheta * H';
            Sigmaztildatheta = Sigma_epsilon + H * tmp1;
            InvSigmaztildatheta = Sigmaztildatheta^(-1);
 
            muz_qandtheta    = ...
                muztheta + tmp1 * ...
                InvSigmaztildatheta * (ztilda - muztildatheta);
            Sigmaz_qandtheta = ...
                reshape(Sigmaztheta - tmp1 * ...
                InvSigmaztildatheta * tmp1',[],1);
            
            muz_theta_temp = muz_theta_temp + ...
                muz_qandtheta * pi_qtGtheta(indq);

            Sigmaz_theta_temp = Sigmaz_theta_temp + ...
                (Sigmaz_qandtheta + VectorSquare(muz_qandtheta))*...
                pi_qtGtheta(indq) ;

        end
        muz_theta(:,indtheta)    = muz_theta_temp;
        Sigmaz_theta(:,indtheta) = ...
            Sigmaz_theta_temp - VectorSquare(muz_theta_temp);
    end
    
    posterior(t).muz = muz_theta * pi_theta';
    posterior(t).varz= ...
        (Sigmaz_theta(1:n+1:end,:)+muz_theta.*muz_theta)* pi_theta' ...
        - posterior(t).muz.* posterior(t).muz;
    posterior(t).mux = mux;
    posterior(t).Sigmax = Sigmax;

    % predict next sampling position
    u = velocity(t);
    phi = phi + turnrate(t);
    stmp = reshape([cos(phi),sin(phi)]',[],1);
    F = sparse(itmp,jtmp,stmp,2*N,N);
    mux_    = mux + F * u;
    Sigma_w = eye(2*size(qt,1))* option.vehicle.ModelUncertanity;
    Sigmax_ = Sigmax + Sigma_w;  
    
    
end
% 
toc

for i=1:463
    DT_x(i) = posterior(i).mux(1);
    DT_y(i) = posterior(i).mux(2);
end
plot(DT_x,DT_y)

