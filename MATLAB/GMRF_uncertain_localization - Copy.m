function GMRF_uncertain_localization()
% =========================================================================
%                   ------------------------------
%                   Approximate Bayesian Inference
%                   ------------------------------
% =========================================================================
M = 0;                                                                      % order of approximation
% profile clear
% profile on

%----- initialize system configuration-------------------------------------
global globaloption allposibleQz allposibleQzInverse
[ globaloption , sys_parameter ] = Configuration(); 
nx = size(globaloption.X_mesh,2) ;                                          % number of X grid
ny = size(globaloption.Y_mesh,2) ;                                          % number of Y grid
ntheta= size(globaloption.hyperparameters_possibilities,1);                 % number of possibilities for $theta$
allposibleQz = zeros((nx*ny)^2,ntheta);
allposibleQzInverse = zeros((nx*ny)^2,ntheta);
for k=1:ntheta
    %---- Construct precision matrix Qz------------------------------------
    allposibleQz(:,k) = reshape(...
        prec_mat_5by5(globaloption.hyperparameters_possibilities(k,1),...
        globaloption.hyperparameters_possibilities(k,1),nx,ny),[],1); 
    allposibleQzInverse(:,k) = reshape(...
        prec_mat_5by5(globaloption.hyperparameters_possibilities(k,1),...
        globaloption.hyperparameters_possibilities(k,1),nx,ny)^(-1),[],1);     
end
option = globaloption; 
gridsize = size(option.grids,1);

%----- generate or load the field------------------------------------------
if isempty(sys_parameter.usesaveddata)
    true_field  =   Generate_new_field();                                   % the data will be saved 
    for t =option.T_mesh 
        true_fieldt = true_field(:,:,t);
        [possible_qt yt] = strobe(true_fieldt);
        AgumentedData(t).possible_q = possible_qt;
        AgumentedData(t).y = yt;
        numberofpossibilities(t) = size(possible_qt.prior_qt,2)^...         % number of possibilities for $q(t)$
             size(possible_qt.prior_qt,1);
    end  
else
    load(sys_parameter.usesaveddata,'AgumentedData','true_fieldt');
end

%----- assign prior --------------------------------------------
prior.beta_mean          = option.beta_prior_mean;
prior.beta_variance      = option.beta_prior_variance;
prior.thetha_probability = option.hyperparameters_possibilities(:,end);     % peior on the $theta$
prior.x_mean = 0;
prior.x_variance = 0;

%----- time loop-----------------------------------------------------------
if (M==0)
    [posteriorTrueQ]  = SPGMRFwithL(...
                prior,...                                                   % prior knowledge
                Data);                                                      % relative measurment at time t-m+1 up to t 
    [posteriorNoisyQ] = SPGMRFwithL(...
                prior,...                                                   % prior knowledge
                Data);                                                      % relative measurment at time t-m+1 up to t                     
else
for t=option.T_mesh
    t
    m = min(M,t);
    Ch = prod(numberofpossibilities(t-m+1:t));                              % number of possibilities for $q(t-m+1:t)$ 
    log_pi_q = -inf * ones(m,Ch);
    do_i_need_to_compute = ones(1,Ch);
    for tt = 1:m                                                            % assign measurment value in Data
        Data(tt).y = AgumentedData(t-m+tt).y;
    end    
    for h=1:Ch                                                              % Some times that the possibilities is zero for some $q$, we don't need to compute estimation for this sampling positions.
        q_choice = quadrotary(h,m*option.agentnumbers);
        for tt = 1:m     
            chosenpoints = (tt-1)*option.agentnumbers+1:...
                tt*option.agentnumbers;
            Data(tt).q = ...
                sum(AgumentedData(t-m+tt).possible_q.support_qt.* ...
                q_choice(chosenpoints,:),2);
            q_prior(:,tt) =...
                sum(AgumentedData(t-m+tt).possible_q.prior_qt.* ...
                q_choice(chosenpoints,:),2);
        end
        if (prod(prod(q_prior))==0)
            do_i_need_to_compute(h) = 0; 
        end
        if do_i_need_to_compute(h)
            [posterior(h)] = SPGMRFwithL(...
                prior(t-m+1),...                                                % prior knowledge
                Data(1:m));                                                 % relative measurment at time t-m+1 up to t 
            log_qtk = sum(log(q_prior(:,1)));
            log_pi_q(1,h) = posterior(h).y_logprobability(1,1) + ...
                log_qtk;                            
            for im = 2:m
                log_qtk = sum(log(q_prior(:,im)));
                log_pi_q(im,h) = log_pi_q(im-1,h) + ...
                    posterior(h).y_logprobability(1,im) + ...
                    log_qtk;    
            end
        end
    end
    
    tmp_c = log(sum(exp(log_pi_q(m,:))));
    
    mux= 0;
    variancex = 0;
    sigmabeta = 0;
    pi_theta = 0*option.hyperparameters_possibilities(:,end);
    for h=1:Ch
        if do_i_need_to_compute(h)
            mux = mux + posterior(h).x_mean(:,m) ...
                * exp(log_pi_q(m,h) - tmp_c);
            pi_theta = pi_theta + posterior(h).thetha_probability...        % update $pi(theta)$
                * exp(log_pi_q(m,h) - tmp_c);
        end
    end
    for h=1:Ch
        if do_i_need_to_compute(h)
            variancex = variancex + (posterior(h).x_variance(:,m) + ...
                (mux-posterior(h).x_mean(:,m)).^2)...
                * exp(log_pi_q(m,h) - tmp_c);
            sigmabeta = sigmabeta + (posterior(h).beta_variance(:,m) + ...
                reshape(diag((mux(gridsize+1:end,1)-...
                posterior(h).x_mean(gridsize+1:end,m)).^2),[],1))...
                * exp(log_pi_q(m,h) - tmp_c);
        end
    end
    
    prior(t+1).beta_mean = mux(gridsize+1:end);
    prior(t+1).beta_variance = sigmabeta; 
    prior(t+1).thetha_probability = pi_theta;
    prior(t+1).x_mean = mux;
    prior(t+1).x_variance = variancex;
   
end
end
%----- save result---------------------------------------------------------
save(sys_parameter.savedata)
% profile viewer


function quad = quadrotary(number,numberofdigits)
%==========================================================================
% It's like Binary but in base 4
%      
number = number-1;
quad = sparse(numberofdigits,4);

for d = numberofdigits:-1:1
    quad(d,floor(mod(number,4^d)/4^(d-1))+1) = 1;
end
% quad = sparse(numberofdigits,4);
% 
% a = dec2bin(number-1);
% binarydigit = size(a,2);
% if mod(binarydigit,2)==1
%     a = ['0' a];
% end
% binarydigit = size(a,2);
% for i = 1 : (binarydigit/2)
%     digit = a(binarydigit - i*2 + 1 : binarydigit - i*2+2);
%     quad(i,:) = [strcmp(digit,'11'),strcmp(digit,'10'),...
%         strcmp(digit,'01'),strcmp(digit,'00')];
%     
% end
% for i = (binarydigit/2)+1 : numberofdigits
%     quad(i,:) = [0 0 0 1];
% end

function [possible_qt yt] = strobe(true_fieldt)
% =========================================================================
%                   ----------------
%                   Sampling process 
%                   ----------------
% In this function "option.agentnumbers" measurement will be measured.
% The relevent positions will be coputed as corner of a rough grids.
% The possiblilty of each corner calculated from the area of the 
%
%
%         q2o---------------oq4
%           |           |   |
%           |    p3     |p1 |
%           |           |   |
%           |-----------*---|
%           |           |   |
%           |    p4     |p2 |
%         q1o---------------oq3
%----- initialize system configuration-------------------------------------
global globaloption
option = globaloption;                                                      % time steps which will be strobed simultanously
%----- sensor network of agents--------------------------------------------
[centerarea,~] = find(...                                                   % put a margin on the possible sampling
    (option.grids(:,1)>(option.X_mesh(1)  +option.margin)) .* ...           %     positions.
    (option.grids(:,1)<(option.X_mesh(end)-option.margin)) .* ... 
    (option.grids(:,2)>(option.Y_mesh(1)  +option.margin)) .* ...
    (option.grids(:,2)<(option.Y_mesh(end)-option.margin)));
selectionsize = size(centerarea,1);
q = centerarea(ceil(selectionsize * rand(option.agentnumbers,1)+1));        % $q$ here is the true position


%----- psossible position using rough grid---------------------------------
step_factor = 5;                                                            % biger step_factor means rough grid
nY = size(option.Y_mesh,2);
ix = floor((q-1)/nY)+1;                                                      
iy = mod(q-1,nY)+1;
ix_minus = floor((ix-1)/step_factor)*step_factor+1;
ix_plus  = (floor((ix-1)/step_factor)+1)*step_factor+1;
iy_minus = floor((iy-1)/step_factor)*step_factor+1;
iy_plus  = (floor((iy-1)/step_factor)+1)*step_factor+1;




possible_qt.support_qt = [(ix_minus - 1) * nY + iy_minus , ...
                          (ix_minus - 1) * nY + iy_plus  , ...
                          (ix_plus  - 1) * nY + iy_minus , ...
                          (ix_plus  - 1) * nY + iy_plus ];
possible_qt.prior_qt   = [(ix_plus - ix) .* (iy_plus - iy) , ...
                          (ix_plus - ix) .* (iy - iy_minus), ...
                          (ix - ix_minus).* (iy_plus - iy) , ...
                          (ix - ix_minus).* (iy - iy_minus)] ...
                          / step_factor^2;
rv = rand(size(q));                      

cs = cumsum(possible_qt.prior_qt,2);
true_q = q;
for i=1:option.agentnumbers
    if ~(i>option.number_uncertain_location)
        if (rv(i)<cs(i,1))
            true_q(i) = possible_qt.support_qt(i,1);
        else if (rv(i)<cs(i,2))
            true_q(i) = possible_qt.support_qt(i,2);
            else if (rv(i)<cs(i,3))
                    true_q(i) = possible_qt.support_qt(i,3);
                else
                    true_q(i) = possible_qt.support_qt(i,4);
                end
            end
        end 
    else
        if (rv(i)<cs(i,1))
            true_q(i) = possible_qt.support_qt(i,1);
            possible_qt.prior_qt(i,:) = [1 0 0 0];
        else if (rv(i)<cs(i,2))
            true_q(i) = possible_qt.support_qt(i,2);
            possible_qt.prior_qt(i,:) = [0 1 0 0];
            else if (rv(i)<cs(i,3))
                    true_q(i) = possible_qt.support_qt(i,3);
                    possible_qt.prior_qt(i,:) = [0 0 1 0];
                else
                    true_q(i) = possible_qt.support_qt(i,4);
                    possible_qt.prior_qt(i,:) = [0 0 0 1];
                end
            end
        end
    end
end
                      
possible_qt.measuredposition   =  q;  
possible_qt.true_q = true_q;

measurement_noise = randn(option.agentnumbers,1) * sqrt(option.s_e2);       % measurement noise
tmp_field_t = reshape(true_fieldt,[],1);                                    % true value of the field
yt = tmp_field_t(true_q)  + measurement_noise;                                   % measurements $y(t)$




