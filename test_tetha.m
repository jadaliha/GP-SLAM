function test_tetha()
% =========================================================================
%       ---------------------------------------------------------
%        Test the effect of Tetha and its convergance properties
%       ---------------------------------------------------------
% =========================================================================
                                                                            % order of approximation
%----- initialize system configuration-------------------------------------
global globaloption allposibleQz allposibleQzInverse sys_parameter
[ globaloption , sys_parameter ] = Configuration();
%----- generate or load the field------------------------------------------
option = globaloption;
[~, VehicleStates, AgumentedData, numberofpossibilities]  =  ...            % [field, VehicleStates, AgumentedData, numberofpossibilities] = Generate_new_field(savetofile)
    Generate_new_field();                                                   % the data will be saved in 'savetofile'

%--------------------------------------------------------------------------
nx = size(option.X_mesh,2) ;                                                % number of X grid
ny = size(option.Y_mesh,2) ;                                                % number of Y grid
nt = size(option.T_mesh,2) ;                                                % number of time steps
ntheta= size(option.hyperparameters_possibilities,1);                       % number of possibilities for $theta$
gridsize = size(option.grids,1);                                            % number of spatial sites

%----- Construct precision matrix Qz---------------------------------------
allposibleQz = zeros((nx*ny)^2,ntheta);
allposibleQzInverse = zeros((nx*ny)^2,ntheta);
for k=1:ntheta 
    tmp1 = prec_mat_5by5(option.hyperparameters_possibilities(k,1)...
        ,option.hyperparameters_possibilities(k,2),nx,ny);
    allposibleQz(:,k) = reshape(tmp1,[],1); 
    allposibleQzInverse(:,k) = reshape(tmp1^(-1),[],1);     
end 
N = option.agentnumbers*2 + nx*ny;
globalmux = zeros(N,ntheta); 
globalQx = zeros(N^2,ntheta);
...
    
%----- assign prior -------------------------------------------------------
prior(1).thetha_probability = option.hyperparameters_possibilities(:,end);  % peior on the $theta$
prior(1).x_mean = 0;
prior(1).x_variance = 10^5;                                                 % without prior knowledge we should use a big variance
    

%----- time loop-----------------------------------------------------------
for tt = option.T_mesh                                                      % assign measurment value in Data
    Data(tt).y = AgumentedData(tt).y;
    Data(tt).q = AgumentedData(tt).possible_q.true_q;  
    Data(tt).u = VehicleStates(tt).u;
    Data(tt).h = VehicleStates(tt).h;
end  
[posteriorTrueQ]  = SPGMRFwithL(...
            prior,...                                                       % prior knowledge
            Data);                                                          % relative measurment at time t-m+1 up to t 

for tt = option.T_mesh                                                      % assign measurment value in Data
    Data(tt).y = AgumentedData(tt).y;
    Data(tt).q = AgumentedData(tt).possible_q.measuredposition;        
end             
[posteriorNoisyQ] = SPGMRFwithL(...
            [],...                                                          % prior knowledge
            Data);                                                          % relative measurment at time t-m+1 up to t                     
   
for t=option.T_mesh 
    t
    m = min(M,t);
    Ch = prod(prod(numberofpossibilities(:,t-m+1:t)));                      % number of possibilities for $q(t-m+1:t)$ 
    log_pi_q = [];

    for h=1:Ch                                                              % Some times that the possibilities is zero for some $q$, we don't need to compute estimation for this sampling positions.
        [Data q_prior] = quadrotary(h,AgumentedData(t-m+1:t));

        [posterior(h)] = SPGMRFwithL(...
            prior(t-m+1),...                                                % prior knowledge
            Data(1:m));                                                     % relative measurment at time t-m+1 up to t 
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
    
    tmp_c = log(sum(exp(log_pi_q(m,:))));
    
    mux= 0;
    variancex = 0;
    sigmabeta = 0;
    pi_theta = 0*option.hyperparameters_possibilities(:,end);
    for h=1:Ch
            mux = mux + posterior(h).x_mean(:,m) ...
                * exp(log_pi_q(m,h) - tmp_c);
            pi_theta = pi_theta + posterior(h).thetha_probability(:,end)... % update $pi(theta)$
                * exp(log_pi_q(m,h) - tmp_c);
    end
    for h=1:Ch
            variancex = variancex + (posterior(h).x_variance(:,m) + ...
                (mux-posterior(h).x_mean(:,m)).^2)...
                * exp(log_pi_q(m,h) - tmp_c);
            sigmabeta = sigmabeta + (posterior(h).beta_variance(:,m) + ...
                reshape(diag((mux(gridsize+1:end,1)-...
                posterior(h).x_mean(gridsize+1:end,m)).^2),[],1))...
                * exp(log_pi_q(m,h) - tmp_c);
    end
    
    prior(t+1).beta_mean = mux(gridsize+1:end);
    prior(t+1).beta_variance = sigmabeta; 
    prior(t+1).thetha_probability = pi_theta;
    prior(t+1).x_mean = mux;
    prior(t+1).x_variance = variancex;
   
end

%----- save result---------------------------------------------------------
save(sys_parameter.savedata, ...
    'M', 'posteriorTrueQ', 'posteriorNoisyQ',...
    'prior','AgumentedData', 'true_field', 'beta', ...
    'option', 'sys_parameter')
% profile viewer
% profile off


