function GMRF_uncertain_localization()
% =========================================================================
%                   ------------------------------
%                   Approximate Bayesian Inference
%                   ------------------------------
% =========================================================================
                                                                            % order of approximation
% profile clear
% profile on
M = 1;

%----- initialize system configuration-------------------------------------
global globaloption allposibleQz allposibleQzInverse sys_parameter
[ globaloption , sys_parameter ] = Configuration();
%----- generate or load the field------------------------------------------
if ~isempty(sys_parameter.usesaveddata)
    load(sys_parameter.usesaveddata,...
        'true_field',...
        'AgumentedData',...
        'beta',...
        'option');
    
    if exist('option','var')
        globaloption = option;
    else
        option = globaloption;
    end
    if ~(exist('beta','var'))                                               % check existance of beta in predefined file
        beta = [];  
    end 
    if ~(exist('AgumentedData','var'))
        for t =option.T_mesh 
            true_fieldt = true_field(:,:,t);
            [possible_qt yt] = strobe(true_fieldt);
            AgumentedData(t).possible_q = possible_qt;
            AgumentedData(t).y = yt;
            numberofpossibilities(:,t) = sum(~(possible_qt.prior_qt==0),2); % number of possibilities for $q(t)$
                 
        end  
    end    
else
    option = globaloption;
    [true_field, beta]  =   Generate_new_field();                           % the data will be saved 
    for t =option.T_mesh 
        true_fieldt = true_field(:,:,t);
        [possible_qt yt] = strobe(true_fieldt);
        AgumentedData(t).possible_q = possible_qt;
        AgumentedData(t).y = yt;
        numberofpossibilities(:,t) = sum(~(possible_qt.prior_qt==0),2);     % number of possibilities for $q(t)$
    end  
end


%--------------------------------------------------------------------------
nx = size(option.X_mesh,2) ;                                                % number of X grid
ny = size(option.Y_mesh,2) ;                                                % number of Y grid
nt = size(option.T_mesh,2) ;                                                % number of time steps
ntheta= size(option.hyperparameters_possibilities,1);                       % number of possibilities for $theta$
allposibleQz = zeros((nx*ny)^2,ntheta);
allposibleQzInverse = zeros((nx*ny)^2,ntheta);
for k=1:ntheta 
    %---- Construct precision matrix Qz------------------------------------
    tmp1 = prec_mat_5by5(option.hyperparameters_possibilities(k,1)...
        ,option.hyperparameters_possibilities(k,2),nx,ny);
    allposibleQz(:,k) = reshape(tmp1,[],1); 
    allposibleQzInverse(:,k) = reshape(tmp1^(-1),[],1);     
end 

gridsize = size(option.grids,1);
if ~exist('numberofpossibilities','var')
    numberofpossibilities = zeros (option.agentnumbers, nt);
    for tt=1:nt
        numberofpossibilities(:,tt) = ...
            sum(~(AgumentedData(tt).possible_q.prior_qt==0),2);             % number of possibilities for $q(t)$
    end
end



%----- assign prior --------------------------------------------
prior(1).beta_mean          = option.beta_prior_mean;
prior(1).beta_variance      = option.beta_prior_variance;
prior(1).thetha_probability = option.hyperparameters_possibilities(:,end);  % peior on the $theta$
prior(1).x_mean = 0;
prior(1).x_variance = 0;
    
%----- time loop-----------------------------------------------------------
% if (M==0)   
    for tt = option.T_mesh                                                  % assign measurment value in Data
        Data(tt).y = AgumentedData(tt).y;
        Data(tt).q = AgumentedData(tt).possible_q.true_q;        
    end  
    [posteriorTrueQ]  = SPGMRFwithL(...
                [],...                                                      % prior knowledge
                Data);                                                      % relative measurment at time t-m+1 up to t 
                
    for tt = option.T_mesh                                                  % assign measurment value in Data
        Data(tt).y = AgumentedData(tt).y;
        Data(tt).q = AgumentedData(tt).possible_q.measuredposition;        
    end             
    [posteriorNoisyQ] = SPGMRFwithL(...
                [],...                                                      % prior knowledge
                Data);                                                      % relative measurment at time t-m+1 up to t                     
% else     
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
% end
%----- save result---------------------------------------------------------
save(sys_parameter.savedata, ...
    'M', 'posteriorTrueQ', 'posteriorNoisyQ',...
    'prior','AgumentedData', 'true_field', 'beta', ...
    'option', 'sys_parameter')
% profile viewer
% profile off

function [Data q_prior] = quadrotary(number,parameters)
%==========================================================================
% It's like Binary but in base 4
    
m = size(parameters,2);
number = number - 1;
for tt = 1:m
    
    Data(tt).y = parameters(tt).y;  % assign measurment value in Data
    
    possible_qt = ~(parameters(tt).possible_q.prior_qt==0);
    npqti = sum(possible_qt,2);
    numberofposs = prod(npqti);
    nnn = mod(number,numberofposs);
    number = (number - nnn)/numberofposs;
    
    for ind = 1:size(npqti,1)
        nii = mod(nnn,npqti(ind));
        nnn = (nnn - nii)/npqti(ind);
        tmp1 = find(cumsum(possible_qt(ind,:))==(nii+1));
        selected = tmp1(1);
        Data(tt).q(ind,1) = ...
            parameters(tt).possible_q.support_qt(ind,selected);
        q_prior(ind,tt) = ...
            parameters(tt).possible_q.prior_qt(ind,selected);
    end
    
end


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
step_factor = option.roughgrid;                                             % biger step_factor means rough grid
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
yt = tmp_field_t(true_q)  + measurement_noise;                              % measurements $y(t)$


function [possible_qt yt] = strobe2(true_fieldt)
global globaloption 
option = globaloption; 
selectionsize = size(option.roughgrids,1);
true_q = option.roughgrids(ceil(selectionsize* ...
    rand(option.agentnumbers,1)+1));                                        % $q$ here is the true position
possible_qt.measuredposition = true_q;

n = option.number_uncertain_location;                                       % number of uncertain agents
nx = size(option.X_mesh,2) ; % number of X grid
ny = size(option.Y_mesh,2) ; % number of Y grid


which_quarter_point = rand(n,2);                                            % The first column determine the quarter and second colum determine where in the quarter


x1 = ones(option.roughgrid,1)* (1:option.roughgrid); 
y1 = x1';
% prob_table = zeros(option.roughgrid);
%
%
%         o |
%       ----|----
%           |
%           
%

prob_table = (x1-1).*(y1-1);
prob_table2 = reshape(prob_table,[],1);
prob_table3 = cumsum(prob_table2);
continus_deci = which_quarter_point(:,2) * prob_table3(end);
tmp1 = (   ones(option.roughgrid^2,1)*continus_deci') > ...
          (prob_table3 * ones(1,n));
tmp2 = -diff([ones(1,option.number_uncertain_location);tmp1]); %?????
[row,~] = find(tmp2);
tmp3 = [mod(row,option.roughgrid)+1,ceil(row/option.roughgrid)];
neighborhood = option.finegridsize * option.roughgrid;
for ind=1:n
    true_cordinate = option.grids(true_q(ind),:);
    
    boundary_condition = ...
        [~(true_cordinate(:,1)+ neighborhood > option.X_mesh(end)),...      % right side is possible?
         ~(true_cordinate(:,1)- neighborhood <-option.X_mesh(end)),...      % left  side is possible?
         ~(true_cordinate(:,2)+ neighborhood > option.Y_mesh(end)),...      % up    side is possible?
         ~(true_cordinate(:,2)- neighborhood <-option.Y_mesh(end))];        % down  side is possible?
    possible_quarter = [boundary_condition(1)*boundary_condition(3),...     % 1st quarter is possible?
                        boundary_condition(2)*boundary_condition(3),...     % 2nd quarter is possible?
                        boundary_condition(2)*boundary_condition(4),...     % 3th quarter is possible?
                        boundary_condition(1)*boundary_condition(4)] ;      % 4th quarter is possible?
                    
    tmp4 = cumsum(  possible_quarter );               
    
    tmp5 = which_quarter_point(ind,1)*tmp4(end);
    tmp3(ind,3) = find(diff(tmp5<[-1,tmp4]));
    
    tmp3(:,1:2) = 10 - tmp3(:,1:2);
    switch tmp3(ind,3)
        case 1 %1st quarter
            possible_qt.measuredposition(ind) = ...
                tmp3(ind,1) * ny + tmp3(ind,2);
%             possible_qt.support_qt(i,2)
%             possible_qt.prior_qt(i,:)
        case 2 %2nd quarter 
        case 3
        case 4
        otherwise
            disp('Something is wrong!');
    end



end

