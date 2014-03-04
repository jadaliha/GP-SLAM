function [posterior]=SPGMRFwithL(prior,data)
% =========================================================================
%                   ------------------------------
%                   Approximate Bayesian Inference
%                   ------------------------------
%
% Latent Gaussian model:
%    y_i = z_i + e_i, e ~ N(0,sigma_w^2), sigma_w = 1/sqrt(gamma*kappa),
%    z_i = F(s_i)*beta + eta(s_i), beta ~ N(0,T/kappa),
%    eta ~ GP(0,inv(Q))
%    
% Structure of Q:
%                 1 
%             2  -2a   2
%         1 -2a 4+a^2 -2a 1
%             2  -2a   2
%                 1 
%
% Latent variables:
%    z = (z_1,...,z_n)^T
%    beta = (beta_1,...,beta_p)^T
%
% Hyperparamters:
%    log(kappa)
%    log(alpha)
%
% Known parameters:
%    F, T, gamma
%
% =========================================================================


global globaloption
option = globaloption;                                                      % initialize system configuration
nt = size(data,2) ;                                                         % number of time grid
nbeta = size(option.A,1);                                                   % size of beta vector
ntheta= size(option.hyperparameters_possibilities,1);                       % number of possibilities for $theta$



%------ initilize for given theta----------------------------------------
if ~isempty(prior)                                                             % if prior is not defied, the prior in configuration will be used.
    option.beta_prior_mean = prior.beta_mean;
    option.beta_prior_variance = prior.beta_variance;
    option.hyperparameters_possibilities(:,end) = prior.thetha_probability;
end

gridsize = size(option.grids,1);
% global allposibleQz;
global allposibleQzInverse ;
marginal_mux = zeros(gridsize+nbeta,nt);
marginal_sigmax = zeros(gridsize+nbeta,nt);
log_pi_y = zeros(1,ntheta);
mux_theta = zeros(gridsize+nbeta,ntheta);
sigmax_theta = zeros(gridsize+nbeta,ntheta);
y = zeros(option.agentnumbers,nt);
q = zeros(option.agentnumbers,nt);
log_pi_y_without_theta = zeros(1,nt);
marginal_mu_beta = zeros(nbeta,nt);
marginal_sigma_beta = zeros(nbeta^2,nt);

%----- initialize beta(o) -------------------------------------------------
% beta statistics for given theta for time t-1
previous_mu_beta_theta    = option.beta_prior_mean * ones(1,ntheta);
previous_Sigma_beta_theta = ...
    reshape(option.beta_prior_variance,[],1) * ones(1,ntheta);



%------ main loop----------------------------------------------------------
Fs = option.Fs; A = option.A; B = option.B;
hyperparameter_posterior = zeros(ntheta,nt);
for t=1:nt 
    
    %step 0 : get measurements y(t)
    y(:,t) = data(t).y;
    q(:,t) = data(t).q;
    Mq = sparse(gridsize+nbeta,option.agentnumbers);
    Fq = zeros(option.agentnumbers,nbeta);
    for i = 1: option.agentnumbers
        Mq(q(i,t),i) = 1;
        Fq(i,:) = option.Fs(q(i,t),:);
    end
        
    for k = 1:ntheta
        %step 4
        mu_beta_ = A * previous_mu_beta_theta (:,k);
        Sigma_beta = reshape(previous_Sigma_beta_theta(:,k),nbeta,nbeta);
        Sigma_beta_ = A*Sigma_beta*A'+B*B';


        %step 5
%         Qz = reshape(allposibleQz(:,k),gridsize,gridsize);
        mux_ = [Fs*mu_beta_;mu_beta_];
%         QF = Qz * Fs;
%         Qx_ = [Qz -QF;
%         -QF' Fs'*QF+Sigma_beta_^(-1)];
        QzInverese = reshape(allposibleQzInverse(:,k),gridsize,gridsize);     
        sigmax11_= QzInverese + Fs*( Sigma_beta_ )*Fs';                     % from block inverse lemma
        sigmax12_= Fs*Sigma_beta_ ;
        %sigmax21_=sigmax12_';
        %sigmax22_= Sigma_beta_ ;
        Sigmax_ = [sigmax11_,sigmax12_;sigmax12_',Sigma_beta_];     
        Sigmax_Mq = Sigmax_ * Mq;
        %step 7
        muy_ = Fq*mu_beta_;
        Sigmay_ = Mq'*Sigmax_Mq + option.s_e2 * eye(option.agentnumbers);

        %step 6
%         Qx = Qx_ + option.s_e2^(-1) * (Mq * Mq'); 
        
        %Sigmax = Qx^(-1);
        Sigmax = Sigmax_ - Sigmax_Mq * (Sigmay_)^(-1)* Sigmax_Mq';          % using woodbury lemma 
            
                         
        mux = mux_ + option.s_e2^(-1) * Sigmax * Mq *(y(:,t)-muy_);
        previous_mu_beta_theta (:,k)   = mux(gridsize+1:end);
        previous_Sigma_beta_theta(:,k) = ....
            reshape(Sigmax(gridsize+1:end,gridsize+1:end),[],1);
        
        
        
        %step 9
        log_pi_y(1,k) = -0.5*logdet(Sigmay_)- ...
            0.5*(y(:,t)-muy_)'*(Sigmay_\(y(:,t)-muy_));
        mux_theta(:,k) = mux;
        sigmax_theta(:,k) = diag(Sigmax);
         
    end
    
    log_pi_y_without_theta(1,t) = log(exp(log_pi_y)...                      %$ pi(yt|D_{1:t-1},qt) = sum( pi(yt| theta, d_{1:t-1},qt) *
        * option.hyperparameters_possibilities(:,end));                     %      pi(theta|D_{1:t-1})  )
    
    %step 8    
    min_log_pi_y = min(log_pi_y);
    for k = 1:ntheta
        option.hyperparameters_possibilities(k,end) = ...
            exp(-min_log_pi_y + log_pi_y(1,k) + ...
            log(option.hyperparameters_possibilities(k,end)));        
    end
    tmp_c = sum(option.hyperparameters_possibilities(:,end));
    option.hyperparameters_possibilities(:,end) = ...
        option.hyperparameters_possibilities(:,end)/tmp_c;
    
    
    %step 10
    marginal_mux(:,t) = ...
        mux_theta * option.hyperparameters_possibilities(:,end);
    marginal_sigmax(:,t) = ...
        (sigmax_theta +(mux_theta-marginal_mux(:,t)*ones(1,ntheta)).^2)*...
        option.hyperparameters_possibilities(:,end);        
    marginal_mu_beta(:,t) = previous_mu_beta_theta * ...
        option.hyperparameters_possibilities(:,end);
    tmp1 = (previous_mu_beta_theta-marginal_mu_beta(:,t)*ones(1,ntheta));
    tmp2 = zeros(nbeta^2,ntheta);
    for k=1:ntheta
        tmp2(:,k)= reshape(tmp1(:,k)*tmp1(:,k)',[],1);  
    end
    marginal_sigma_beta(:,t) = (previous_Sigma_beta_theta + tmp2)...
        * option.hyperparameters_possibilities(:,end);  
    hyperparameter_posterior(:,t) = option.hyperparameters_possibilities(:,end);
    %------------------
    
end
    %----generate output----
    posterior.beta_mean = marginal_mu_beta;
    posterior.beta_variance = marginal_sigma_beta;
    posterior.thetha_probability = hyperparameter_posterior;
    posterior.x_mean = marginal_mux;
    posterior.x_variance = marginal_sigmax;
    posterior.y_logprobability = log_pi_y_without_theta;
