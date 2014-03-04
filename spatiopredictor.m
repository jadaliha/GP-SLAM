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
clear all;
close all; pause(0.5);
clc;

[ option , sys_parameter ] = Configuration(); % initialize system configuration


nx = size(option.X_mesh,2) ; % number of X grid
ny = size(option.Y_mesh,2) ; % number of Y grid
nt = size(option.T_mesh,2) ; % number of time grid
nbeta = size(option.A,1);    % size of beta vector
ntheta= size(option.hyperparameters_possibilities,1);

%----- generate or load the field--------
true_field          = Generate_new_field();
estimated_field     = zeros(ny,nx,nt);
estimated_variance  = zeros(ny,nx,nt);
estimation_error    = zeros(ny,nx,nt);




%---- sensor network of agents------
y = zeros(option.agentnumbers,nt); %measurements
gridsize = size(option.grids,1);
q = ceil(gridsize * rand(option.agentnumbers,nt)+1);
measurement_noise = randn(option.agentnumbers,nt) * sqrt(option.s_e2);
% % with movment constrain
% for t = 2:nt
%     
% end

% main loop
Fs = option.Fs; A = option.A; B = option.B;

fig_field = figure;
set(fig_field,'position',[100 50 600 800]);
subplot 411
plt_truefield     = imagesc(option.X_mesh,option.Y_mesh,true_field(:,:,1));
hold on
plt_agents = plot(option.grids(q(:,1),1),option.grids(q(:,1),2),'w','marker','o','linestyle','none');
colorbar; title('True field');
subplot 412
plt_estimatedfield= imagesc(option.X_mesh,option.Y_mesh,estimated_field(:,:,1));
colorbar; title('Estimated');
subplot 413
plt_variancefield = imagesc(option.X_mesh,option.Y_mesh,estimated_variance(:,:,1));
colorbar; title('Variance');
subplot 414
plt_errorfield    = imagesc(option.X_mesh,option.Y_mesh,estimation_error(:,:,1));
colorbar; title('Error^2');

%-------- Initilize for given theta -------
previous_mu_beta_theta = zeros(nbeta,ntheta);
previous_Sigma_beta_theta = zeros(nbeta^2,ntheta);
allposibleQz = sparse(gridsize^2,ntheta);
allposibleQzInverse = zeros(gridsize^2,ntheta);
for k=1:ntheta
    %---- Construct precision matrix Qz
    allposibleQz(:,k) = reshape(prec_mat_5by5(option.hyperparameters_possibilities(k,1),...
        option.hyperparameters_possibilities(k,1),nx,ny),[],1); 
    allposibleQzInverse(:,k) = reshape(prec_mat_5by5(option.hyperparameters_possibilities(k,1),...
        option.hyperparameters_possibilities(k,1),nx,ny)^(-1),[],1);     
    %----- initialize beta(o) --------------------
    % beta statistics for given theta for time t-1
    previous_mu_beta_theta (:,k) = ones(nbeta,1).*option.beta_prior_mean;
    previous_Sigma_beta_theta(:,k) =  reshape(eye(nbeta,nbeta)*option.beta_prior_variance,[],1);
end
%----------------------------------    
marginal_mux = zeros(gridsize+nbeta,nt);
marginal_sigmax = zeros(gridsize+nbeta,nt);
log_pi_y = zeros(1,ntheta);
mux_theta = zeros(gridsize+nbeta,ntheta);
sigmax_theta = zeros(gridsize+nbeta,ntheta);
for t=1:nt 
    
    %step 0 : get measurements y(t)
    tmp_field_t = reshape(true_field(:,:,t),[],1);
    y(:,t) = tmp_field_t(q(:,t))  + measurement_noise(:,t);
    
    for k = 1:ntheta
        %step 4
        mu_beta_ = A * previous_mu_beta_theta (:,k);
        Sigma_beta = reshape(previous_Sigma_beta_theta(:,k),nbeta,nbeta);
        Sigma_beta_ = A*Sigma_beta*A'+B*B';


        %step 5
        Qz = reshape(allposibleQz(:,k),nx*ny,nx*ny);
        mux_ = [Fs*mu_beta_;mu_beta_];
        QF = Qz * Fs;
        Qx_ = [Qz -QF;
        -QF' Fs'*QF+Sigma_beta_];

        %step 6
        Mq = sparse(gridsize+nbeta,option.agentnumbers);
        Fq = zeros(option.agentnumbers,nbeta);
        for i = 1: option.agentnumbers
            Mq(q(i,t),i) = 1;
            Fq(i,:) = option.Fs(q(i,t),:);
        end
        Qx = Qx_ + option.s_e2^(-1) * Mq * Mq'; 
%         Qz_inv = reshape(allposibleQzInverse(:,k),nx*ny,nx*ny);
%         U = option.s_e2^(-0.5) * Mq(1:gridsize,:);
%         A_inv = Qz_inv - Qz_inv *
%         U*((eye(option.agentnumbers)+U'*Qz_inv*U)\U')*Qz_inv;
        Sigmax = Qx^(-1);
        mux = mux_ + option.s_e2^(-1) * (Qx\Mq) *(y(:,t)-Fq*mu_beta_);
        previous_mu_beta_theta (:,k)   = mux(gridsize+1:end);
        previous_Sigma_beta_theta(:,k) = reshape(Sigmax(gridsize+1:end,gridsize+1:end),[],1);
        

        %step 7
        muy_ = Fq*mu_beta_;
        Sigmay_ = Mq'*(Qx_\Mq) + option.s_e2;
        Sigmay_ = 0.5*(Sigmay_+Sigmay_');
        
        %step 8
        log_pi_y(1,k) = -0.5*logdet(Sigmay_)-0.5*(y(:,t)-muy_)'*(Sigmay_\(y(:,t)-muy_));
        mux_theta(:,k) = mux;
        sigmax_theta(:,k) = diag(Sigmax);
         
    end
    
    %step 10    
    min_log_pi_y = min(log_pi_y);
    for k = 1:ntheta
        option.hyperparameters_possibilities(k,end) = exp(-min_log_pi_y + ...
            log_pi_y(1,k) + log(option.hyperparameters_possibilities(k,end)));        
    end
    tmp_c = sum(option.hyperparameters_possibilities(:,end));
    option.hyperparameters_possibilities(:,end) = option.hyperparameters_possibilities(:,end)/tmp_c;
    marginal_mux(:,t) = mux_theta * option.hyperparameters_possibilities(:,end);
    tmp_diagonal = sparse(diag(option.hyperparameters_possibilities(:,end)));
    
    marginal_sigmax(:,t) = (sigmax_theta * option.hyperparameters_possibilities(:,end)+...
        diag((mux_theta-marginal_mux(:,t)*ones(1,ntheta)) *tmp_diagonal* (mux_theta-marginal_mux(:,t)*ones(1,ntheta))'));        
    
    %----plot field----
    estimated_field(:,:,t)    = reshape(marginal_mux(1:gridsize,t),ny,nx);
    estimated_variance(:,:,t) = reshape(marginal_sigmax(1:gridsize,t),ny,nx);
    estimation_error(:,:,t)   = (true_field(:,:,t)-estimated_field(:,:,t)).^2;
    set(plt_truefield,'CData',true_field(:,:,t));
    set(plt_agents,'XData',option.grids(q(:,t),1));
    set(plt_agents,'YData',option.grids(q(:,t),2));
    set(plt_estimatedfield,'CData',estimated_field(:,:,t));
    set(plt_variancefield,'CData',estimated_variance(:,:,t));
    set(plt_errorfield,'CData',estimation_error(:,:,t));
    refreshdata(fig_field);
    pause(1);
    saveas(fig_field,['./fig/' num2str(t) 'frame.jpg']);

    %------------------
    
end

% 
% % -------------------------------------------------------------------------
% % Bayesian Inference
% % -------------------------------------------------------------------------
% % Initialize display and plots
% fprintf(1,'==========================================\n');
% fprintf(1,'%-5s %-10s %-10s %-10s\n',...
%     'time','beta','variance','rms error');
% h1 = figure;
% set(gcf,'position',[600 600 600 200]);
% h2 = figure;
% set(gcf,'position',[600 200 600 200]);
% 
% % Initialize Q, b and C_diag
% Qz = prec_mat_5by5(kappa,alpha);
% QF = Qz * F;
% Q = [Qz -QF;
%     -QF' F'*QF+T];
% b = zeros(n+p,1);
% % C_diag = diag(Q_star\speye(n+p));
% temp = Qz\speye(n,1);
% C_diag(1:n,1) = temp(1)*ones(n,1)+F*(inv(T)*F(1,:)');
% C_diag(n+1:n+p,1) = diag(inv(T));
% 
% % Stop time
% tf = 1;
% % Plot time
% tplot = [1 5 20 tf];
% % Initialize sampling locations
% q = zeros(tf,1);
% % Initialize results
% est_beta = zeros(tf,p);
% var_beta = zeros(tf,p);
% rms = zeros(tf,1);
% aver_var = zeros(tf,1);
% 
% % Start iteration
% for t = 1:tf
%     
%     % Select sampling locations
%     if t == 1
% %         q(t) = index(r+1,r+1);
%         q(t) = index(30,20);
%     else
%         q(t) = planning(q(t-1),C_diag);
%     end
%     
%     % Noisy observations
%     X(t,:) = S(q(t),:);
%     y(t,1) = z_true(q(t)) + sig_w*randn(1);
%     
%     % Update C_diag
%     u = zeros(n+p,1);
%     u(q(t)) = 1/sig_w;
%     h = Q\u;
%     C_diag = C_diag - h.^2/(1+u'*h);
%     
%     % Updated Q and b
%     Q(q(t),q(t)) = Q(q(t),q(t)) + 1/sig_w^2;
%     b(q(t)) = b(q(t)) + y(t)/sig_w^2;
%     x = Q\b;
%     
%     % Compute estimate of beta and RMS error
%     est_beta(t,:) = x(end-p+1:end);
%     var_beta(t,:) = C_diag(end-p+1:end);
%     rms(t) = sqrt(sum((x(arr)-z_true(arr)).^2)/length(arr));
%     aver_var(t) = sqrt(sum(C_diag(arr))/length(arr));
%     
%     % Display results
%     fprintf(1,'%-5i %-10.2f %-10.2f %-10.2f \n',...
%         t,est_beta(t,:),var_beta(t,:),rms(t));
%     
%     % Plot predicted field
%     figure(h1);
%     imagesc(s1,s2,reshape(x(1:n),n2,n1));
%     colorbar;
%     caxis(clim);
% %     hold on;
% %     line([r+1 r+1 n1-r n1-r r+1],[r+1 n2-r n2-r r+1 r+1],'color','k');
% %     hold off;
%     axis([r+1-0.5 n1-r+0.5 r+1-0.5 n2-r+0.5]);
%     pause(0.5);
%     
%     % Plot prediction error variance
%     figure(h2);
%     imagesc(s1,s2,reshape(C_diag(1:n),n2,n1));
%     colorbar;
%     hold on;
%     plot(X(:,1),X(:,2),'w');
%     plot(X(:,1),X(:,2),'wo');
%     plot(X(end,1),X(end,2),'o',...
%         'markeredgecolor','w','markerfacecolor','w');
% %     line([r+1 r+1 n1-r n1-r r+1],[r+1 n2-r n2-r r+1 r+1],'color','k');
%     hold off;
%     axis([r+1-0.5 n1-r+0.5 r+1-0.5 n2-r+0.5]);
%     pause(0.5);
%     
%     % Save plots
%     if plot_flag == 1 && ismember(t,tplot)
%         exportfig(h1,strcat('./fig/prediction-',num2str(t),'.eps'),...
%             'width',6.5,'height',2,'color','rgb',...
%             'fontmode','fixed','fontsize',14,...
%             'linemode','fixed','linewidth',2);
%         exportfig(h2,strcat('./fig/variance-',num2str(t),'.eps'),...
%             'width',6.5,'height',2,'color','rgb',...
%             'fontmode','fixed','fontsize',14,...
%             'linemode','fixed','linewidth',2);
%     end
%     
% end
% 
% % -------------------------------------------------------------------------
% % PLOT RESULTS
% % -------------------------------------------------------------------------
% % Plot estimated beta
% figure
% hold on;
% beta_grid = beta-20:0.1:beta+20;
% str_legend = cell(length(tplot),1);
% for i = 1:length(tplot)
%     t = tplot(i);
%     plot(beta_grid,normpdf(beta_grid,est_beta(t),var_beta(t)),color{i});
%     str_legend{i} = strcat('t=',num2str(t));
% end
% hold off;
% box on;
% axis square;
% xlabel('$\beta$','interpreter','latex');
% legend(str_legend);
% if plot_flag == 1
%     exportfig(gcf,'./fig/beta.eps',...
%         'width',6.5,'height',6.5,'color','rgb',...
%         'fontmode','fixed','fontsize',14,...
%         'linemode','fixed','linewidth',2);
% end
% 
% % Plot RMS error
% figure
% plot(rms,'marker','.');
% axis square;
% xlabel('$t$','interpreter','latex');
% if plot_flag == 1
%     exportfig(gcf,'./fig/rms.eps',...
%         'width',6.5,'height',6.5,'color','rgb',...
%         'fontmode','fixed','fontsize',14,...
%         'linemode','fixed','linewidth',2);
% end
% 
% % Plot average variance
% figure
% plot(aver_var,'marker','.');
% axis square;
% xlabel('$t$','interpreter','latex');
% if plot_flag == 1
%     exportfig(gcf,'./fig/aver_var.eps',...
%         'width',6.5,'height',6.5,'color','rgb',...
%         'fontmode','fixed','fontsize',14,...
%         'linemode','fixed','linewidth',2);
% end