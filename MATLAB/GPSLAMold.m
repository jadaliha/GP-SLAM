% =========================================================================
%       ---------------------------------------------------------
%        Test the effect of Tetha and its convergance properties
%       ---------------------------------------------------------
% =========================================================================
tic                                                                            
clear all
clc
% profile clear    
% profile on
%----- initialize system configuration-------------------------------------
global globaloption  sys_parameter
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
n = size(option.grids,1);                                            % number of spatial sites
N = option.agentnumbers;


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
    t
    xtilda = reshape(option.grids(...
        AgumentedData(1,t).possible_q.measuredposition,:)',[],1);
    ztilda = AgumentedData(t).y;

    
    qt= AgumentedData(1,t).possible_q.support_qt{1,1}';
    for tmpindex = 2:N
        qtN = AgumentedData(1,t).possible_q.support_qt{1,tmpindex}';
        qt= [kron(ones(1,size(qtN,2)),qt);...
            kron(qtN,ones(1,size(qt,2)))];
    end
    Sigma_e = eye(2*size(qt,1))* option.vehicle.ObservationNoise;


    
    numberofpossibleqt = size(qt,2);
    Sigma_epsilon = eye(size(qt,1))* option.s_e2;
    f_qtandtheta = zeros(numberofpossibleqt,ntheta);
    muz_qandtheta =    zeros(n, 1);
    Sigmaz_qandtheta = zeros(n^2,1);
    for qindex = 1:numberofpossibleqt
        q = qt(:,qindex);
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
            f_qtandtheta(qindex,indtheta) = ...
                f_qtildaGqt  + f_theta(indtheta) + f_qt + f_ztilda;                % pi(q{t},theta|D{t})
            

        end
    end
    tmp_c = log(sum(sum(exp(f_qtandtheta))));
    if (abs(tmp_c)>100)
        pause
    end
    f_qtandtheta = f_qtandtheta - tmp_c;
    pi_qtandtheta = exp(f_qtandtheta);
    
    pi_theta = sum(pi_qtandtheta,1);
    f_theta = log(pi_theta);
    pi_q = sum(pi_qtandtheta,2);
    
    for indj = 1:N
        mux(indj*2-1:indj*2,1) = ...
            sum((pi_q * ones(1,2)) .* option.grids(qt(indj,:),:));
    end
    mux2 = zeros(2*N,2*N);
    for indq = 1:numberofpossibleqt
        tmp1 = reshape(option.grids(qt(:,indq),:)',[],1);
        mux2 = mux2 + tmp1*tmp1'*pi_q(indq,1);
    end
    Sigmax = mux2 - mux * mux';


    for indtheta = 1:ntheta
        pi_qtGtheta = (pi_qtandtheta(:,indtheta)/pi_theta(indtheta));
        muz_theta_temp = parallel.gpu.GPUArray.zeros(n,1);
        Sigmaz_theta_temp = parallel.gpu.GPUArray.zeros(n^2,1);
        for qindex = 1:numberofpossibleqt
            q = qt(:,qindex);
            H = sparse(1:N,q',ones(1,N),N,n);                               % find the map $H_{qt}$ from $q{t}$ to spacial sites S.
            muztheta = muz_theta(:,indtheta);
            Sigmaztheta      = reshape(Sigmaz_theta(:,indtheta),n,n);
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
                muz_qandtheta * pi_qtGtheta(qindex);

            Sigmaz_theta_temp = Sigmaz_theta_temp + ...
                (Sigmaz_qandtheta + VectorSquare(muz_qandtheta))*...
                pi_qtGtheta(qindex) ;

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
    u = VehicleStates(t).u;
    phi = VehicleStates(1,t).h';
    stmp = reshape([cos(phi),sin(phi)]',[],1);
    F = sparse(itmp,jtmp,stmp,2*N,N);
    mux_    = mux + F * u;
    Sigma_w = eye(2*size(qt,1))* option.vehicle.ModelUncertanity;
    Sigmax_ = Sigmax + Sigma_w;  
    
if (sys_parameter.plot_enabled)
    if ~exist('fig_field','var')
        fig_field = figure;
        set(fig_field,'position',[10 50 400 300]);
        hold on
        t = 1;
        plt_truefield = imagesc(option.X_mesh,option.Y_mesh,reshape(posterior(t).muz,nx,ny));
        for jind= 1:N
            plt_robots_pos(jind) = plot(posterior(t).mux (2*jind-1),...
                posterior(t).mux (2*jind),'m','marker','.');
            plt_robots_unc(jind) = ...
                error_ellipse(posterior(t).Sigmax(2*jind-1:2*jind,2*jind-1:2*jind),...
                posterior(t).mux (2*jind-1:2*jind));
        end
        xlim([option.X_mesh(1), option.X_mesh(end)]);
        ylim([option.Y_mesh(1), option.Y_mesh(end)]);
        colorbar; title('GMRF');    
        colormap('gray');
    end
    
        set(plt_truefield,'CData',reshape(posterior(t).muz,nx,ny));
        for jind= 1:N
            set(plt_robots_pos(jind),'XData',posterior(t).mux (2*jind-1),...
                'YData',posterior(t).mux (2*jind));
      
            delete(plt_robots_unc(jind))
            plt_robots_unc(jind) = ...
                error_ellipse(posterior(t).Sigmax(2*jind-1:2*jind,2*jind-1:2*jind),...
                posterior(t).mux (2*jind-1:2*jind));
        end        
        refreshdata(fig_field);
%         pause(0.2);
    
end    
    
end
% 
toc

% 
% profile viewer
% profile off
