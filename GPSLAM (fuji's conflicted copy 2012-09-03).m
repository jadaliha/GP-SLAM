% =========================================================================
%       ---------------------------------------------------------
%        Test the effect of Tetha and its convergance properties
%       ---------------------------------------------------------
% =========================================================================
                                                                            
                                                                            
profile clear    
profile on
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
gridsize = size(option.grids,1);                                            % number of spatial sites
N = option.agentnumbers;
n = gridsize;

%----- Construct precision matrix Qz---------------------------------------
muz_theta = zeros(n,ntheta);
Sigmaz_theta = zeros(n^2,ntheta);
for k=1:ntheta 
    tmp1 = prec_mat_5by5(option.hyperparameters_possibilities(k,1)...
        ,option.hyperparameters_possibilities(k,2),nx,ny);
    Sigmaz_theta(:,k) = reshape(tmp1^(-1),[],1);
end 



f_theta = log(option.hyperparameters_possibilities(:,end));
itmp = (1:2*N);
jtmp = kron((1:N)',[1;1]);
mux = zeros(2*N,1);
Sigmax = eye(2*N)*1000;
for t=option.T_mesh 
    xtilda = reshape(option.grids(...
        AgumentedData(1,t).possible_q.measuredposition,:)',[],1);
    ztilda = AgumentedData(t).y;
    u = VehicleStates(t).u;
    
    qt= AgumentedData(1,t).possible_q.support_qt{1,1}';
    for tmpindex = 2:N
        qtN = AgumentedData(1,t).possible_q.support_qt{1,tmpindex}';
        qt= [kron(ones(1,size(qtN,2)),qt);...
            kron(qtN,ones(1,size(qt,2)))];
    end

    phi = VehicleStates(1,t).h';
    stmp = reshape([sin(phi),cos(phi)]',[],1);
    F = sparse(itmp,jtmp,stmp,2*N,N);
    Sigma_e = eye(2*size(qt,1))* option.vehicle.ObservationNoise;
    Sigma_w = eye(2*size(qt,1))* option.vehicle.ModelUncertanity;
    mux_    = mux + F * u;
    Sigmax_ = Sigmax + Sigma_w;
    
    numberofpossibleqt = size(qt,2);
    Sigma_epsilon = eye(size(qt,1))* option.s_e2;
    f_qtandtheta = zeros(numberofpossibleqt,ntheta);
    muz_qandtheta =    zeros(n, numberofpossibleqt*ntheta);
    Sigmaz_qandtheta = zeros(n^2,numberofpossibleqt*ntheta);
    for qindex = 1:numberofpossibleqt
        q = qt(:,qindex);
        H = sparse(1:N,q',ones(1,N),N,gridsize);                            % find the map $H_{qt}$ from $q{t}$ to spacial sites S.
        x = reshape(option.grids(q,:)',[],1);
        f_qt  = logmvnpdf(x,mux_,Sigmax_);
        f_qtildaGqt  = logmvnpdf(xtilda,x,Sigma_e);
        for k = 1:ntheta
            muztheta = muz_theta(:,k);
            Sigmaztheta      = reshape(Sigmaz_theta(:,k),n,n);
            muztildatheta    = H * muztheta;                                % $mu_{tilde{zt}|theta,D{t-1},\q{t}}  =  H_{qt} mu_{z| theta,D{t-1}}$            
            Sigmaztildatheta = Sigma_epsilon + H * Sigmaztheta * H';
            f_ztilda = logmvnpdf(ztilda,muztildatheta,Sigmaztildatheta);
            % approximate the distribution of $\q{t},\theta |\D{0}{t}$           
            f_qtandtheta(qindex,k) = ...
                f_qtildaGqt  + f_theta(k) + f_qt + f_ztilda;                % pi(q{t},theta|D{t})
            
            
            

            InvSigmaztildatheta = Sigmaztildatheta^(-1);
            
            muz_qandtheta(:,(k-1)*numberofpossibleqt + qindex)    = ...
                muztheta + Sigmaztheta * H' * ...
                InvSigmaztildatheta * (ztilda - muztildatheta);
            Sigmaz_qandtheta(:,(k-1)*numberofpossibleqt + qindex) = ...
                reshape(Sigmaztheta - Sigmaztheta * H' * ...
                InvSigmaztildatheta * H * Sigmaztheta,[],1);
            
            
         
        end
    end
    tmp_c = log(sum(sum(exp(f_qtandtheta))));
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
    
    Sigmaz_theta = zeros(n^2,ntheta);
    for indtheta = 1:ntheta
        pi_qtGtheta = (pi_qtandtheta(:,indtheta)/pi_theta(indtheta));
        index = (indtheta-1)*numberofpossibleqt + (1:numberofpossibleqt);
        muz_theta(:,indtheta) = muz_qandtheta(:,index) ...
            * pi_qtGtheta;

        Sigmaz_theta(:,indtheta) = ...
            (Sigmaz_qandtheta(:,index)+...
            VectorSquare(muz_qandtheta(:,index)))*pi_qtGtheta - ...
            VectorSquare(muz_theta(:,indtheta));
    end
    
    posterior(t).muz        = muz_theta * pi_theta';
    posterior(t).Varz       = (Sigmaz_theta(1:n+1:end,:)+ ...
        muz_theta .* muz_theta)* pi_theta' - ...
        posterior(t).muz .* posterior(t).muz;
    posterior(t).mux        = mux;
    posterior(t).Sigmax = zeros(N,4);
    for indj = 1:N
        posterior(t).Sigmax(N,:)     = reshape(...
            Sigmax(indj*2-1:indj*2,indj*2-1:indj*2)',1,[]);
    end
    
    
    
end
% 
% if (sys_parameter.plot_enabled)
%     fig_field = figure;
%     set(fig_field,'position',[10 50 400 300]);
%     hold on
%     t = 1;
%     plt_truefield = imagesc(option.X_mesh,option.Y_mesh,reshape(muxtime(2*N+1:end,t),nx,ny));
%     for n= 1:N
%         plt_robots_pos(n) = plot(muxtime(2*n-1,t),...
%             muxtime(2*n,t),'m','marker','o');
%     end
%     xlim([option.X_mesh(1), option.X_mesh(end)]);
%     ylim([option.Y_mesh(1), option.Y_mesh(end)]);
%     colorbar; title('GMRF');    
%     colormap('gray');
%     for t = option.T_mesh
%         set(plt_truefield,'CData',reshape(muxtime(2*N+1:end,t),nx,ny));
%         for n= 1:option.agentnumbers
%             set(plt_robots_pos(n),'XData',muxtime(2*n-1,t),...
%                 'YData',muxtime(2*n,t));
%       
% 
%         end        
%         refreshdata(fig_field);
%         pause(0.02);
%     end
% end
% 
profile viewer
profile off
