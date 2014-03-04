% function plotting2()
close all
clc

nt = size(true_field,3);
nx = size(true_field,2); % number of X grid
ny = size(true_field,1); % number of Y grid

% ----- Existance of parameters -------------------------------------------
existbeta = exist('beta','var');
if existbeta
    existbeta = ~isempty(beta);
end


% ----- Main code ---------------------------------------------------------
if existbeta % if the true beta is known
    for t=1:nt
        fig_field = figure('name',['time = ', num2str(t)]);
        set(fig_field,'position',[10 50 900 220]);
        subplot 131
        plt_truefield = imagesc(option.X_mesh,option.Y_mesh,true_field(:,:,t));
        colorbar; title('(a)'); axis off;
        tmp1 = get(get(plt_truefield,'parent'),'position');
        tmp1(1) = 0.02; tmp1(2) = 0.08; tmp1(3) = 0.20;
        set(get(plt_truefield,'parent'),'position',tmp1);
        subplot 132
        plt_rb = imagesc(option.X_mesh,option.Y_mesh,reshape(Fs*beta(:,t),ny,nx));
        colorbar; title('(b)');  axis off;
        tmp1(1) = 0.36; tmp1(2) = 0.08; tmp1(3) = 0.20;
        set(get(plt_rb,'parent'),'position',tmp1);    

        subplot 133
        plt_gmrf = imagesc(option.X_mesh,option.Y_mesh,true_field(:,:,t)-...
            reshape(Fs*beta(:,t),ny,nx));
        colorbar; title('(c)');  axis off;
        tmp1(1) = 0.70; tmp1(2) = 0.08; tmp1(3) = 0.20;
        set(get(plt_gmrf,'parent'),'position',tmp1);  
    end


for t=1:(nt+1)
    B1(:,t) = prior(t).beta_mean;
end
% for n = 1: size(beta,1)
%     figure('name',['radial basis function' num2str(n)])
%     h = bar([beta(n,:)' , posteriorTrueQ.beta_mean(n,:)',posteriorNoisyQ.beta_mean(n,:)', B1(n,:)']);
%     %set(get(h(1),'BaseLine'),'LineWidth',2,'LineStyle',':')
% end
figure('name','radial basis function RMS Error & Hyperparameter vector posterior probability')
subplot(2,1,1)
plt_beta = plot((0:nt)',sqrt([...
    sum((beta - [option.beta_prior_mean, posteriorTrueQ.beta_mean]).^2)',...
    sum((beta - [option.beta_prior_mean, posteriorNoisyQ.beta_mean]).^2)',...
    sum((beta - B1).^2)']));
title('(a)');
set(plt_beta(1),'markersize',10,'linewidth',2,...
    'Marker','o','LineStyle',':');
set(plt_beta(2),'markersize',10,'linewidth',2,...
    'Marker','pentagram','LineStyle','--');
set(plt_beta(3),'markersize',10,'linewidth',2,...
    'Marker','square','LineStyle','-.');
for n = 1:(nt+1)
    theta_probability(n) = prior(n).thetha_probability(1,end);
end
subplot(2,1,2)
plt_theta = plot((0:nt)',...
           [...
           [option.hyperparameters_possibilities(1,end), ...
             posteriorTrueQ.thetha_probability(1,:)];...
             [option.hyperparameters_possibilities(1,end), ...
             posteriorNoisyQ.thetha_probability(1,:)];
             theta_probability]');
title('(b)')         
set(plt_theta(1),'markersize',10,'linewidth',2,...
    'Marker','o','LineStyle',':');
set(plt_theta(2),'markersize',10,'linewidth',2,...
    'Marker','pentagram','LineStyle','--');
set(plt_theta(3),'markersize',10,'linewidth',2,...
    'Marker','square','LineStyle','-.');
end

for t=1:nt
fig_field = figure('name',['time = ', num2str(t)]);
set(fig_field,'position',[100 50 900 800]);
fieldsize = size(true_field(:,:,1));
max_val = -inf;
min_val =  inf;
max_var = -inf;
max_err = -inf;
gridsize = size(option.grids,1);

%% True Position
        %----plot estimated field & support of positions-------------------
subplot(3,3,1)
estimated_field = reshape(posteriorTrueQ.x_mean(1:gridsize,t),fieldsize);
max_val = max(max(max(estimated_field)),max_val);
min_val = min(min(min(estimated_field)),min_val);
plt_trueposition_estimatedfield= ...
    imagesc(option.X_mesh,option.Y_mesh,estimated_field);
axis off; title('(a)');

        %----plot estimation variance--------------------------------------
subplot(3,3,4)
estimated_variance = ...
    reshape(posteriorTrueQ.x_variance(1:gridsize,t),fieldsize);
max_var = max(max(max(estimated_variance)),max_var);
plt_trueposition_variancefield = ...
    imagesc(option.X_mesh,option.Y_mesh,estimated_variance);
axis off; title('(d)');
MeanSquareVariance(t,1) = mean(mean( estimated_variance));
        %----plot estimation error-----------------------------------------
subplot(3,3,7)
estimation_error = (true_field(:,:,t) - estimated_field).^2;
plt_trueposition_errorfield = ...
    imagesc(option.X_mesh,option.Y_mesh,estimation_error);
max_err = max(max(max(estimation_error)),max_err);
axis off; title('(g)');
hold on
q  = AgumentedData(t).possible_q.true_q;
plot(option.grids(q,1),option.grids(q,2),....
    'm','marker','o','linestyle','none');

MeanSquareError(t,1) = mean(mean( estimation_error));

%% Noisy Position
        %----plot estimated field & support of positions-------------------
subplot(3,3,2)
estimated_field = reshape(posteriorNoisyQ.x_mean(1:gridsize,t),fieldsize);
max_val = max(max(max(estimated_field)),max_val);
min_val = min(min(min(estimated_field)),min_val);
plt_noisyposition_estimatedfield= ...
    imagesc(option.X_mesh,option.Y_mesh,estimated_field);
axis off; 
title('(b)');

        %----plot estimation variance--------------------------------------
subplot(3,3,5)
estimated_variance = ...
    reshape(posteriorNoisyQ.x_variance(1:gridsize,t),fieldsize);
max_var = max(max(max(estimated_variance)),max_var);
plt_noisyposition_variancefield = ...
    imagesc(option.X_mesh,option.Y_mesh,estimated_variance);
axis off; title('(e)');
MeanSquareVariance(t,2) = mean(mean( estimated_variance));
        %----plot estimation error-----------------------------------------
subplot(3,3,8)
estimation_error = (true_field(:,:,t) - estimated_field).^2;
max_err = max(max(max(estimation_error)),max_err);
plt_noisyposition_errorfield = ...
    imagesc(option.X_mesh,option.Y_mesh,estimation_error);
axis off; title('(h)');
hold on
for i=1:option.agentnumbers
    if ~(AgumentedData(t).possible_q.prior_qt(i,:)==1)
        q = AgumentedData(t).possible_q.measuredposition(i);
        plot(option.grids(q,1),option.grids(q,2),...
            'm','marker','*');
    else
        q  = AgumentedData(t).possible_q.true_q(i);
        plot(option.grids(q,1),option.grids(q,2),....
            'm','marker','o','linestyle','none');
    end
end
MeanSquareError(t,2) = mean(mean( estimation_error));

%% Uncertain Position
        %----plot estimated field & support of positions-------------------
subplot(3,3,3)
estimated_field = reshape(prior(t+1).x_mean(1:gridsize),fieldsize);
max_val = max(max(max(estimated_field)),max_val);
min_val = min(min(min(estimated_field)),min_val);
plt_uncertainposition_estimatedfield= ...
    imagesc(option.X_mesh,option.Y_mesh,estimated_field);
axis off;  caxis([-1 1]); colorbar; title('(c)'); 

        %----plot estimation variance--------------------------------------
subplot(3,3,6)
estimated_variance = reshape(prior(t+1).x_variance(1:gridsize),fieldsize);
max_var = max(max(max(estimated_variance)),max_var);
plt_uncertainposition_variancefield = ...
    imagesc(option.X_mesh,option.Y_mesh,estimated_variance);
colorbar; title('(f)'); axis off;
caxis([0 0.4]);
MeanSquareVariance(t,3) = mean(mean( estimated_variance));
        %----plot estimation error-----------------------------------------
subplot(3,3,9)
estimation_error = (true_field(:,:,t) - estimated_field).^2;
max_err = max(max(max(estimation_error)),max_err);
plt_uncertainposition_errorfield = ...
    imagesc(option.X_mesh,option.Y_mesh,estimation_error);
hold on
caxis([0 0.4]); colorbar; title('(i)'); axis off;
for i=1:option.agentnumbers
    if ~(AgumentedData(t).possible_q.prior_qt(i,:)==1)
        q = AgumentedData(t).possible_q.support_qt(i,:);
        q = [q(1) q(2) q(4) q(3) q(1)];
        plot(option.grids(q,1),option.grids(q,2),'m');
    else
        q  = AgumentedData(t).possible_q.true_q(i);
        plot(option.grids(q,1),option.grids(q,2),....
            'm','marker','o','linestyle','none');
    end
end
MeanSquareError(t,3) = mean(mean( estimation_error));

%% Edit plots
est_colorbar_step = 0.25;
est_caxis = [floor(min_val/est_colorbar_step)*est_colorbar_step, ...
                ceil(max_val/est_colorbar_step)*est_colorbar_step];
set(get(plt_trueposition_estimatedfield,'parent'),'clim',est_caxis);
set(get(plt_noisyposition_estimatedfield,'parent'),'clim',est_caxis);
set(get(plt_uncertainposition_estimatedfield,'parent'),'clim',est_caxis);

var_colorbar_step = 0.001;
var_caxis = [0, ceil(max_var/var_colorbar_step)*var_colorbar_step];
set(get(plt_trueposition_variancefield,'parent'),'clim',var_caxis);
set(get(plt_noisyposition_variancefield,'parent'),'clim',var_caxis);
set(get(plt_uncertainposition_variancefield,'parent'),'clim',var_caxis);

err_colorbar_step = 0.05;
err_caxis = [0, ceil(max_err/err_colorbar_step)*err_colorbar_step];
set(get(plt_trueposition_errorfield,'parent'),'clim',err_caxis);
set(get(plt_noisyposition_errorfield,'parent'),'clim',err_caxis);
set(get(plt_uncertainposition_errorfield,'parent'),'clim',err_caxis);

tmp1 = get(get(plt_uncertainposition_estimatedfield,'parent'),'position');
tmp1(3) = tmp1(3)+0.07;
set(get(plt_uncertainposition_estimatedfield,'parent'),'position',tmp1);

tmp1 = get(get(plt_uncertainposition_variancefield,'parent'),'position');
tmp1(3) = tmp1(3)+0.07;
set(get(plt_uncertainposition_variancefield,'parent'),'position',tmp1);

tmp1 = get(get(plt_uncertainposition_errorfield,'parent'),'position');
tmp1(3) = tmp1(3)+0.07;
set(get(plt_uncertainposition_errorfield,'parent'),'position',tmp1);

pause(1)
saveas(fig_field,[ num2str(t) 'frame.jpg']);
end



mean(MeanSquareError)
mean(MeanSquareVariance)
