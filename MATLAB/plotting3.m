for t=1:10
    B1(:,t) = prior(t+1).beta_mean;
end
for n = 1: size(beta,1)
    figure('name',['radial basis function' num2str(n)])
    h = bar([beta(n,:)' , posteriorTrueQ.beta_mean(n,:)',posteriorNoisyQ.beta_mean(n,:)', B1(n,:)']);
    %set(get(h(1),'BaseLine'),'LineWidth',2,'LineStyle',':')
end


for t=1:10
fig_field = figure;
set(fig_field,'position',[100 50 900 800]);
fieldsize = size(true_field(:,:,1));
max_val = -inf;
min_val =  inf;
max_var = -inf;
max_err = -inf;

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


%% Edit plots
est_colorbar_step = 0.025;
est_caxis = [floor(min_val/est_colorbar_step)*est_colorbar_step, ...
                ceil(max_val/est_colorbar_step)*est_colorbar_step];
set(get(plt_trueposition_estimatedfield,'parent'),'clim',est_caxis);
set(get(plt_noisyposition_estimatedfield,'parent'),'clim',est_caxis);

var_colorbar_step = 0.005;
var_caxis = [0, ceil(max_var/var_colorbar_step)*var_colorbar_step];
set(get(plt_trueposition_variancefield,'parent'),'clim',var_caxis);
set(get(plt_noisyposition_variancefield,'parent'),'clim',var_caxis);

err_colorbar_step = 0.005;
err_caxis = [0, ceil(max_err/err_colorbar_step)*err_colorbar_step];
set(get(plt_trueposition_errorfield,'parent'),'clim',err_caxis);
set(get(plt_noisyposition_errorfield,'parent'),'clim',err_caxis);

pause(1)
saveas(fig_field,[ num2str(t) 'frame.jpg']);
end


mean(MeanSquareError)
mean(MeanSquareVariance)