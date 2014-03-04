% function plotting2()
close all
clear all
clc
load('test13alien.mat')



trajectory_time = (10:30);

x_width = (option.X_mesh(end) - option.X_mesh(1) + option.finegridsize);
x_widthDIV2 = x_width/2;
y_width = (option.Y_mesh(end) - option.Y_mesh(1) + option.finegridsize);
y_widthDIV2 = y_width/2;
fig_field = figure;
set(fig_field,'position',[10 50 400 300]);
hold on
plt_truefield = imagesc(option.X_mesh,option.Y_mesh,true_field(:,:,1));

% set(fig_field,'YDir','reverse');

x = option.grids(AgumentedData(trajectory_time(1)).possible_q.true_q,1)';
y = option.grids(AgumentedData(trajectory_time(1)).possible_q.true_q,2)';
x2 = x;
y2 = y;
x3 = posterior(trajectory_time(1)).mux(1);
y3 = posterior(trajectory_time(1)).mux(2);
x4 = option.grids(AgumentedData(trajectory_time(1)).possible_q.measuredposition,1)';
y4 = option.grids(AgumentedData(trajectory_time(1)).possible_q.measuredposition,2)';
error_ellipse('mu',posterior(trajectory_time(1)).mux,...
    'C',posterior(trajectory_time(1)).Sigmax)
for t=trajectory_time(2:end)
    xtmp = option.grids(AgumentedData(t).possible_q.true_q,1)';
    ytmp = option.grids(AgumentedData(t).possible_q.true_q,2)';
    x2 = [x2;xtmp];
    y2 = [y2;ytmp];
    muxtmp = posterior(t).mux(1);
    muytmp = posterior(t).mux(2);
    if ((muxtmp-xtmp)>x_widthDIV2)
        muxtmp = muxtmp - x_width;
    else if ((muxtmp-xtmp)<-x_widthDIV2)
            muxtmp = muxtmp + x_width;
        end
    end
    if ((muytmp-ytmp)>y_widthDIV2)
        muytmp = muytmp - y_width;
    else if ((muytmp-ytmp)<-y_widthDIV2)
            muytmp = muytmp + y_width;
        end
    end    
    x3 = [x3;muxtmp];
    y3 = [y3;muytmp];

    
    noisyxtmp = option.grids(AgumentedData(t).possible_q.measuredposition,1);
    noisyytmp = option.grids(AgumentedData(t).possible_q.measuredposition,2);
    if ((noisyxtmp-xtmp)>x_widthDIV2)
        noisyxtmp = noisyxtmp - x_width;
    else if ((noisyxtmp-xtmp)<-x_widthDIV2)
            noisyxtmp = noisyxtmp + x_width;
        end
    end
    if ((noisyytmp-ytmp)>y_widthDIV2)
        noisyytmp = noisyytmp - y_width;
    else if ((noisyytmp-ytmp)<-y_widthDIV2)
            noisyytmp = noisyytmp + y_width;
        end
    end    
    x4 = [x4;noisyxtmp];
    y4 = [y4;noisyytmp];    
    
    error_ellipse('mu',[muxtmp,muytmp],...
    'C',posterior(t).Sigmax)
    for n= 1:option.agentnumbers
        if (xtmp(n) - x(end,n))> x_widthDIV2
             xdiff(n) = - x_width;
        else if (xtmp(n) - x(end,n))< -x_widthDIV2
                xdiff(n) = x_width;
            else
                xdiff(n) = 0;
            end
        end
        if (ytmp(n) - y(end,n))> y_widthDIV2
             ydiff(n) = - y_width;
        else if (ytmp(n) - y(end,n))< -y_widthDIV2
                ydiff(n) = y_width;
            else
                ydiff(n) = 0;
            end
        end            
    end
    if ((norm(xdiff)==0)&&(norm(ydiff)==0))
        x = [x;xtmp];
        y = [y;ytmp];
    else
        x = [x;xtmp+xdiff;ones(1,option.agentnumbers)*inf;x(end,:)-xdiff;xtmp];
        y = [y;ytmp+ydiff;ones(1,option.agentnumbers)*inf;y(end,:)-ydiff;ytmp];
        x3= [x3(1:end-1);ones(1,option.agentnumbers)*inf;muxtmp];
        y3= [y3(1:end-1);ones(1,option.agentnumbers)*inf;muytmp];
    end
end
plot(x,y,'r','marker','.')

plot(x3,y3,'b','marker','.')
plot(x4,y4,'r','marker','*','linestyle','none')


xlim([option.X_mesh(1), option.X_mesh(end)]);
ylim([option.Y_mesh(1), option.Y_mesh(end)]);
colorbar; title('GMRF');    
colormap('gray');





nt = size(true_field,3);
nx = size(true_field,2); % number of X grid
ny = size(true_field,1); % number of Y grid

for t=101:100
fig_field = figure('name',['Compare estimated field Case 1,2 and 3      time = ', num2str(t)]);
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
estimated_field = reshape(posteriorTrueQ(t).muz,fieldsize);
max_val = max(max(max(estimated_field)),max_val);
min_val = min(min(min(estimated_field)),min_val);
plt_trueposition_estimatedfield= ...
    imagesc(option.X_mesh,option.Y_mesh,estimated_field);
axis on; title('(a)');

        %----plot estimation variance--------------------------------------
subplot(3,3,4)
estimated_variance = ...
    reshape(posteriorTrueQ(t).varz,fieldsize);
max_var = max(max(max(estimated_variance)),max_var);
plt_trueposition_variancefield = ...
    imagesc(option.X_mesh,option.Y_mesh,estimated_variance);
axis on; title('(d)');
MeanSquareVariance(t,1) = mean(mean( estimated_variance));
        %----plot estimation error-----------------------------------------
subplot(3,3,7)
estimation_error = (true_field(:,:,t) - estimated_field).^2;
plt_trueposition_errorfield = ...
    imagesc(option.X_mesh,option.Y_mesh,estimation_error);
max_err = max(max(max(estimation_error)),max_err);
axis on; title('(g)');
hold on
q  = AgumentedData(t).possible_q.true_q;
plot(option.grids(q,1),option.grids(q,2),....
    'm','marker','o','linestyle','none');

MeanSquareError(t,1) = mean(mean( estimation_error));

%% Noisy Position
        %----plot estimated field & support of positions-------------------
subplot(3,3,2)
estimated_field = reshape(posteriorNoisyQ(t).muz,fieldsize);
max_val = max(max(max(estimated_field)),max_val);
min_val = min(min(min(estimated_field)),min_val);
plt_noisyposition_estimatedfield= ...
    imagesc(option.X_mesh,option.Y_mesh,estimated_field);
axis on; 
title('(b)');

        %----plot estimation variance--------------------------------------
subplot(3,3,5)
estimated_variance = ...
    reshape(posteriorNoisyQ(t).varz,fieldsize);
max_var = max(max(max(estimated_variance)),max_var);
plt_noisyposition_variancefield = ...
    imagesc(option.X_mesh,option.Y_mesh,estimated_variance);
axis on; title('(e)');
MeanSquareVariance(t,2) = mean(mean( estimated_variance));
        %----plot estimation error-----------------------------------------
subplot(3,3,8)
estimation_error = (true_field(:,:,t) - estimated_field).^2;
max_err = max(max(max(estimation_error)),max_err);
plt_noisyposition_errorfield = ...
    imagesc(option.X_mesh,option.Y_mesh,estimation_error);
axis on; title('(h)');
hold on
for i=1:option.agentnumbers

        q = AgumentedData(t).possible_q.measuredposition(i);
        plot(option.grids(q,1),option.grids(q,2),...
            'm','marker','*');


end
MeanSquareError(t,2) = mean(mean( estimation_error));

%% Uncertain Position
        %----plot estimated field & support of positions-------------------
subplot(3,3,3)
estimated_field = reshape(posterior(t).muz,fieldsize);
max_val = max(max(max(estimated_field)),max_val);
min_val = min(min(min(estimated_field)),min_val);
plt_uncertainposition_estimatedfield= ...
    imagesc(option.X_mesh,option.Y_mesh,estimated_field);
axis on;  caxis([-1 1]); colorbar; title('(c)'); 

        %----plot estimation variance--------------------------------------
subplot(3,3,6)
estimated_variance = reshape(posterior(t).varz,fieldsize);
max_var = max(max(max(estimated_variance)),max_var);
plt_uncertainposition_variancefield = ...
    imagesc(option.X_mesh,option.Y_mesh,estimated_variance);
colorbar; title('(f)'); axis on;
caxis([0 0.4]);
MeanSquareVariance(t,3) = mean(mean( estimated_variance));
        %----plot estimation error-----------------------------------------
subplot(3,3,9)
estimation_error = (true_field(:,:,t) - estimated_field).^2;
max_err = max(max(max(estimation_error)),max_err);
plt_uncertainposition_errorfield = ...
    imagesc(option.X_mesh,option.Y_mesh,estimation_error);
hold on
caxis([0 0.4]); colorbar; title('(i)'); axis on;
for i=1:option.agentnumbers
        q = AgumentedData(t).possible_q.support_qt{1,i};
        plot(option.grids(q,1),option.grids(q,2),'m','marker','.','linestyle','none');
end
MeanSquareError(t,3) = mean(mean( estimation_error));

%% Edit plots
est_colorbar_step = 0.001;
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
% saveas(fig_field,[ num2str(t) 'frame.jpg']);
end



mean(MeanSquareError)
mean(MeanSquareVariance)

% for n = 1: size(realY,1)
%     plot(Xbar(n*2-1),Xbar(n*2),'m','marker','o','markersize',10,'MarkerFaceColor',[1 1 1]);
%     text(Xbar(n*2-1),Xbar(n*2),num2str(n),'fontsize'...
%         ,7,'horizontalalignment',    'center','fontweight','bold');
% end
