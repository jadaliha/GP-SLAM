function [field, VehicleStates, AgumentedData, numberofpossibilities] = Generate_new_field(savetofile)
% Generate_new_field creates a spatio-temporal field 
% Authour : Mahdi Jadaliha
% e-mail  : jadaliha@gmail.com
% =========================================================================
global globaloption sys_parameter                                           
if isempty(globaloption);
    [option , sys_parameter ] = Configuration();                            % initialize system configuration
else
    option = globaloption;                                                  % use global parameters if exist 
end


nx = size(option.X_mesh,2) ;                                                % number of X grid
ny = size(option.Y_mesh,2) ;                                                % number of Y grid
gridsize = size(option.grids,1);                                            % number of spacial sites 
                                                                            % e.g. for regular latice number of $gridsize = nx * ny$
nt = size(option.T_mesh,2) ;                                                % number of time grid

centerarea = option.grids;                                                  % In this case sampling postion include all the interested area
selectionsize = size(centerarea,1);

possible_radious = sqrt(option.vehicle.ObservationNoise);                                                     % determine sampling position uncertainity support

field = zeros(ny,nx,nt);                                                    % initializing memory to keep true field
eta = zeros(nx*ny);                                                         % field in time t

rand_num = randn(nx*ny,nt);                                                 % a unique random field generaator variable
cs = rand(option.agentnumbers,nt);                                          % generate random variavle to choose uncertain samplin positions
rand_u = rand(option.agentnumbers * 2, nt);                                 % control input
measurement_noise = randn(option.agentnumbers,nt) * sqrt(option.s_e2);      % measurement noise
wt = sqrt(option.vehicle.ModelUncertanity) * randn(2 * option.agentnumbers , nt);

if ~isfield(option, 'StartPositions')
    t = 1;
    for n=1:option.agentnumbers
                                                                            % Assume the heading can be measured accuratly trough the compass
        q(n)  = ceil(gridsize * rand(1,1)+1);                               % q0 here is the initial true position
        h(n) = rand(1,1) * 2 * pi;
        x(n) = option.grids(q(n),1);
        y(n) = option.grids(q(n),2);
    end
else
    q = option.StartPositions;
end
%---- Construct precision matrix Qz
Qz = prec_mat_5by5(option.kappa,option.alpha,nx,ny); 
L = chol(Qz,'lower');


%% Plotting 
if (sys_parameter.plot_enabled)
    fig_field = figure;
    set(fig_field,'position',[10 50 400 300]);
    hold on
    plt_truefield = imagesc(option.X_mesh,option.Y_mesh,field(:,:,1));
    for n= 1:option.agentnumbers
        plt_robots_pos(n) = plot(x(n),...
            y(n),'marker','o');
        tmp_c = 0.05;
        plt_robots_h(n) = plot(x(n) +...
            tmp_c * cos(h(n)),...
            y(n)+ tmp_c * sin(h(n)) ,'marker','.');
        plt_robots_qsupport(n) = plot(x(n) ,...
            y(n),'m','marker','.','linestyle','none');
    end
    xlim([option.X_mesh(1), option.X_mesh(end)]);
    ylim([option.Y_mesh(1), option.Y_mesh(end)]);
    colorbar; title('GMRF');    
    colormap('gray');
end

x_width = option.X_mesh(end) - option.X_mesh(1) + option.finegridsize;
y_width = option.Y_mesh(end) - option.Y_mesh(1) + option.finegridsize;
%% main loop
for t=1:nt
    if (t==1)
        eta = L'\rand_num(:,1);                                             % generate initial value of field
    else
        eta = eta + option.c_randomwallk * (L'\rand_num(:,t));                % time varyation of the field
    end
    field(:,:,t)= reshape(eta,ny,nx);                                       % reshape vector to 2D space
    
    if sys_parameter.torus                                                  % Tures case
        for n= 1:option.agentnumbers
 
            %------------------
            VehicleStates(t).x(n) = x(n);
            VehicleStates(t).y(n) = y(n);  
            VehicleStates(t).h(n) = h(n);      
            r2 = min([(option.grids(:,1) - x(n)).^2,...
                (option.grids(:,1) - x(n) + x_width).^2,...
                (option.grids(:,1) - x(n) - x_width).^2],[],2)+  ...
                min([(option.grids(:,2) - y(n)).^2,...
                (option.grids(:,2) - y(n) + y_width).^2,...
                (option.grids(:,2) - y(n) - y_width).^2],[],2);
            [~,VehicleStates(t).q(n)] = min(r2);
            
%----- psossible position using Gaussian distribuation---------------------
            Ix = find(r2<possible_radious^2);
            tmp3 = exp(-0.5 * r2(Ix)/option.vehicle.ObservationNoise);
            tmp3 = tmp3/(sum(tmp3));
            tmp4 = cumsum(tmp3);
            tmp5 = find(tmp4>cs(n,t));
             
            possible_qt.measuredposition(n)  = Ix(tmp5(1));
            possible_qt.true_q(n)            = VehicleStates(t).q(n);
            
            noisy_x = centerarea(possible_qt.measuredposition(n),1);
            noisy_y = centerarea(possible_qt.measuredposition(n),2);
            r2 = min([(option.grids(:,1) - noisy_x).^2,...
                (option.grids(:,1) - noisy_x + x_width).^2,...
                (option.grids(:,1) - noisy_x - x_width).^2],[],2)+  ...
                min([(option.grids(:,2) - noisy_y).^2,...
                (option.grids(:,2) - noisy_y + y_width).^2,...
                (option.grids(:,2) - noisy_y - y_width).^2],[],2);  
            
            Ix = find(r2<possible_radious^2);
            tmp3 = exp(- 0.5 * r2(Ix)/option.vehicle.ObservationNoise);
            tmp3 = tmp3/(sum(tmp3));    
            possible_qt.support_qt(n)        = {Ix};
            possible_qt.prior_qt(n)          = {tmp3};
            possible_qt.N_possible_qt(n)     = size(tmp3,1);         
            if (sys_parameter.plot_enabled)
                set(plt_robots_pos(n),'XData',x(n),...
                    'YData',y(n));
                set(plt_robots_h(n),'XData',x(n) +...
                    tmp_c * cos(h(n)),...
                    'YData', y(n) +...
                    tmp_c * sin(h(n)));
                support_qt = possible_qt.support_qt(n);
                set(plt_robots_qsupport(n),...
                    'XData',centerarea(support_qt{1},1),...
                    'YData',centerarea(support_qt{1},2));
            end          
            
            u = rand_u(2*n-1:2*n,t);                            
            v = 2.5*sum(u);
            VehicleStates(t).u(n,1) = v;
            w = diff(u);
            x(n) = circulate (x(n) + v * cos(h(n) + wt(2*n-1,t)),...
                [option.X_mesh(1), option.X_mesh(end)]);
            y(n) = circulate (y(n) + v * sin(h(n) + wt(2*n  ,t)),...
                [option.Y_mesh(1), option.Y_mesh(end)]);
            h(n) = circulate (h(n) + w ,[0, 2*pi]);               
        end
    else
    end


                      
    true_q = possible_qt.true_q ;

    yt = eta(true_q)  + measurement_noise(:,t);                              % measurements $y(t)$    
    AgumentedData(t).possible_q = possible_qt;
    AgumentedData(t).y = yt;
    numberofpossibilities(:,t) = prod(possible_qt.N_possible_qt);           % number of possibilities for $q(t)$    
    
    %----plot field----
    if (sys_parameter.plot_enabled)
        set(plt_truefield,'CData',field(:,:,t));  
        refreshdata(fig_field);
        pause(1);
    end

    
end

if (nargin == 1)
    save(savetofile, 'option', 'field');
end
end

function vout = circulate(vin,limits)
a = limits(2) - limits(1);
while ((vin>limits(2)) || (vin<limits(1)))
    if (vin> limits(2))
        vin = vin - a;
    else if (vin < limits(1))
        vin = vin + a;
        end
    end
end
vout = vin;
end
