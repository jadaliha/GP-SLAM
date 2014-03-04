function experimentaldeta()

clc
clear
filename = 'EX1';
open([filename '.jpg']);
raw_data = eval(filename);
raw_data2 = ((0.2*raw_data(:,:,1)+0.5*raw_data(:,:,2)+0.3*raw_data(:,:,3)));

imagesc(raw_data2)

disp([min(min(raw_data2)) max(max(raw_data2))]);
nx1 = size(raw_data2,2);
ny1 = size(raw_data2,1);

%-----reduce resulotion----------------------------------------------------
target_resolution = [61, 61];
data = zeros(target_resolution(2), target_resolution(1));
x_step = floor(nx1/target_resolution(1));
y_step = floor(ny1/target_resolution(2));

for indy = (1:y_step:ny1-y_step)
    for indx = (1:x_step:nx1-x_step)
        data((indy-1)/y_step+1,(indx-1)/x_step+1) = ...
            mean(mean(raw_data2(indy:indy+y_step-1,indx:indx+x_step-1)));
    end
end
true_fieldt = data(1:target_resolution(2), 1:target_resolution(1))/max(max(data));
imagesc(data)

for n = 1:10
    true_field(:,:,n) = true_fieldt;
end

save('exdata3.mat','true_field')