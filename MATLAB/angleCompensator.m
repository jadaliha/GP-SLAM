close all
clear
clc
[~, ~, raw] = xlsread('C:\Users\dohuan\Documents\GitHub\GP-SLAM\MATLAB\FeaturesAndLocations.xlsx','Sheet1','A2:AL464');
data = reshape([raw{:}],size(raw));
jpgfile = data(:,6);
time = data(:,5);
turnrate = data(:,2);
nt=size(data,1);
currAng=0; %in radian
for i=1:nt
    if (i==1)
        currAng=currAng+turnrate(i)*(time(i)); 
    else
    currAng=currAng+turnrate(i)*(time(i)-time(i-1)); 
    end
    %angTrack(i)=currAng;
    s=sprintf('./image_d/Doughnut_%d.jpg',jpgfile(i));
    img=imread(s);
    unwrap_img=imageUnwrapper(img,currAng);
    s=sprintf('./image_p_angleFixed/Panoramic_%d.jpg',jpgfile(i));
    imwrite(unwrap_img,s);
    %imshow(imresize(unwrap_img,[128 128]));
end
