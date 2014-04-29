close all
clc
clear
nf=16;
%nt=558;
[~, ~, raw] = xlsread('C:\Users\dohuan\Documents\GitHub\GP-SLAM\MATLAB\FeaturesAndLocations.xlsx','Sheet1','A2:AL464');
data = reshape([raw{:}],size(raw));
jpgfile = data(:,6);
nt=size(data,1);
% s=sprintf('./image_p/Panoramic_1.jpg');
% img=imread(s);
% img_gr=imresize(rgb2gray(img),[128,128]);
% F=fft2(img_gr);
% FFT=reshape(reshape(abs(F(1,2:1+nf)),1,[]),[],1);
%C=reshape(fftshift(F),[],1);
feature=zeros(nt,1);
for i=1:nt
    s=sprintf('./image_p_angleFixed/Panoramic_%d.jpg',jpgfile(i));
    img=imread(s);
    img_gr=imresize(rgb2gray(img),[128,128]);
    %B=imhist(img_gr);
    F=fft2(img_gr);
    FFT=reshape(reshape(abs(F(1,2:1+nf)),1,[]),[],1);
    C=abs(F(1,:));
    C=C';
    fit=polyfit((1:128)',log(C),1);
    feature(i)=fit(1);
    figure(1);
    movegui(figure(1),'southwest');
    subplot(1,2,1);
    imshow(img_gr);
    subplot(1,2,2);
    hold on;
    plot(i,fit(1),'o','MarkerFaceColor','red');
    hold off;
    %feature(i)=mean(FFT)/var(FFT);
    %plot(log(abs(F(1,2:1+nf))));
    figure(2);
    movegui(figure(2),'northwest');
    %surf(log(abs(fftshift(F))));
    %contour(log(abs(F)));
    %plot(C);
    angleCurr=mean(angle(F(1,2:1+nf)));
    if i==1
       oldAngle=angleCurr; 
    end
    angleRate=angleCurr-oldAngle;
    hold on;
    plot(i,angleRate,'*r');
    hold off;
    oldAngle=angleCurr;
    M(i)=getframe;
    %pause(0.1);
end
%movie(M,1,5);
%close all
%plot(feature)