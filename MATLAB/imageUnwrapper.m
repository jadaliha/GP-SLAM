function output=imageUnwrapper(B,ang) %ang in radian

VD = struct('globe_radious',[40 79],...
    'currentimage',B,...
    'globe_center',[80 80]);
output=zeros(40,503,3,class(VD.currentimage)); %panoramic image
output_size=size(output);
nf=16; % histogram with 16 bins

for index1 = 1:output_size(1)
    for index2 = 1:output_size(2)
        tmp_th = 2*pi*index2/output_size(2)-ang;
        tmp_r = VD.globe_radious(1) + index1/output_size(1) * ...
            (VD.globe_radious(2)-VD.globe_radious(1));
        x1 = (tmp_r * sin(tmp_th)+VD.globe_center(1));
        x2 = (tmp_r * cos(tmp_th)+VD.globe_center(2));
        
        output(index1,index2,:) = ...
            VD.currentimage(round(x1),round(x2),:);
    end
end 
% gray_image = imresize(rgb2gray(output),[128,128]);
% F=fft2(gray_image);
% FFT=reshape(reshape(abs(F(1,2:1+nf)),1,[]),[],1);
% HIST=imhist(gray_image,nf);
% NFFT=zeros(sizeFFT(1),sizeFFT(2));
end
% 
% 
% %% Absolute FFT features
% 
% Gray_image = imresize(rgb2gray(RGB_image),[128,128]);
% F = fft2(Gray_image);
% y = reshape(reshape(abs(F(1,2:1+nf)),1,[]),[],1);
% 
% 
% %% Absolute HISTOGRAM features
% 
% Gray_image = rgb2gray(RGB_image);
% [y,~] = imhist(Gray_image,nf);

