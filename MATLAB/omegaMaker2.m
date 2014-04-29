function [omega,true_point]=omegaMaker2(grids,IX_true,IX_noise,noise_threshold)
    noise_threshold=noise_threshold/3.5; %compensate for scaling x and y axis
    sizeGrid=size(grids,1);
    true_point=grids(IX_true(1),:);
    noise_point=grids(IX_noise(1),:);
    
    dist=sqrt((grids(:,1)-ones(sizeGrid,1)*noise_point(1)).^2+...
        (grids(:,2)-ones(sizeGrid,1)*noise_point(2)).^2);
    [IC,IX]=sort(dist);
    thres_IX= IC<=noise_threshold;
    omega=IX(thres_IX);
    %omega=reshape(omega,[],1);
end