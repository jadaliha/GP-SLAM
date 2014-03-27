clear
close all
clc
Xmax=10;
Ymax=10;
X_mesh = (-Xmax:1:Xmax);
Y_mesh = (-Ymax:1:Ymax);
sizeX=size(X_mesh,2);
sizeY=size(Y_mesh,2);
[tmp_S1, tmp_S2] = meshgrid(X_mesh,Y_mesh);
grids = [reshape(tmp_S1,[],1) reshape(tmp_S2,[],1)];

point1=[-5 5];
point2=[0 0];
dist1=sqrt((grids(:,1) - point1(1)).^2 + ...
        (grids(:,2) - point1(2)).^2);
dist2=sqrt((grids(:,1) - point2(1)).^2 + ...
        (grids(:,2) - point2(2)).^2);
[IC_1,IX_1] = sort(dist1);
[IC_2,IX_2] = sort(dist2);

[Omg,rec]=omegaMaker(grids,IX_1,IX_2,sizeX,sizeY);
OmgPlot=grids(Omg,:);
hold on;
plot(tmp_S1,tmp_S2,'bo');
plot(OmgPlot(:,1),OmgPlot(:,2),'r*');
plot([point1(1) point2(1)],[point1(2) point2(2)],'y*');
plot(grids(gridIX(box.point(2),box.point(1)),1),...
    grids(gridIX(box.point(2),box.point(1)),2),'b*');
rectangle('Position',[rec.point(1) (rec.point(2)-rec.height) rec.width rec.height]);
hold off;