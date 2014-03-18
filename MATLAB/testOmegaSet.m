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

point1=[-1 2];
point2=[-6 9];
dist1=sqrt((grids(:,1) - point1(1)).^2 + ...
        (grids(:,2) - point1(2)).^2);
dist2=sqrt((grids(:,1) - point2(1)).^2 + ...
        (grids(:,2) - point2(2)).^2);
[IC_1,IX_1] = sort(dist1);
[IC_2,IX_2] = sort(dist2);
%pointIX.IX=sort([IX_1(1);IX_2(1)]);%index of two points on grids
% gridIX=1:1:size(grids,1);
% gridIX=reshape(gridIX,sizeY,[]);
% gridIX=flipud(gridIX);
% 
% pointIX.point=cell(2,1);
% [pointIX.point{1}(2), pointIX.point{1}(1)]=find(gridIX==pointIX.IX(1)); % 1 is x-axis, 2 is y-axis
% [pointIX.point{2}(2), pointIX.point{2}(1)]=find(gridIX==pointIX.IX(2));
% %box.point=zeros(2,1) box.height box.width
% size_x=abs(pointIX.point{2}(1)-pointIX.point{1}(1));
% size_y=abs(pointIX.point{2}(2)-pointIX.point{1}(2));
% 
% if size_x==0
%     if (size_y==0)
%         box.point(1)=pointIX.point{1}(1)-1;
%         box.point(2)=pointIX.point{1}(2)-1;
%         box.height=2;
%         box.width=2;
%     else
%         box.point(1)=round(min(pointIX.point{1}(1),pointIX.point{2}(1))-size_y/2);
%         box.point(2)=min(pointIX.point{1}(2),pointIX.point{2}(2));
%         box.height=size_y;
%         box.width=size_y; 
%     end
% else if size_y==0
%         box.point(1)=min(pointIX.point{1}(1),pointIX.point{2}(1));
%         box.point(2)=round(min(pointIX.point{1}(2),pointIX.point{2}(2))-size_x/2);
%         box.height=size_x;
%         box.width=size_x;
%     else
%         box.point(1)=min(pointIX.point{1}(1),pointIX.point{2}(1));
%         box.point(2)=min(pointIX.point{1}(2),pointIX.point{2}(2));
%         box.height=size_y;
%         box.width=size_x;
%     end
% end
% Omg=gridIX(box.point(2):box.point(2)+box.height,...
%            box.point(1):box.point(1)+box.width);
% Omg=reshape(Omg,[],1);
Omg=omegaMaker(grids,IX_1,IX_2,sizeX,sizeY);
OmgPlot=grids(Omg,:);
hold on;
plot(tmp_S1,tmp_S2,'bo');
plot(OmgPlot(:,1),OmgPlot(:,2),'r*');
plot([point1(1);point2(1)],[point1(2);point2(2)],'b*');
hold off;