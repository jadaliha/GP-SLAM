
function [omega,rec,true_point]=omegaMaker(grids,IX_true,IX_noise,sizeX,sizeY)
    %persistent pointIX size_x size_y box gridIX
    true_point=grids(IX_true(1),:);
    pointIX.IX=sort([IX_true(1);IX_noise(1)]);
    gridIX=1:1:size(grids,1);
    gridIX=reshape(gridIX,sizeY,[]);
    gridIX=flipud(gridIX);

    pointIX.point=cell(2,1);
    [pointIX.point{1}(2), pointIX.point{1}(1)]=find(gridIX==pointIX.IX(1)); % 1 is x-axis, 2 is y-axis
    [pointIX.point{2}(2), pointIX.point{2}(1)]=find(gridIX==pointIX.IX(2));
    % box.point=zeros(2,1) box.height box.width
    %size_x size_y : size of the box
    size_x=abs(pointIX.point{2}(1)-pointIX.point{1}(1));
    size_y=abs(pointIX.point{2}(2)-pointIX.point{1}(2));

    if size_x==0
        if (size_y==0)
            box.point(1)=pointIX.point{1}(1)-1;
            box.point(2)=pointIX.point{1}(2)-1;
            if box.point(1)<1 
                box.point(1)=pointIX.point{1}(1);
            end
            if box.point(2)<1 
                box.point(2)=pointIX.point{1}(2);
            end
            box.height=2;
            box.width=2;
            if (box.point(1)+box.width)>sizeX
                box.width=sizeX-box.point(1);
            end
            if (box.point(2)+box.height)>sizeY
                box.height=sizeY-box.point(2);
            end
        else
            box.point(1)=round(min(pointIX.point{1}(1),pointIX.point{2}(1))-size_y/2);
            box.point(2)=min(pointIX.point{1}(2),pointIX.point{2}(2));
            if box.point(1)<1
               box.point(1)=1; 
            end
            box.height=size_y;
            box.width=size_y;
            if (box.point(1)+box.width)>sizeX
                box.width=sizeX-box.point(1);
            end
            if (box.point(2)+box.height)>sizeY
                box.height=sizeY-box.point(2);
            end
        end
    else if size_y==0
            box.point(1)=min(pointIX.point{1}(1),pointIX.point{2}(1));
            box.point(2)=round(min(pointIX.point{1}(2),pointIX.point{2}(2))-size_x/2);
            if box.point(2)<1
               box.point(2)=1; 
            end
            box.height=size_x;
            box.width=size_x;
            if (box.point(1)+box.width)>sizeX
                box.width=sizeX-box.point(1);
            end
            if (box.point(2)+box.height)>sizeY
                box.height=sizeY-box.point(2);
            end
        else
            box.point(1)=min(pointIX.point{1}(1),pointIX.point{2}(1))-1;
            box.point(2)=min(pointIX.point{1}(2),pointIX.point{2}(2))-1;
            box.height=size_y+2;
            box.width=size_x+2;
            if box.point(1)<1
               box.point(1)=1; 
            end
            if box.point(2)<1
               box.point(2)=1; 
            end
            if (box.point(1)+box.width)>sizeX
                box.width=sizeX-box.point(1);
            end
            if (box.point(2)+box.height)>sizeY
                box.height=sizeY-box.point(2);
            end
        end
    end
    assignin('base','pointIX',pointIX);
    assignin('base','gridIX',gridIX);
    assignin('base','box',box);
    assignin('base','size_x',size_x);
    assignin('base','size_y',size_y);
    omega=gridIX(box.point(2):box.point(2)+box.height,...
               box.point(1):box.point(1)+box.width);
    omega=reshape(omega,[],1);
    rec.point=grids(gridIX(box.point(2),box.point(1)),:);
    heightTmp=grids(gridIX(box.point(2)+box.height),:)-grids(gridIX(box.point(2)),:);
    widthTmp=grids(gridIX(box.point(1)+box.width),:)-grids(gridIX(box.point(1)),:);
    rec.height=abs(nonzeros(heightTmp));
    rec.width=abs(nonzeros(widthTmp));
    %clear pointIX size_x size_y box gridIX
end