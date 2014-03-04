function [ Q ] = prec_mat(n1,n2)
%GENERATEQ Generate Q Matrix
%   m = rows of original matrix
%   n = columns of original matrix


n = n1;
m = n2;

% Q = zeros(m*n);

size = m*n;
%Q = spalloc(size, size, size*13);
%Q = sparse(size, size);
I = zeros(size*13, 1);
J = zeros(size*13, 1);
K = zeros(size*13, 1);
curr = 1;

wait = waitbar(0,'Processing Q Matrix...');

for x = 1:n
    waitbar(x/n,wait, ['Processing Q Matrix (' num2str(x) ' of ' num2str(n) ' columns)']);
    for y = 1:m
        index = (x-1)*m + y;
        
        %waitbar(index/size,wait, ['Processing Q Matrix (' num2str(index) ' of ' num2str(size)]);
        
        % Start the Process
        
        % Centers (No Missing Items)
        if (x >= 3) && (y >= 3) && (x <= n-2) && (y <= m-2)
            %Q = Center(Q, m, index);
            [ I, J, K, curr ] = Center( I, J, K, curr, m, index );
        else
            HCut = 0;
            VCut = 0;
            if (x <= 1)
                HCut = -2;
            elseif (x == 2)
                HCut = -1;
            elseif (x >= n)
                HCut = 2;
            elseif (x == n-1)
                HCut = 1;
            end
            if (y <= 1)
                VCut = 2;
            elseif (y == 2)
                VCut = 1;
            elseif (y >= m)
                VCut = -2;
            elseif (y == m-1)
                VCut = -1;
            end
            %Q = Cut(Q, m, index, HCut, VCut);
            [ I, J, K, curr ] = Cut(I, J, K, curr, m, index, HCut, VCut);
        end
        
    end
end

waitbar(1,wait,'Generating Sparse Structure...');

Q = sparse(I(1:curr-1), J(1:curr-1), K(1:curr-1), size, size);

close(wait);
pause(0.2);

end

% COMPLETE!
function [ I, J, K, curr ] = Center( I, J, K, curr, m, index )

%Q(index,index) = 20;
I(curr) = index; J(curr) = index; K(curr) = 20; curr = curr+1;
negeights = [-1 1 -m m];
%Q(index,index+negeights) = -8;
for pos = 1:4
    I(curr) = index; J(curr) = index+negeights(pos); K(curr) = -8; curr = curr+1;
end
twos = [(-m-1) (-m+1) (m+1) (m-1)];
%Q(index,index+twos) = 2;
for pos = 1:4
    I(curr) = index; J(curr) = index+twos(pos); K(curr) = 2; curr = curr+1;
end
ones = 2*negeights;
%Q(index,index+ones) = 1;
for pos = 1:4
    I(curr) = index; J(curr) = index+ones(pos); K(curr) = 1; curr = curr+1;
end

end

% COMPLETE!
function [ I, J, K, curr ] = Cut( I, J, K, curr, m, index, HCut, VCut)
% Cut is positive if top or right is cut off
% Cut is negative if bottom or left is cut off

pos = [0 0 (index-2) 0 0;
    0 ((index-m)-1) (index-1) ((index+m)-1) 0;
    (index-(2*m)) (index-m) index (index+m) (index+(2*m));
    0 ((index-m)+1) (index+1) ((index+m)+1) 0;
    0 0 (index+2) 0 0;];

M = zeros(5,5); %#ok<NASGU>

if abs(HCut) == 2 && abs(VCut) == 2
    M = [0 0 1 0 0; 0 0 -4 2 0; 0 0 4 -4 1; 0 0 0 0 0; 0 0 0 0 0]; %Case 5
elseif abs(HCut) == 1 && abs(VCut) == 1
    M = [0 0 1 0 0; 0 2 -8 2 0; 0 -6 18 -8 1; 0 2 -6 2 0; 0 0 0 0 0]; %Case 4
else
    if HCut == 0 || VCut == 0
        cut = max(abs(HCut), abs(VCut));
        if cut == 1
            M = [0 0 1 0 0; 0 2 -8 2 0; 1 -8 19 -8 1; 0 2 -6 2 0; 0 0 0 0 0]; %Case 2
        else
            M = [0 0 1 0 0; 0 2 -6 2 0; 1 -6 11 -6 1; 0 0 0 0 0; 0 0 0 0 0]; %Case 3
        end
    else
        M = [0 0 1 0 0; 0 2 -6 2 0; 0 -4 10 -6 1; 0 0 0 0 0; 0 0 0 0 0]; %Case 1
    end
    if abs(HCut) > abs(VCut)
        M = rot90(M, -1);
        M = flipud(M);
    end
end

if HCut > 0
    M = fliplr(M);
end
if VCut > 0
    M = flipud(M);
end

for i = 1:5
    for j = 1:5
        if M(i,j) ~= 0
            %Q(index,pos(i,j)) = M(i,j);
            I(curr) = index; J(curr) = pos(i, j); K(curr) = M(i, j); curr = curr+1;
        end
    end
end
            
end

