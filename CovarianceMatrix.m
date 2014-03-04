function K = CovarianceMatrix(xstar,X,theta)
% updated code for optimized matrix calculation 2/21/2013
if (nargin==3)
    sigma2f = theta(1);
    C1_inv = -0.5 * diag(theta(2:1 + size(xstar,2))).^(-2);
    diagonalized_C1 = 1;
else
    global system
    sigma2f=system.sigma2f;
    C1_inv = -0.5*(system.Sigma2X)^(-1);
    diagonalized_C1 = ...
        (sum(abs(diag(C1_inv) - diag(-0.5*(system.Sigma2X)^(-1))))<0.0001); % check if C1_inv is diagonal
end
m = size(xstar,1);

if ~diagonalized_C1           
    if(nargin>1)
        n = size(X,1);
        K = zeros(m,n);

        for j=1:n
            for i=1:m
                K(i,j) = exp((xstar(i,:)-X(j,:))*C1_inv*(xstar(i,:)-X(j,:))');

            end
        end
    else
        K = eye(m,m);
        for j=1:m-1
            for i=j+1:m
                tmp = exp((xstar(i,:)-xstar(j,:))*C1_inv*(xstar(i,:)-xstar(j,:))');
                K(i,j) = tmp;
                K(j,i) = tmp;
            end
        end
    end    
else
    % optimized matrix calculation
    if(nargin==1)
        X = xstar;
    end
    n = size(X,1);
    K = exp((xstar(:,1) * ones(1,n) - ones(m,1) * X(:,1)').^2*C1_inv(1,1) + ...
        (xstar(:,2) * ones(1,n) - ones(m,1) * X(:,2)').^2*C1_inv(2,2));
end

K = sigma2f* K;