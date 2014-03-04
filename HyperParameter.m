function [p fval] = HyperParameter(f,x,p0)

if (nargin < 3)
    global system
    theta0 = [var(f) 1 1 system.sigma2w]'; % sig_f^2, sig_x sig_y
else 
    theta0 = p0;
end

optimset('algorithm','trust-region-reflective','gradobj','on');
exitflag = 0;
counter = 0;
while ~exitflag
    [p,fval,exitflag] = fminsearch(@(theta) -likelihood(x,f,theta),... %changed "comp_L"-->"likelihood" in 4/24/2013
        theta0);
    counter = counter + 1;
    exitflag = exitflag || (counter==3);
end
p = abs(p);
 
end
function L = likelihood(X,f,theta)
    theta = abs(theta);
    n = length(f);
    K = CovarianceMatrix(X,X,theta(1:end-1));
    sigma2w = abs(theta(end));
    C = K + sigma2w * eye(n);
    v = C\f;
    L = -1/2*f'*v - 1/2*logdet(C) - n/2*log(2*pi);
end
