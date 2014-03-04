function [ F ] = logmvnpdf2( x,mu,Sigma,option )
Sigma = 0.5*(Sigma + Sigma');
x = x- mu;
mu = 0*mu;
x_width = option.X_mesh(end) - option.X_mesh(1) + option.finegridsize;
y_width = option.Y_mesh(end) - option.Y_mesh(1) + option.finegridsize;
tmpx = kron(ones(option.agentnumbers,1),[x_width;y_width]);
tmpx = [x - tmpx, x, x + tmpx];
for jind=1:option.agentnumbers
    [~,index] =        min(abs(tmpx(jind*2-1,:)));
    x(jind*2-1) = tmpx(jind*2-1,index);
    [~,index] =        min(abs(tmpx(jind*2,:)));
    x(jind*2) = tmpx(jind*2,index);
end
F = log(mvnpdf(x,mu,Sigma));
end

