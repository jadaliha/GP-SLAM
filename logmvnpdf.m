function [ F ] = logmvnpdf( x,mu,Sigma )
F = log(mvnpdf(x,mu,Sigma));
% log_pi_y(k) = -0.5*logdet(Sigmay)- 0.5*(y-muy)'*InvSigmay*(y-muy);
end

