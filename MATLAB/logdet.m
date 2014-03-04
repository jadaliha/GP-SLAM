function y = logdet(A)

L = chol(A,'lower');
y = 2*sum(log(diag(L)));