function [ K ] = VectorSquare( A )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   y = kron(x,ones(size(x,1),1)).*kron(ones(size(x,1),1),x);
[ma,na] = size(A);
[ia,~] = meshgrid(1:ma);
ja = (1:na);
K = A(ia,ja).*A(ia',ja);
end
