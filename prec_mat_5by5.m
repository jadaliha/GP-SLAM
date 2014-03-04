function Q = prec_mat_5by5(kappa,alpha,n1,n2)

% a = 4 + alpha;
% Q = sparse(n1*n2,n1*n2);
% 
% for i = 1:n1
%     for j = 1:n2
%         k = (i-1)*n2+j;
%         
%         Q(k,k) = 4 + a^2;
%         
%         Q(k,mod(i,n1)*n2+j) = -2*a;
%         Q(k,mod(i-2,n1)*n2+j) = -2*a;
%         Q(k,(i-1)*n2+mod(j,n2)+1) = -2*a;
%         Q(k,(i-1)*n2+mod(j-2,n2)+1) = -2*a;
%         
%         Q(k,mod(i,n1)*n2+mod(j,n2)+1) = 2;
%         Q(k,mod(i-2,n1)*n2+mod(j,n2)+1) = 2;
%         Q(k,mod(i,n1)*n2+mod(j-2,n2)+1) = 2;
%         Q(k,mod(i-2,n1)*n2+mod(j-2,n2)+1) = 2;
%         
%         Q(k,mod(i+1,n1)*n2+j) = 1;
%         Q(k,mod(i-3,n1)*n2+j) = 1;
%         Q(k,(i-1)*n2+mod(j+1,n2)+1) = 1;
%         Q(k,(i-1)*n2+mod(j-3,n2)+1) = 1;
% 
%         
%     end
% end

Q = prec_mat_3by3(1,alpha,n1,n2);
Q = Q^2;
Q = kappa*Q;

