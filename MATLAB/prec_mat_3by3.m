function Q = prec_mat_3by3(kappa,alpha,n1,n2)

a = 4 + alpha;
Q = sparse(n1*n2,n1*n2);

global sys_parameter

if sys_parameter.torus
    for i = 1:n1
        for j = 1:n2
            k = (i-1)*n2+j;

            Q(k,k) = a;

            Q(k,mod(i,n1)*n2+j) = -1;
            Q(k,mod(i-2,n1)*n2+j) = -1;
            Q(k,(i-1)*n2+mod(j,n2)+1) = -1;
            Q(k,(i-1)*n2+mod(j-2,n2)+1) = -1;


        end
    end
else
    for i = 1:n1
        for j = 1:n2
            k = (i-1)*n2+j;

            Q(k,k) = a;
            if (mod(i,n1)==i)
                Q(k,i*n2+j) = -1;
            else
                Q(k,k) = Q(k,k) - 1;
            end
            if (mod(i-2,n1)==i-2)
                Q(k,(i-2)*n2+j) = -1;
            else
                Q(k,k) = Q(k,k) - 1;                
            end
            if (mod(j,n2)==j)
                Q(k,(i-1)*n2+mod(j,n2)+1) = -1;
            else
                Q(k,k) = Q(k,k) - 1;                
            end
            if (mod(j-2,n2)==j-2)
                Q(k,(i-1)*n2+(j-2)+1) = -1;
            else
                Q(k,k) = Q(k,k) - 1;                
            end

        end
    end    
end


