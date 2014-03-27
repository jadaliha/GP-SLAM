function Y = downSampling(X,ratio)
[m n]=size(X);
if mod(m,ratio)==0
    Y=zeros(m/ratio,n);
    for i=1:size(Y,1)
        Y(i,:)=X(i*ratio,:);
    end
else
   error('myApp:argChk', 'ratio does not match')
end
end