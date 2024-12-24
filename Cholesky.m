function L=Cholesky(A)
m=length(A);
L=zeros(m,m);
% 
for i=1:m
    L(i,i)=sqrt(A(i,i)-sum(L(i,1:i-1).^2));
    for j=i+1:m
            L(j,i)=(A(i,j)-sum(L(i,1:i-1).*L(j,1:i-1)))/L(i,i);
    end
end
