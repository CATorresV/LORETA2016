function [C,sigv] = A_corrent(X)

[n,p] = size(X);
sigv = zeros(p,1);
for i = 1 : p
    sigv(i) = mean(pdist(X(:,i)));
end
C = zeros(n,n);
for i = 1 : n
    for ii = 1 : n
        for j = 1 : p
            C(i,ii) = C(i,ii) + exp(-(X(i,j)-X(ii,j))^2/(2*sigv(j)^2));
        end
        C(i,ii) = (1/p)*C(i,ii);
    end
end