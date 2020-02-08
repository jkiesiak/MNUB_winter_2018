function [p, xK] = macierzP(xn, yn, K)
iks=@(u) (2*((1:u)-1))/(u-1)-1;  
[z,z]=size(xn);
W = zeros(z,K);
xK = zeros(K);
xK = iks(K);
for n = 1:z
    for k = 1:K
        if (abs(xn(n) - xK(k)) < 1)
            W(n,k) = 2*(abs(xn(n)-xK(k)))^3 - 3*((xn(n)-xK(k))^2)+1;            
        else
            W(n,k) = 0;
        end 
    end
end
p = (W'*W) \ (W'*(yn)');
end