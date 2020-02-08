function [xn, yn] = wyliczWezly(N)
F =@(u) (u+1/3).*(u+1/3)+exp(-u-2);
iks=@(w) (2*((1:w)-1))/(w-1)-1;  
x1000 = linspace(-1,1,1000);
for i = 1:N
    xn = iks(i);
    yn = F(xn);
end 

end
