function [bladSra, bladMaxa] = wyznaczBledy(yZ) 
x = linspace(-1,1,100);
Fun =@(u) (u+1/3).*(u+1/3)+exp(-u-2);
bladSra=norm(yZ-Fun(x))/norm(Fun(x));
bladMaxa=norm(yZ-Fun(x),Inf)/norm(Fun(x),Inf);
end



