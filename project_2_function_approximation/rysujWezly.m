function [xn, yn] = rysujWezly(N)
F =@(u) (u+1/3).*(u+1/3)+exp(-u-2);
iks=@(w) (2*((1:w)-1))/(w-1)-1;  
x1000 = linspace(-1,1,1000);
for i = 1:N
    xn = iks(i);
    yn = F(xn);
end 
% figure;
% plot(xn,yn,'o',x1000,F(x1000));
% title(['Wykres funkcji f(x) na podstawie ci¹gu próbek N=',num2str(N)]);
% legend({'ci¹g próbek N','funkcja y(x)'},'Location','Best');
% xlabel("x");  %nazwa osi x
% ylabel("f(x)");   %nazwa osi y
% grid on;
% hold on;
end
