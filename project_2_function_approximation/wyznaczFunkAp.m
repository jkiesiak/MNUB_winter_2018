function [yZ] = wyznaczFunkAp(xK, p, K, xn)
z=100;
Fz=@(u) (u+1/3).*(u+1/3)+exp(-u-2);
x100 = linspace(-1,1,z);
for n = 1:z
    for k = 1:K
        if (abs(x100(n) - xK(k)) < 1)
            Wp(n,k) = 2*(abs(x100(n)-xK(k)))^3 - 3*((x100(n)-xK(k))^2)+1;            
        else
            Wp(n,k) = 0;
        end
    end
end
%mnzenie ka¿dego elementu do wykresu
for n=1:z
        yZ(n)=sum(Wp(n,:).*p');
end
% %%%Generacja wykresów
% figure;
% plot(x100,yZ,'r',xn,Fz(xn),'*',x1000,Fz(x1000),'g' );
% title(['Wykres funkcji f(x) po aproksymacji K=',num2str(K)]);
% legend({'aproksymacja funkcji','ci¹g próbek f(x)','funkcja y'},'Location','Best');
% xlabel("x");  %nazwa osi x
% ylabel("f(x)");   %nazwa osi y
% grid on;
% hold on;
end
