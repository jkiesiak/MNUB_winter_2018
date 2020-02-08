clear all
close all

N = [10 20 30];
K = 8;
x100= linspace(-1,1,100);
x1000 = linspace(-1,1,1000);


%wz�r funkcji - funkcja anonimowa
F=@(u) (u+1/3).*(u+1/3)+exp(-u-2);
iks=@(w) (2*((1:w)-1))/(w-1)-1;  


%%%Zadanie 1 - sprz�dzanie wykresu funkcji
for i=1:length(N)
    [xn, yn] = wyliczWezly(N(i));   
    
    figure;
    plot(xn,yn,'o',x1000,F(x1000));
    title(['Wykres funkcji f(x) na podstawie ci�gu pr�bek N=',num2str(N(i))]);
    legend({'ci�g pr�bek N','funkcja y(x)'},'Location','Best');
    xlabel("x");  %nazwa osi x
    ylabel("f(x)");   %nazwa osi y
    grid on;
    hold on;
end


%%%Zadanie 2 - aproksymacja funkcji dla kilku par N i K
N2 = [10 20 30];
K2 = [3 6 15];
for i = 1:length(N2)
    for j = 1:length(K2)
        [xn, yn] = wyliczWezly(N2(i));     
        [p, xK] = macierzP(xn, yn, K2(j));
        [yZ] = wyznaczFunkAp(xK, p, K2(j), xn);
        
        figure;
        plot(x100,yZ,'r',xn,F(xn),'*',x1000,F(x1000),'g' );
        title(['Wykres funkcji f(x) po aproksymacji K= ',num2str(K2(j)),' dla N=',num2str(N2(i))]);
        legend({'aproksymacja funkcji','ci�g pr�bek f(x)','funkcja y'},'Location','Best');
        xlabel("x");  %nazwa osi x
        ylabel("f(x)");   %nazwa osi y
        grid on;
        hold on;
    end

end

N=5:50;
%%%Zadanie 3 - wykresy b��d�w
for i = N 
[xn, yn] = wyliczWezly(i);    
    for j = 1: i-1
    [p, xK] = macierzP(xn, yn, j);
    [yZ] = wyznaczFunkAp(xK, p, j);
    [bladSra, bladMaxa] = wyznaczBledy(yZ);    
    bladSr(i,j) = bladSra;
    bladMax(i,j) = bladMaxa;    
    end    
end
figure;
surf(log10(bladSr));     %rysuje wykres
title("Wykres wska�niks b��du �redniokwadratowego w skali logarytmicznej");
xlabel("N");  %nazwa osi x
ylabel("K");    %nazwa osi y
zlabel("log10(z b��du �redniokwadratowego)"); %nazwa osi z
grid on;     %tworzenie linie na wykresie
hold on;
figure;
surf(log10(bladMax));     %rysuje wykres
% set(bladMax,'ZScale','log');
title("Wykres wska�nika b��du maksymalnego w skali logarytmicznej");
xlabel("N");  %nazwa osi x
ylabel("K");    %nazwa osi y
zlabel("log10(wska�nik b��du maksymalnego)"); %nazwa osi z
grid on;       %tworzenie linie na wykresie
hold on;



%%%Zadanie 4 - aproksymacja przy zaburzonych danych
bladSred1=[];
for z = logspace(-5,-1)
    for i = N   %N = 5:50
    [xn, yn] = wyliczWezly(i);    
    yn = yn + randn(size(yn))*z; %addytywne zaburzenie danych
        for j = 3: i-1         %K<N
        [p, xK] = macierzP(xn, yn, j);
        [yZ] = wyznaczFunkAp(xK, p, j);
        [bladSra, bladMaxa] = wyznaczBledy(yZ);  
        bladSred(i-4,j-2) = bladSra;
        end
    end
    bladSred1 = [bladSred1;min(min(bladSred(bladSred>0)))];
end   
figure;
loglog(logspace(-5,-1), bladSred1, 'o');
xlabel("delta [10-5;10-1]");
ylabel("Warto�ci b��du �redniokwadratowego");
title("Wykres warto�ci b��du �redniokwadratowego");
hold on;
p=polyfit(logspace(-5,-1),bladSred1',7);
semilogx(logspace(-5,-1),polyval(p, logspace(-5,-1)),'r');
hold on;
grid on;

%% Funkcja wyliczaj�ca w�z�y
function [xn, yn] = wyliczWezly(N)
F =@(u) (u+1/3).*(u+1/3)+exp(-u-2);
iks=@(w) (2*((1:w)-1))/(w-1)-1;  
x1000 = linspace(-1,1,1000);
for i = 1:N
    xn = iks(i);
    yn = F(xn);
end 
end


%% Funkcja wyznaczaj�ca macierz-fi
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


%% Funkcja aproksymuj�ca w�z�y
function [yZ] = wyznaczFunkAp(xK, p, K, xn)
z=100;
Fz=@(u) (u+1/3).*(u+1/3)+exp(-u-2);
x100 = linspace(-1,1,z);
x1000 = linspace(-1,1,1000);
for n = 1:z
    for k = 1:K
        if (abs(x100(n) - xK(k)) < 1)
            Wp(n,k) = 2*(abs(x100(n)-xK(k)))^3 - 3*((x100(n)-xK(k))^2)+1;            
        else
            Wp(n,k) = 0;
        end
    end
end
%mnzenie ka�dego elementu do wykresu
for n=1:z
        yZ(n)=sum(Wp(n,:).*p');
end
end


%% Funkcja zwracaj�ca b��d �redniokwadratowy i maksymalny
function [bladSra, bladMaxa] = wyznaczBledy(yZ) 
x = linspace(-1,1,100);
Fun =@(u) (u+1/3).*(u+1/3)+exp(-u-2);
bladSra=norm(yZ-Fun(x))/norm(Fun(x));
bladMaxa=norm(yZ-Fun(x),Inf)/norm(Fun(x),Inf);
end

