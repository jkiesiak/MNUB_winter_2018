clear all


tt = 0:0.01:10;
d = length(tt);
%rozwi¹zanie równania funkcj¹ ode45
opts = odeset('RelTol',2e-2,'AbsTol',1e-6);
[T, Y] = ode45(@funk, tt, [0 1], opts);
%generacja wykresu za pomoc¹ funkcji ode45
figure(1)
plot(T,Y(:,1),'-',T,Y(:,2))
title('12y''=-4y''-3y');
xlabel('Time t');
ylabel('Solution y');
grid on 
hold on

h = 0.01;

%rozwi¹zanie równania za pomoc¹ metody Radau
[Tr, Yr] = metRadau(h, d);
figure(1)
plot(Tr,Yr(:,1),'-.',Tr,Yr(:,2),'-')
legend('y_1', 'y_2','y_1metRad', 'y_2metRad');
hold on
% BLEDY
[det2r, detInfr] = wyznaczBledy(Y, Yr);
figure(5)
plot(tt,log10(det2r(:,1)),tt,log10(det2r(:,2)));     %rysuje wykres
% plot(tt,(det2r(:,1)),tt,(det2r(:,2)));     %rysuje wykres
title("Wykres wskaŸnika b³êdu œredniokwadratowego w skali logarytmicznej");
xlabel('Time t');
ylabel('Value of error');
grid on;     %tworzenie linie na wykresie
hold on;
figure(6)
plot(tt,log10(detInfr(:,1)),tt,log10(detInfr(:,2)));     %rysuje wykres
% set(bladMax,'ZScale','log');
title("Wykres wskaŸnika b³êdu maksymalnego w skali logarytmicznej");
xlabel('Time t');
ylabel('Value of error');
grid on;       %tworzenie linie na wykresie
hold on;


%rozwi¹zanie równania za pomoc¹ otwartej metody Eulera
[Te, Ye] = euler(h, d);
figure(1)
plot(Te,Ye(:,1),'-',Te,Ye(:,2),'-.')
legend('y_1', 'y_2','y_1metRad', 'y_2metRad','y_1Eul', 'y_2Eul');
hold on
% BLEDY 
[det2e, detInfe] = wyznaczBledy(Y, Ye);
figure(5)
plot(tt,log10(det2e(:,1)),tt,log10(det2e(:,2)));     %rysuje wykres
title("Wykres wskaŸnika b³êdu œredniokwadratowego w skali logarytmicznej");
xlabel('Time t');
ylabel('Value of error');
legend('y sr Rad','y sr Eul');
grid on;     %tworzenie linie na wykresie
hold on;
figure(6)
plot(tt,log10(detInfe(:,1)),tt,log10(detInfe(:,2)));     %rysuje wykres
% set(bladMax,'ZScale','log');
title("Wykres wskaŸnika b³êdu maksymalnego w skali logarytmicznej");
xlabel('Time t');
ylabel('Value of error');
legend('y max Rad','y max Eul');
grid on;       %tworzenie linie na wykresie
hold on;



%% functions
function dydt = funk(t, y) %funkcja zadana
A = [0, 1; -(1/4), -(1/3)];
dydt = A*[ y(1); y(2)];
end

function [t, y] = euler(h, d) % otwarta metoda Eulera
y = zeros(2,d);
t = zeros(d,1);
%wartoœci pocz¹tkowe
y(1,1)=0;
y(2,1)=1;
for i=1:d-1
    t(i+1) = t(i)+h;
    y(:,i+1)=y(:,i) + h*funk(t(i), y(:,i));
end
y=y';
end


function [t, y] = metRadau(h, d) %metoda Radau
t = 0:h:10;
I = eye(2,2);
Z = zeros(2,2);
%nasza macierzA z funkcji funk
A = [0, 1; -(1/4), -(1/3)];
sqrt6=(6)^(1/2);
%wspó³czynniki z tablicy Butchera
B = [I-(1/9)*A*h, -((-1 - sqrt6)/18)*A*h, -((-1 + sqrt6)/18)*A*h;
     (-1/9)*A*h, I-(11/45 + 7*sqrt6/360)*A*h, -(11/45 - 43*sqrt6/360)*A*h;
     (-1/9)*A*h, -(11/45 + 43*sqrt6/360)*A*h, I-(11/45 - 7*sqrt6/360)*A*h];
y = zeros(2,d);
%wartoœci pocz¹tkowe
y(1,1) = 0;
y(2,1) = 1;
for i=1:d-1
    t(i+1)=t(i)+h;
    Ay = [A*y(:,i);A*y(:,i);A*y(:,i)];
    F = B\Ay;
    y(:,i+1)=y(:,i) + h*((1/6)*F(1:2)+(2/3)*F(3:4)+(1/6)*F(5:6));
end
t=t';
y=y';

end

function [det2, detInf] = wyznaczBledy(Y, Ywe)
[x1, x2] =size(Y);
for k =1:x2
    for i = 1:x1
        det2(i,k)=norm(Ywe(i,k)-Y(i,k)) / norm(Y(i,k));
        detInf(i,k)=norm(Ywe(i,k)-Y(i,k),Inf) / norm(Y(i,k));
    end
end
end