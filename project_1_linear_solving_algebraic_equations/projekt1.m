clear all
close all
  %Zadanie 1
A = MacierzA(alfa,N);

% %Zadanie 2
alfa = detMacierzA(A, N);

%zadanie 3
Y1 = wyznaczOdwrotnaLLt(A);
Y = wyznaczOdwrotnaLU(A);

%zadanie 4 i 5
B = MacierzB(k,N);

%% Zadanie 1
%Wygenerowanie macierzy o zadanym rozmiarze n i zadanej wartosci x
function A=MacierzA(alfa,N)
%A-zwracana wygenerowana zgodnie ze wzorem macierz
A=zeros(N);
x = sqrt((alfa)^2 + 0.5) - 1; %zadana w poleceniu funkcja
A(1,1) = (x)^2; %macierz min. 1x1, wiÍc wyrzucam przed pÍtlÍ
for i=N:-1:1   %Iteracja wierszy (od koÒca)
    for j=N:-1:1   %Iteracja kolumn (od koÒca)
        if (i>=2 && j>=2)
            A(i,j) = min(i,j)*4/9;
        elseif ((i==1 && j>1) || (j==1 && i>1 ))
            A(i,j)=(2/3)*x;
        end
    end
end
end

%% Zadanie 2

%Znalezienie najmniejszej dodatniej wartosci alfa, dla ktorej wartosc
%wyznacznika macierzy o danym rozmiarze wynosi 0
function alfa = detMacierzA(A, N)
[N, N] = size(A);
A = @(z) det(MacierzA(z, N));
for i = linspace(0, 5 , 100000)      
   if (A(i) < 10^(-14))             %W praktyce wartoúÊ wyznacznika musi dπøyÊ do zera
       alfa = i;
       break;
   end
end
end

%blok wylicza wartoúci z zaleønoúci wyznacznika macierzy oraz 
%wskaünika uwarunkowania kaødej z macierzy od alfa
function [OSy, OSyy] = wykresyZal(N, A)
alfan = detMacierzA(A, N);
% alfan = sqrt(2*2)/2;
OSx = linspace(alfan-0.01,alfan+0.01,1000); % x - tworzy wektor punktÛw 
                                          %dla ktÛrych rysowane sπ wykresy
OSy = zeros(1,1000);
OSyy = zeros(1,1000);
A = @(z) det(MacierzA(z,N));
 % pÍtla przypisuje wartosÊ wyznacznika macierzy dla kaødego x
    for i=1:1000              
          OSy(i) = A(OSx(i));
    end

 A = @(z) cond(MacierzA(z,N));
% pÍtla przypisuje wartoúÊ wskaünika uwarunkowania macierzy dla kaødego x 
    for i=1:1000                         
          OSyy(i) = A(OSx(i));
    end
end

%blok rysuje wykresy z zaleønoúci wyznacznika macierzy oraz 
%wskaünika uwarunkowania kaødej z macierzy od alfa
% alfa = 0;
N = [3 10 20];
z = 1000;
alfan = sqrt(2*2)/2;

OSx = linspace(alfan-0.01,alfan+0.01,1000);
OSy = zeros(3,z);
OSyy = zeros(3,z);
OSp = zeros(3,z);
OSpp = zeros(3,z);
j = 0;

for i = 1 : length(N)
    A = MacierzA(alfan,N(i));
    alfa = detMacierzA(A, N(i));
    [OSp, OSpp] = wykresyZal(N(i), A);    
    j = j+1;
    OSy(j,:) = OSp;
    OSyy(j,:) = OSpp;
end

figure;
semilogy(OSx, OSy(1,:), OSx, OSy(2,:), OSx, OSy(3,:));
axis([alfan-0.01 alfan+0.01 0 inf])     
xlabel('Alfa w przedziale alfa-0.01,alfa+0.01');
ylabel('Wartoúci wyznacznika macierzy');
title('Wykres wyznacznika macierzy od alfa');
legend('wyznacznik macierzy N=3','wyznacznik macierzy N=10','wyznacznik macierzy N=20', 'Best', 'Best');
grid on;
hold on;

figure;
semilogy(OSx, OSyy(1,:), OSx, OSyy(2,:), OSx, OSyy(3,:));
axis([alfan-0.01 alfan+0.01 0 inf])     
xlabel('Alfa w przedziale alfa-0.01,alfa+0.01');
ylabel('Wartoúci wskaünika uwarunkowania');
title('Wykres wskaünika uwarunkowania od alfa');
legend('wskaünik uwarunkowania N=3','wskaünik uwarunkowania N=10','wskaünik uwarunkowania N=20', 'Best', 'Best');
grid on;
hold on;



%% Zadanie 3
function Y1 = wyznaczOdwrotnaLLt(A)
%rozk≥ad Cholewskiego-Banachiewicza
L = chol(A, 'lower');
Lt=L';          %transpozycja dolnej macierzy trÛjkπtnej
[N, N] = size(L);
[N, N] = size(Lt);
Y1 = zeros(N);
E1 = eye(N);         %macierz jednostkowa
X1 = zeros(N);
% L*X = E              ;wyszukujemy X
for ii = 1:N
    for jj = 1:N
            s1 = 0;
            for kk = 1 : jj-1
                s1 = s1 + X1(kk,ii)*L(jj,kk);     
            end
            X1(jj,ii) = ( E1(jj,ii) - s1 ) * ( 1 / L(jj,jj) ); 
    end
end
% checkEye1 = L*X1 ;
% Lt*Y=X               ;wyszukujemy Y
for i = 1 : N
    for j = N : -1 : 1
            s2 = 0;
            for v =  N: -1 : j+1
                s2 = s2 + Y1(v,i)*Lt(j,v);     
            end
            Y1(j,i) = ( X1(j,i) - s2 ) * ( 1 / Lt(j,j) ); 
    end
end
% check2shouldbeX1 = Lt*Y1 ;
% X1 = X1 ;
% checkbeEye1 = A*Y1 ;
% matrixA1 = L*Lt; 
% Eye1 = L*Lt*Y1 ;
end

function Y = wyznaczOdwrotnaLU(A)
[L1, U, P] = lu(A);         %P - macierz permutacji
[N, N] = size(L1);
[N, N] = size(U);
E = eye(N);
X = zeros(N);
Y = zeros(N);
% L*X = E              ;wyszukujemy X
for ii = 1:N
    for jj = 1:N
        s1 = 0;
            for kk = 1 : jj-1
                s1 = s1 + X(kk,ii)*L1(jj,kk);     
            end
        X(jj,ii) = ( E(jj,ii) - s1 ) * ( 1 / L1(jj,jj) ); 
    end
end
checkShouldBeEye2 = L1*X ;
% U*Y=X                ;wyszukujemy Y
for i = 1 : N
    for j = N : -1 : 1
            s2 = 0;
            for v =  N: -1 : j+1
                s2 = s2 + Y(v,i)*U(j,v);     
            end
            Y(j,i) = ( X(j,i) - s2 ) * ( 1 / U(j,j) ); 
    end
end
% check2shouldbeX2 = U*Y ;
% X2 = X ;
% checkbeEye2 = L1*U*Y ;
end

%% Zadanie 4 i 5
function B = MacierzB(k, N)
%B-zwracana wygenerowana zgodnie ze wzorem macierz      ;analogicznie do
%MacierzA
B=zeros(N);
x = 2^k / 300; %zadana w poleceniu funkcja
B(1,1) = (x)^2; %macierz min. 1x1
for i=N:-1:1   %Iteracja wierszy (od koÒca)
    for j=N:-1:1   %Iteracja kolumn (od koÒca)
        if (i>=2 && j>=2)
            B(i,j) = min(i,j)*4/9;
        elseif ((i==1 && j>1) || (j==1 && i>1 ))
            B(i,j)=(2/3)*x;
        end
    end
end
end

%Norma druga i nieskonczona
function [deltaSr, deltaMx ]= wskBlad(B)
    deltaSr = max(sqrt(eigs(B' * B)));
    deltaMx = max(sum(abs(B), 2));
end

%Wykresy wzkaünikÛw
%inicjalizacja zmiennych
z = 1000;
N = [3 10 20];
k = 0:1:21;
K = 22;
B = zeros(z);
X = 2.^k / 300;

det2_inv = zeros(2, K);
det2_LLt = zeros(2, K);
det2_LU = zeros(2, K);
detInf_inv = zeros(2, K);
detInf_LLt = zeros(2, K);
detInf_LU = zeros(2, K);

for i = 1 : length(N)
    z=0;
    for j = 1 : length(k)
        z=z+1;
        %generacja macierzy jednostkowej i permutacji w zaleønoúci od
        %rozmiaru N
        I = eye(N(i));
        P = zeros(N(i));
        B = MacierzB(k(j),N(i));     
        
% dla polecenia wbudowanego inv
        B_inv = inv(B);
        BB_inv = B*B_inv - I;
        [deltaSr, deltaMx ]= wskBlad(BB_inv);
        det2_inv(1, z) = deltaSr;
        detInf_inv(1, z) = deltaMx;
        det2_inv(2, z) = norm(BB_inv, 2);
        detInf_inv(2, z) = norm(BB_inv, inf);

        
% dla metody LLt
        B_LLt = wyznaczOdwrotnaLLt(B);
        BB_LLt = B*B_LLt - I;
        [deltaSr, deltaMx ]= wskBlad(BB_LLt);
        det2_LLt(1, z) = deltaSr;
        detInf_LLt(1,z) = deltaMx;        
        det2_LLt(2, z) = norm(BB_LLt, 2);
        detInf_LLt(2, z) = norm(BB_LLt, inf);

% dla metody LU
        [B_LU, P] = wyznaczOdwrotnaLU(B);
        BB_LU = B*B_LU*P - I;
        [deltaSr, deltaMx ]= wskBlad(BB_LU);     
        det2_LU(1, z) = deltaSr;
        detInf_LU(1, z) = deltaMx;
        det2_LU(2, z) = norm(BB_LU, 2);
        detInf_LU(2,z) = norm(BB_LU, inf);
    end

    figure;
%     loglog(X, det2_inv(1,:),'-o', X, det2_LLt(1,:),'-o', X, det2_LU(1,:),'-o');
    loglog(X, det2_inv(1,:)./X,'-o', X, det2_LLt(1,:)./X,'-o', X, det2_LU(1,:)./X,'-o');
        %dla bledu wzglednego
    xlabel('Wartoúci X(k)');
    ylabel('Wartoúci b≥Ídu úredniokwadratowego wzglÍdnego');
    title('Wykres wskaünika b≥Ídu úredniokwadratowego wzglÍdnego');
    legend('b≥πd úred. inv','b≥πd úred. LLt','b≥πd úred. LU','Location', 'Best');
%         legend('b≥πd úred. inv','b≥πd úred. LLt','b≥πd úred. LU','Location', 'Best');
    grid on;
    hold on;

    figure;
%     loglog(X, detInf_inv(1,:),'-o', X, detInf_LLt(1,:),'-o', X, detInf_LU(1,:),'-o');
    loglog(X, detInf_inv(1,:)./X,'-o', X, detInf_LLt(1,:)./X,'-o', X, detInf_LU(1,:)./X,'-o');
    xlabel('Wartoúci X(k)');
    ylabel('Wartoúci b≥Ídu maksymalnego wzglÍdnego');
    title('Wykres wskaünika b≥Ídu maksymalnego wzglÍdnego');
    legend('b≥πd max. inv','b≥πd max. LLt','b≥πd max. LU', 'Location', 'Best');
%         legend('b≥πd max. inv','b≥πd max. LLt','b≥πd max. LU', 'Location', 'Best');
    grid on;
    hold on;

%         figure;
%     semilogy(X, det2_inv(2,:), X, det2_LLt(2,:), X, det2_LU(2,:));
%     xlabel('Alfa w przedziale alfa-0.01,alfa+0.01');
%     ylabel('Wartoúci b≥Ídu');
%     title('Wykres wskaünika b≥Ídu úredniokwadratowego');
%     legend('wyznacznik macierzy N=3','wyznacznik macierzy N=10','wyznacznik macierzy N=20', 'Best', 'Best');
%     grid on;
%     hold on;
% 
%     figure;
%     semilogy(X, detInf_inv(2,:), X, detInf_LLt(2,:), X, detInf_LU(2,:));
%     xlabel('Alfa w przedziale alfa-0.01,alfa+0.01');
%     ylabel('Wartoúci b≥Ídu maksymalnego');
%     title('Wykres wskaünika b≥Ídu maksymalnego');
%     legend('wskaünik uwarunkowania N=3','wskaünik uwarunkowania N=10','wskaünik uwarunkowania N=20', 'Best', 'Best');
%     grid on;
%     hold on;
end
