
function [OSy, OSyy] = wykresyZal(N, A)
alfan = wyznaczAlfe(A, N);
% alfan = sqrt(2*2)/2;
OSx = linspace(alfan-0.01,alfan+0.01,1000); % x - tworzy wektor punkt�w 
                                          %dla kt�rych rysowane s� wykresy
OSy = zeros(1,1000);
OSyy = zeros(1,1000);
A = @(z) det(MacierzA(z,N));
 % p�tla przypisuje wartos� wyznacznika macierzy dla ka�dego x
    for i=1:1000              
          OSy(i) = A(OSx(i));
    end

 A = @(z) cond(MacierzA(z,N));
% p�tla przypisuje warto�� wska�nika uwarunkowania macierzy dla ka�dego x 
    for i=1:1000                         
          OSyy(i) = A(OSx(i));
    end
end
