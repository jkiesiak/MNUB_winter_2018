%Znalezienie najmniejszej dodatniej wartosci alfa, dla ktorej wartosc
%wyznacznika macierzy o danym rozmiarze wynosi 0

function alfa = detMacierzA(A, N)
[N, N] = size(A);
W = zeros(N);
W = @(z) det(MacierzA(z, N));

for i = linspace(0, 5 , 100000)      
%    if (A(i) < 0)          
   if (W(i) < 10^(-14))             %W praktyce wartoœæ wyznacznika musi d¹¿yæ do zera
       alfa = i;
       break;
   end
end
end