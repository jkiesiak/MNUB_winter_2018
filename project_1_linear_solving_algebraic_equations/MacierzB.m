%Wygenerowanie macierzy o zadanym rozmiarze n i zadanej wartosci x
function B = MacierzB(k,N)
%A-zwracana wygenerowana zgodnie ze wzorem macierz
B=zeros(N);
% x = @(z) sqrt((z)^2 + 0.5) - 1; %zadana w poleceniu funkcja
x = 2^k / 300; %zadana w poleceniu funkcja
% A(1,1) = (x(alfa))^2; %macierz min. 1x1, wiêc wyrzucam przed pêtlê
B(1,1) = (x)^2; %macierz min. 1x1, wiêc wyrzucam przed pêtlê
for i=N:-1:1   %Iteracja wierszy (od koñca)
    for j=N:-1:1   %Iteracja kolumn (od koñca)
        if (i>=2 && j>=2)
            B(i,j) = min(i,j)*4/9;
        elseif ((i==1 && j>1) || (j==1 && i>1 ))
            B(i,j)=(2/3)*x;
        end
    end
end


end