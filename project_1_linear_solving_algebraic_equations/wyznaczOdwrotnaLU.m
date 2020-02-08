function [Y, P] = wyznaczOdwrotnaLU(B)
% B = [1 2 3; 
%     2 5 8; 
%     3 8 14];
% L1 = chol(B, 'lower');
[L1, U, P] = lu(B);
[N, N] = size(B);
[N, N] = size(L1);
[N, N] = size(U);
E = eye(N);
X = zeros(N);
Y = zeros(N);
for i = 1:N
    for j = 1:N
        s1 = 0;
            for k = 1 : j-1
                s1 = s1 + X(k,i)*L1(j,k);     
            end
        X(j,i) = ( E(j,i) - s1 ) * ( 1 / L1(j,j) ); 
    end
end
% checkShouldBeEye = L1*X

% U*Y=X
% wyszukujemy Y
for i = 1 : N
    for j = N : -1 : 1
            s2 = 0;
            for v =  N: -1 : j+1
                s2 = s2 + Y(v,i)*U(j,v);     
            end
            Y(j,i) = ( X(j,i) - s2 ) * ( 1 / U(j,j) ); 
    end
end
Y = Y
% check2shouldbeX = U*Y
% X = X
% checkbeEye = L1*U*Y
% P'*B*Y
% matrixInput = P'*L1*U
% check3 = A*check2
% check4 = P'*L1*U*Y
end
