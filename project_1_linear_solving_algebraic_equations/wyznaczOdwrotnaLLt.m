function Y = wyznaczOdwrotnaLLt(B)
% A = [1 2 3; 
%     2 5 8; 
%     3 8 14];
L = chol(B, 'lower');
Lt=L';
[N, N] = size(B);
[N, N] = size(L);
[N, N] = size(Lt);
Y = zeros(N);
E = eye(N);
X = zeros(N);
% L*X = E
for i = 1:N
    for j = 1:N
            s1 = 0;
            for k = 1 : j-1
                s1 = s1 + X(k,i)*L(j,k);     
            end
            X(j,i) = ( E(j,i) - s1 ) * ( 1 / L(j,j) ); 
    end
end
% checkEye1 = L*X

% Lt*Y=X
% wyszukujemy Y

for i = 1 : N
    for j = N : -1 : 1
            s2 = 0;
            for v =  N: -1 : j+1
                s2 = s2 + Y(v,i)*Lt(j,v);     
            end
            Y(j,i) = ( X(j,i) - s2 ) * ( 1 / Lt(j,j) ); 
    end
end
Y = Y
% check2shouldbeX = Lt*Y
% X = X
% checkbeEye = B*Y
% matrixInput = L*Lt
% Eye = L*Lt*Y
% check3 = A*check2
% check4 = P'*L1*U*Y
end