function [kt, et, vct, Dkt] = LSE (A, y, Dy);

P = inv(Dy);
N = A'*P*A;
c = A'*P*y;
iN = inv(N);
kt = iN * c;
et = y - A*kt;
vct = et'*P*et / (length(y) - rank(A));
Dkt = vct * iN;
