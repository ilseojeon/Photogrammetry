% input
% m: number of unknowns, n: number of observations
% A; Design matrix %nxm
% y; observations  %nx1
% Dy; dispersion of observation y % mxm

% N; square matrix %nxn
% c; the right back side of least square equation ksihat=inv(A'PA)*(c=A'Py)

% output
% kt; ksihat
% et; residuals
% vct; variance component estimate
% Dkt; cofactor estimate

function [kt, et, vct, Dkt] = LSE (A, y, Dy);

P = inv(Dy);
N = A'*P*A;
c = A'*P*y;
iN = inv(N);
kt = iN * c;
et = y - A*kt;
vct = et'*P*et / (length(y) - rank(A));
Dkt = vct * iN;
