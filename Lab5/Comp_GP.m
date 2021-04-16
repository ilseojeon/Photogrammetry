function GP = Comp_GP(TP, IO, EO)

for m = 1:2,
    R{m} = Rot3D(EO(m,4), EO(m,5), EO(m,6));
end

for n=1:size(TP,1),
    for m = 1:2,
        p{m} = [TP(n,2*m)-IO(1), TP(n,2*m+1)-IO(2), -IO(3)]';
        C{m} = EO(m,1:3)';
    end
    A = [R{1}'*p{1}, -R{2}'*p{2}];
    y = C{2}-C{1};
    kt = inv(A'*A)*A'*y;
    for m = 1:2,
        P{m} = kt(m) * R{m}'*p{m} + C{m};
    end
    GP(n,:) = (P{1}+P{2})'/2;
end
