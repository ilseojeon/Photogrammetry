function GP = Comp_GP_with_Visualization(TP, IO, EO, BP)

Vis_IEO(IO,EO,BP)

for m = 1:2,
    R{m} = Rot3D(EO(m,4), EO(m,5), EO(m,6));
end

for n=1:size(TP,1),
    for m = 1:2,
        p{m} = [TP(n,2*m)-IO(1), TP(n,2*m+1)-IO(2), -IO(3)]';
        C{m} = EO(m,1:3)';
        p_G{m}(n,:) = (R{m}'*p{m} + C{m})';
    end
    A = [R{1}'*p{1}, -R{2}'*p{2}];
    y = C{2}-C{1};
    kt = inv(A'*A)*A'*y;
    for m = 1:2,
        P{m} = kt(m) * R{m}'*p{m} + C{m};
    end
    GP(n,:) = (P{1}+P{2})'/2;
    plot3([C{1}(1) P{1}(1)], [C{1}(2) P{1}(2)], [C{1}(3) P{1}(3)], 'm-');
    plot3([C{2}(1) P{2}(1)], [C{2}(2) P{2}(2)], [C{2}(3) P{2}(3)], 'b-');
    plot3(GP(n,1), GP(n,2), GP(n,3), 'r*');
end
plot3(p_G{1}(:,1), p_G{1}(:,2), p_G{1}(:,3), 'm.');
plot3(p_G{2}(:,1), p_G{2}(:,2), p_G{2}(:,3), 'b.');

