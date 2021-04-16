function Vis_IEO(IO,EO,BP)

no_IM = size(EO,1);
for m = 1:no_IM-1,
    base(m,:) = EO(m+1,:) - EO(m,:);
    lenBase(m) = norm(base(m,:));
end
avgLenBase = mean(lenBase);
lenArr = avgLenBase / 5;
for m = 1:no_IM,
    plot3(EO(m,1), EO(m,2), EO(m,3), 'rx');
    text(EO(m,1), EO(m,2), EO(m,3), sprintf('C%d',m) );
    R{m} = Rot3D(EO(m,4), EO(m,5), EO(m,6));
    for k = 1:3,
        h = quiver3(EO(m,1), EO(m,2), EO(m,3), ...
            lenArr*R{m}(k,1), lenArr*R{m}(k,2), lenArr*R{m}(k,3));
        set( h, 'LineWidth', 2 );
        text(EO(m,1)+lenArr*R{m}(k,1), EO(m,2)+lenArr*R{m}(k,2), ...
            EO(m,3)+lenArr*R{m}(k,3), sprintf('%c_c', 'W'+k) );
    end
end
no_BP = size(BP,1);
for m = 1:no_IM,
    for n = 1:no_BP,
        BP_G{m}(n,:) = R{m}'*[BP(n,2)-IO(1), BP(n,3)-IO(2), -IO(3)]' + EO(m,1:3)';
    end
    PP_G{m} = R{m}'*[0, 0, -IO(3)]' + EO(m,1:3)';
    idx = [1 2 3 4 1];
    plot3(BP_G{m}(idx,1), BP_G{m}(idx,2), BP_G{m}(idx,3), 'r-');
    plot3(PP_G{m}(1), PP_G{m}(2), PP_G{m}(3), 'rx');
    plot3([EO(m,1),PP_G{m}(1)], [EO(m,2),PP_G{m}(2)], [EO(m,3),PP_G{m}(3)], 'r:');
    text(PP_G{m}(1), PP_G{m}(2), PP_G{m}(3), 'PP');
end

