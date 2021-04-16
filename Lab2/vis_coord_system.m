function vis_coord_system (pos, rot, len, coord_name, cstr)

cs = {'r', 'g', 'b'};
hold on
for n = 1:3 %%축이 3개라서 3까지임
    cstr = cs{n};
    pt = pos + len * rot(:,n);
    hold on
    h = plot3([pos(1), pt(1)], [pos(2), pt(2)], [pos(3), pt(3)]); %%한 축의 x, y, z 값
    set(h, 'LineWidth', 2, 'Color', cstr);
    h = text( pt(1),  pt(2),  pt(3), sprintf('%c', n+'x'-1) );
    set(h, 'Color', cstr, 'FontWeight', 'bold', 'FontSize', 9, ...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
end
cstr = 'k';
h = text( pos(1),  pos(2),  pos(3), coord_name );
set(h, 'Color', cstr, 'FontWeight', 'bold', 'FontSize', 9, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
