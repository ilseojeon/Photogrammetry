function stat = comp_stat(va)

nv = length(va);
stat(1,:) = min(va);   % min
stat(2,:) = max(va);   % max
stat(3,:) = mean(va);   % average
stat(4,:) = sqrt(sum((va-stat(3,:)).^2)/(nv-1));  % std. dv.
stat(5,:) = sqrt(sum((va).^2)/(nv-1));      % RMSE