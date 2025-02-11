function i = find_nearest(surf, x)
    [~,i] = min(abs(surf-x));
end