function yi = interp1nosort(x, Y, xi)

if any(diff(x) == 0)
    yi = nan;
    return;
end
    
[x, ind] = sort(x);
yi = interp1(x, Y(ind), xi);
