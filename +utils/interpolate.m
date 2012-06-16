function xResult = interpolate(x, n, method) 

xResult = interp1(1:length(x), x, 1 : (length(x)-1) / (n - 1) : length(x), method);

end