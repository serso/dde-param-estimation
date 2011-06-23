function xResult = interpolate(x, N, method) 

xResult = interp1(1:1:length(x), x, 1 : (length(x)-1) / (N - 1) : length(x), method);

end