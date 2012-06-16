function [xd, pd] = getDerivatives(derivatives, n, np)
    out = deal( derivatives );
    xd = out(1:n);
    pd = out(n+1:n+np);
end