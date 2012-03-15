function [xd, thetad] = getDerivatives(derivatives, Nx, p)
    out = deal( derivatives );
    xd = out(1:Nx);
    thetad = out (Nx+1:Nx+p);
end