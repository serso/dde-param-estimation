function [xResult, thetaResult] = getDerivatives(derivatives, Nx, p)
out = deal( derivatives );
xResult = out(1:Nx);
thetaResult = out (Nx+1:Nx+p);
end