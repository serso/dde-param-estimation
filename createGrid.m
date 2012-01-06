% Creates grid of N elements iteratig through min(t) = t(1) to max(t) =
% t(length(t)) (assuming t is natural sorted).
% Elements of t closest to grid knots are used instead of exact knots 
% (i.e. the grid is not uniformly distributed)
%
% @return tResult - grid of N elements (from min(t) to max(t) with use of some elelemts of t)
% @return tUsedResult - array, i-th element shows how many points of t are closest to current (i-th) knot 
function [tResult, tUsedResult] = createGrid(t, N)

eps = 0.0000001;

tMin = t(1);
tMax = t(length(t));

if ( tMax <= tMin )
    throw (MException ('ArgumentCheck:OutOfRange', 'tMax has to be more than tMin'));
end

tGrid = tMin : ( tMax - tMin ) / ( N - 1 ) : tMax;
tResult = zeros(N, 1);

tUsedResult = zeros(N, 1);

% index for t array
tIndex = 1;

%index for result array
i = 1;

for tPos = tGrid
    
    if ( t(tIndex) < tPos + eps)
        tUsedResult(i) = 1;
        
        tResult(i) = t(tIndex);
        tIndex = tIndex + 1;
        
        if ( tIndex <= length(t) && t(tIndex) <= tPos )
            
            warning ('ArgumentCheck:IllegalArgument', 'More that one element found in one step. Increase number of elements of grid (N) to prevent this warning!');
            
            while (tIndex <= length(t) && t(tIndex) <= tPos)
                tUsedResult(i) = tUsedResult(i) + 1;
                tIndex = tIndex + 1;
            end
        end
    else
        tUsedResult(i) = 0;
        tResult(i) = tPos;
    end;
    
    i = i + 1;
end

end