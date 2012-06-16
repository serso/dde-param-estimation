% Creates grid of N elements iteratig through min(t) = t(1) to max(t) =
% t(length(t)) (assuming t is natural sorted).
% Elements of t closest to grid knots are used instead of exact knots 
% (i.e. the grid is not uniformly distributed)
%
% @return tResult - grid of N elements (from min(t) to max(t) with use of some elelemts of t)
% @return tUsedResult - array, i-th element shows how many points of t are closest to current (i-th) knot 
function [tResult, tUsedResult, tEven] = ddeParamEst_createGrid(t, n)

tError = 10^-12;

tMin = t(1);
tMax = t(length(t));

% check if first element is less than last element
if ( tMax <= tMin )
    throw (MException ('ArgumentCheck:OutOfRange', 'tMax has to be more than tMin'));
end

% create uniformly distributed t grid with n elements
tGrid = tMin : ( tMax - tMin ) / (n - 1) : tMax;
tEven = true;

% prepare result values
tResult = zeros(n, 1);
tUsedResult = zeros(n, 1);

% index for t array
tIndex = 1;

tLength = length(t);

% index for result array
i = 1;

for tPos = tGrid
    
    if ( t(tIndex) < tPos + tError)
%         if ( abs(t(tIndex) - tPos ) >= tError )
%             tEven = false;
%         end 
        
        tUsedResult(i) = 1;
        tResult(i) = t(tIndex);
        tIndex = tIndex + 1;
        
        if ( tIndex <= tLength && t(tIndex) <= tPos )
            
            warning ('ArgumentCheck:IllegalArgument', 'More that one element found in one step. Increase number of elements of grid (n) to prevent this warning!');
            fprintf('tIndex = %i, t(tIndex) = %f, tPos = %f\n', tIndex, t(tIndex), tPos);
            
            % we need to skip ALL elements of t array which belong to
            % current step
            while (tIndex <= length(t) && t(tIndex) <= tPos)
                % next element in still in current step => increase indices
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

tEven = true;

% if ( ~tEven )
%     disp('Warning: supplied t grid is not even!');
% end

end