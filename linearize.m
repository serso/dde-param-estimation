% Function linearize x grid.
%
% @return grid of [length(x) - 1] * (numberOfDotes + 1) + 1 size
% e.g. if initial x grid has 3 points, and number of added dotes is 3
% then we will have result grid of 10 points where lineResult(1) = x(1)
% lineResult(5) = x(2) and lineresult(10) = x(3). Other points will be set
% as linear approximation of specified points
%
% @param x - intial grid to be linerized
% @param numberOfDotes - number of dotes for linear approximation between
% two next points of specified grid x
function result = linearize (x, numberOfDotes)

size = (length(x) - 1) * (numberOfDotes + 1) + 1;

result = zeros(size, 1);

for k = 1:1:(length(x) - 1)
    
    for m = 1:1:(numberOfDotes + 1)
        
        result( (numberOfDotes + 1)  * ( k - 1) + m ) = x(k) + (m - 1) * (x(k + 1) - x(k)) / (numberOfDotes + 1);
    
    end
    
end

% setting right border for result
result(size) = x(length(x));

end