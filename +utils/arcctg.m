function [ res ] = arcctg( x )
%ARCCTG calculates inverse contangent of x

if (x >= 0) 
    res = asin(1/sqrt(1+x^2)); 
else
    res = pi - asin(1/sqrt(1+x^2)); 
end;

end

