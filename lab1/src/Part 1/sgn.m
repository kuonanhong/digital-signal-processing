function y = sgn(x)
% Signum function:
% For each element of x, sgn(x) returns 1 if the element
% is greater than or equal to zero, and -1 if it is
% less than zero.
y = (x>=0) + (-1)*(x<0);