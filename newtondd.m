function [y, c] = newtondd(x, x_i, y_i)
% Computes coefficients (c) of interpolating polynomial for interpolation 
% points (x_i, y_i) using Newton basis and then evaluates at specified 
% inputs (x) using nested form with base points (nest_base.m)
% Input: x_i and y_i are vectors containing the x and y coordinates
%        of the n+1 interpolation points
%        x is a vector containing the x values to evaluate
% Output: coefficients (c) of interpolating polynomial using Newton basis
%         in nested form (leading coefficient first)
%         y values (y) obtained by evaluating interpolating polynomial at x
%         values
% Dependency: Uses nest_base.m to evaluate interpolating polynomial in
%             nested form with base points (Newton basis)

% F[x_i]
c = y_i;

% Newton divided differences to calculate F[x_0, x_1], F[x_0, x_1, x_2],...
for i = 1:(length(c)-1);
    c(i+1:end) = (c(i+1:end) - c(i:end-1)) ./ (x_i(i+1:end) - x_i(1:end-i));
end

% Performing Newton divided differences with vectorized operations ^^ is the
% most efficient in MATLAB, but it can also be done with a double for loop,
% which might be more intuitive (see alternative code below)
% for i = 1:(length(c)-1)
%     for j = (length(c)-1):-1:i
%         c(j+1) = (c(j+1) - c(j)) / (x_i(j+1) - x_i(j+1-i));
%     end
% end

% Reverse order of coefficient vector (c) to prepare for use with
% nest_base.m
c = c(end:-1:1);

% Evaluate the polynomial we've built at the specified x values using
% nested form (Newton basis)
y = nest_base(c, x, x_i((end-1):-1:1));

end
