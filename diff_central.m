function [v, a] = diff_central(t, s)
% diff_central
% Computes first and second derivatives of s(t) using finite differences.
% Input:
%   t : 1D vector of time values (assumed equally spaced)
%   s : 1D vector of function values s(t)
% Output:
%   v : approximate first derivative (velocity)
%   a : approximate second derivative (acceleration)

dt = t(2) - t(1);
N  = length(s);

v = zeros(size(s));
a = zeros(size(s));

% --- First derivative v(t) ---

% forward diff at first point
v(1) = (s(2) - s(1)) / dt;

% backward diff at last point
v(N) = (s(N) - s(N-1)) / dt;

% centered diff in the interior
v(2:N-1) = (s(3:N) - s(1:N-2)) / (2*dt);

% --- Second derivative a(t) ---

% use standard 3-point stencil
a(2:N-1) = (s(3:N) - 2*s(2:N-1) + s(1:N-2)) / dt^2;

% simple one-sided estimates at the endpoints
a(1) = a(2);
a(N) = a(N-1);

end
