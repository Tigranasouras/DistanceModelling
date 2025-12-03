% Program 3.5 (Sauer) – Calculation of spline coefficients
% Natural cubic spline: S''(x_1) = S''(x_n) = 0
%
% INPUT:
%   x, y : data vectors (x(1) < ... < x(n))
% OUTPUT:
%   coeff : (n-1 x 3) matrix of [b_i, c_i, d_i] for i = 1..n-1
%           S_i(x) = y_i + b_i*dx + c_i*dx^2 + d_i*dx^3,  dx = x - x_i

function coeff = splinecoeff_natural(x, y)

    n = length(x);
    v1 = 0; vn = 0;  %#ok<NASGU>  % not used for natural, but kept for book similarity

    A = zeros(n,n);
    r = zeros(n,1);
    dx = zeros(n-1,1);
    dy = zeros(n-1,1);

    % define the deltas
    for i = 1:n-1
        dx(i) = x(i+1) - x(i);
        dy(i) = y(i+1) - y(i);
    end

    % load interior rows of A and RHS r
    for i = 2:n-1
        A(i, i-1:i+1) = [dx(i-1), 2*(dx(i-1)+dx(i)), dx(i)];
        r(i) = 3*(dy(i)/dx(i) - dy(i-1)/dx(i-1));
    end

    % --- Natural spline endpoint conditions ---
    A(1,1) = 1;      % c_1 = 0
    A(n,n) = 1;      % c_n = 0

    coeff = zeros(n,3);      % columns: b, c, d
    coeff(:,2) = A \ r;      % solve for c_i (second-derivative related)

    % solve for b_i and d_i for i = 1..n-1
    for i = 1:n-1
        coeff(i,3) = (coeff(i+1,2) - coeff(i,2)) / (3*dx(i));          % d_i
        coeff(i,1) = dy(i)/dx(i) - dx(i)*(2*coeff(i,2)+coeff(i+1,2))/3; % b_i
    end

    coeff = coeff(1:n-1,1:3);   % only intervals 1..n-1
end
