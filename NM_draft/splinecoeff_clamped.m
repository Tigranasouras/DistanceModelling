function coeff = splinecoeff_clamped(x, y, v1, vn)
% splinecoeff_clamped
% Sauer Program 3.5 with *clamped* endpoint conditions:
%   S'(x_1) = v1,  S'(x_n) = vn

    n = length(x);
    A = zeros(n,n);
    r = zeros(n,1);
    dx = zeros(n-1,1);
    dy = zeros(n-1,1);

    % define the deltas
    for i = 1:n-1
        dx(i) = x(i+1) - x(i);
        dy(i) = y(i+1) - y(i);
    end

    % interior rows
    for i = 2:n-1
        A(i, i-1:i+1) = [dx(i-1)  2*(dx(i-1)+dx(i))  dx(i)];
        r(i) = 3*( dy(i)/dx(i) - dy(i-1)/dx(i-1) );
    end

    % CLAMPED endpoint conditions (from Sauer)
    A(1,1:2)   = [2*dx(1)  dx(1)];
    r(1)       = 3*( dy(1)/dx(1) - v1 );

    A(n,n-1:n) = [dx(n-1)  2*dx(n-1)];
    r(n)       = 3*( vn - dy(n-1)/dx(n-1) );

    coeff = zeros(n,3);
    coeff(:,2) = A \ r;                    % c_i

    % solve for b_i and d_i
    for i = 1:n-1
        coeff(i,3) = (coeff(i+1,2) - coeff(i,2)) / (3*dx(i));           % d_i
        coeff(i,1) = dy(i)/dx(i) - dx(i)*(2*coeff(i,2)+coeff(i+1,2))/3; % b_i
    end

    coeff = coeff(1:n-1,1:3);
end
