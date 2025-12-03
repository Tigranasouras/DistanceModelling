function [s_vals, v_vals, a_vals] = eval_spline_nat(x, y, coeff, t_query)
% eval_spline_nat
% Evaluate a cubic spline (built from splinecoeff_natural or splinecoeff_clamped)
% and its first two derivatives at a set of query times.
%
% INPUT:
%   x       : knot locations (same x_i used to build coeff),  size n x 1 or 1 x n
%   y       : data values at knots (y_i), same size as x
%   coeff   : (n-1 x 3) matrix [b_i, c_i, d_i] for each interval [x_i, x_{i+1}]
%             so that on interval i,
%                  S_i(t) = y_i + b_i*dx + c_i*dx^2 + d_i*dx^3,
%             where dx = t - x_i
%   t_query : vector of times t where we want to evaluate the spline
%
% OUTPUT:
%   s_vals  : S(t_query)     – position values
%   v_vals  : S'(t_query)    – velocity values
%   a_vals  : S''(t_query)   – acceleration values

    % t_query to a column vector.
    % Simplifies indexing with t_query(k).
    t_query = t_query(:);

    % Number of spline intervals is n-1 if we have n knots.
    n = length(x) - 1;

    % Preallocate output vectors with same size as t_query
    s_vals = zeros(size(t_query));
    v_vals = zeros(size(t_query));
    a_vals = zeros(size(t_query));

    % Loop over each query point t_query(k)
    for k = 1:length(t_query)
        t = t_query(k);   % current time at which we evaluate S, S', S''

        % This is to Choose which interval [x(i), x(i+1)] contains t
        % If t lies to the left of the first knot, clamp to first interval.
        if t <= x(1)
            i = 1;
        % If t lies to the right of the last knot, clamp to last interval.
        elseif t >= x(end)
            i = n;
        else
            % Otherwise, find the last x(i) such that x(i) <= t.
            % This gives us the correct interval index i.
            i = find(x <= t, 1, 'last');

            % Safety check: if we accidentally land on the last knot,
            % move back one so that we still have [x(i), x(i+1)].
            if i == length(x)
                i = i - 1;
            end
        end

        % Local coordinate on interval i: shift by left endpoint x(i)
        dx = t - x(i);

        % Extract coefficients for this interval:
        %   S_i(t) = y(i) + b*dx + c*dx^2 + d*dx^3
        b = coeff(i,1);
        c = coeff(i,2);
        d = coeff(i,3);

        % Evaluate S(t) using nested multiplication (like Program 3.6)
        % S_i(t) = y_i + dx*(b + dx*(c + dx*d))
        s = ((d*dx + c)*dx + b)*dx + y(i);

        % Evaluate S'(t) and S''(t) analytically
        % From S_i(t) = y_i + b*dx + c*dx^2 + d*dx^3,
        %   S'_i(t)  = b + 2*c*dx + 3*d*dx^2
        %   S''_i(t) = 2*c + 6*d*dx
        v = b + dx*(2*c + 3*d*dx);
        a = 2*c + 6*d*dx;

        % Store results at the k-th query point
        s_vals(k) = s;
        v_vals(k) = v;
        a_vals(k) = a;
    end

    % Reshape outputs to match the original shape of t_query input.
    % (If t_query was a row vector when passed in, outputs will be row vectors too.)
    s_vals = reshape(s_vals, size(t_query));
    v_vals = reshape(v_vals, size(t_query));
    a_vals = reshape(a_vals, size(t_query));
end
