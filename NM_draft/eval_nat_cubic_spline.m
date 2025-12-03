function [s_vals, v_vals, a_vals] = eval_nat_cubic_spline(x, coeffs, t_query)
% eval_nat_cubic_spline
% Evaluates S(t), S'(t), S''(t) for the spline built by
% build_natural_cubic_spline.
%
% INPUT:
%   x       : knots x_0..x_n (same as used in builder)
%   coeffs  : n x 4 matrix [a_i, b_i, c_i, d_i]
%   t_query : vector of times to evaluate
%
% OUTPUT:
%   s_vals : S(t)
%   v_vals : S'(t)
%   a_vals : S''(t)

    t_query = t_query(:);
    n = length(x) - 1;

    s_vals = zeros(size(t_query));
    v_vals = zeros(size(t_query));
    a_vals = zeros(size(t_query));

    for k = 1:length(t_query)
        t = t_query(k);

        % pick the right interval [x_i, x_{i+1}]
        if t <= x(1)
            i = 1;
        elseif t >= x(end)
            i = n;
        else
            i = find(x <= t, 1, 'last');
            if i == length(x)
                i = i-1;
            end
        end

        dt = t - x(i);
        a_i = coeffs(i,1);
        b_i = coeffs(i,2);
        c_i = coeffs(i,3);
        d_i = coeffs(i,4);

        s_vals(k) = a_i + b_i*dt + c_i*dt^2 + d_i*dt^3;
        v_vals(k) =        b_i   + 2*c_i*dt + 3*d_i*dt^2;
        a_vals(k) =                 2*c_i    + 6*d_i*dt;
    end

    % reshape back
    s_vals = reshape(s_vals, size(t_query));
    v_vals = reshape(v_vals, size(t_query));
    a_vals = reshape(a_vals, size(t_query));
end
