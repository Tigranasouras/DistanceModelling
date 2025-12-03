function [s_vals, v_vals, a_vals] = eval_spline_nat(x, y, coeff, t_query)
% eval_spline_nat
% Evaluate natural cubic spline and its first two derivatives.
%
% INPUT:
%   x, y   : original data points (same used in splinecoeff_natural)
%   coeff  : (n-1 x 3) matrix [b_i, c_i, d_i]
%   t_query: vector of times to evaluate
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

        % choose interval i so that t in [x(i), x(i+1)]
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

        dx = t - x(i);          % local coordinate in this interval
        b = coeff(i,1);
        c = coeff(i,2);
        d = coeff(i,3);

        % nested form like Program 3.6:
        % S_i(x) = y_i + dx*(b + dx*(c + dx*d))
        s = ((d*dx + c)*dx + b)*dx + y(i);

        % derivatives:
        % S'_i(x) = b + 2*c*dx + 3*d*dx^2
        % S''_i(x) = 2*c + 6*d*dx
        v = b + dx*(2*c + 3*d*dx);
        a = 2*c + 6*d*dx;

        s_vals(k) = s;
        v_vals(k) = v;
        a_vals(k) = a;
    end

    % reshape back to same shape as t_query input
    s_vals = reshape(s_vals, size(t_query));
    v_vals = reshape(v_vals, size(t_query));
    a_vals = reshape(a_vals, size(t_query));
end


