% race_interp_driver.m
% MTH/CSC 4150 - Distance runner interpolation project
% Uses Newton interpolating polynomial (newtondd.m) and natural cubic spline
% plus simple finite differences for velocity and acceleration.

clear; clc; close all;

%% --- Data: EDIT THESE VECTORS TO YOUR ACTUAL RACE SPLITS ---

% Example dummy data (times in seconds, distances in meters)
t_i = [   0   60  120  180  240  300 ];   % split times
s_i = [   0 1200 2400 3600 4800 5000 ];   % cumulative distances

% Make sure row vectors
t_i = t_i(:).'; 
s_i = s_i(:).';

%% --- Set up evaluation grid ---

t_min  = t_i(1);
t_max  = t_i(end);
dt     = 0.5;                       % time step for plots & diffs (seconds)
t_fine = t_min:dt:t_max;            % fine time grid

%% --- Polynomial interpolation (Newton form via newtondd.m) ---

% newtondd returns p(t_fine) and the coefficient vector c
[s_poly, c_poly] = newtondd(t_fine, t_i, s_i);

%% --- Natural cubic spline interpolation ---

s_spline = spline(t_i, s_i, t_fine);

%% --- Numerical differentiation (finite differences from §5.1) ---

% First and second derivatives for both models
[v_poly,  a_poly]  = diff_central(t_fine, s_poly);
[v_spl,   a_spl]   = diff_central(t_fine, s_spline);

%% --- Plot 1: Distance vs time ---

figure(1); clf; hold on;
colors = lines;

plot(t_fine, s_poly,  'color', colors(1,:), 'linewidth', 3);
plot(t_fine, s_spline,'--',   'color', colors(2,:), 'linewidth', 3);
plot(t_i,    s_i, 'o', 'color', colors(3,:), 'linewidth', 2, 'markersize', 8);

set(gca, 'fontsize', 14);
xlabel('time (s)');
ylabel('distance (m)');
title('Position vs. Time: Polynomial vs. Natural Cubic Spline');
legend('Newton interpolating polynomial', 'Natural cubic spline', 'split data', ...
       'location', 'northwest');
grid on;
hold off;

%% --- Plot 2: Velocity vs time ---

figure(2); clf; hold on;

plot(t_fine, v_poly,  'color', colors(1,:), 'linewidth', 3);
plot(t_fine, v_spl,   '--',    'color', colors(2,:), 'linewidth', 3);

set(gca, 'fontsize', 14);
xlabel('time (s)');
ylabel('velocity (m/s)');
title('Velocity vs. Time from Interpolants');
legend('poly-based v(t)', 'spline-based v(t)', 'location', 'best');
grid on;
hold off;

%% --- Plot 3: Acceleration vs time ---

figure(3); clf; hold on;

plot(t_fine, a_poly,  'color', colors(1,:), 'linewidth', 3);
plot(t_fine, a_spl,   '--',    'color', colors(2,:), 'linewidth', 3);

set(gca, 'fontsize', 14);
xlabel('time (s)');
ylabel('acceleration (m/s^2)');
title('Acceleration vs. Time from Interpolants');
legend('poly-based a(t)', 'spline-based a(t)', 'location', 'best');
grid on;
hold off;

%% --- Very simple "accuracy" check (optional, qualitative) ---
% Just compare how smooth the spline is vs. the polynomial in the plots.
% You can discuss these differences in the Discussion section of your report.

