% race_interp_driver.m
% MTH/CSC 4150 - Distance runner interpolation project
% Uses Newton interpolating polynomial (newtondd.m), a natural cubic
% spline, and a clamped cubic spline built with Sauer's Program 3.5 style.
% Velocities/accelerations:
%   - polynomial: finite differences (diff_central.m)
%   - splines: analytic derivatives of cubic pieces (eval_spline_nat.m)

clear; clc; close all;

% Units:
%   time in seconds (from 0:00:00 gun time)
%   distance in meters (mile & km markers converted to meters)
%
% Splits used (3M, 5K, 4M, ..., MAR) come from the
%  marathon times of professor

t_i = [ ...
      0, ...          % start
   1761, ...          %  3M  0:29:21
   1828, ...          %  5K  0:30:28
   2352, ...          %  4M  0:39:12
   2934, ...          %  5M  0:48:54
   3602, ...          %  6M  1:00:02
   3726, ...          % 10K  1:02:06
   4183, ...          %  7M  1:09:43
   4787, ...          %  8M  1:19:47
   5406, ...          %  9M  1:30:06
   5600, ...          % 15K  1:33:20
   5995, ...          % 10M  1:39:55
   6577, ...          % 11M  1:49:37
   7196, ...          % 12M  1:59:56
   7442, ...          % 20K  2:04:02
   7865, ...          % HALF ~21.1 km 2:11:05
   8388, ...          % 14M  2:19:48
   9083, ...          % 15M  2:31:23
   9411, ...          % 25K  2:36:51
   9671, ...          % 16M  2:41:11
  10258, ...          % 17M  2:50:58
  10823, ...          % 18M  3:00:23
  11195, ...          % 30K  3:06:35
  11399, ...          % 19M  3:09:59
  11989, ...          % 20M  3:19:49
  12585, ...          % 21M  3:29:45
  13043, ...          % 35K  3:37:23
  13194, ...          % 22M  3:39:54
  13787, ...          % 23M  3:49:47
  14417, ...          % 24M  4:00:17
  14928, ...          % 40K  4:08:48
  15017, ...          % 25M  4:10:17
  15603, ...          % 26M  4:20:03
  15723  ];           % MAR ~42.2 km 4:22:03

s_i = [ ...
      0, ...          % start
   4828, ...          %  3M  ~3*1609
   5000, ...          %  5K
   6437, ...          %  4M
   8047, ...          %  5M
   9656, ...          %  6M
  10000, ...          % 10K
  11265, ...          %  7M
  12875, ...          %  8M
  14484, ...          %  9M
  15000, ...          % 15K
  16093, ...          % 10M
  17703, ...          % 11M
  19312, ...          % 12M
  20000, ...          % 20K
  21082, ...          % HALF ~21.1 km
  22531, ...          % 14M
  24140, ...          % 15M
  25000, ...          % 25K
  25749, ...          % 16M
  27359, ...          % 17M
  28968, ...          % 18M
  30000, ...          % 30K
  30577, ...          % 19M
  32187, ...          % 20M
  33796, ...          % 21M
  35000, ...          % 35K
  35405, ...          % 22M
  37015, ...          % 23M
  38624, ...          % 24M
  40000, ...          % 40K
  40234, ...          % 25M
  41843, ...          % 26M
  42165  ];           % MAR ~42.2 km

% Make sure row vectors
t_i = t_i(:).';
s_i = s_i(:).';

%% Set up evaluation grid

t_min  = t_i(1);
t_max  = t_i(end);
dt     = 0.5;                   % time step for plots & diffs (seconds)
t_fine = t_min:dt:t_max;        % fine time grid

%% Polynomial interpolation (Newton form via newtondd.m)

% newtondd returns p(t_fine) and the coefficient vector c
[s_poly, c_poly] = newtondd(t_fine, t_i, s_i); 

%% Natural cubic spline via Sauer Program 3.5 style

coeff_nat = splinecoeff_natural(t_i, s_i);
[s_nat, v_nat, a_nat] = eval_spline_nat(t_i, s_i, coeff_nat, t_fine);

%% Clamped cubic spline via Sauer Program 3.5 (clamped endpoints)

% Estimate endpoint velocities from one-sided differences
v1 = (s_i(2)       - s_i(1))     / (t_i(2)       - t_i(1));       % left slope
vn = (s_i(end)     - s_i(end-1)) / (t_i(end)     - t_i(end-1));   % right slope

coeff_clamp = splinecoeff_clamped(t_i, s_i, v1, vn);
[s_clamp, v_clamp, a_clamp] = eval_spline_nat(t_i, s_i, coeff_clamp, t_fine);

%% Numerical differentiation (finite differences) for polynomial ONLY

[v_poly, a_poly] = diff_central(t_fine, s_poly);

%% Plot 1: Distance vs time

figure(1); clf; hold on;
colors = lines;

plot(t_fine, s_poly,    'color', colors(1,:), 'linewidth', 3);
plot(t_fine, s_nat,    '--',     'color', colors(2,:), 'linewidth', 3);
plot(t_fine, s_clamp,  ':',      'color', colors(4,:), 'linewidth', 3);
plot(t_i,    s_i, 'o', 'color', colors(3,:), 'linewidth', 2, 'markersize', 8);

set(gca, 'fontsize', 14);
xlabel('time (s)');
ylabel('distance (m)');
title('Position vs. Time: Polynomial vs. Natural vs. Clamped Splines');
legend('Newton interpolating polynomial', ...
       'Natural cubic spline', ...
       'Clamped cubic spline', ...
       'split data', ...
       'location', 'northwest');
grid on;
hold off;

%% Plot 2: Velocity vs time

figure(2); clf; hold on;

plot(t_fine, v_poly,    'color', colors(1,:), 'linewidth', 3);
plot(t_fine, v_nat,    '--',     'color', colors(2,:), 'linewidth', 3);
plot(t_fine, v_clamp,  ':',      'color', colors(4,:), 'linewidth', 3);

set(gca, 'fontsize', 14);
xlabel('time (s)');
ylabel('velocity (m/s)');
title('Velocity vs. Time from Interpolants');
legend('poly-based v(t) (finite diff)', ...
       'natural spline v(t) (analytic)', ...
       'clamped spline v(t) (analytic)', ...
       'location', 'best');
grid on;
hold off;

%% Plot 3: Acceleration vs time

figure(3); clf; hold on;

plot(t_fine, a_poly,    'color', colors(1,:), 'linewidth', 3);
plot(t_fine, a_nat,    '--',     'color', colors(2,:), 'linewidth', 3);
plot(t_fine, a_clamp,  ':',      'color', colors(4,:), 'linewidth', 3);

set(gca, 'fontsize', 14);
xlabel('time (s)');
ylabel('acceleration (m/s^2)');
title('Acceleration vs. Time from Interpolants');
legend('poly-based a(t) (finite diff)', ...
       'natural spline a(t) (analytic)', ...
       'clamped spline a(t) (analytic)', ...
       'location', 'best');
grid on;
hold off;

%% Simple "accuracy" / behavior check (optional)
% Compare how smooth the spline versions are vs. the polynomial in the plots.
% Natural spline forces S''(t_0)=S''(t_n)=0, so its acceleration is ~0 at
% the endpoints. The clamped spline matches the observed endpoint velocities
% v1 and vn, so its acceleration is allowed to be nonzero at the start/finish.
