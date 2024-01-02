
%{
MATLAB Codes for explicit solution of Saint-Venant equations
© Hamed Valae 
Dec 2011
Advanced Hydraulics Programming Project
%}


% Clear memory, command window, and close all figures
clear all;
clc;
close all;

% Input parameters
L = input('Enter the Length of the Channel or River (ft): ');
b = input('Enter the Width of the Channel or River (ft): ');
s0 = input('Enter the Slope (s0): ');
y0 = input('Enter the Normal Depth (y0): ');
n = input('Enter Manning''s roughness coefficient (n): ');
time = input('Enter the Total Computation Time (min): ');
N = input('Enter the Number of Sections: ');
time = time * 60; % Convert time to seconds

g = 32.172; % Acceleration due to gravity (ft/s^2)
A0 = b * y0; % Initial Area
R0 = b * y0 / (b + 2 * y0); % Initial hydraulic radius
Q0 = 1.486 / n * sqrt(s0) * R0^(2/3) * b * y0; % Initial discharge
V0 = Q0 / A0; % Initial Velocity
a = sqrt(g * y0) + V0;
T = L / a;
Dx = L / (N - 1); % Delta x

Dt = floor((Dx / a) / 10) * 10 - 20; % Delta t

FT = ceil(time / Dt); % Calculate cycles of loop

% Variable Matrix for Q, Y, R, V
Q = zeros(FT, 1);
Y = zeros(N, FT);
V = zeros(N, FT);
R = zeros(N, FT);
courant = zeros(FT, 1);

% Calculate initial Courant number
courant(1) = Dx / a;

X = linspace(0, L, N);

% Initialize Y, R, Q, V at time = 0
Y(:, 1) = y0;
R(:, 1) = R0;
Q(1) = Q0;
V(:, 1) = V0;

% Loop for times of calculations
for j = 2:FT
    t = t + Dt;

    % Check and calculate Q for any time
    if t <= 1200
        Q(j) = Q0 + (Q0 / 1200) * t;
    elseif t > 1200 && t <= 1800
        Q(j) = -Q0 / 400 * (t - 1200) + 2 * Q0;
    else
        Q(j) = Q0 / 2;
    end

    % Calculate Y, V, R at upstream condition
    ys = Y(2, j - 1);
    A = b * ys;
    P = (b + 2 * ys);
    Vs = V(2, j - 1);
    
    % Define the Saint Venant equation and solve it
    M = solve(Q(j) / (b * y) - Vs - sqrt(g / ys) * (y - ys) + g * ((Q(j))^2 * n^2 / (A^2 * (A / P)^(4/3)) - s0) * Dt);

    if double(M(1, 1)) > 0
        Y(1, j) = M(1, 1);
    else
        Y(1, j) = M(2, 1);
    end

    R(1, j) = (b * Y(1, j) / (b + 2 * Y(1, j))); % Calculate R according to Y
    V(1, j) = Q(j) / (b * Y(1, j)); % Calculate V according to Y

    % Loop for calculating data at points between upstream and downstream conditions
    for i = 2:N - 1
        Y(i, j) = Y(i, j - 1) - (Dt / (2 * Dx)) * (V(i, j - 1) * (Y(i + 1, j - 1) - Y(i - 1, j - 1)) + Y(i, j - 1) * (V(i + 1, j - 1) - V(i - 1, j - 1)));
        R(i, j) = (b * Y(i, j) / (b + 2 * Y(i, j)));
        sf = V(i, j - 1)^2 * n^2 / (R(i, j)^(4 / 3));
        V(i, j) = V(i, j - 1) - V(i, j - 1) * ((V(i + 1, j - 1) - V(i - 1, j - 1)) * Dt / (2 * Dx)) - g * Dt * (Y(i + 1, j - 1) - Y(i - 1, j - 1)) / (2 * Dx) + g * Dt * (s0 - sf);
    end

    % Calculate Y, V, R at downstream condition
    Y(N, j) = (V(N - 1, j) * Y(N - 1, j) * b / 264)^(2 / 3) + 2.32;
    V(N, j) = (264 * (Y(N, j) - 2.32)^(1.5)) / (b * Y(N, j));
    R(N, j) = (b * Y(N, j) / (b + 2 * Y(N, j)));

    % Calculate Courant number for any time
    courant(j) = max(Dx / (sqrt(g * Y(:, j)) + V(:, j)));
end

% Plotting
xmax = 18;
figure;
for i = 1:3:20
    plot(X(:), Y(:, i), '-b', 'LineWidth', 1);
    title('X-Y Graph for Special Times of Computations');
    xlabel('X = 0 ~ 12000 (ft)');
    ylabel('Y (ft)');
    box on;
    axis([0 L 0 xmax]);
    hold on;
    grid on;
end

ff = floor(FT / 2);
ee = floor(FT / 2 + 30);
figure;
for i = ff:3:ee
    plot(X(:), Y(:, i), ':k', 'LineWidth', 1);
    title('X-Y Graph for Special Times of Computations');
    xlabel('X = 0 ~ 12000 (ft)');
    ylabel('Y (ft)');
    box on;
    axis([0 L 0 xmax]);
    hold on;
end

ff = floor(FT - 30);
ee = floor(FT);
figure;
for i = ff:3:ee
    plot(X(:), Y(:, i), '--r', 'LineWidth', 1);
    title('X-Y Graph for Special Times of Computations');
    xlabel('X = 0 ~ 12000 (ft)');
    ylabel('Y (ft)');
    box on;
    axis([0 L 0 xmax]);
    hold on;
end

% Display Results
disp('|/////////////////////////////////////////////////////////////////|');
disp('|////  Results of Explicit Solution of Saint-Venant Equations /////|');
disp('|////////////////////////© Hamed Valae ///////////////////////////|');
disp('   ');
disp('X = ');
disp(X);
disp('Y = ');
disp(Y);
disp('V = ');
disp(V);
disp('R = ');
disp(R);
