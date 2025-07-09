alpha = 1.2; % similarity parameter for porous medium wall effects
% i.e. du/dy = alpha / sqrt(k) * u_wall - Q/A
k = 1e-10; % m^2
h = 1e-4; % m
mu = 100; % Pa.s
Q = 1e-12 / h; % volume flow rate per unit depth
dpdx = -Q * mu / k; % darcy's law

% channel flow (above the porous material)
y = linspace(0, h, 100); % discretised y coordinates
s = h / sqrt(k); % dimensionless height
u_wall = - k / (2 * mu) * (s^2 + 2 * alpha * s) / (1 + alpha * s) * dpdx;
u = u_wall * (1 + y .* alpha / sqrt(k))...
    + 1/(2 * mu) * (y.^2 + 2 * alpha * y * sqrt(k)) * dpdx;
max_u = max(u);

% composite flow (within the porous material)
u_darcy = -k / mu * dpdx;

% BL approximation
slope = (u(2)-u(1)) / (y(2)-y(1));
du = u(1) - u_darcy;
dy = 2 * du / slope;
coeff = du * dy ^ -2;
BL_y = linspace(-dy, 0, 100);
BL_u = coeff * (BL_y + dy) .^ 2 + u_darcy;

figure
hold on
grid on

% plot background
yline(h, color='black', LineWidth=2)
yline(0, color='black', LineWidth=2)
yline(-h, color='black', LineWidth=2)
rectangle('Position', [0, -h, max_u*1.1, h], 'FaceColor', 'black', 'FaceAlpha', 0.3);

% plot parabolic flow
plot(u, y, color='black', LineWidth=2)

% plot darcy flow
plot([u_darcy u_darcy], [-h min(BL_y)], color='black', LineWidth=2)

% plot BL
plot(BL_u, BL_y, color='black', LineWidth=2, LineStyle="--")

% formatting
xlim([0, max_u*1.1]);

