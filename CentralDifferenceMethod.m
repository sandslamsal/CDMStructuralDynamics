clc
clear

% modal parameters
wn = 6.283; % eigen-frequency (rad/s)
eta = 0.05; % damping ratio
wd = wn.*sqrt(1-eta.^2); % damped eigen frequency (rad/s)
M = 0.2533; % mass
K = wn^2.*M; % stiffness
C = 2*eta.*M.*wn; % damping

% time domain
t = 0:0.05:10;  % time from 0 to 1 with a step of 0.1
dt = 0.05;  % time step

%t = linspace(0,300,2000); % time
%dt = median(diff(t)); % time step

F0 = 10; % amplitude of force
w = wd; % pulsation of the harmonic force
F = F0.*sin(pi()*t/0.6); % expression of the harmonic force

% Initial calculations
a0 = (F(1) - C*0 - K*0)/M;  % initial acceleration
u_1 = 0 - dt*0 + (dt^2)/2*a0;  % displacement at time t = -dt
k_hat = M/(dt^2) + C/(2*dt);
a = M/(dt^2) - C/(2*dt);
b = K - 2*M/(dt^2);

% initialize
u_minus1 = 0;  % displacement at time t = -dt
u = zeros(size(t));
P_hat = zeros(size(t));
u_next = zeros(size(t));
u_theoretical = zeros(size(t));  % initialize theoretical u

% Calculations for each time step
for i = 1:length(t)-1
    P_hat(i) = F(i) - a*u_minus1 - b*u(i);
    u_next(i+1) = P_hat(i) / k_hat;
    
    % update for next iteration
    u_minus1 = u(i);
    u(i+1) = u_next(i+1);
end



% tabulate the output
T = table(t(1:end-1)', F(1:end-1)', u(1:end-1)', u(1:end-1)', P_hat(1:end-1)', u_next(2:end)',  ...
    'VariableNames', {'ti', 'Pi', 'Ui_minus1', 'Ui', 'Phat_i', 'Ui_plus1'});
disp(T);


plot(t,u);
hold on
plot(t,u_theoretical)
title('Central Difference Method Displacement');
ylabel('Displacement (m)');
xlabel('Time (s)');
hold off

