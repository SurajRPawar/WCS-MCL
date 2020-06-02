%--------------------------------------------------------------------------
% Simulate Voice Coil model
%--------------------------------------------------------------------------

close all; clear; clc;

Lvclv = 2.4e-3;
Rvclv = 3;
rvc = 14.5;
m = 0.7;
Rplv = 0.01;

parameters = [Lvclv; Rvclv; rvc; m; Rplv];

x0 = [0; 0; 0];

t = [0 : 0.001 : 1];
steps = numel(t);
F = 0;
num_states = numel(x0);
x = zeros(steps,num_states); x(1,:) = x0.';
u = zeros(1,steps);

for i = 2:steps
    u(i) = 10*sin(2*pi*2*t(i-1));
    inputs = [u(i); F];
    xp = x((i-1),:).';
    tp = t(i-1);
    tc = t(i);
    [tout,xout] = ode45(@(t,x)model_vc(tp, xp,parameters, inputs), [tp tc],xp);
    xc = xout(end,:);
    x(i,:) = xc;
end

figure;
subplot(3,1,1);
plot(t,x(:,1));
title('Current (A)');

subplot(3,1,2);
plot(t,x(:,2));
title('Velocity (m/s)');

subplot(3,1,3);
plot(t,100*x(:,3));
title('Position (cm)');