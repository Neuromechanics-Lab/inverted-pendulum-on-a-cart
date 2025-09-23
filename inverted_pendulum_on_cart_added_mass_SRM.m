% Autumn Routt
% 08/27/2025
% Inverted pendulum on a cart with added mass
% Torque at pivot point modelled with SRM-style controller as:
% T = kp*theta + kv*theta_dot + ka*theta_ddot
% And common time delay, lambda

clc; clear; 
% close all;

% Defining variables
M = 73.5; % body mass (kg) (73.5kg = 50th percentile for women in US)
m = 0; % added mass (kg) (9kg = CDC recommended weight gain for 30 weeks pregant w normal starting BMI)
h = 1.704; % overall height (m) (1.612m = 50th percentile for women in US)
l = 0.543*h; % body COM height (m) (avg COM height in women is 0.543*overall height)
x_a = 0.87; % added mass height (m)
y_a = 0.15; % added mass horizontal offset from pendulum arm (m)

theta_a = atan((m*y_a)/(M*l+m*x_a));
l_lumped = sqrt(((M*l+m*x_a)/(M+m))^2+((m*y_a)/(M+m))^2);
I_lumped = M*l^2+m*(x_a^2+y_a^2);

kp = 680; % angle gain
kv = 10; % angular velocity gain
ka = 0; % angular acceleration gain
delay = 0; % common time delay (ms), must be <2s and must be an integer

simTime = 2; % how much time is simulated (seconds)
timestep = 0.001;
pertDuration = 10; % number of timesteps cart takes to accelerate and decelerate
cart_acc_time = 500; % number of time steps before cart begins accelerating
cart_dec_time = 1000; % number of time steps before cart begins decelerating

% Defining cart acceleration profile 
temp_t = 0:timestep:simTime;
temp_acc = zeros(size(temp_t));
temp_acc((0:pertDuration)+cart_acc_time) = ...
    -cos((0:pertDuration)*2*pi/pertDuration)+1; % acceleration
temp_acc((0:pertDuration)+cart_dec_time) = ...
    cos((0:pertDuration)*2*pi/pertDuration)-1; % deceleration 
cart_acc_spline = spline(temp_t,temp_acc*50);


% use the forward Euler method to find solution with the time delay
x_sim = zeros(2000,2); % x_sim = [angle, angular velocity]
t_sim = zeros(2000,1);
ang_acc = zeros(2000,1);
cartTrq = []; gravTrq = []; muscTrq = [];
for iter = 2000:2000+size(temp_t,2)
    [dX,cart_trq,gravity_trq,musc_trq] = dPendulumStatesAndTrqs(t_sim, x_sim, ang_acc, cart_acc_spline, M, m, l_lumped, theta_a, I_lumped, kp, kv, ka, iter, delay);
    new_x1 = x_sim(iter,1)+timestep*dX(1,:);
    if new_x1>=deg2rad(90)
        new_x2 = 0;
    else 
        new_x2 = x_sim(iter,2)+timestep*dX(2,:);
    end 
    x_sim = [x_sim;new_x1,new_x2];
    t_sim = [t_sim;(iter-2000)*timestep];
    ang_acc = [ang_acc;dX(2,:)];
    cartTrq = [cartTrq;cart_trq];
    gravTrq = [gravTrq;gravity_trq];
    muscTrq = [muscTrq;musc_trq];
end 
x_sim = x_sim(2001:size(x_sim,1),:);
t_sim = t_sim(2001:size(t_sim,1),:);
ang_acc = ang_acc(2001:size(ang_acc,1),:);


% plot result (ang, ang vel, ang acc, trqs)
% figure
numplots = 22;
subplot(numplots,1,1:2)
plot(temp_t, temp_acc);
ylabel('cart acceleration')
titleString = ['k_p=',num2str(kp),', k_v=',num2str(kv),', k_a=',num2str(ka),...
    ', delay=',num2str(delay),', M=',num2str(M),', l=',num2str(l),...
    ', m=',num2str(m),', x_a=',num2str(x_a),', y_a=',num2str(y_a)];
% title(titleString)

subplot(numplots,1,4:7)
plot(t_sim, x_sim(:,1));
ylabel('angle (rad)')
hold on

subplot(numplots,1,9:12)
plot(t_sim, x_sim(:,2))
% xlabel('time (s)')
ylabel('angular velocity (rad/s)')
hold on

subplot(numplots,1,14:17)
plot(t_sim,ang_acc(:,1))
% xlabel('time (s)')
ylabel('angular acceleration (rad/s^2)')
hold on

subplot(numplots,1,19:22)
% plot(t_sim,gravTrq)
% hold on
% plot(t_sim,cartTrq)
plot(t_sim,muscTrq)
xlabel('time (s)')
ylabel('torque (N*m)')
% legend('Gravity','Cart','Muscles')
hold on

fileString = ['kp',num2str(kp),'_kv',num2str(kv),'_ka',num2str(ka),'_delay'...
    ,num2str(delay),'_M',num2str(M),'_l',num2str(l),'_m',num2str(m)...
    ,'_xa',num2str(x_a),'_ya',num2str(y_a)];
figName = ['Output/ModelOutput_',fileString,'.fig'];
savefig(figName)