% Autumn Routt
% 08/27/2025
% Inverted pendulum on a cart with added mass
% Torque at pivot point modelled as:
% T = kp*theta + kv*theta_dot + ka*theta_ddot
% And common time delay, lambda
% Optimizing gain values to minimize ending angle

clc; clear; close all;

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

%%
prev_best_ending_angle = deg2rad(90);
kp_best = 0; kv_best = 0; ka_best = 0; delay_best = 0;
% kp_alt = []; kv_alt = []; ka_alt = []; delay_alt = [];
best_vals = [];
for kp = 600:20:700 % angle gain
    for kv = 1100:20:1200 % angular velocity gain
        % for ka = 0:0.25:1 % angular acceleration gain
            ka = 0;
            % current = [kp,kv,ka]
            % for delay = 0:50:500 % common time delay (ms), must be <2s and must be an integer
                delay = 0;
                % use the forward Euler method to find solution with the time delay
                x_sim = zeros(2000,2); % x_sim = [angle, angular velocity]
                t_sim = zeros(2000,1);
                ang_acc = zeros(2000,1);
                for iter = 2000:2000+size(temp_t,2)
                    dX = dPendulumStates(t_sim, x_sim, ang_acc, cart_acc_spline, M, m, l_lumped, theta_a, I_lumped, kp, kv, ka, iter, delay);
                    new_x1 = x_sim(iter,1)+timestep*dX(1,:);
                    if new_x1>=deg2rad(90)
                        new_x2 = 0;
                    else 
                        new_x2 = x_sim(iter,2)+timestep*dX(2,:);
                    end 
                    x_sim = [x_sim;new_x1,new_x2];
                    t_sim = [t_sim;(iter-2000)*timestep];
                    ang_acc = [ang_acc;dX(2,:)];
                end 
                x_sim = x_sim(2001:size(x_sim,1),:);
                t_sim = t_sim(2001:size(t_sim,1),:);
                ang_acc = ang_acc(2001:size(ang_acc,1),:);

                if abs(x_sim(end,1))<=prev_best_ending_angle
                    prev_best_ending_angle = abs(x_sim(end,1));
                    kp_best = kp;
                    kv_best = kv;
                    ka_best = ka;
                    delay_best = delay;
                    best_vals = [best_vals;kp,kv,ka,delay,x_sim(end,1)];
                %     kp_alt = [kp_alt,kp];
                %     kv_alt = [kv_alt,kv];
                %     ka_alt = [ka_alt,ka];
                %     delay_alt = [delay_alt,delay];
                % elseif abs(x_sim(end,1))==prev_best_ending_angle
                %     kp_alt = [kp_alt,kp];
                %     kv_alt = [kv_alt,kv];
                %     ka_alt = [ka_alt,ka];
                %     delay_alt = [delay_alt,delay];
                end 
            % end
        % end
    end
end

%%
for vals = 1:size(best_vals,1)
    kp_best = best_vals(vals,1); kv_best = best_vals(vals,2);
    ka_best = best_vals(vals,3); delay_best = best_vals(vals,4);
    % use the forward Euler method to find solution with the time delay
    x_sim = zeros(2000,2); % x_sim = [angle, angular velocity]
    t_sim = zeros(2000,1);
    ang_acc = zeros(2000,1);
    for iter = 2000:2000+size(temp_t,2)
        dX = dPendulumStates(t_sim, x_sim, ang_acc, cart_acc_spline, M, m, l_lumped, theta_a, I_lumped, kp_best, kv_best, ka_best, timestep, iter, delay_best);
        new_x1 = x_sim(iter,1)+timestep*dX(1,:);
        new_x2 = x_sim(iter,2)+timestep*dX(2,:);
        x_sim = [x_sim;new_x1,new_x2];
        t_sim = [t_sim;(iter-2000)*timestep];
        ang_acc = [ang_acc;dX(2,:)];
    end 
    x_sim = x_sim(2001:size(x_sim,1),:);
    t_sim = t_sim(2001:size(t_sim,1),:);
    ang_acc = ang_acc(2001:size(ang_acc,1),:);
    
    % plot result 
    subplot(15,1,1)
    plot(temp_t, temp_acc);
    ylabel('cart acceleration')
    hold on
    subplot(15,1,3:6)
    plot(t_sim, x_sim(:,1));
    ylabel('angle (rad)')
    hold on
    subplot(15,1,8:11)
    plot(t_sim, x_sim(:,2))
    xlabel('time (s)')
    ylabel('angular velocity (rad/s)')
    hold on
    subplot(15,1,12:15)
    plot(t_sim,ang_acc(:,1))
    xlabel('time (s)')
    ylabel('angular acceleration (rad/s^2)')
    linkaxes(get(gcf, 'Children'),'x')
    hold on
end 
