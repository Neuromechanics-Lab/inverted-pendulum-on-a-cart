% Autumn Routt
% 08/27/2025
% Inverted pendulum on a cart with added mass
% Torque at pivot point modelled as:
% T = kp*theta + kv*theta_dot + ka*theta_ddot
% And common time delay, lambda
% Optimizing gain values and time delay to stabilize pendulum quickly 

clc; clear; close all;

% Defining variables
M = 77; % body mass (kg) (77kg = avg for women in US)
m = 0; % added mass (kg) (9kg = CDC recommended weight gain for 30 weeks pregant w normal starting BMI)
l = 0.87; % body COM height (m) (0.87m = avg for women in US)
x_a = 1; % added mass height (m)
y_a = 0.5; % added mass horizontal offset from pendulum arm (m)
% k = 20; % torsional spring constant, need to add as a param in func if you want to uncomment this

alpha = atan((m*y_a)/(M*l+m*x_a));
l_lumped = sqrt(((M*l+m*x_a)/(M+m))^2+((m*y_a)/(M+m))^2);
I_lumped = M*l^2+m*(x_a^2+y_a^2);

% Define acceleration profile. 
timestep = 0.001;
temp_t = 0:timestep:5;
temp_acc = zeros(size(temp_t));
perturb_num = 10;
temp_acc( (0:perturb_num)+500) = ...
    -cos((0:perturb_num)*2*pi/perturb_num)+1; % acceleration
temp_acc( (0:perturb_num)+1000) = ...
    cos((0:perturb_num)*2*pi/perturb_num)-1; % deceleration 
cart_acc_spline = spline(temp_t,temp_acc*50);

%%
prev_best_ending_angle = deg2rad(90);
kp_best = 0; kv_best = 0; ka_best = 0; delay_best = 0;
% kp_alt = []; kv_alt = []; ka_alt = []; delay_alt = [];
best_vals = [];
for kp = 0:250:1000 % angle gain
    current=kp
    for kv = 0:250:1000 % angular velocity gain
        for ka = 0:250:1000 % angular acceleration gain
            for delay = 0:50:500 % common time delay (ms), must be <2s and must be an integer
                
                % use the forward Euler method to find solution with the time delay
                x_sim = zeros(2000,2); % x_sim = [angle, angular velocity]
                t_sim = zeros(2000,1);
                ang_acc = zeros(2000,1);
                for iter = 2000:2000+size(temp_t,2)
                    dX = dPendulumStates(t_sim, x_sim, ang_acc, cart_acc_spline, M, m, l_lumped, alpha, I_lumped, kp, kv, ka, timestep, iter, delay);
                    new_x1 = x_sim(iter,1)+timestep*dX(1,:);
                    new_x2 = x_sim(iter,2)+timestep*dX(2,:);
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
            end
        end
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
        dX = dPendulumStates(t_sim, x_sim, ang_acc, cart_acc_spline, M, m, l_lumped, alpha, I_lumped, kp_best, kv_best, ka_best, timestep, iter, delay_best);
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
legend('1','2','3','4','5','6','7','8','9','10')