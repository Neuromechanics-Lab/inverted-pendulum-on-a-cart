% Autumn Routt
% 09/10/2025
% Inverted pendulum on a cart with added mass
% Torque at pivot point modelled with SRM-style controller as:
% T = kp*theta + kv*theta_dot + ka*theta_ddot
% And common time delay, lambda
% Sweeping gain values and plotting influence on outcome variables

clc; clear; 
close all;

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

% kp = 0; % angle gain
% kv = 0; % angular velocity gain
ka = 0; % angular acceleration gain (infinitely expanding sin wave when ka>=59)
delay = 0; % common time delay (ms), must be <2s and must be an integer

simTime = 2; % how much time is simulated (seconds)
timestep = 0.001;
pertDuration = 10; % number of time steps cart takes to accelerate and decelerate
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
dec_end = cart_dec_time+pertDuration;

% Looping through potential gain values to make surf and contour plots
settlingTimes = []; muscGrossWork = []; muscNetWork = []; muscImpulse = [];
kp_array = 0:10:1000; kv_array = 0:10:1000;
kp_iter = 0;
for kp = kp_array
    kp_iter = kp_iter + 1;
    kv_iter = 0;

    for kv = kv_array
        kv_iter = kv_iter + 1;

        % use the forward Euler method to find solution with the time delay
        x_sim = zeros(2000,2); % x_sim = [angle, angular velocity]
        t_sim = zeros(2000,1);
        ang_acc = zeros(2000,1);
        cartTrq = []; gravTrq = []; muscTrq = [];
        settlingTime = simTime;
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
            if settlingTime==simTime && (iter-2000)>dec_end && abs(x_sim(iter))<deg2rad(90)
                recentAngVelVals = x_sim(iter-(0:500),2);
                angVelNearZero = abs(recentAngVelVals) < 0.01;
                if sum(angVelNearZero)==501
                    settlingTime = (iter-2500)*timestep;
                end 
            end 
        end 
        x_sim = x_sim(2001:size(x_sim,1),:);
        t_sim = t_sim(2001:size(t_sim,1),:);
        ang_acc = ang_acc(2001:size(ang_acc,1),:);

        x_simInfo = stepinfo(x_sim(:,1),t_sim);
        settlingTimes(kp_iter,kv_iter) = settlingTime;
        % settlingTimes(kp_iter,kv_iter) = x_simInfo.SettlingTime;

        muscPower = muscTrq.*x_sim(:,2);
        muscGrossWork(kp_iter,kv_iter) = trapz(abs(muscPower));
        muscNetWork(kp_iter,kv_iter) = trapz(muscPower);
        muscImpulse(kp_iter,kv_iter) = trapz(muscTrq);
    end
end 

%%
% plot gains sweep with outcome variables

fig=figure;
surf(kp_array,kv_array,settlingTimes')
xlabel('kp')
ylabel('kv')
zlabel('Settling Time (s)')
figName = ['Output/SurfSettlingTime_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.fig'];
savefig(figName)
pngName = ['Output/SurfSettlingTime_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.png'];
saveas(fig,pngName)

fig=figure;
surf(kp_array,kv_array,muscGrossWork')
xlabel('kp')
ylabel('kv')
zlabel('Muscle Gross Work')
figName = ['Output/SurfGrossWork_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.fig'];
savefig(figName)
pngName = ['Output/SurfGrossWork_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.png'];
saveas(fig,pngName)

fig=figure;
surf(kp_array,kv_array,muscNetWork')
xlabel('kp')
ylabel('kv')
zlabel('Muscle Net Work')
figName = ['Output/SurfNetWork_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.fig'];
savefig(figName)
pngName = ['Output/SurfNetWork_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.png'];
saveas(fig,pngName)

fig=figure;
surf(kp_array,kv_array,muscImpulse')
xlabel('kp')
ylabel('kv')
zlabel('Muscle Impulse')
figName = ['Output/SurfImpulse_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.fig'];
savefig(figName)
pngName = ['Output/SurfImpulse_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.png'];
saveas(fig,pngName)

fig=figure;
contour(kp_array,kv_array,settlingTimes')
xlabel('kp')
ylabel('kv')
title('Settling Time (s)')
colorbar;
figName = ['Output/ContourSettlingTime_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.fig'];
savefig(figName)
pngName = ['Output/ContourSettlingTime_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.png'];
saveas(fig,pngName)

fig=figure;
contour(kp_array,kv_array,muscGrossWork')
xlabel('kp')
ylabel('kv')
title('Muscle Gross Work')
colorbar;
figName = ['Output/ContourGrossWork_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.fig'];
savefig(figName)
pngName = ['Output/ContourGrossWork_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.png'];
saveas(fig,pngName)

fig=figure;
contour(kp_array,kv_array,muscNetWork')
xlabel('kp')
ylabel('kv')
title('Muscle Net Work')
colorbar;
figName = ['Output/ContourNetWork_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.fig'];
savefig(figName)
pngName = ['Output/ContourNetWork_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.png'];
saveas(fig,pngName)

fig=figure;
contour(kp_array,kv_array,muscImpulse')
xlabel('kp')
ylabel('kv')
title('Muscle Impulse')
colorbar;
figName = ['Output/ContourImpulse_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.fig'];
savefig(figName)
pngName = ['Output/ContourImpulse_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.png'];
saveas(fig,pngName)


% figure % plotting surf plots on top of each other
% hold on;
% s1 = surf(kp_array,kv_array,settlingTimes');
% s2 = surf(kp_array,kv_array,muscPowerInt');
% set(s1, 'FaceAlpha', 0.6, 'FaceColor', [0, 0, 1]);
% set(s1, 'EdgeColor', 'none'); 
% set(s2, 'FaceAlpha', 0.6, 'FaceColor', [1, 0, 0]); 
% set(s2, 'EdgeColor', 'none'); 
% xlabel('kp')
% ylabel('kv')
% zlabel('Combined Settling Time (s) and Area Under Power Curve')
% figName = ['Output/SurfCombined_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.fig'];
% savefig(figName)
% hold off;