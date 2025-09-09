% Autumn Routt
% 08/27/2025
% Inverted pendulum on a cart with added mass
% Torque at pivot point modelled as:
% T = kp*theta + kv*theta_dot + ka*theta_ddot
% And common time delay, lambda

clc; clear; 
close all;

% for y_a=[-0.5 0 0.5] % for when I want to vary one component

% Defining variables
M = 77; % body mass (kg) (77kg = avg for women in US)
m = 0; % added mass (kg) (9kg = CDC recommended weight gain for 30 weeks pregant w normal starting BMI)
l = 0.87; % body COM height (m) (0.87m = avg for women in US)
x_a = 0.87; % added mass height (m)
y_a = 0.15; % added mass horizontal offset from pendulum arm (m)
% k = 20; % torsional spring constant, need to add as a param in func if you want to uncomment this

alpha = atan((m*y_a)/(M*l+m*x_a));
l_lumped = sqrt(((M*l+m*x_a)/(M+m))^2+((m*y_a)/(M+m))^2);
I_lumped = M*l^2+m*(x_a^2+y_a^2);

% kp = 660; % angle gain
% kv = 1200; % angular velocity gain
ka = 0; % angular acceleration gain (infinitely expanding sin wave when ka>=59)
delay = 0; % common time delay (ms), must be <2s and must be an integer

% Define acceleration profile. 
timestep = 0.001;
temp_t = 0:timestep:2;
temp_acc = zeros(size(temp_t));
perturb_num = 10;
temp_acc( (0:perturb_num)+500) = ...
    -cos((0:perturb_num)*2*pi/perturb_num)+1; % acceleration
temp_acc( (0:perturb_num)+1000) = ...
    cos((0:perturb_num)*2*pi/perturb_num)-1; % deceleration 
cart_acc_spline = spline(temp_t,temp_acc*50);

% clear dPendulumStates % trying to reset pastX, which is a persistent variable

% looping through potential gain values to make a 3d plot
settlingTimes = []; muscPowerInt = [];
kp_array = 500:20:1000; kv_array = 1000:100:3000;
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
        for iter = 2000:2000+size(temp_t,2)
            [dX,cart_trq,gravity_trq,musc_trq] = dPendulumStatesAndTrqs(t_sim, x_sim, ang_acc, cart_acc_spline, M, m, l_lumped, alpha, I_lumped, kp, kv, ka, timestep, iter, delay);
            new_x1 = x_sim(iter,1)+timestep*dX(1,:);
            new_x2 = x_sim(iter,2)+timestep*dX(2,:);
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

        x_simInfo = stepinfo(x_sim(:,1),t_sim);
        settlingTimes(kp_iter,kv_iter) = x_simInfo.SettlingTime;

        muscPower = muscTrq.*x_sim(:,2);
        muscPowerInt(kp_iter,kv_iter) = trapz(abs(muscPower));

    end
end 

%% plot gains sweep with optimizing options
figure(1)
surf(kp_array,kv_array,settlingTimes')
xlabel('kp')
ylabel('kv')
zlabel('Settling Time (s)')
figName = ['SettlingTime_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.fig'];
savefig(figName)

figure(2)
surf(kp_array,kv_array,muscPowerInt')
xlabel('kp')
ylabel('kv')
zlabel('Area Under Power Curve')
figName = ['Power_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.fig'];
savefig(figName)

% clear dPendulumStates % trying to reset pastX, which is a persistent variable
% dX = [;];
% for iter = 1:size(t_sim,1)
%     new_dX = dPendulumStates(t_sim(iter,1), x_sim(iter,:), cart_acc_spline, M, m, l_lumped, alpha, I_lumped, kp, kv, ka, timestep);
%     dX = [dX,new_dX];
% end 

%% plot result (ang, ang vel, ang acc, trqs)
figure(1)
numplots = 21;
subplot(numplots,1,1)
plot(temp_t, temp_acc);
ylabel('cart acceleration')
subplot(numplots,1,3:6)
plot(t_sim, x_sim(:,1));
% hold on
% ylim([-0.1 0.4])
ylabel('angle (rad)')
subplot(numplots,1,8:11)
plot(t_sim, x_sim(:,2))
xlabel('time (s)')
ylabel('angular velocity (rad/s)')
linkaxes(get(gcf, 'Children'),'x')
% hold on
subplot(numplots,1,13:16)
plot(t_sim,ang_acc(:,1))
xlabel('time (s)')
ylabel('angular acceleration (rad/s^2)')
% hold on
subplot(numplots,1,18:21)
plot(t_sim,gravTrq)
hold on
plot(t_sim,cartTrq)
plot(t_sim,muscTrq)
xlabel('time (s)')
ylabel('torque (N*m)')
legend('Gravity','Cart','Muscles')

% legend('kp=0,kv=0,ka=0,delay=0,added mass=0','kp=660,kv=1200,ka=0,delay=0,added mass=0')

% legend('0kg added','9kg added','18kg added','Location','northwest')
% legend('xa=0.75m','xa=1m','xa=1.25m','Location','northwest')
% legend('ya=-0.5m','ya=0m','ya=0.5m','Location','northwest')
% title('m=9kg, xa=1m, M=77kg, l=0.87m')
% print(gcf,'AddedMassOffset','-dsvg','-r300');
% exportgraphics(gcf,'MyFigure.svg','ContentType','vector');

% end % end for loop




% %% ODE for inverted pendulum 
% function dX = dPendulumStates(t, X, cart_acc_spline, M, m, l_lumped, alpha, I_lumped, kp, kv, ka, timestep)
% % X = [angle; angular velocity]
% 
% persistent pastX
% % if exist('pastX','var') == 0
% %     pastX = zeros(2,2/timestep);
% % end 
% 
% if abs(X(1))>=deg2rad(90) % making the pendulum stop falling once it hits the cart
%     X(2)=0;
%     gravity_ang_acc = 0;
%     cart_ang_acc = 0;
%     spring_ang_acc = 0;
% else 
%     cart_ang_acc = ((M+m)*ppval(cart_acc_spline,t)*cos(X(1)+alpha)*l_lumped)/I_lumped;                              
%     gravity_ang_acc = ((M+m)*9.81*sin(X(1)+alpha)*l_lumped)/I_lumped;
% 
%     % spring_ang_acc = -k*X(1)/I_lumped; % adding a torsional spring at pivot point
% 
%     % % what if the spring was exactly as stiff as it needs to be to resist gravity?
%     % if X(1)~= 0
%     %     k = (gravity_ang_acc*I_lumped)/X(1);
%     % else
%     %     k = 0;
%     % end 
% 
%     % now what if the ankle had a gain for angle, ang vel, and ang acc?
%     if isempty(pastX)
%         spring_ang_acc = (kv*X(2) + kp*X(1))/I_lumped;
%     else
%         spring_ang_acc = (ka*((X(2)-pastX(2,size(pastX,2)))/timestep) + kv*X(2) + kp*X(1))/I_lumped;
%     end 
% 
% end
% 
% % % and a time delay?
% % if t>delay
% %     length(pastX)
% %     delayedIndex = length(pastX)-delay/timestep
% %     spring_ang_acc = (ka*((pastX(2,delayedIndex)-pastX(2,delayedIndex-1))/timestep)...
% %         + kv*pastX(2,delayedIndex) + kp*pastX(1,delayedIndex))/I_lumped;
% % else
% %     spring_ang_acc = 0;
% % end 
% 
% % spring_ang_acc = 0;
% 
% % dX = [angular velocity; angular acceleration] 
% dX = [X(2); gravity_ang_acc + cart_ang_acc - spring_ang_acc];
% 
% pastX = [pastX,X];
% end