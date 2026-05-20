% OLD AND OUTDATED, USE inverted_pendulum_on_cart_added_mass_SRM_func
% Autumn Routt
% 10/23/2025
% Inverted pendulum on a cart with added mass
% Torque at pivot point modelled with SRM-style controller as:
% T = kp*theta + kv*theta_dot + ka*theta_ddot
% And common time delay, lambda
% Sweeping gain values and plotting influence on outcome variables
% Tried to implement parfor on the outer loop, but haven't gotten it to work yet

clc; clear; 
close all;

% Defining variables
M = 73.5; % body mass (kg) (73.5kg = 50th percentile for women in US)
% m = 0; % added mass (kg) (9kg = CDC recommended weight gain for 30 weeks pregant w normal starting BMI)
h = 1.612; % overall height (m) (1.612m = 50th percentile for women in US)
l = 0.543*h; % body COM height (m) (avg COM height in women is 0.543*overall height)
x_a = l; % added mass height (m)
y_a = 0.15; % added mass horizontal offset from pendulum arm (m)
leaningStart = 1; % binary variable to decide if the pendulum will start at a lean to compensate for added mass (1) or at theta zero (0)

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
kp_array = 500:100:9000;
kv_array = 500:100:9000;
m_array = 0:10:70;

minGains = zeros(2, length(m_array));
minCFvals = zeros(size(m_array));
OVsAtMin = zeros(3, length(m_array));
peakTrqs = zeros(size(m_array));

parfor i = 1:8
    m = m_array(i);
    theta_a = atan((m*y_a)/(M*l+m*x_a));
    l_lumped = sqrt(((M*l+m*x_a)/(M+m))^2+((m*y_a)/(M+m))^2);
    I_lumped = M*l^2+m*(x_a^2+y_a^2);
    kp_iter = 0;
    progress = 0;
    loading = waitbar(0,'Sweeping gain values...');

    settlingTimes = zeros(length(kp_array),length(kv_array));
    muscGrossWork = zeros(length(kp_array),length(kv_array));
    muscNetWork = zeros(length(kp_array),length(kv_array));
    muscImpulse = zeros(length(kp_array),length(kv_array));
    muscAvgTrq = zeros(length(kp_array),length(kv_array));
    muscSumTrqSqrd = zeros(length(kp_array),length(kv_array));
    thetaSumSqrd = zeros(length(kp_array),length(kv_array));
    thetaDotSumSqrd = zeros(length(kp_array),length(kv_array));

    kp_n = length(kp_array);
    kv_n = length(kv_array);

    % for kp = kp_array
    %     kp_iter = kp_iter + 1;
    %     kv_iter = 0;
    for j = 1:kp_n
        kp_iter = j;
        kp = kp_array(kp_iter);

        % for kv = kv_array
        %     kv_iter = kv_iter + 1;
        for k = 1:kv_n
            kv_iter = k;
            kv = kv_array(kv_iter);
            numVals = length(kp_array)*length(kv_array);
            progress = progress + 1;
            waitbar(progress/numVals,loading)
            % use the forward Euler method to find solution with the time delay
            % x_sim = zeros(2000,2); % x_sim = [angle, angular velocity]
            if leaningStart==1
                x_sim(1:2000,1) = -1*theta_a;
            end 
            % t_sim = zeros(2000,1);
            % ang_acc = zeros(2000,1);
            % cartTrq = []; gravTrq = []; muscTrq = [];
            settlingTime = simTime;

            x_sim = zeros(2000+size(temp_t,2),2);
            t_sim = zeros(2000+size(temp_t,2),1);
            ang_acc = zeros(2000+size(temp_t,2),1);
            cartTrq = zeros(size(temp_t,2),1);
            gravTrq = zeros(size(temp_t,2),1);
            muscTrq = zeros(size(temp_t,2),1);

            for iter = 2000:2000+size(temp_t,2)
                [dX,cart_trq,gravity_trq,musc_trq] = dPendulumStatesAndTrqs(t_sim, x_sim, ang_acc, cart_acc_spline, M, m, l_lumped, theta_a, I_lumped, kp, kv, ka, iter, delay);
                new_x1 = x_sim(iter,1)+timestep*dX(1,:);
                if new_x1>=deg2rad(90)
                    new_x2 = 0;
                else 
                    new_x2 = x_sim(iter,2)+timestep*dX(2,:);
                end 
                % x_sim = [x_sim;new_x1,new_x2];
                % t_sim = [t_sim;(iter-2000)*timestep];
                % ang_acc = [ang_acc;dX(2,:)];
                % cartTrq = [cartTrq;cart_trq];
                % gravTrq = [gravTrq;gravity_trq];
                % muscTrq = [muscTrq;musc_trq];

                x_sim(iter,:) = [new_x1, new_x2];
                t_sim(iter,1) = (iter-2000)*timestep;
                ang_acc(iter,1) = dX(2,:);
                cartTrq(iter-1999,1) = cart_trq;
                gravTrq(iter-1999,1) = gravity_trq;
                muscTrq(iter-1999,1) = musc_trq;

                if settlingTime==simTime && (iter-2000)>dec_end && abs(x_sim(iter,1))<deg2rad(90)
                    recentAngVelVals = x_sim(iter-(0:50),2);
                    angVelNearZero = abs(recentAngVelVals) < 0.02;
                    % slopeAngVel = (x_sim(iter,2)-x_sim(iter-50,2))/(50*timestep);
                    recentAngAccVals = ang_acc(iter-(0:50),1);
                    angAccNearZero = abs(recentAngAccVals) < 0.02;
                    % if sum(angVelNearZero)==51 && slopeAngVel < 0.1
                    if sum(angVelNearZero)==51 && sum(angAccNearZero)==51
                        settlingTime = (iter-2050)*timestep;
                    end 
                end 
            end 
            x_sim = x_sim(2000:size(x_sim,1),:);
            t_sim = t_sim(2000:size(t_sim,1),:);
            ang_acc = ang_acc(2000:size(ang_acc,1),:);
            % cartTrq = cartTrq(2:size(cartTrq,1),:);
            % gravTrq = gravTrq(2:size(gravTrq,1),:);
            % muscTrq = muscTrq(2:size(muscTrq,1),:);
    
            % x_simInfo = stepinfo(x_sim(:,1),t_sim);
            settlingTimes(kp_iter,kv_iter) = settlingTime;
            % settlingTimes(kp_iter,kv_iter) = x_simInfo.SettlingTime;
    
            muscPower = muscTrq.*x_sim(:,2);
            muscGrossWork(kp_iter,kv_iter) = trapz(abs(muscPower));
            muscNetWork(kp_iter,kv_iter) = trapz(muscPower);
            muscImpulse(kp_iter,kv_iter) = trapz(muscTrq);
            muscAvgTrq(kp_iter,kv_iter) = mean(abs(muscTrq));
    
            muscSumTrqSqrd(kp_iter,kv_iter) = sum(muscTrq.^2);
            thetaSumSqrd(kp_iter,kv_iter) = sum(x_sim(:,1).^2);
            thetaDotSumSqrd(kp_iter,kv_iter) = sum(x_sim(:,2).^2);
        end
    end 

close(loading);

% plot combined outcome variables using cost functions

% Second attempt at a cost function, making adjustable weights for outcome variables

% Assigning weights
C1 = 1; % Torque sum
C2 = 10^7; % Theta 
C3 = 0.5*10^7; % Theta dot 
% C1 = 1; % Torque sum
% C2 = 1; % Theta
% C3 = 1; % Theta dot 

OV1 = C1*muscSumTrqSqrd;
OV2 = C2*thetaSumSqrd;
OV3 = C3*thetaDotSumSqrd;
CF = C1*muscSumTrqSqrd+C2*thetaSumSqrd+C3*thetaDotSumSqrd;
minCFval = min(CF,[],"all");
[minCFx, minCFy] = find(CF==minCFval);
% minGains = [minGains; kp_array(minCFx), kv_array(minCFy)];
minGains(:,i) = [kp_array(minCFx); kv_array(minCFy)];
% minGains(2,i) = kv_array(minCFy);
% minCFvals = [minCFvals; minCFval];
minCFvals(i) = minCFval;
OVsAtMin(:,i) = [OV1(minCFx,minCFy); OV2(minCFx,minCFy); OV3(minCFx,minCFy)];
% OVsAtMin(2,i) = OV2(minCFx,minCFy);
% OVsAtMin(3,i) = OV3(minCFx,minCFy);
peakTrqs(i) = max(abs(muscTrq));

CFstr = ['CF = ',num2str(C1),'*sum(muscTrq^2)+',num2str(C2),'*sum(theta^2)+',num2str(C3),'*sum(theta dot^2)'];
fig=figure;
surf(kp_array,kv_array,CF')
xlabel('kp')
ylabel('kv')
zlabel(CFstr)
figName = ['Output/SurfWeightedCF_C1',num2str(C1),'_C2',num2str(C2),'_C3',num2str(C3),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.fig'];
savefig(figName)
pngName = ['Output/SurfWeightedCF_C1',num2str(C1),'_C2',num2str(C2),'_C3',num2str(C3),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.png'];
saveas(fig,pngName)

fig=figure;
contour(kp_array,kv_array,CF')
xlabel('kp')
ylabel('kv')
title(CFstr)
colorbar;
figName = ['Output/ContourWeightedCF_C1',num2str(C1),'_C2',num2str(C2),'_C3',num2str(C3),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.fig'];
savefig(figName)
pngName = ['Output/ContourWeightedCF_C1',num2str(C1),'_C2',num2str(C2),'_C3',num2str(C3),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.png'];
saveas(fig,pngName)

OV1str = [num2str(C1),'*sum(muscTrq^2)'];
fig=figure;
surf(kp_array,kv_array,OV1')
xlabel('kp')
ylabel('kv')
zlabel(OV1str)
figName = ['Output/Surf_C1',num2str(C1),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.fig'];
savefig(figName)
pngName = ['Output/Surf_C1',num2str(C1),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.png'];
saveas(fig,pngName)

fig=figure;
contour(kp_array,kv_array,OV1')
xlabel('kp')
ylabel('kv')
title(OV1str)
colorbar;
figName = ['Output/Contour_C1',num2str(C1),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.fig'];
savefig(figName)
pngName = ['Output/Contour_C1',num2str(C1),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.png'];
saveas(fig,pngName)

OV2str = [num2str(C2),'*sum(theta^2)'];
fig=figure;
surf(kp_array,kv_array,OV2')
xlabel('kp')
ylabel('kv')
zlabel(OV2str)
figName = ['Output/Surf_C2',num2str(C2),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.fig'];
savefig(figName)
pngName = ['Output/Surf_C2',num2str(C2),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.png'];
saveas(fig,pngName)

fig=figure;
contour(kp_array,kv_array,OV2')
xlabel('kp')
ylabel('kv')
title(OV2str)
colorbar;
figName = ['Output/Contour_C2',num2str(C2),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.fig'];
savefig(figName)
pngName = ['Output/Contour_C2',num2str(C2),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.png'];
saveas(fig,pngName)

OV3str = [num2str(C3),'*sum(theta dot^2)'];
fig=figure;
surf(kp_array,kv_array,OV3')
xlabel('kp')
ylabel('kv')
zlabel(OV3str)
figName = ['Output/Surf_C3',num2str(C3),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.fig'];
savefig(figName)
pngName = ['Output/Surf_C3',num2str(C3),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.png'];
saveas(fig,pngName)

fig=figure;
contour(kp_array,kv_array,OV3')
xlabel('kp')
ylabel('kv')
title(OV3str)
colorbar;
figName = ['Output/Contour_C3',num2str(C3),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.fig'];
savefig(figName)
pngName = ['Output/Contour_C3',num2str(C3),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.png'];
saveas(fig,pngName)

end 

% %% First attempts at a cost function, tried normalizing outcomes and punishing big gains
% 
% STnorm = settlingTimes./max(settlingTimes,[],'all');
% GWnorm = muscGrossWork./max(muscGrossWork,[],'all');
% NWnorm = abs(muscNetWork)./max(abs(muscNetWork),[],'all');
% Inorm = abs(muscImpulse)./max(abs(muscImpulse),[],'all');
% ATnorm = abs(muscAvgTrq)./max(abs(muscAvgTrq),[],'all');
% kpnorm = kp_array./max(kp_array,[],'all');
% kvnorm = kv_array./max(kv_array,[],'all');
% 
% exp = 2;
% 
% CF1 = (STnorm.^exp+GWnorm.^exp+kpnorm.^exp+kvnorm.^exp)/4;
% CF2 = (STnorm.^exp+NWnorm.^exp+kpnorm.^exp+kvnorm.^exp)/4;
% CF3 = (STnorm.^exp+Inorm.^exp+kpnorm.^exp+kvnorm.^exp)/4;
% CF4 = (STnorm.^exp+ATnorm.^exp+kpnorm.^exp+kvnorm.^exp)/4;
% 
% % CF1 = (STnorm.^exp+GWnorm.^exp)/2;
% % CF2 = (STnorm.^exp+NWnorm.^exp)/2;
% % CF3 = (STnorm.^exp+Inorm.^exp)/2;
% 
% fig=figure;
% surf(kp_array,kv_array,CF1')
% xlabel('kp')
% ylabel('kv')
% zlabel('Settling Time (s) and Gross Work')
% figName = ['Output/SurfSTnormGWnormExp',int2str(exp),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.fig'];
% savefig(figName)
% pngName = ['Output/SurfSTnormGWnormExp',int2str(exp),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.png'];
% saveas(fig,pngName)
% 
% fig=figure;
% contour(kp_array,kv_array,CF1')
% xlabel('kp')
% ylabel('kv')
% title('Settling Time (s) and Gross Work')
% colorbar;
% figName = ['Output/ContourSTnormGWnormExp',int2str(exp),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.fig'];
% savefig(figName)
% pngName = ['Output/ContourSTnormGWnormExp',int2str(exp),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.png'];
% saveas(fig,pngName)
% 
% fig=figure;
% surf(kp_array,kv_array,CF2')
% xlabel('kp')
% ylabel('kv')
% zlabel('Settling Time (s) and Net Work')
% figName = ['Output/SurfSTnormNWnormExp',int2str(exp),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.fig'];
% savefig(figName)
% pngName = ['Output/SurfSTnormNWnormExp',int2str(exp),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.png'];
% saveas(fig,pngName)
% 
% fig=figure;
% contour(kp_array,kv_array,CF2')
% xlabel('kp')
% ylabel('kv')
% title('Settling Time (s) and Net Work')
% colorbar;
% figName = ['Output/ContourSTnormNWnormExp',int2str(exp),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.fig'];
% savefig(figName)
% pngName = ['Output/ContourSTnormNWnormExp',int2str(exp),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.png'];
% saveas(fig,pngName)
% 
% fig=figure;
% surf(kp_array,kv_array,CF3')
% xlabel('kp')
% ylabel('kv')
% zlabel('Settling Time (s) and Impulse')
% figName = ['Output/SurfSTnormInormExp',int2str(exp),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.fig'];
% savefig(figName)
% pngName = ['Output/SurfSTnormInormExp',int2str(exp),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.png'];
% saveas(fig,pngName)
% 
% fig=figure;
% contour(kp_array,kv_array,CF3')
% xlabel('kp')
% ylabel('kv')
% title('Settling Time (s) and Impulse')
% colorbar;
% figName = ['Output/ContourSTnormInormExp',int2str(exp),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.fig'];
% savefig(figName)
% pngName = ['Output/ContourSTnormInormExp',int2str(exp),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.png'];
% saveas(fig,pngName)
% 
% fig=figure;
% surf(kp_array,kv_array,CF4')
% xlabel('kp')
% ylabel('kv')
% zlabel('Settling Time (s) and Avg Torque')
% figName = ['Output/SurfSTnormATnormExp',int2str(exp),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.fig'];
% savefig(figName)
% pngName = ['Output/SurfSTnormATnormExp',int2str(exp),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.png'];
% saveas(fig,pngName)
% 
% fig=figure;
% contour(kp_array,kv_array,CF4')
% xlabel('kp')
% ylabel('kv')
% title('Settling Time (s) and Avg Torque')
% colorbar;
% figName = ['Output/ContourSTnormATnormExp',int2str(exp),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.fig'];
% savefig(figName)
% pngName = ['Output/ContourSTnormATnormExp',int2str(exp),'_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'_m',int2str(m),'.png'];
% saveas(fig,pngName)



%%
% plot gains sweep with outcome variables

% fig=figure;
% surf(kp_array,kv_array,settlingTimes')
% xlabel('kp')
% ylabel('kv')
% zlabel('Settling Time (s)')
% figName = ['Output/SurfSettlingTime_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.fig'];
% savefig(figName)
% pngName = ['Output/SurfSettlingTime_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.png'];
% saveas(fig,pngName)
% 
% fig=figure;
% surf(kp_array,kv_array,muscGrossWork')
% xlabel('kp')
% ylabel('kv')
% zlabel('Muscle Gross Work')
% figName = ['Output/SurfGrossWork_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.fig'];
% savefig(figName)
% pngName = ['Output/SurfGrossWork_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.png'];
% saveas(fig,pngName)
% 
% fig=figure;
% surf(kp_array,kv_array,muscNetWork')
% xlabel('kp')
% ylabel('kv')
% zlabel('Muscle Net Work')
% figName = ['Output/SurfNetWork_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.fig'];
% savefig(figName)
% pngName = ['Output/SurfNetWork_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.png'];
% saveas(fig,pngName)
% 
% fig=figure;
% surf(kp_array,kv_array,muscImpulse')
% xlabel('kp')
% ylabel('kv')
% zlabel('Muscle Impulse')
% figName = ['Output/SurfImpulse_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.fig'];
% savefig(figName)
% pngName = ['Output/SurfImpulse_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.png'];
% saveas(fig,pngName)
% 
% fig=figure;
% contour(kp_array,kv_array,settlingTimes')
% xlabel('kp')
% ylabel('kv')
% title('Settling Time (s)')
% colorbar;
% figName = ['Output/ContourSettlingTime_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.fig'];
% savefig(figName)
% pngName = ['Output/ContourSettlingTime_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.png'];
% saveas(fig,pngName)
% 
% fig=figure;
% contour(kp_array,kv_array,muscGrossWork')
% xlabel('kp')
% ylabel('kv')
% title('Muscle Gross Work')
% colorbar;
% figName = ['Output/ContourGrossWork_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.fig'];
% savefig(figName)
% pngName = ['Output/ContourGrossWork_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.png'];
% saveas(fig,pngName)
% 
% fig=figure;
% contour(kp_array,kv_array,muscNetWork')
% xlabel('kp')
% ylabel('kv')
% title('Muscle Net Work')
% colorbar;
% figName = ['Output/ContourNetWork_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.fig'];
% savefig(figName)
% pngName = ['Output/ContourNetWork_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.png'];
% saveas(fig,pngName)
% 
% fig=figure;
% contour(kp_array,kv_array,muscImpulse')
% xlabel('kp')
% ylabel('kv')
% title('Muscle Impulse')
% colorbar;
% figName = ['Output/ContourImpulse_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.fig'];
% savefig(figName)
% pngName = ['Output/ContourImpulse_kp',int2str(kp_array(1,1)),'to',int2str(kp_array(1,end)),'_kv',int2str(kv_array(1,1)),'to',int2str(kv_array(1,end)),'.png'];
% saveas(fig,pngName)


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