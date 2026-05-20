% Autumn Routt
% 11/6/2025
% Evaluating cost function term trends to see if weighting will affect
% presence of energy/stability tradeoff 

clc; clear;
close all;

% str = '_leanGoal_withFoot_5cmPert_withMusc_Backward';
% load(['addedMassSweepVars',str,'.mat'])
% 
% figure;
% plot(m_array,OVsAtMin(1,:))
% hold on
% plot(m_array,OVsAtMin(2,:))
% plot(m_array,OVsAtMin(3,:))

m_array = 0:10:50;
muscSumTrqSqrd = zeros(1,length(m_array));
thetaSumSqrd = zeros(1,length(m_array));
thetaDotSumSqrd = zeros(1,length(m_array));

for i = 1:length(m_array)
    output = 1;

    kp = 550; % angle gain
    kv = 350; % angular velocity gain
    ka = 0; % angular acceleration gain
    delay = 0; % common time delay (ms), must be <2s and must be an integer
    
    m = m_array(i); % added mass (kg) (9kg = CDC recommended weight gain for 30 weeks pregant w normal starting BMI)
    y_a = 0.15; % added mass horizontal offset from pendulum arm (m), >0 is front-loaded, <0 is back-loaded
    M = 73.5; % body mass (kg) (73.5kg = 50th percentile for women in US)
    h = 1.612; % overall height (m) (1.612m = 50th percentile for women in US)
    l = 0.543*h; % body COM height (m) (avg COM height in women is 0.543*overall height)
    x_a = l; % added mass height (m)
    leaningStart = 1; % binary variable to decide if the pendulum will start at a lean to compensate for added mass (1) or at theta zero (0)

    simTime = 2; % how much time is simulated (seconds)
    timestep = 0.001;
    pertDuration = 150; % number of timesteps cart takes to accelerate and decelerate
    cart_acc_time = 5; % number of time steps before cart begins accelerating
    cart_dec_time = 255; % number of time steps before cart begins decelerating
    pertMag = 2.5; % max magnitude of cart acceleration (m/s^2)
    pertDir = 1; % 1 = cart moves right, -1 = cart moves left

    [settlingTime, muscPower, muscGrossWork, muscNetWork, muscImpulse, muscAvgTrq,muscTrq,x_sim,t_sim]...
        = inverted_pendulum_on_cart_added_mass_SRM_func(output,kp,kv,ka,delay,m,y_a,M,h,l,x_a,leaningStart,...
    simTime,timestep,pertDuration,cart_acc_time,cart_dec_time,pertMag,pertDir);
    
    muscSumTrqSqrd(i) = sum(muscTrq.^2);
    adjustedTheta = x_sim(:,1)-x_sim(1,1);
    thetaSumSqrd(i) = sum(adjustedTheta.^2);
    thetaDotSumSqrd(i) = sum(x_sim(:,2).^2);

end 
legend('0kg','10kg','20kg','30kg','40kg','50kg')

% Assigning weights
C1 = 1; % Torque sum
C2 = 10^6; % Theta
C3 = 10^5; % Theta dot
% C1 = 1; % Torque sum
% C2 = 9*10^5; % Theta 
% C3 = 0.9*10^5; % Theta dot 
% C1 = 1; % Torque sum
% C2 = 1.1592e+06; % Theta 
% C3 =  9.1053e+04; % Theta dot 

figure;
plot(m_array,muscSumTrqSqrd.*C1)
hold on
plot(m_array,thetaSumSqrd.*C2)
plot(m_array,thetaDotSumSqrd.*C3)
plot(m_array,(thetaSumSqrd.*C2+thetaDotSumSqrd.*C3)./2)
% plot(m_array,(thetaSumSqrd.*C2+thetaDotSumSqrd.*C3))
ylabel('Cost Function Terms')
xlabel('Added Mass (kg)')
legend('Musc Trq','Theta','Theta Dot', '(Theta+Theta Dot)/2')