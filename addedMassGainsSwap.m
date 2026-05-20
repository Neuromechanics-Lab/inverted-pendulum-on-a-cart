% Autumn Routt
% 10/23/2025
% Generating figures to show the outcome variables 
% for different added mass conditions at their optimal gain values
% and at each other's optimal gain values.

clc; clear;
close all;

% addpath('C:\Users\aroutt3\OneDrive - Georgia Institute of Technology\Modelling & Simulation\Old Inverted Pendulum Outputs\Dec 2025')
addpath('/Users/autumnroutt/Library/CloudStorage/OneDrive-GeorgiaInstituteofTechnology/Modelling & Simulation/inverted-pendulum-on-a-cart-MEW06724D/Output')
str = '_NP01pert_muscTrq';
% load(['addedMassSweepVars',str,'.mat'])
load(['fullOutput',str,'.mat'])
m_array = 0:10:50;

output = 1;

ka = 0; % angular acceleration gain
delay = 0; % common time delay (ms), must be <2s and must be an integer

y_a = 0.15; % added mass horizontal offset from pendulum arm (m), >0 is front-loaded, <0 is back-loaded
M = 73.5; % body mass (kg) (73.5kg = 50th percentile for women in US)
h = 1.612; % overall height (m) (1.612m = 50th percentile for women in US)
l = 0.543*h; % body COM height (m) (avg COM height in women is 0.543*overall height)
x_a = l; % added mass height (m)
leaningStart = 1; % binary variable to decide if the pendulum will start at a lean to compensate for added mass (1) or at theta zero (0)

% simTime = 2; % how much time is simulated (seconds)
% timestep = 0.001;
% pertDuration = 150; % number of timesteps cart takes to accelerate and decelerate
% cart_acc_time = 5; % number of time steps before cart begins accelerating
% cart_dec_time = 255; % number of time steps before cart begins decelerating
% pertMag = 2.5; % max magnitude of cart acceleration (m/s^2)
% pertDir = -1; % 1 = cart moves right, -1 = cart moves left

peakCOMangDisp_noAdapt = zeros(1,length(m_array));
peakTrqs_noAdapt = zeros(1,length(m_array));
OVsAtMin_noAdapt = zeros(3,length(m_array));

x_sim_noAdapt = zeros(2*length(m_array),1000*simTime+2);

for i = 1:length(m_array)
    m=m_array(i);

    kp = minGains(1,1); % angle gain
    kv = minGains(2,1); % angular velocity gain

    [settlingTime, muscPower, muscGrossWork, muscNetWork, muscImpulse, muscAvgTrq, muscTrq, x_sim, t_sim]...
    = inverted_pendulum_on_cart_added_mass_SRM_func(output,kp,kv,ka,delay,m,y_a,M,h,l,x_a,leaningStart);

    peakCOMangDisp_noAdapt(1,i) = rad2deg(max(abs(x_sim(:,1)-x_sim(1,1))));
    peakTrqs_noAdapt(1,i) = max(abs(muscTrq));

    muscSumTrqSqrd = sum(muscTrq.^2);
    % thetaSumSqrd = sum(x_sim(:,1).^2);
    adjustedTheta = x_sim(:,1)-x_sim(1,1);
    thetaSumSqrd = sum(adjustedTheta.^2);
    thetaDotSumSqrd = sum(x_sim(:,2).^2);

    x_sim_noAdapt(i,:) = rad2deg(x_sim(:,1)');
    x_sim_noAdapt(i+length(m_array),:) = rad2deg(x_sim(:,2)');

    % C1 = 1; % Torque sum
    % C2 = 10^7; % Theta 
    % C3 = 0.5*10^7; % Theta dot 
    % C1 = 1; % Torque sum
    % C2 = 10^5; % Theta 
    % C3 = 10^5; % Theta dot 
    C1 = 1; % Torque sum
    C2 = 9*10^5; % Theta 
    C3 = 0.9*10^5; % Theta dot 
    % C1 = 1; % Torque sum
    % C2 = 10^7; % Theta 
    % C3 = 10^5; % Theta dot 

    OV1 = C1*muscSumTrqSqrd;
    OV2 = C2*thetaSumSqrd;
    OV3 = C3*thetaDotSumSqrd;

    OVsAtMin_noAdapt(1,i) = OV1;
    OVsAtMin_noAdapt(2,i) = OV2;
    OVsAtMin_noAdapt(3,i) = OV3;
end 

legend('0kg - Unadapted','10kg - Unadapted','20kg - Unadapted','30kg - Unadapted','40kg - Unadapted','50kg - Unadapted','60kg - Unadapted','70kg - Unadapted')
savefig(['UnadaptedTraces',str])
saveas(gcf,['UnadaptedTraces',str,'.png'])

x_sim_Adapt = zeros(2*length(m_array),1000*simTime+2);

figure;
for i = 1:length(m_array)
    m = m_array(i);

    kp = minGains(1,i); % angle gain
    kv = minGains(2,i); % angular velocity gain
    
    [settlingTime, muscPower, muscGrossWork, muscNetWork, muscImpulse, muscAvgTrq, muscTrq, x_sim, t_sim]...
    = inverted_pendulum_on_cart_added_mass_SRM_func(output,kp,kv,ka,delay,m,y_a,M,h,l,x_a,leaningStart,...
    simTime,timestep,pertDuration,cart_acc_time,cart_dec_time,pertMag,pertDir);

    x_sim_Adapt(i,:) = rad2deg(x_sim(:,1)');
    x_sim_Adapt(i+length(m_array),:) = rad2deg(x_sim(:,2)');
end 

legend('0kg - Adapted','10kg - Adapted','20kg - Adapted','30kg - Adapted','40kg - Adapted','50kg - Adapted','60kg - Adapted','70kg - Adapted')
savefig(['AdaptedTraces',str])
saveas(gcf,['AdaptedTraces',str,'.png'])

figure;
plot(m_array,peakCOMangDispAtMin(:,(1:length(m_array))));
hold on
plot(m_array,peakCOMangDisp_noAdapt);
xlabel('Added Mass (kg)')
ylabel('Peak COM Angle Disp (deg)')
legend('Adapted','Unadapted')
savefig(['AngleGainsComp',str])
saveas(gcf,['AngleGainsComp',str,'.png'])

figure;
plot(m_array,peakTrqsAtMin(:,(1:length(m_array)))./M);
hold on
plot(m_array,peakTrqs_noAdapt./M)
xlabel('Added Mass (kg)')
ylabel('Peak Torque (N*m/kg)')
legend('Adapted','Unadapted')
savefig(['TorqueGainsComp',str])
saveas(gcf,['TorqueGainsComp',str,'.png'])

figure;
plot(m_array,OVsAtMin(1,(1:length(m_array))));
hold on
plot(m_array,OVsAtMin_noAdapt(1,:));
xlabel('Added Mass (kg)')
ylabel('J_{energy}')
legend('Adapted','Unadapted')
savefig(['EnergyGainsComp',str])
saveas(gcf,['EnergyGainsComp',str,'.png'])

% figure;
% plot(m_array,OVsAtMin(2,:)+OVsAtMin(3,:));
% hold on
% plot(m_array,OVsAtMin_noAdapt(2,:)+OVsAtMin_noAdapt(3,:));
% xlabel('Added Mass (kg)')
% ylabel('J_{instability}')
% legend('Adapted','Unadapted')
% savefig(['InstabilityGainsComp',str])
% saveas(gcf,['InstabilityGainsComp',str,'.png'])

% Averaging instability terms, new CF
figure;
plot(m_array,(OVsAtMin(2,(1:length(m_array)))+OVsAtMin(3,(1:length(m_array))))./2);
hold on
plot(m_array,(OVsAtMin_noAdapt(2,(1:length(m_array)))+OVsAtMin_noAdapt(3,(1:length(m_array))))./2);
xlabel('Added Mass (kg)')
ylabel('J_{instability}')
legend('Adapted','Unadapted')
savefig(['InstabilityGainsComp',str])
saveas(gcf,['InstabilityGainsComp',str,'.png'])

figure;
plot(m_array,OVsAtMin(1,(1:length(m_array))));
hold on
plot(m_array,OVsAtMin_noAdapt(1,:));
plot(m_array,(OVsAtMin(2,(1:length(m_array)))+OVsAtMin(3,(1:length(m_array))))./2);
plot(m_array,(OVsAtMin_noAdapt(2,(1:length(m_array)))+OVsAtMin_noAdapt(3,(1:length(m_array))))./2);
plot(m_array, minCFvals(:,(1:length(m_array))))
plot(m_array,OVsAtMin_noAdapt(1,(1:length(m_array)))+(OVsAtMin_noAdapt(2,(1:length(m_array)))+OVsAtMin_noAdapt(3,(1:length(m_array))))./2)
xlabel('Added Mass (kg)')
ylabel('J')
legend('Adapted J_{energy}','Unadapted J_{energy}','Adapted J_{instability}','Unadapted J_{instability}','Adapted J','Unadapted J')
savefig(['EnergyInstabilityGainsComp',str])
saveas(gcf,['EnergyInstabilityGainsComp',str,'.png'])

figure;
plot(m_array,OVsAtMin(1,(1:length(m_array)))./minCFvals(:,1));
hold on
plot(m_array,OVsAtMin_noAdapt(1,:)./minCFvals(:,1));
plot(m_array,((OVsAtMin(2,(1:length(m_array)))+OVsAtMin(3,(1:length(m_array))))./2)./minCFvals(:,1));
plot(m_array,((OVsAtMin_noAdapt(2,(1:length(m_array)))+OVsAtMin_noAdapt(3,(1:length(m_array))))./2)./minCFvals(:,1));
plot(m_array, minCFvals(:,(1:length(m_array)))./minCFvals(:,1))
plot(m_array,(OVsAtMin_noAdapt(1,(1:length(m_array)))+(OVsAtMin_noAdapt(2,(1:length(m_array)))+OVsAtMin_noAdapt(3,(1:length(m_array))))./2)./minCFvals(:,1))
xlabel('Added Mass (kg)')
ylabel('J')
legend('Adapted J_{energy}','Unadapted J_{energy}','Adapted J_{instability}','Unadapted J_{instability}','Adapted J','Unadapted J')
savefig(['NormalizedEnergyInstabilityGainsComp',str])
saveas(gcf,['NormalizedEnergyInstabilityGainsComp',str,'.png'])

% figure;
% plot(m_array,OVsAtMin(1,:)-OVsAtMin(1,1));
% hold on
% plot(m_array,OVsAtMin_noAdapt(1,:)-OVsAtMin_noAdapt(1,1));
% xlabel('Added Mass (kg)')
% ylabel('Delta J_{energy}')
% legend('Adapted','Unadapted')
% savefig(['EnergyGainsComp',str])
% saveas(gcf,['EnergyGainsComp',str,'.png'])
% 
% figure;
% plot(m_array,(OVsAtMin(2,:)+OVsAtMin(3,:))-(OVsAtMin(2,1)+OVsAtMin(3,1)));
% hold on
% plot(m_array,(OVsAtMin_noAdapt(2,:)+OVsAtMin_noAdapt(3,:))-(OVsAtMin_noAdapt(2,1)+OVsAtMin_noAdapt(3,1)));
% xlabel('Added Mass (kg)')
% ylabel('Delta J_{instability}')
% legend('Adapted','Unadapted')
% savefig(['InstabilityGainsComp',str])
% saveas(gcf,['InstabilityGainsComp',str,'.png'])
% 
% figure;
% plot(m_array,(OVsAtMin(2,:)+OVsAtMin(3,:))-(OVsAtMin(2,1)+OVsAtMin(3,1)));
% hold on
% plot(m_array,(OVsAtMin_noAdapt(2,:)+OVsAtMin_noAdapt(3,:))-(OVsAtMin_noAdapt(2,1)+OVsAtMin_noAdapt(3,1)));
% plot(m_array,OVsAtMin(1,:)-OVsAtMin(1,1));
% plot(m_array,OVsAtMin_noAdapt(1,:)-OVsAtMin_noAdapt(1,1));
% xlabel('Added Mass (kg)')
% ylabel('Delta J')
% legend('Adapted Instability','Unadapted Instability', 'Adapted Energy', 'Unadapted Energy')
% savefig(['InstabilityGainsComp',str])
% saveas(gcf,['InstabilityGainsComp',str,'.png'])


% figure;
% plot(m_array,(OVsAtMin(2,:)+OVsAtMin(3,:))-(OVsAtMin_noAdapt(2,:)+OVsAtMin_noAdapt(3,:)));
% hold on;
% plot(m_array,OVsAtMin(1,:)-OVsAtMin_noAdapt(1,:));
% xlabel('Added Mass (kg)')
% ylabel('Delta J (Unadapted-Adapted)')
% legend('Instability','Energy')
% savefig(['DeltaJGainsComp',str])
% saveas(gcf,['DeltaJGainsComp',str,'.png'])

% Averaging instability terms, new CF
figure;
plot(m_array,(OVsAtMin(2,(1:length(m_array)))+OVsAtMin(3,(1:length(m_array))))./2-(OVsAtMin_noAdapt(2,(1:length(m_array)))+OVsAtMin_noAdapt(3,(1:length(m_array))))./2);
hold on;
plot(m_array,OVsAtMin(1,(1:length(m_array)))-OVsAtMin_noAdapt(1,(1:length(m_array))));
xlabel('Added Mass (kg)')
ylabel('Delta J (Adapted-Unadapted)')
legend('Instability','Energy')
savefig(['DeltaJGainsComp',str])
saveas(gcf,['DeltaJGainsComp',str,'.png'])

% Averaging instability terms, new CF
figure;
plot(m_array,((OVsAtMin_noAdapt(2,(1:length(m_array)))+OVsAtMin_noAdapt(3,(1:length(m_array))))./2-(OVsAtMin(2,(1:length(m_array)))+OVsAtMin(3,(1:length(m_array))))./2)./abs(((OVsAtMin_noAdapt(2,(1:length(m_array)))+OVsAtMin_noAdapt(3,(1:length(m_array))))./2))*100);
hold on;
plot(m_array,(OVsAtMin_noAdapt(1,(1:length(m_array)))-OVsAtMin(1,(1:length(m_array))))./abs(OVsAtMin_noAdapt(1,(1:length(m_array))))*100);
xlabel('Added Mass (kg)')
ylabel('Delta J ((Unadapted-Adapted)/Unadapted)')
legend('Instability','Energy')
savefig(['NormDeltaJGainsComp',str])
saveas(gcf,['NormDeltaJGainsComp',str,'.png'])

% Percent of J for energy and instability
minCFvals_noAdapt = OVsAtMin_noAdapt(1,(1:length(m_array))) + (OVsAtMin_noAdapt(2,(1:length(m_array)))+OVsAtMin_noAdapt(3,(1:length(m_array))))./2;
percInstabilityAdapt = 100*((OVsAtMin(2,(1:length(m_array)))+OVsAtMin(3,(1:length(m_array))))./2)./minCFvals(1,1:length(m_array));
percInstabilityUnadapt = 100*((OVsAtMin_noAdapt(2,(1:length(m_array)))+OVsAtMin_noAdapt(3,(1:length(m_array))))./2)./minCFvals_noAdapt(1,1:length(m_array));
percEnergyAdapt = 100*OVsAtMin(1,(1:length(m_array)))./minCFvals(1,1:length(m_array));
percEnergyUnadapt = 100*OVsAtMin_noAdapt(1,(1:length(m_array)))./minCFvals_noAdapt(1,1:length(m_array));
figure;
plot(m_array,percInstabilityAdapt);
hold on;
plot(m_array,percInstabilityUnadapt);
plot(m_array,percEnergyAdapt);
plot(m_array,percEnergyUnadapt);
xlabel('Added Mass (kg)')
ylabel('% of J')
legend('Instability % Adapted','Instability % Unadapted','Energy % Adapted','Energy % Unadapted')
savefig(['PercentGainsComp',str])
saveas(gcf,['PercentGainsComp',str,'.png'])

% Change in percent of J for energy and instability
figure;
plot(m_array,percInstabilityAdapt-percInstabilityUnadapt);
hold on;
plot(m_array,percEnergyAdapt-percEnergyUnadapt);
xlabel('Added Mass (kg)')
ylabel('Delta % of J')
legend('Instability','Energy')
savefig(['DeltaPercentGainsComp',str])
saveas(gcf,['DeltaPercentGainsComp',str,'.png'])

% figure;
% for i = 1:length(m_array)
%     plot(x_sim_noAdapt(i+length(m_array),:),x_sim_noAdapt(i,:));
%     hold on
% end 
% xlabel('Theta Dot')
% ylabel('Theta')
% legend('0kg - Unadapted','10kg - Unadapted','20kg - Unadapted','30kg - Unadapted','40kg - Unadapted','50kg - Unadapted','60kg - Unadapted','70kg - Unadapted')
% savefig(['ThetaVsThetaDotUnadapted',str])
% saveas(gcf,['ThetaVsThetaDotUnadapted',str,'.png'])
% 
% figure;
% for i = 1:length(m_array)
%     plot(x_sim_Adapt(i+length(m_array),:),x_sim_Adapt(i,:));
%     hold on
% end 
% xlabel('Theta Dot')
% ylabel('Theta')
% legend('0kg - Adapted','10kg - Adapted','20kg - Adapted','30kg - Adapted','40kg - Adapted','50kg - Adapted','60kg - Adapted','70kg - Adapted')
% savefig(['ThetaVsThetaDotAdapted',str])
% saveas(gcf,['ThetaVsThetaDotAdapted',str,'.png'])

figure;
plot(m_array, minGains(1,:))
hold on;
plot(m_array, minGains(2,:))
xlabel('Added Mass (kg)')
ylabel('Gain Values')
legend('Kp','Kv')
savefig(['GainVals',str])
saveas(gcf,['GainVals',str,'.png'])