% Autumn Routt
% 10/23/2025
% Generating figues to show the outcome variables 
% for different added mass conditions at their optimal gain values
% and at each other's optimal gain values.

clc; clear;
close all;

str = '5cmPert';
load(['addedMassSweepVars',str,'.mat'])

y_a = 0.15; % added mass horizontal offset from pendulum arm (m)
M = 73.5; % body mass (kg) (73.5kg = 50th percentile for women in US)
h = 1.612; % overall height (m) (1.612m = 50th percentile for women in US)
l = 0.543*h; % body COM height (m) (avg COM height in women is 0.543*overall height)
x_a = l; % added mass height (m)
leaningStart = 1; % binary variable to decide if the pendulum will start at a lean to compensate for added mass (1) or at theta zero (0)

peakCOMangDisp_noAdapt = zeros(1,length(m_array));
peakTrqs_noAdapt = zeros(1,length(m_array));
OVsAtMin_noAdapt = zeros(3,length(m_array));

x_sim_noAdapt = zeros(2*length(m_array),2002);

for i = 1:length(m_array)
    m=m_array(i);

    [settlingTime, muscPower, muscGrossWork, muscNetWork, muscImpulse, muscAvgTrq, muscTrq, x_sim, t_sim]...
    = inverted_pendulum_on_cart_added_mass_SRM_func(1,minGains(1,1),minGains(2,1),0,0,m,y_a,M,h,l,x_a,leaningStart);

    peakCOMangDisp_noAdapt(1,i) = abs(max(x_sim(:,1))-min(x_sim(1,1)));
    peakTrqs_noAdapt(1,i) = max(abs(muscTrq));

    muscSumTrqSqrd = sum(muscTrq.^2);
    % thetaSumSqrd = sum(x_sim(:,1).^2);
    adjustedTheta = x_sim(:,1)-x_sim(1,1);
    thetaSumSqrd = sum(adjustedTheta.^2);
    thetaDotSumSqrd = sum(x_sim(:,2).^2);

    x_sim_noAdapt(i,:) = rad2deg(x_sim(:,1)');
    x_sim_noAdapt(i+length(m_array),:) = rad2deg(x_sim(:,2)');

    C1 = 1; % Torque sum
    C2 = 10^7; % Theta 
    C3 = 0.5*10^7; % Theta dot 

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

x_sim_Adapt = zeros(2*length(m_array),2002);

figure;
for i = 1:length(m_array)
    m = m_array(i);

    [settlingTime, muscPower, muscGrossWork, muscNetWork, muscImpulse, muscAvgTrq, muscTrq, x_sim, t_sim]...
    = inverted_pendulum_on_cart_added_mass_SRM_func(1,minGains(1,i),minGains(2,i),0,0,m,y_a,M,h,l,x_a,leaningStart);

    x_sim_Adapt(i,:) = rad2deg(x_sim(:,1)');
    x_sim_Adapt(i+length(m_array),:) = rad2deg(x_sim(:,2)');
end 

legend('0kg - Adapted','10kg - Adapted','20kg - Adapted','30kg - Adapted','40kg - Adapted','50kg - Adapted','60kg - Adapted','70kg - Adapted')
savefig(['AdaptedTraces',str])
saveas(gcf,['AdaptedTraces',str,'.png'])

figure;
plot(m_array,rad2deg(peakCOMangDispAtMin));
hold on
plot(m_array,rad2deg(peakCOMangDisp_noAdapt));
xlabel('Added Mass (kg)')
ylabel('Peak COM Angle Disp (deg)')
legend('Adapted','Unadapted')
savefig(['AngleGainsComp',str])
saveas(gcf,['AngleGainsComp',str,'.png'])

figure;
plot(m_array,peakTrqsAtMin./M);
hold on
plot(m_array,peakTrqs_noAdapt./M)
xlabel('Added Mass (kg)')
ylabel('Peak Torque (N*m/kg)')
legend('Adapted','Unadapted')
savefig(['TorqueGainsComp',str])
saveas(gcf,['TorqueGainsComp',str,'.png'])

figure;
plot(m_array,OVsAtMin(1,:));
hold on
plot(m_array,OVsAtMin_noAdapt(1,:));
xlabel('Added Mass (kg)')
ylabel('J_{energy}')
legend('Adapted','Unadapted')
savefig(['EnergyGainsComp',str])
saveas(gcf,['EnergyGainsComp',str,'.png'])

figure;
plot(m_array,OVsAtMin(2,:)+OVsAtMin(3,:));
hold on
plot(m_array,OVsAtMin_noAdapt(2,:)+OVsAtMin_noAdapt(3,:));
xlabel('Added Mass (kg)')
ylabel('J_{instability}')
legend('Adapted','Unadapted')
savefig(['InstabilityGainsComp',str])
saveas(gcf,['InstabilityGainsComp',str,'.png'])

figure;
for i = 1:length(m_array)
    plot(x_sim_noAdapt(i+length(m_array),:),x_sim_noAdapt(i,:));
    hold on
end 
xlabel('Theta Dot')
ylabel('Theta')
legend('0kg - Unadapted','10kg - Unadapted','20kg - Unadapted','30kg - Unadapted','40kg - Unadapted','50kg - Unadapted','60kg - Unadapted','70kg - Unadapted')
savefig(['ThetaVsThetaDotUnadapted',str])
saveas(gcf,['ThetaVsThetaDotUnadapted',str,'.png'])

figure;
for i = 1:length(m_array)
    plot(x_sim_Adapt(i+length(m_array),:),x_sim_Adapt(i,:));
    hold on
end 
xlabel('Theta Dot')
ylabel('Theta')
legend('0kg - Adapted','10kg - Adapted','20kg - Adapted','30kg - Adapted','40kg - Adapted','50kg - Adapted','60kg - Adapted','70kg - Adapted')
savefig(['ThetaVsThetaDotAdapted',str])
saveas(gcf,['ThetaVsThetaDotAdapted',str,'.png'])