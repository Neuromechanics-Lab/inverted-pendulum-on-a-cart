% Autumn Routt
% 08/25/2025
% Inverted pendulum on a cart with added mass

clc; clear; close all;

% for y_a=[-0.5 0 0.5] % for when I want to vary one component

% Defining variables
M = 77; % body mass (kg) (77kg = avg for women in US)
m = 9; % added mass (kg) (9kg = CDC recommended weight gain for 30 weeks pregant w normal starting BMI)
l = 0.87; % body COM height (m) (0.87m = avg for women in US)
x_a = 1; % added mass height (m)
y_a = 0.5; % added mass horizontal offset from pendulum arm (m)
% k = 20; % torsional spring constant, need to add as a param in func if you want to uncomment this

theta_a = atan((m*y_a)/(M*l+m*x_a));
l_lumped = sqrt(((M*l+m*x_a)/(M+m))^2+((m*y_a)/(M+m))^2);
I_lumped = M*l^2+m*(x_a^2+y_a^2);

% Code from here down is modified from Hansol Ryu's PendulumOnCart.m

% Define acceleration profile. 
% I defined as a time series and passed it to ODE as a spline. 
% It accelerates at time 0, decelerates at time 0.5. 
temp_t = 0:0.001:2;
temp_acc = zeros(size(temp_t));
perturb_num = 10;
temp_acc( (0:perturb_num)+1000) = ...
    -cos((0:perturb_num)*2*pi/perturb_num)+1; % acceleration
temp_acc( (0:perturb_num)+1500) = ...
    cos((0:perturb_num)*2*pi/perturb_num)-1; % decelataion 
cart_acc_spline = spline(temp_t,temp_acc*50);

% run ODE solver
odeopt = odeset('maxstep',1e-2);
[t_sim,x_sim] = ode45(@(t,x)dPendulumStates(t, x, cart_acc_spline, M, m, l_lumped, theta_a, I_lumped), ...
    temp_t, ... simulation time range 
    [0, 0], ... initial condition [angle, angular velocity] 
    odeopt);

% plot result 
figure
subplot(5,1,1)
plot(temp_t, temp_acc);
ylabel('cart acceleration')
subplot(5,1,2:5)
plot(t_sim, x_sim(:,1));
% ylim([-0.1 0.4])
xlabel('time (s)')
ylabel('angle (rad)')
linkaxes(get(gcf, 'Children'))
% hold on
% legend('0kg added','9kg added','18kg added','Location','northwest')
% legend('xa=0.75m','xa=1m','xa=1.25m','Location','northwest')
% legend('ya=-0.5m','ya=0m','ya=0.5m','Location','northwest')
% title('m=9kg, xa=1m, M=77kg, l=0.87m')
% print(gcf,'AddedMassOffset','-dsvg','-r300');
% exportgraphics(gcf,'MyFigure.svg','ContentType','vector');

% end % end for loop

%% ODE for inverted pendulum 
function dX = dPendulumStates(t, X, cart_acc_spline, M, m, l_lumped, alpha, I_lumped)
% X = [angle; angular velocity]

if abs(X(1))>=deg2rad(90) % making the pendulum stop falling once it hits the cart
    X(2)=0;
    gravity_ang_acc = 0;
    cart_ang_acc = 0;
    spring_ang_acc = 0;
else 
    cart_ang_acc = ((M+m)*ppval(cart_acc_spline,t)*cos(X(1)+alpha)*l_lumped)/I_lumped;                              
    gravity_ang_acc = ((M+m)*9.81*sin(X(1)+alpha)*l_lumped)/I_lumped;

    % what if the spring was exactly as stiff as it needs to be to resist gravity?
    if X(1)~= 0
        k = (gravity_ang_acc*I_lumped)/X(1);
    else
        k = 0;
    end 
    spring_ang_acc = -k*X(1)/I_lumped; % adding a torsional spring at pivot point
end

% dX = [angular velocity; angular acceleration] 
dX = [X(2); gravity_ang_acc + cart_ang_acc + spring_ang_acc];
end