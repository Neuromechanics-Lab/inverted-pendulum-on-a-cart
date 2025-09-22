% From Hansol Ryu 2025

% Define acceleration profile. 
% I defined as a time series and passed it to ODE as a spline. 
% It accelerates at time 0, decelerates at time 0.5. 
temp_t = -1:0.001:1;
temp_acc = zeros(size(temp_t));
perturb_num = 10;
temp_acc( (0:perturb_num)+1000) = ...
    -cos((0:perturb_num)*2*pi/perturb_num)+1; % acceleration
temp_acc( (0:perturb_num)+1500) = ...
    cos((0:perturb_num)*2*pi/perturb_num)-1; % decelataion 
cart_acc_spline = spline(temp_t,temp_acc*50);

% run ODE solver
odeopt = odeset('maxstep',1e-2);
[t_sim,x_sim] = ode45(@(t,x)dPendulumStates(t,x,cart_acc_spline), ...
    -1:0.001:1, ... simulation time range 
    [0, 0], ... initial condition [angle, angular velocity] 
    odeopt);

% plot result 
figure
subplot(5,1,1)
plot(temp_t, temp_acc);
ylabel('cart acceleration')
subplot(5,1,2:5)
plot(t_sim, x_sim(:,1));
xlabel('time')
ylabel('angle')
linkaxes(get(gcf, 'Children'))

%% ODE for inverted pendulum 
function dX = dPendulumStates(t,X, cart_acc_spline)
% X = [angle; angular velocity]
cart_ang_acc = ppval(cart_acc_spline,t)*cos(X(1)); % /l is omitted (normalization)                              
gravity_ang_acc = +9.81*sin(X(1)); % /l

% dX = [angular verlocity; angular acceleration] 
dX = [X(2); gravity_ang_acc + cart_ang_acc];
end