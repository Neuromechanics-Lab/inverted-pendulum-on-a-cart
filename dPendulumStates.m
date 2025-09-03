%% ODE for inverted pendulum 
function dX = dPendulumStates(t, X, ang_acc, cart_acc_spline, M, m, l_lumped, alpha, I_lumped, kp, kv, ka, timestep, iter, delay)
% X = [angle; angular velocity]

% persistent pastX % persistent variable to check past timestep doesn't
%                  % really work here because ode45 jumps all around in time
% if exist('pastX','var') == 0
%     pastX = zeros(2,2/timestep);
% end 

if abs(X(iter-1,1))>=deg2rad(90) % making the pendulum stop falling once it hits the cart
    X(iter,2)=0;
    % gravity_ang_acc = 0;
    % cart_ang_acc = 0;
    % spring_ang_acc = 0;
    total_ang_acc = 0;
else 
    % cart_ang_acc = ((M+m)*ppval(cart_acc_spline,t)*cos(X(1)+alpha)*l_lumped)/I_lumped;                              
    % gravity_ang_acc = ((M+m)*9.81*sin(X(1)+alpha)*l_lumped)/I_lumped;
    cart_trq = ((M+m)*ppval(cart_acc_spline,t(iter))*cos(X(iter,1)+alpha)*l_lumped);                              
    gravity_trq = ((M+m)*9.81*sin(X(iter,1)+alpha)*l_lumped);

    % spring_ang_acc = -k*X(1)/I_lumped; % adding a torsional spring at pivot point

    % % what if the spring was exactly as stiff as it needs to be to resist gravity?
    % if X(1)~= 0
    %     k = (gravity_ang_acc*I_lumped)/X(1);
    % else
    %     k = 0;
    % end 

    % now what if the ankle had a gain for angle, ang vel, and ang acc?
    % total_ang_acc = (cart_ang_acc+gravity_ang_acc-kp*X(iter-1,1)-kv*X(iter-1,2))/(I_lumped+ka);
    total_ang_acc = (cart_trq+gravity_trq-kp*X(iter-delay,1)-kv*X(iter-delay,2)-ka*ang_acc(iter-delay,1))/I_lumped;

    % if isempty(pastX)
    %     spring_ang_acc = (kv*X(2) + kp*X(1))/I_lumped;
    % else
    %     spring_ang_acc = (ka*((X(2)-pastX(2,size(pastX,2)))/timestep) + kv*X(2) + kp*X(1))/I_lumped;
    % end 

end

% % and a time delay? 
% if t>delay
%     length(pastX)
%     delayedIndex = length(pastX)-delay/timestep
%     spring_ang_acc = (ka*((pastX(2,delayedIndex)-pastX(2,delayedIndex-1))/timestep)...
%         + kv*pastX(2,delayedIndex) + kp*pastX(1,delayedIndex))/I_lumped;
% else
%     spring_ang_acc = 0;
% end 

% spring_ang_acc = 0;

% dX = [angular velocity; angular acceleration] 
% dX = [X(2); gravity_ang_acc + cart_ang_acc - spring_ang_acc];
dX = [X(iter,2); total_ang_acc];

% pastX = [pastX,X];
end