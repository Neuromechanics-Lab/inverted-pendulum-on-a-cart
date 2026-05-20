% Autumn Routt
% 5/19/2026
% Calculates angle, angular velocity, angular acceleration, and torques
% for inverted pendulum on a cart with SRM-style controller. Torques are
% found using a muscle model that takes in activation.

function [dX,cart_trq,gravity_trq,musc_trq,cart_acc,cart_vel,cart_pos] = dPendulumStatesAndTrqsWActivation(t, X, ang_acc, cart_acc_spline, M, m, l_lumped, theta_a, I_lumped, kp, kv, ka, iter, delay, prevMuscTrq, timestep, cartAcc, cartVel, cartPos)
% X = [angle; angular velocity]
% if iter>4000
%         disp(iter)
% end 

cart_acc = ppval(cart_acc_spline,t(iter));

cartVel = fnint(cart_acc_spline);
cart_vel = ppval(cartVel,t(iter));

cartPos = fnint(cartVel);
cart_pos = ppval(cartPos,t(iter));

d1 = 0.1473; % distance from malleolus to toes (Wunderlich 2001)
d2 = 0.0965; % distance from back of heel to malleolus (Wunderlich 2001)

if abs(X(iter-1,1))>=deg2rad(90) % making the pendulum stop falling once it hits the cart
    % display('fell down')
    X(iter,2) = 0;
    cart_trq = 0;
    gravity_trq = 0;
    musc_trq = 0;
    total_ang_acc = 0;
else
    % if iter>2002
    %     disp(iter)
    % end 
    cart_trq = ((M+m)*ppval(cart_acc_spline,t(iter))*cos(X(iter,1)+theta_a)*l_lumped);                              
    gravity_trq = ((M+m)*9.81*sin(X(iter,1)+theta_a)*l_lumped);

    musc_act = -kp*(X(iter-delay,1)+theta_a)-kv*X(iter-delay,2)-ka*ang_acc(iter-delay,1);

    % musc_trq = musc_act*exp((-1*iter)/40); % Lockhart & Ting 2007, need to change this if timestep is no longer 1ms

    tau = 40*timestep; % Lockhart & Ting 2007, time for muscle to generate torque from activation
    dt = timestep; % timestep = 1 ms
    dT = (musc_act-prevMuscTrq)*(dt/tau);
    musc_trq = prevMuscTrq + dT;

    total_ang_acc = (cart_trq+gravity_trq+musc_trq)/I_lumped;
end

% dX = [angular velocity; angular acceleration] 
dX = [X(iter,2); total_ang_acc];

end