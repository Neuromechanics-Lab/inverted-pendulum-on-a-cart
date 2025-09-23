% Autumn Routt
% 9/4/2025
% Calculates angle, angular velocity, angular acceleration, and torques
% for inverted pendulum on a cart with SRM-style controller

function [dX,cart_trq,gravity_trq,musc_trq] = dPendulumStatesAndTrqs(t, X, ang_acc, cart_acc_spline, M, m, l_lumped, theta_a, I_lumped, kp, kv, ka, iter, delay)
% X = [angle; angular velocity]

if abs(X(iter-1,1))>=deg2rad(90) % making the pendulum stop falling once it hits the cart
    X(iter,2)=0;
    cart_trq = 0;
    gravity_trq = 0;
    musc_trq = 0;
    total_ang_acc = 0;
else 
    cart_trq = ((M+m)*ppval(cart_acc_spline,t(iter))*cos(X(iter,1)+theta_a)*l_lumped);                              
    gravity_trq = ((M+m)*9.81*sin(X(iter,1)+theta_a)*l_lumped);

    % SRM-style controller for torque about the pivot point
    musc_trq = -kp*(X(iter-delay,1)+theta_a)-kv*X(iter-delay,2)-ka*ang_acc(iter-delay,1);
    total_ang_acc = (cart_trq+gravity_trq+musc_trq)/I_lumped;
end

% dX = [angular velocity; angular acceleration] 
dX = [X(iter,2); total_ang_acc];

end