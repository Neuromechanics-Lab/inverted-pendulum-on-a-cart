% Autumn Routt
% 9/19/25
% Function to find torques for inverted pendulum on a cart with added mass
% Uses SRM-style controller for muscle torque

function [dX,cart_trq,gravity_trq,musc_trq] = dPendulumDemo(t, X, ang_acc, cart_acc_spline, M, m, l_lumped, theta_a, I_lumped, kp, kv, ka, iter, delay)

if abs(X(iter-1,1))>=deg2rad(90) % making the pendulum stop falling once it hits the cart
    X(iter,2)=0;
    cart_trq = 0;
    gravity_trq = 0;
    musc_trq = 0;
    total_ang_acc = 0;
else 
    cart_trq = ((M+m)*ppval(cart_acc_spline,t(iter))*cos(X(iter,1)+theta_a)*l_lumped);                              
    gravity_trq = ((M+m)*9.81*sin(X(iter,1)+theta_a)*l_lumped);

    musc_trq = -kp*(X(iter-delay,1)+theta_a)-kv*X(iter-delay,2)-ka*ang_acc(iter-delay,1); % SRM-style controller
    total_ang_acc = (cart_trq+gravity_trq+musc_trq)/I_lumped;

end

dX = [X(iter,2); total_ang_acc];

end