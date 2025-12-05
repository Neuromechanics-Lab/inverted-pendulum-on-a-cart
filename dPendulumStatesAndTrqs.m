% Autumn Routt
% 9/4/2025
% Calculates angle, angular velocity, angular acceleration, and torques
% for inverted pendulum on a cart with SRM-style controller

function [dX,cart_trq,gravity_trq,musc_trq,cart_acc,cart_vel,cart_pos] = dPendulumStatesAndTrqs(t, X, ang_acc, cart_acc_spline, M, m, l_lumped, theta_a, I_lumped, kp, kv, ka, iter, delay, prevMuscTrq, timestep, cartAcc, cartVel, cartPos)
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

    % SRM-style controller for torque about the pivot point
    musc_trq = -kp*(X(iter-delay,1)+theta_a)-kv*X(iter-delay,2)-ka*ang_acc(iter-delay,1);

    % A "foot" which dictates max torque before a step would be necessary
    if musc_trq > (((M+m)*9.81)/2)*d2
        % display('before')
        % musc_trq
        % display('after')
        musc_trq = (((M+m)*9.81)/2)*d2;
        % display('capped CCW')
    elseif musc_trq < -1*(((M+m)*9.81)/2)*d1
        % display('before')
        % musc_trq
        % display('after')
        musc_trq = -1*(((M+m)*9.81)/2)*d1;
        % display('capped CW')
    end 

    % musc_act = -kp*(X(iter-delay,1)+theta_a)-kv*X(iter-delay,2)-ka*ang_acc(iter-delay,1);
    % 
    % % A "foot" which dictates max torque before a step would be necessary
    % if musc_act > (((M+m)*9.81)/2)*d2
    %     % display('before')
    %     % musc_trq
    %     % display('after')
    %     musc_act = (((M+m)*9.81)/2)*d2;
    %     % display('capped CCW')
    % elseif musc_act < -1*(((M+m)*9.81)/2)*d1
    %     % display('before')
    %     % musc_trq
    %     % display('after')
    %     musc_act = -1*(((M+m)*9.81)/2)*d1;
    %     % display('capped CW')
    % end 
    % 
    % musc_trq = musc_act*exp((-1*iter)/40); % Lockhart & Ting 2007, need to change this if timestep is no longer 1ms

    % capping muscle torque derivative 
    % (using ankle torque ramp rate from Jakubowski 2024 as a reference)
    if abs(musc_trq-prevMuscTrq)/timestep > 250
        if musc_trq-prevMuscTrq>0
            musc_trq = 250*timestep+prevMuscTrq;
        else
            musc_trq = -250*timestep+prevMuscTrq;
        end 
    end 

    total_ang_acc = (cart_trq+gravity_trq+musc_trq)/I_lumped;
end

% dX = [angular velocity; angular acceleration] 
dX = [X(iter,2); total_ang_acc];

end