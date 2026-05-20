% Autumn Routt
% 10/23/2025
% Inverted pendulum on a cart with added mass
% Torque at pivot point modelled with SRM-style controller as:
% T = kp*theta + kv*theta_dot + ka*theta_ddot
% And common time delay, lambda
% Turning inverted_pendulum_on_cart_added_mass_SRM into a function


function [settlingTime, muscPower, muscGrossWork, muscNetWork, muscImpulse, muscAvgTrq,muscTrq,x_sim,t_sim]...
    = inverted_pendulum_on_cart_added_mass_SRM_func(output,kp,kv,ka,delay,m,y_a,M,h,l,x_a,leaningStart,...
    simTime,timestep,pertDuration,cart_acc_time,cart_dec_time,pertMag,pertDir)

    arguments
        % Defining variables  y y
        output = 1;

        kp = 450; % angle gain 450
        kv = 350; % angular velocity gain 350
        ka = 0; % angular acceleration gain 0
        delay = 0; % common time delay (ms), must be <2s and must be an integer
        
        m = 0; % added mass (kg) (9kg = CDC recommended weight gain for 30 weeks pregant w normal starting BMI)
        y_a = 0.15; % added mass horizontal offset from pendulum arm (m), >0 is front-loaded, <0 is back-loaded
        M = 73.5; % body mass (kg) (73.5kg = 50th percentile for women in US)
        h = 1.612; % overall height (m) (1.612m = 50th percentile for women in US)
        l = 0.543*h; % body COM height (m) (avg COM height in women is 0.543*overall height)
        x_a = l; % added mass height (m)
        leaningStart = 1; % binary variable to decide if the pendulum will start at a lean to compensate for added mass (1) or at theta zero (0)

        simTime = 2; % how much time is simulated (seconds)
        timestep = 0.001;
        pertDuration = 145; % number of timesteps cart takes to accelerate and decelerate
        cart_acc_time = 0; % number of time steps before cart begins accelerating
        cart_dec_time = 480; % number of time steps before cart begins decelerating
        pertMag = 2.7; % max magnitude of cart acceleration (m/s^2)
        pertDir = 1; % 1 = cart moves right, -1 = cart moves left
    end

    theta_a = atan((m*y_a)/(M*l+m*x_a));
    l_lumped = sqrt(((M*l+m*x_a)/(M+m))^2+((m*y_a)/(M+m))^2);
    I_lumped = M*l^2+m*(x_a^2+y_a^2);
    
    % Defining cart acceleration profile 
    temp_t = 0:timestep:simTime;
    temp_acc = zeros(size(temp_t));
    temp_acc((1:pertDuration)+cart_acc_time) = ...
        -cos((1:pertDuration)*2*pi/pertDuration)+1; % acceleration
    temp_acc((1:pertDuration)+cart_dec_time) = ...
        cos((1:pertDuration)*2*pi/pertDuration)-1; % deceleration 
    cart_acc_spline = spline(temp_t,temp_acc*0.5*(pertMag)*pertDir);
    dec_end = cart_dec_time+pertDuration;
    
    % use the forward Euler method to find solution with the time delay
    x_sim = zeros(2000,2); % x_sim = [angle, angular velocity]
    if leaningStart==1
        x_sim(1:2000,1) = -1*theta_a;
    end 
    t_sim = zeros(2000,1);
    ang_acc = zeros(2000,1);
    cartTrq = []; gravTrq = []; muscTrq = []; 
    cartAcc = []; cartVel = []; cartPos = [];
    settlingTime = simTime;
    prevMuscTrq = 0;
    for iter = 2000:2000+size(temp_t,2)
        % [dX,cart_trq,gravity_trq,musc_trq,cart_acc,cart_vel,cart_pos] = dPendulumStatesAndTrqs(t_sim, x_sim, ang_acc, cart_acc_spline, M, m, l_lumped, theta_a, I_lumped, kp, kv, ka, iter, delay, prevMuscTrq, timestep, cartAcc, cartVel, cartPos);
        [dX,cart_trq,gravity_trq,musc_trq,cart_acc,cart_vel,cart_pos] = dPendulumStatesAndTrqsWActivation(t_sim, x_sim, ang_acc, cart_acc_spline, M, m, l_lumped, theta_a, I_lumped, kp, kv, ka, iter, delay, prevMuscTrq, timestep, cartAcc, cartVel, cartPos);
        new_x1 = x_sim(iter,1)+timestep*dX(1,:);
        if abs(new_x1)>=deg2rad(90)
            new_x2 = 0;
        else 
            new_x2 = x_sim(iter,2)+timestep*dX(2,:);
        end 
        x_sim = [x_sim;new_x1,new_x2];
        t_sim = [t_sim;(iter-2000)*timestep];
        ang_acc = [ang_acc;dX(2,:)];
        cartTrq = [cartTrq;cart_trq];
        gravTrq = [gravTrq;gravity_trq];
        muscTrq = [muscTrq;musc_trq];
        prevMuscTrq = musc_trq;
        cartAcc = [cartAcc;cart_acc];
        cartVel = [cartVel;cart_vel];
        cartPos = [cartPos;cart_pos];
        if settlingTime==simTime && (iter-2000)>dec_end && abs(x_sim(iter,1))<deg2rad(90)
            recentAngVelVals = x_sim(iter-(0:50),2);
            angVelNearZero = abs(recentAngVelVals) < 0.02;
            % slopeAngVel = (x_sim(iter,2)-x_sim(iter-50,2))/(50*timestep);
            recentAngAccVals = ang_acc(iter-(0:50),1);
            angAccNearZero = abs(recentAngAccVals) < 0.02;
            % if sum(angVelNearZero)==51 && slopeAngVel < 0.1
            if sum(angVelNearZero)==51 && sum(angAccNearZero)==51
                settlingTime = (iter-2050)*timestep;
            end 
        end 
    end 
    x_sim = x_sim(2001:size(x_sim,1),:);
    t_sim = t_sim(2001:size(t_sim,1),:);
    ang_acc = ang_acc(2001:size(ang_acc,1),:);
    
    muscPower = muscTrq.*x_sim(:,2);
    muscGrossWork = trapz(abs(muscPower));
    muscNetWork = trapz(muscPower);
    muscImpulse = trapz(muscTrq);
    muscAvgTrq = mean(abs(muscTrq));

    if output==1
        % fprintf('%d kg vals \n', m)
        fprintf('Settling time = %d \n',settlingTime);
        fprintf('Gross work = %d \n',muscGrossWork);
        fprintf('Net work = %d \n',muscNetWork);
        fprintf('Impulse = %d \n',muscImpulse);
        fprintf('Avg torque = %d \n \n',muscAvgTrq);
        
        % plot result (ang, ang vel, ang acc, trqs)
        % fig=figure;
        numplots = 22;
        subplot(numplots,1,1:2)
        plot(temp_t, temp_acc*0.5*(pertMag));
        ylabel('cart acc (m/s^2)')
        % plot(t_sim, cartPos*100)
        % ylabel('cart pos (m)')
        % ylim([-0.5,6])
        titleString = ['k_p=',num2str(kp),', k_v=',num2str(kv),', k_a=',num2str(ka),...
            ', delay=',num2str(delay),', M=',num2str(M),', l=',num2str(l),...
            ', m=',num2str(m),', x_a=',num2str(x_a),', y_a=',num2str(y_a)];
        % title(titleString)
        hold on
        
        subplot(numplots,1,4:7)
        plot(t_sim, rad2deg(x_sim(:,1)));
        ylabel('ang (deg)')
        xlim([0 length(t_sim)/1000])
        hold on
        
        subplot(numplots,1,9:12)
        plot(t_sim, rad2deg(x_sim(:,2)))
        % xlabel('time (s)')
        ylabel('ang vel (deg/s)')
        xlim([0 length(t_sim)/1000])
        hold on
        
        subplot(numplots,1,14:17)
        plot(t_sim,rad2deg(ang_acc(:,1)))
        % xlabel('time (s)')
        ylabel('ang acc (deg/s^2)')
        xlim([0 length(t_sim)/1000])
        hold on
        
        subplot(numplots,1,19:22)
        % plot(t_sim,gravTrq)
        % hold on
        % plot(t_sim,cartTrq)
        % plot(t_sim,muscTrq./M)
        plot(t_sim,muscTrq)
        xlabel('time (s)')
        % ylabel('musc trq (N*m/kg)')
        ylabel('musc trq (N*m)')
        % legend('Gravity','Cart','Muscles')
        % ylim([-60 20])
        xlim([0 length(t_sim)/1000])
        hold on
    
        % fileString = ['kp',num2str(kp),'_kv',num2str(kv),'_ka',num2str(ka),'_delay'...
        %     ,num2str(delay),'_M',num2str(M),'_l',num2str(l),'_m',num2str(m)...
        %     ,'_xa',num2str(x_a),'_ya',num2str(y_a)];
        % figName = ['Output/ModelOutput_',fileString,'.fig'];
        % savefig(figName)
        % saveas(fig,['Output/ModelOutput_',fileString,'.png']);
    
        % figure;
        % plot(t_sim,cartVel);
        % 
        % figure;
        % subplot(3,1,1)
        % plot(t_sim,cartAcc);
        % ylabel('Cart Acc (m/s^2)')
        % 
        % subplot(3,1,2)
        % plot(t_sim,cartVel)
        % ylabel('Cart Vel (m/s)')
        % 
        % subplot(3,1,3)
        % plot(t_sim,cartPos)
        % xlabel('Time (s)')
        % ylabel('Cart Pos (m)')
        
        % figure;
        % numplots=3;
        % subplot(numplots,1,1)
        % plot(temp_t, temp_acc*0.5*(pertMag));
        % ylabel('cart acc (m/s^2)')
        % xlim([0 length(t_sim)/1000])
        % subplot(numplots,1,2)
        % plot(t_sim,cartVel)
        % ylabel('cart vel (m/s)')
        % xlim([0 length(t_sim)/1000])
        % subplot(numplots,1,3)
        % plot(t_sim,cartPos)
        % ylabel('cart pos (m)')
        % xlim([0 length(t_sim)/1000])

    end 
end 