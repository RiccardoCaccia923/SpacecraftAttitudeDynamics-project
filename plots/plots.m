function plots(out,data)
    
    % plot options
    set(groot,'defaultTextInterpreter','latex',...
              'defaultAxesTickLabelInterpreter','latex',...
              'defaultLegendInterpreter','latex',...
              'defaultColorbarTickLabelInterpreter','latex');


    % define parameters
    t           = out.tout;
    w           = out.w;
    r           = out.r_N;
    err_norm    = out.errOrtNormA;
    err_rht     = out.rightHandleCheck;
    err_att     = out.A_err;
    err_w       = out.wError;
    err_theta   = out.thetaErr;
    ggPert      = out.ggPert;
    % dragPert    = out.dragPert;
    magFPert    = out.magFPert;
    B           = out.b_N;
    % dragF       = out.dragForce;
    Bmeas       = out.magnetometer;
    q_est_Dav   = out.q_est_DavqMethod;
    q_corr      = out.q_corr;
    w_corr      = out.w_corr;
    P_corr      = out.P_corr;
    attContrErr = out.alphaState(:,1:3);
    wContrErr   = out.alphaState(:,4:6);
    err_statAtt = out.thetaErrQ;
    err_dynAtt  = out.thetaErrMekf; 
    magDipole   = out.magDipole;
    magTorque   = out.magEffTorque;
    hDotWheel   = out.hDotWheel;
    hWheel      = out.hWheel;
    

    
    %% 1: angular velocity w
    for i = 1 : length(t)
        wNorm(i) = norm(w(i,:));
    end
    tP = t/data.T;
    figure('Name','1. Angular Velocity w','NumberTitle','off')
    plot(tP,w(:,1))
    hold on
    plot(tP,w(:,2))
    plot(tP,w(:,3))
    plot(tP,wNorm(:))
    % title('angular velocity w')
    grid on
    xlabel('$ T_{orb}$')
    ylabel('$w [rad/s]$')
    legend('$w_x$','$w_y$','$w_z$','$|w|$')

    %% 2: errors
    % 2.1 ERROR: orto-norm DCM 
    figure('Name','2.1 ortoNormERROR','NumberTitle','off')
    plot(t,err_norm)
    grid on
    xlabel('$time [s]$')
    ylabel('ortoNormalization error')
    legend('$err$')

    % 2.2 ERROR: right-handle vectors
    figure('Name','2.2 rightHandERROR','NumberTitle','off')
    plot(t,err_rht)
    grid on
    xlabel('$time [s]$')
    ylabel('right-handedness error for DCM')
    legend('$err$')

    % 2.3 ERROR: attitude
    % figure('Name','2.3 ortoNormERROR','NumberTitle','off')
    % subplot(2,1,1)
    % plot(t,err_att.e_deg(:,1))
    % hold on
    % plot(t,err_att.e_deg(:,2))
    % plot(t,err_att.e_deg(:,3))
    % grid on
    % xlabel('$time [s]$')
    % ylabel('$err attitude$')
    % legend('$err_11$','$err_21$','$err_31$')
    % 
    % % 2.3 ERROR: angular velocity 
    % subplot(2,1,2)
    % plot(t,err_w(:,1))
    % hold on
    % plot(t,err_w(:,2))
    % plot(t,err_w(:,3))
    % grid on
    % xlabel('$time [s]$')
    % ylabel('$err ang vel$')
    % legend('$err_x$','$err_y$','$err_z$')

    %% 3: environment
     figure('Name','3. Magnetic Field vs Position','NumberTitle','off')
    yyaxis left
    plot(t, r(:,1)); 
    hold on;
    plot(t, r(:,2));
    plot(t, r(:,3));
    ylabel('Position [km]');
    hold on
    yyaxis right
    plot(t, B(:,1));
    plot(t, B(:,2));
    plot(t, B(:,3));
    ylabel('Magnetic field [T]');
    xlabel('t [s]');
    legend('$r_x$','$r_y$','$r_z$','$B_x$','$B_y$','$B_z$','Location','northwest')
    for i = 1: length(t)
    Bnorm(i) = norm(B(i,:));
    end
    subplot(2,1,1)
    plot(t,r(:,1))
    hold on
    plot(t,r(:,2))
    plot(t,r(:,3))
    grid on
    xlabel('time [s]')
    ylabel('position [km]')
    legend('$r_x$','$r_y$','$r_z$')
        eq_idx = find(r(:,3).*circshift(r(:,3),-1) < 0);
    for k = 1:length(eq_idx)
        xline(t(eq_idx(k)), '--k','HandleVisibility','off');
    end
    [~, north_idx] = max(r(:,3));
    [~, south_idx] = min(r(:,3)); 
    xline(t(north_idx), '--b','HandleVisibility','off');
    xline(t(south_idx), '--r','HandleVisibility','off');

    subplot(2,1,2)
    plot(t,B(:,1))
    hold on
    plot(t,B(:,2))
    plot(t,B(:,3))
    plot(t,Bnorm)
    grid on
    xlabel('time [s]')
    ylabel('magnetic field [T]')
    legend('$B_x$','$B_y$','$B_z$','$|B|$') %'Location','northoutside','Orientation','horizontal'
    for k = 1:length(eq_idx)
        xline(t(eq_idx(k)), '--k','HandleVisibility','off');
    end
    xline(t(north_idx), '--b','HandleVisibility','off');
    xline(t(south_idx), '--r','HandleVisibility','off');


    %% 4: perturbations
    % 4.1: Gravity Gradient Perturbation
    figure('Name','4.1 Gravity Gradient','NumberTitle','off')
    subplot(2,1,1)
    plot(t,ggPert(:,1))
    hold on
    plot(t,ggPert(:,2))
    plot(t,ggPert(:,3))
    grid on
    xlabel('$time [s]$')
    ylabel('$Torque [Nm]$')
    legend('$GG_x$','$GG_y$','$GG_z$')
    ggPertNorm = zeros(length(t),1);
    for i = 1 : length(t)
        ggPertNorm(i) = norm(ggPert(i,:));
    end
    subplot(2,1,2)
    plot(t,ggPertNorm)
    grid on
    xlabel('$time [s]$')
    ylabel('$Torque [Nm]$')
    legend('$|GG|$')

    % 4.2: Amtmospheric Drag Perturbation
    % figure('Name','4.2a Atmospheric ForceVsDrag ','NumberTitle','off')
    % subplot(2,1,1)
    % dragFNorm = zeros(length(t),1);
    % for i = 1 : length(t)
    %     dragFNorm(i) = norm(dragF(i,:));
    % end
    % plot(t,dragFNorm)
    % grid on
    % xlabel('$time [s]$')
    % ylabel('$Force [Nm]$')
    % legend('$|dragForce|$')
    % subplot(2,1,2)
    % dragNorm = zeros(length(t),1);
    % for i = 1 : length(t)
    %     dragNorm(i) = norm(dragPert(i,:));
    % end
    % plot(t,dragNorm)
    % grid on
    % xlabel('$time [s]$')
    % ylabel('$Torque [Nm]$')
    % legend('$|dragTorque|$')
    % 
    % figure('Name','4.2b Atmospheric Drag ','NumberTitle','off')
    % subplot(2,1,1)
    % plot(t,dragPert(:,1))
    % hold on
    % plot(t,dragPert(:,2))
    % plot(t,dragPert(:,3))
    % grid on
    % xlabel('$time [s]$')
    % ylabel('$Torque [Nm]$')
    % legend('$drag_x$','$drag_y$','$drag_z$')
    % dragNorm = zeros(length(t),1);
    % for i = 1 : length(t)
    %     dragNorm(i) = norm(dragPert(i,:));
    % end
    % subplot(2,1,2)
    % plot(t,dragNorm)
    % grid on
    % xlabel('$time [s]$')
    % ylabel('$Torque [Nm]$')
    % legend('$|dragTorque|$')

    % 4.3: Magnetic Field Perturbation
    figure('Name','4.3 Magnetic Field ','NumberTitle','off')
    subplot(2,1,1)
    plot(t,magFPert(:,1))
    hold on
    plot(t,magFPert(:,2))
    plot(t,magFPert(:,3))
    grid on
    xlabel('$time [s]$')
    ylabel('$Torque [Nm]$')
    legend('$magF_x$','$magF_y$','$magF_z$')
    magFNorm = zeros(length(t),1);
    for i = 1 : length(t)
        magFNorm(i) = norm(magFPert(i,:));
    end
    subplot(2,1,2)
    plot(t,magFNorm)
    grid on
    xlabel('$time [s]$')
    ylabel('$Torque [Nm]$')
    legend('$|magF|$')

    % %% 5: sensors
    % % 5.1 magnetometer error on each axis
    % for i = 1:length(t)
    %     BmeasNorm(i) = norm(Bmeas.signals.values(i,:));
    %     Breal_b(i,:) = out.A_BN(:,:,i)*B(i,:)';
    % end
    % figure('Name','5.1 Magnetometer ','NumberTitle','off')
    % plot(Bmeas.time,Bmeas.signals.values(:,1))
    % hold on
    % plot(Bmeas.time,Bmeas.signals.values(:,2))
    % plot(Bmeas.time,Bmeas.signals.values(:,3))
    % plot(Bmeas.time,BmeasNorm)
    % grid on
    % xlabel('$time [s]$')
    % ylabel('$B measured [T]$')
    % legend('$Bmeas_x$','$Bmeas_y$','$Bmeas_z$','$|B|$')
    % 
    % figure('Name','5.2 Magnetometer Measure Error ','NumberTitle','off')
    % magMeasErr = Bmeas.signals.values - Breal_b;
    % plot(Bmeas.time(2:end),magMeasErr((2:end),:))
    % grid on
    % xlabel('$time [s]$')
    % ylabel('$measurement error [T]$')
    % legend('$err_x$','$err_y$','$err_z$')
    % 
    % % 5.1 magnetometer Allan Variance diagram
    % xErr = magMeasErr((2:end),1);   
    % yErr = magMeasErr((2:end),1);
    % zErr = magMeasErr((2:end),1);
    % Ts = Bmeas.time(2:end);
    % N = length(xErr);
    % m_vec = round(logspace(0, log10(N/10), 50));
    % [taux, adevx] = allan_deviation(xErr, Ts, m_vec);
    % [tauy, adevy] = allan_deviation(yErr, Ts, m_vec);
    % [tauz, adevz] = allan_deviation(zErr, Ts, m_vec);
    % 
    % figure('Name','5.3 Allan Variance Diagram','NumberTitle','off')
    % colors = lines(3);
    % loglog(taux, adevx, 'LineWidth', 1.5);
    % hold on
    % loglog(tauy, adevy, 'LineWidth', 1.5);
    % loglog(tauz, adevz, 'LineWidth', 1.5);
    % grid on;
    % t1 = 0.6;    % fine white-noise
    % t2 = 1.2;    % fine transition
    % t3 = 10;     % inizio bias-instability
    % t4 = 100;    % inizio brown    
    % xline(t1, '--k', 'LineWidth', 1);
    % xline(t2, '--k', 'LineWidth', 1);
    % xline(t3, '--k', 'LineWidth', 1);
    % yl = ylim;
    % text(sqrt(1e-1*t1),  yl(2)/2, 'White noise',       'HorizontalAlignment','center');
    % text(sqrt(t1*t2),    yl(2)/3, 'Transition',        'HorizontalAlignment','center');
    % text(sqrt(t2*t3),    yl(2)/2, 'Bias instability',  'HorizontalAlignment','center');
    % text(sqrt(t3*t4),    yl(2)/3, 'Random walk',       'HorizontalAlignment','center');
    % xlabel('$\tau$ [s]');
    % ylabel('$\sigma_A/tau$ [T]');
    % 
    % % 6: STATIC attitude determination
    % % dqStatic = zeros(length(t), 4);
    % % err_qStatic_deg = zeros(length(t), 3);
    % % q_est_DavSF = [q_est_Dav(:,4),q_est_Dav(:,1),q_est_Dav(:,2),q_est_Dav(:,3)];
    % % for i = 1 : length(t)
    % %     q_true = dcm2quat(out.A_BN(:,:,i));
    % %     q_estStatic = q_est_DavSF(i,:);
    % %     dqStatic(i,:) = quatmultiply(quatconj(q_estStatic), q_true);
    % %     if dqStatic(i,1) < 0
    % %     dqStatic(i,:) = -dqStatic(i,:);
    % %     end
    % %     err_qStatic_deg(i, :) = rad2deg(2*dqStatic(i,2:4));
    % % end
    % % 
    % % figure('Name','6.1 Static Attitude Estimation Error','NumberTitle','off')
    % % subplot(2,1,1)
    % % plot(t,err_qStatic_deg(:,1))
    % % hold on
    % % plot(t,err_qStatic_deg(:,2))
    % % plot(t,err_qStatic_deg(:,3))
    % % xlabel('$time [s]$')
    % % ylabel('$error [deg]$')
    % % ylim([-1 1]); 
    % % legend('$roll$','$pitch$','$yaw$')
    % % 
    % % % 7: DYNAMIC attitude determination
    % % dq = zeros(length(t), 4);
    % % err_q_B_deg = zeros(length(t), 3);
    % % err_w_B = zeros(length(t), 3);
    % % for i = 1 : length(t)
    % %     q_true = dcm2quat(out.A_BN(:,:,i));
    % %     q_estCorr = q_corr(i,:);
    % %     dq(i,:) = quatmultiply(quatconj(q_estCorr), q_true);
    % %     if dq(i,1) < 0
    % %     dq(i,:) = -dq(i,:);
    % %     end
    % %     err_q_B_deg(i, :) = rad2deg(2*dq(i,2:4));
    % %     err_w_B(i,:) = out.w(i,:) - w_corr(:,:,i);
    % %     err_cov_att(i,:) = 3 * sqrt(diag(P_corr(1:3,1:3,i))); % 3-sigma
    % %     err_cov_w(i,:) = 3 * sqrt(diag(P_corr(4:6,4:6,i))); % 3-sigma
    % % end
    % 
    % %figure('Name','7.1 Dynamic Attitude Estimation Error','NumberTitle','off')
    % % subplot(2,1,2)
    % % plot(t,err_q_B_deg(:,1))
    % % hold on
    % % plot(t,err_q_B_deg(:,2))
    % % plot(t,err_q_B_deg(:,3))
    % % xlabel('$time [s]$')
    % % ylabel('$error [deg]$')
    % % ylim([-1 1]); 
    % % legend('$roll$','$pitch$','$yaw$')
    % % 
    % % figure('Name','7.2 Dynamic Angular Velocity Estimation Error','NumberTitle','off')
    % % plot(t,err_w_B(:,1))
    % % hold on
    % % plot(t,err_w_B(:,2))
    % % plot(t,err_w_B(:,3))
    % % xlabel('$time [s]$')
    % % ylabel('$error [rad/s]$')
    % % ylim([-0.01 0.01]);
    % % legend('$w_{roll}$','$w_{pitch}$','$w_{yaw}$')
    % % 
    % % figure('Name','7.3 Covariance Estimation Analysis','NumberTitle','off')
    % % subplot(2,1,1)
    % % plot(t,rad2deg(err_cov_att))
    % % grid on
    % % xlabel('$time [s]$')
    % % ylabel('$3\sigma$ Uncertainty [deg]')
    % % title('\textbf{Attitude Estimation Uncertainty}')
    % % legend('roll','pitch','yaw')
    % % 
    % % subplot(2,1,2)
    % % plot(t,err_cov_w(:,1))
    % % hold on
    % % plot(t,err_cov_w(:,2))
    % % plot(t,err_cov_w(:,3))
    % % grid on
    % % xlabel('$time [s]$')
    % % ylabel('$3\sigma$ Uncertainty [rad/s]')
    % % title('\textbf{Angular Velocity Estimation Uncertainty}')
    % % legend('$3\sigma_{\omega_x}$', '$3\sigma_{\omega_y}$', '$3\sigma_{\omega_z}$')
    % 
    % 
    % % %% 8: POINTING PHASE ERRORS
    % enableIdx = 0;
    % for i = 1 : length(out.tout)
    %     if out.enablePoint(i) == 0
    %         enableIdx = i;
    %         break
    %     end
    % end
    % 
    % % 8.1: Attitude Error
    % figure('Name','Attitude Control','NumberTitle','off')
    % plot(t(enableIdx:end),rad2deg(attContrErr))
    % grid on
    % title('Attitude Estimation Error - Pointing Phase')
    % xlabel('Time [s]')
    % ylabel(' Error [deg]$')
    % ylim([-70 70]);
    % legend('Roll','Pitch','Yaw')
    % 
    % % 8.2: Angular Velocity Error
    % figure('Name','Angular Velocity Control','NumberTitle','off')
    % plot(t(enableIdx:end),rad2deg(wContrErr))
    % grid on
    % title('Angular Velocity Estimation Error - Pointing Phase')
    % xlabel('Time [s]')
    % ylabel('Error[deg/s]$')
    % ylim([-0.1 0.1]);
    % legend('Roll','Pitch','Yaw')
    % 
    % %% 9: ACTUATORS
    % figure('Name','Magnetorquers','NumberTitle','off')
    % subplot(2,1,1)
    % plot(t./data.T,magDipole(:,1))
    % hold on
    % plot(t./data.T,magDipole(:,2))
    % plot(t./data.T,magDipole(:,3))
    % yline(30,'r--')
    % yline(-30,'r--')
    % grid on
    % xlabel('$ T_{orb}$')
    % ylabel('$m [Am^2]$')
    % legend('$m_x$','$m_y$','$m_z$')
    % 
    % subplot(2,1,2)
    % plot(t./data.T,magTorque(:,1))
    % hold on
    % plot(t./data.T,magTorque(:,2))
    % plot(t./data.T,magTorque(:,3))
    % grid on
    % xlabel('$ T_{orb}$')
    % ylabel('$Torque [Nm]$')
    % legend('$\tau_x$','$\tau_y$','$\tau_z$')
    % 
    % figure('Name','Reaction Wheel','NumberTitle','off')
    % subplot(2,1,1)
    % plot(t(enableIdx:end)./data.T,hDotWheel(:,3))
    % % yline(0.265,'r--')
    % % yline(-0.265,'r--')
    % grid on
    % xlabel('$T_{orb}$')
    % ylabel('Torque [Nms]')
    % legend('$hDot_{wheel}$')
    % 
    % subplot(2,1,2)
    % plot(t(enableIdx:end)./data.T,hWheel(:,3))
    % % yline(18,'r--')
    % % yline(-18,'r--')
    % grid on
    % xlabel('$T_{orb}$')
    % ylabel('Momentum [Nms]')
    % legend('$h_{wheel}$')
    % 
    % 
    % %% 10: RESULTS
    % tP = t(enableIdx:end)./data.T;
    % figure('Name','10.1 Detumbling action','NumberTitle','off')
    % plot(t./data.T,out.w(:,1))
    % hold on
    % plot(t./data.T,out.w(:,2))
    % plot(t./data.T,out.w(:,3))
    % plot(t./data.T,out.w_norm(:))
    % yline(0.005,'r--')
    % % yline(-0.005,'r--')
    % xline(t(enableIdx)/data.T,'g--')
    % grid on
    % grid minor
    % xlim([0,3])
    % xlabel('$ T_{orb}$')
    % ylabel('$w [rad/s]$')
    % legend('$w_x$','$w_y$','$w_z$','$|w|$','$DET_{treshold}$','$DET_{END}$')
    % 
    % figure('Name','10.2 RPE','NumberTitle','off')
    % % plot(tP,w((enableIdx:end),1))
    % % hold on
    % % plot(tP,w((enableIdx:end),2))
    % % plot(tP,w((enableIdx:end),3))
    % plot(tP,wNorm(enableIdx:end)*180/pi)
    % yline(0,'k--')
    % yline(5e-3,'r--')
    % % title('angular velocity w')
    % grid on
    % grid minor
    % axis tight
    % ylim([-1e-3,6e-3])
    % xlabel('$ T_{orb}$')
    % ylabel('$w_{err} [rad/s]$')
    % % legend('$|w_{err}|$')%('$Roll$','$Pitch$','$Yaw$'),
    % 
    % figure('Name','10.3 APE','NumberTitle','off')
    % plot(tP,err_theta(enableIdx:end))
    % yline(0,'k--')
    % yline(5,'r--')
    % grid on
    % grid minor
    % axis tight
    % ylim([-5,20])
    % xlabel('$ T_{orb}$')
    % ylabel('$ \Theta_{err} [deg]$')
    % 
    % figure('Name','10.4 static att theta err','NumberTitle','off')
    % plot(tP,err_statAtt(enableIdx:end))
    % yline(0,'k--')
    % yline(5e-1,'r--')
    % grid on
    % grid minor
    % axis tight
    % ylim([-2,2])
    % xlabel('$ T_{orb}$')
    % ylabel('$ \Theta_{err} [deg]$')
    % 
    % figure('Name','10.5 dynamic att theta err','NumberTitle','off')
    % plot(tP,err_dynAtt(enableIdx:end))
    % yline(0,'k--')
    % yline(5e-1,'r--')
    % grid on
    % grid minor
    % axis tight
    % ylim([-2,2])
    % xlabel('$ T_{orb}$')
    % ylabel('$ \Theta_{err} [deg]$')
    % 
    % figure('Name','10.6 DCM components Error','NumberTitle','off')
    % subplot(2,1,1)
    % yline(1,'k--')
    % hold on
    % plot(tP,squeeze(err_att(1,1,(enableIdx:end))))
    % plot(tP,squeeze(err_att(2,2,(enableIdx:end))))
    % plot(tP,squeeze(err_att(3,3,(enableIdx:end))))
    % grid on
    % grid minor
    % axis tight
    % ylim([0.95,1.05])
    % xlabel('$ T_{orb}$')
    % legend('$A_{err}11$','$A_{err}22$','$A_{err}33$')
    % 
    % subplot(2,1,2)
    % yline(0,'k--')
    % hold on
    % plot(tP,squeeze(err_att(1,2,(enableIdx:end))))
    % plot(tP,squeeze(err_att(1,3,(enableIdx:end))))
    % plot(tP,squeeze(err_att(2,1,(enableIdx:end))))
    % plot(tP,squeeze(err_att(2,3,(enableIdx:end))))
    % plot(tP,squeeze(err_att(3,1,(enableIdx:end))))
    % plot(tP,squeeze(err_att(3,2,(enableIdx:end))))
    % 
    % grid on
    % grid minor
    % axis tight
    % ylim([-0.3,0.3])
    % xlabel('$ T_{orb}$')
    % legend('$A_{err}12$','$A_{err}13$','$A_{err}21$','$A_{err}23$','$A_{err}31$','$A_{err}32$')
end