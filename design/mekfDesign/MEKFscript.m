%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MULTIPLICATIVE EXTENDED KALMAN FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% DEFINIZIONE STATO DI ERRORE
% dTheta = [alphaX , alphaY , alphaZ]';               % errore di assetto
% dOmega = [alphaDotX , alphaDotY , alphaDotZ]';      % errore di velocita angolare
% dx = [dTheta ; dOmega]';                            % stato di errore

%% INIZIALIZZAZIONE MEKF
N = length(out.tout);
dt = out.tout(2) - out.tout(1);
fprintf('Analisi dati Simulink: %d campioni, dt = %.4f s\n', N, dt);

q_corr = [1, 0, 0, 0];   
w_corr = [0, 0, 0];

P_corr = eye(6)*1e-8;       

ssNoiseVariance = sensors.starSensor.noiseDensity/sqrt(2*sensors.starSensor.sampleTime)*1e-2;
magNoiseVariance = sensors.magnetometer.noiseDensity/sqrt(2*sensors.magnetometer.sampleTime)*1e3;      

J = data.Imatrix;

hist_q_est = zeros(4, N);
hist_w_est = zeros(3, N);
hist_err_cov = zeros(3, N);

%% LOOP FILTRO
fprintf('Esecuzione MEKF sui dati registrati...\n');

for k = 1 : N
    % PREDICTION STEP
    % INPUT:
    % q_corr (1x4): quaternione corretto al passo k-1 [w, x, y, z] SCALAR FIRST
    % w_corr (3x1): velocità angolare corretta al passo k-1
    % P_corr (6x6): matrice di covarianza corretta al passo k-1
    % J (3x3): matrice d'inerzia
    % h_RW (3x1): momento angolare della wheel
    % hDot_RW (3x1): torque della wheel
    % dist_Tq (3x1): torque delle disturbances
    % control_Tq (3x1): torque di controllo
    % dt: sampleTime simlazione (0.01)

    % Misure
    q_meas_in = out.q_est_DavqMethod(:,:,k)'; 
    q_raw = q_meas_in; 
    if size(q_raw,1) == 1 
        q_raw = q_raw'; 
    end % Forziamo a colonna per comodità
    q_meas = [q_raw(4); q_raw(1); q_raw(2); q_raw(3)]';

    b_meas_B = out.magnetometer.signals.values(k,:);   
    
    % Riferimenti (Modello magnetico all'istante k)
    b_meas_N = out.magReal_N(k,:);                   
    
    % Controllo 
    control_Tq = [0,0,0];  % Per ora zero 
    h_RW       = [0,0,0];
    hDot_RW    = [0,0,0];
    dist_Tq    = [0,0,0]'; % out.A_BN(:,:,k)*out.dist_Tq(k,:)';
    
    %% propagazione dello stato nominale
    % cinematica quaernione
    w_corr_quat = [0, w_corr];
    qDot_pred = 0.5 * quatmultiply(q_corr , w_corr_quat);
    q_pred = q_corr + qDot_pred * dt;
    q_pred = q_pred/norm(q_pred);

    % dinamica sistema
    wDot_pred = J \ (-cross(w_corr,J*w_corr') - cross(w_corr,h_RW) - hDot_RW + dist_Tq' + control_Tq)';
    w_pred = w_corr + wDot_pred'*dt;
    
    %% propagazione matrice di covarianza
    % matrice di sistema continua (F)
    F = [zeros(3) , eye(3) ; 
         zeros(3) , zeros(3)];      % SEMPLIFICAZIONE PER PUNTAMENTO INERZIALE -> piccole rotazioni
    
    % matrice di transizione discreta (PHI)
    PHI = eye(6) + F*dt;
    
    % matrice di covarianza del rumore di processo discreta (Q)
    qTheta = 1e-9;                                      % noise density assetto -> molto piccolo
    qOmega = 1e-7;                                      % noise density velocita angolare -> grande -> non ho i gyro
    Q = [eye(3,3)*qTheta , zeros(3,3) ; 
            zeros(3,3) , eye(3,3)*qOmega] *dt;
    
    % MATRICE DI COVARIANZA
    P_pred = PHI * P_corr * PHI' + Q;
    P_pred = 0.5 * (P_pred + P_pred');                  % mantenimento simmetria

    %% MEASURE STEP
    if size(q_pred,1) > 1
        q_pred = q_pred';                               % metti in riga
    end
    dq = quatmultiply(quatconj(q_pred),q_meas);
    if dq(1) < 0
            dq = -dq;                                   % evita ambiguita
    end
    res_ss = dq(2:4) * 2;                               % rotazioni piccole -> err vettoriale ~ err angolare/2
    b_pred = quatrotate(q_pred,b_meas_N);
    res_mag = b_meas_B - b_pred;
    b_pred_skew = [0 , -b_pred(3) , b_pred(2) ;
                   b_pred(3) , 0 , -b_pred(1) ; 
                   -b_pred(2) , b_pred(1) , 0 ];
    
    y = [res_ss , res_mag]';
    

    %% matrice di misura (H)
    H_ss = [eye(3) , zeros(3)];
    H_mag = [b_pred_skew , zeros(3,3)];
    H = [H_ss ; H_mag];

    %% matrice di rumore di misura (R)
    R_ss_mat = eye(3) * ssNoiseVariance; 
    R_mag_mat = eye(3) * magNoiseVariance;
    R = blkdiag(R_ss_mat, R_mag_mat);

    %% kalman gain (K)
    % S = H * P_pred * H' + R;
    % epsilon = 1e-12; 
    % S_reg = S + eye(size(S)) * epsilon; 
    % K = (P_pred * H') / S_reg;

    % 1. Calcoliamo S (Denominatore) e Num_K (Numeratore) separati
    S = H * P_pred * H' + R;
    Num_K = P_pred * H';

    % 2. CACCIA ALL'INFINITO (Il vero colpevole probabile)
    % isinf() trova numeri troppo grandi (> 1.7e308) che isnan() ignora
    if any(isinf(P_pred(:)))
        error('CRASH al tempo %f: P_pred è diventata INFINITA! Il filtro sta divergendo. Riduci Q o P0.', k*dt);
    end

    if any(isinf(S(:)))
        error('CRASH al tempo %f: La matrice S contiene valori INFINITI. P_pred è esplosa.', k*dt);
    end

    % 3. CACCIA AL NAN MATEMATICO (0 * Inf)
    if any(isnan(S(:)))
        % Se arriviamo qui, P o H contenevano Inf, e moltiplicando per 0 hanno generato NaN
        error('CRASH al tempo %f: S è diventata NaN (probabile operazione 0 * Inf).', k*dt);
    end

    % 4. Fix di Regolarizzazione (Anti-singolarità)
    % Aggiungiamo un epsilon minuscolo per evitare divisioni per zero perfette
    S = S + eye(size(S)) * 1e-12;

    % 5. Calcolo K Protetto
    % Usiamo try-catch per vedere se l'errore RCOND avviene proprio qui
    try
        K = Num_K / S;
    catch ME
        fprintf('--- ERRORE MATEMATICO RILEVATO ---\n');
        disp('Ecco la matrice S che non si riesce a invertire:');
        disp(S);
        rethrow(ME); % Rilancia l'errore per fermare Simulink
    end
    
    %% CORRECTION STEP
    % aggiornamento stato
    dx_corr = K * y;

    % correzione assetto (multiplicative)
    dTheta = 0.5*+dx_corr(1:3);
    dq = [1 ; dTheta]';
    dq = dq / norm(dq);
    q_corr = quatmultiply(q_pred,dq);
    q_corr = q_corr/norm(q_corr);

    % correzione velocità angolari
    dw = dx_corr(4:6)';
    w_corr = w_pred + dw;

    % aggiornamento matrice covarianza: formula Joseph
    P_corr = (eye(6) - K * H) * P_pred * (eye(6) - K * H)' + K * R * K';
    P_corr = 0.5 * (P_corr + P_corr');      % Simmetria
    min_cov = 1e-10;                        % Valore minimo di incertezza
    for i = 1:6
        if P_corr(i,i) < min_cov
            P_corr(i,i) = min_cov;
        end
    end

    %% SALVATAGGIO RISULTATI
    hist_q_est(:, k) = q_corr;
    hist_w_est(:, k) = w_corr;
    hist_err_cov(:, k) = 3 * sqrt(diag(P_corr(1:3,1:3))); % 3-sigma

end

%% 5. ANALISI RISULTATI (Grafici)
% Calcolo errore rispetto alla verità Simulink
sim_q_true = dcm2quat(out.A_BN); 
sim_w_true = out.w;
sim_time = out.tout;

err_deg = zeros(3, N);
for k=1:N
    dq_err = quatmultiply(quatconj(hist_q_est(:,k)'), sim_q_true(k,:));
    if dq_err(1) < 0
        dq_err = -dq_err; 
    end
    err_deg(:,k) = rad2deg(dq_err(2:4) * 2);
end

figure(1); 
subplot(2,1,1); 
plot(sim_time, err_deg); 
grid on;
title('Errore Assetto EKF vs Simulink Truth');
ylabel('Errore (deg)'); legend('Roll', 'Pitch', 'Yaw');
% ylim([-1 1]); 

subplot(2,1,2);
plot(sim_time, sim_w_true - hist_w_est'); 
grid on;
title('Errore Stima Velocità Angolare');
ylabel('Errore (rad/s)'); legend('Wx', 'Wy', 'Wz');
% ylim([-0.05 0.05]); 