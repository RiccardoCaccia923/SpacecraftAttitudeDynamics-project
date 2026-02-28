function vec = attitudeErrComp(A_BL,A_BL_ref)
% A_BL      : 3x3xN   (DCM Body->L al tempo k)
% A_BL_ref  : 3x3xN OPPURE 3x3 (fisso) OPPURE vuoto
% mat, vec  : strutture con risultati per tutti gli N

N = size(A_BL, 3);

% prealloc
R_err      = zeros(3,3,N);
diag_cos   = zeros(3,N);
phi_deg    = zeros(3,N);
theta_deg1 = zeros(1,N);   % dal modo "mat"
e_rad      = zeros(3,N);
e_deg      = zeros(3,N);
theta_deg2 = zeros(1,N);   % dal modo "vec"

hasRefSeq  = (nargin >= 2) && ~isempty(A_BL_ref) && (size(A_BL_ref,3) == N);
hasRefFix  = (nargin >= 2) && ~isempty(A_BL_ref) && (size(A_BL_ref,3) == 1);

for k = 1:N

    A = A_BL(:,:,k);

    if hasRefSeq
        R = A_BL_ref(:,:,k).' * A;
    elseif hasRefFix
        R = A_BL_ref.' * A;
    else
        R = A;                       % riferimento = identità
    end

    % --- proiezione su SO(3) ---
    [U,~,V] = svd(R);
    R = U*V.';
    if det(R) < 0
        U(:,3) = -U(:,3);
        R = U*V.';
    end

    R_err(:,:,k) = R;

    % ---------- modo 1: analisi a coseni ----------
    % dc = diag(R);
    % dc = max(-1, min(1, dc));
    % ph = acos(dc);                  % [rad]
    % diag_cos(:,k) = dc;
    % phi_deg(:,k)  = rad2deg(ph);
    % 
    c = (trace(R) - 1)/2;
    c = max(-1, min(1, c));
    theta = acos(c);                % [rad]
    theta_deg1(k) = rad2deg(theta);

    % ---------- modo 2: vettore errore ----------
    if theta < 1e-8
        e = 0.5 * [ R(3,2)-R(2,3);
                    R(1,3)-R(3,1);
                    R(2,1)-R(1,2) ];
    else
        s = 2*sin(theta);
        u = [ (R(3,2)-R(2,3))/s;
              (R(1,3)-R(3,1))/s;
              (R(2,1)-R(1,2))/s ];
        n = norm(u);
        if n > 0, u = u/n; end
        e = theta * u;
    end

    e_rad(:,k)   = e;
    e_deg(:,k)   = rad2deg(e);
    theta_deg2(k)= norm(e_deg(:,k),2);   % stessa cosa di theta_deg1(k)
end

% mat.R_err     = R_err;
% mat.diag_cos  = diag_cos;
% mat.phi_deg   = phi_deg;
% mat.theta_deg = theta_deg1;

vec.e_rad     = e_rad';
vec.e_deg     = e_deg';
vec.theta_deg = theta_deg2';
end
