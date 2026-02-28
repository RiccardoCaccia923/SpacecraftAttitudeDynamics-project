function [tau, adev] = allan_deviation(x, Ts, m_vec)

x = x(:);               
N = length(x);

if isscalar(Ts)
    Ts = Ts;
else
    t  = Ts(:);
    Ts = mean(diff(t));  
end

% --- se m_vec non passato, creo vettore logaritmico ---
if nargin < 3 || isempty(m_vec)
    m_vec = round(logspace(0, log10(N/10), 50));
    m_vec = unique(m_vec);
else
    m_vec = m_vec(:).';  % rendo riga
end

tau  = m_vec * Ts;       % intervallo integrazione
adev = zeros(size(m_vec));

for k = 1:length(m_vec)
    m = m_vec(k);
    K = floor(N / m);    % numero di blocchi
    if K < 2
        adev(k) = NaN;
        continue;
    end

    % media su blocco m-size
    x_block = reshape(x(1:K*m), m, K);   % m × K
    x_bar   = mean(x_block, 1);          % 1 × K

    % Allan deviation
    adev(k) = sqrt( 0.5 * mean( diff(x_bar).^2 ) );
end

end
