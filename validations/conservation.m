function conservation(Imatrix,t,w)

    n = length(t);
    hNorm = zeros(n,1);
    T = zeros(n,1);
    
    for i = 1 : n
    
        % angular momentum conservation
        h = Imatrix*w(i,:)';
        hNorm(i) = norm(h);
    
        % kinetic energy conservation
        T(i) = w(i,:)*(Imatrix*w(i,:)');
    end


    % set a tollerance
    toll = 1e-6;

    figure
    subplot(2,1,1)
    plot(t,hNorm)
    grid on
    hold on
    title('angular momentum')
    hdrift = max(hNorm)-min(hNorm);
    if hdrift < toll
        disp(['angular momentum conserved, max drift : ', num2str(hdrift)])
    else
        disp(['angular momentum NOT conserved, max drift : ', num2str(hdrift)])
    end

    subplot(2,1,2)
    plot(t,T)
    title('kinetic energy')
    Tdrift = max(T)-min(T);
    if Tdrift < toll
        disp(['kinetic energy conserved, max drift : ', num2str(Tdrift)])
    else
        disp(['kinetic energy NOT conserved, max drift : ', num2str(Tdrift)])
    end

end