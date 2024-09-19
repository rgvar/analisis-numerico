function euler_adelante

    dt = 0.25;
    N = 1 / dt;
    f = @(t,y) 2*t*y;

    t(1)=0;
    u(1)=1;
    for k=2: N + 1
        t(k) = t(k-1) + dt;
        u(k) = u(k-1) + dt*f(t(k-1),u(k-1));
        % fprintf("%d: t=%f, u=%f \n",k,t(k),u(k));
    endfor

    figure(1)
    plot(t,u,'ob')
    grid on

endfunction
