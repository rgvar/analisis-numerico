function euler_adelante
    # explícito
    dt = 0.001;
    NDt = 1000;
    f = @(t,y) 2*t*y;

    t = zeros(NDt, 1);
    y = zeros(NDt, 1);

    t(1)=0;
    u(1)=1;
    for k=1: NDt - 1
        k1 = dt*f(t(k),u(k));
        t(k+1) = t(k) + dt;
        u(k+1) = u(k) + k1;
        fprintf("%d: t=%f, u=%f \n",k,t(k),u(k));
    endfor

    figure(1)
    plot(t,u,'b')
    grid on

endfunction
