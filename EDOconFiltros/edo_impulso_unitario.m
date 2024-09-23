function edo_impulso_unitario

    p = 1.2;
    Dt = 0.0001;
    fs = 50000;
    NDt = fs;
    f = @ (t,y) 1 - p * y;

    t = zeros(NDt,1);
    y = zeros(NDt,1);

    t(1) = 0;
    y(1) = 0;
    # aproximación EDO por método euler explícito (adelante)
    for j=1: NDt - 1
        k1 = Dt*f(t(j),y(j));
        y(j+1) = y(j) + k1;
        t(j+1) = t(j) + Dt;
    endfor

    # calculo derivada primera central
    # y' = (-y(s-1) + 0*y(s) + y(s+1)) / (2*Dt);
    h = zeros (NDt,1);
    for j=2: NDt - 1
        h(j)= (-y(j-1) + y(j+1))/(2*Dt);
    endfor

    figure(1)
    plot(t,y,'r')
    hold on
    plot(t,h,'b')
    grid on;

endfunction
