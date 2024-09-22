function filtro_edo

    p = -2;
    Dt = 0.1;
    fs = 5;
    NDt = fs;
    f= @ (t,y) 1 -(p * y);

    t = zeros(NDt,1);
    y = zeros(NDt,1);

    t(1) = 0;
    y(1) = 0;
    for j=1: NDt
        k1=Dt*f(t(j),y(j));
        y(j+1) = y(j) + k1;
        t(j+1) = t(j) + Dt;
    endfor

    figure(1)
    plot(t,y,'or')
    grid on

endfunction
