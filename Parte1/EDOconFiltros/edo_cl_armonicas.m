function edo_cl_armonicas

    p=1.2;
    fs=16000;
    Dt=1/fs;
    NDt=fs*6;

    g = @ (t) 4*sin(5*t) + 0.8*sin(120*t);
    f = @ (t,u) g(t) - p*u;

    t = zeros(NDt, 1);
    y = zeros(NDt, 1);
    h = zeros(NDt, 1);

    t(1) = 0;
    y(1) = 6;
    h(1) = 0;
    for j=1:NDt - 1
        k1 = Dt * f(t(j),y(j));
        t(j+1) = t(j) + Dt;
        y(j+1) = y(j) + k1;
        h(j+1) = g(t(j));
    endfor

    figure(1)
    plot(t,y,'r')
    hold on
    plot(t,h,'b')
    grid on;

endfunction
