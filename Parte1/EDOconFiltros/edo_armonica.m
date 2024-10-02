function edo_armonica

    p = 1.2;
    Dt = 6.25E-5;
    fs = 16000;
    NDt = fs*6;

    #A=4;
    #w=5;
    A = 4/5;
    w = 120;


    f = @ (t,y) A*sin(w*t) - p * y;

    t = zeros(NDt, 1);
    y = zeros(NDt, 1);
    g = zeros(NDt, 1);

    t(1)=0;
    y(1)=6;
    #y(1)=0.1;
    g(1)=0;
    for j=1:NDt - 1
        k1 = Dt * f(t(j),y(j));
        t(j+1) = t(j) + Dt;
        y(j+1) = y(j) + k1;
        g(j+1) = A*sin(w*t(j));
    endfor

    figure(1)
    plot(t,y,'r')
    hold on
    plot(t,g,'b')
    grid on



endfunction
