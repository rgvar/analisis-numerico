

function runge_kutta

    # w = 1 euler modificado
    # w = 1/2 euler mejorado

    w=1/2;
    Dt = 0.01;
    NDt = 1000;

    f = @ (t,y) t/2 - y/2 ;

    y = zeros(NDt,1);
    t = zeros(NDt,1);

    y(1) = 4;
    t(1) = 0;
    for j=1: NDt - 1
        k1= Dt * f(t(j),y(j));
            tg = t(j) + Dt/(2*w);
            yg = y(j) + k1/(2*w);
        k2= Dt * f(tg, yg);
        y(j+1) = y(j) + (1-w)*k1 + w*k2;
        t(j+1) = t(j) + Dt;
        # fprintf("vuelta %d: t=%.4f;y=%f;k1=%f; tg=%f; yg=%f; k2=%f; \n",j,t(j),y(j),k1,tg,yg,k2);
    endfor

    figure(1)
    plot(t,y,'r')
    grid on


endfunction



