

function runge_kutta_2

    Dt=1E-2;
    NDt= 1500;
    w=1;

    t= zeros(1, NDt);
    y= zeros(2, NDt);

    f = @(t,y) [ y(2) ; -4*y(1) + 10 * sin(3*t) ];

    t(1)= 0;
    y(:,1)= 0;

    for j=1: NDt -1
        k1= Dt*f(t(j),y(:,j));
            tg = t(j) + Dt/(2*w);
            yg = y(j) + k1/(2*w);
        k2 = Dt*f(tg,yg);
        t(j+1) = t(j) + Dt;
        y(:,j+1) = y(:,j) + (1-w)*k1 + w*k2;

    endfor

    figure(1)
    plot(y(1,:),y(2,:),'r')
    grid on

endfunction
