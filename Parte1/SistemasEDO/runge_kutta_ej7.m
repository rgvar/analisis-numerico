function runge_kutta_ej7

    Dt=0.02;
    NDt=100;
    w=1;

    f = @(t,y) [ -10*y(1) + 4*y(2) ; -4*y(1)];

    t= zeros(1,NDt);
    y= zeros(2,NDt);

    t(1)=0;
    y(1,1)= 5;
    y(2,1)= 3;

    for j=1:NDt -1
        k1= Dt * f(t(j),y(:,j));
            tg= t(j) + Dt/(2*w);
            yg= y(:,j) + k1/(2*w);
        k2= Dt * f(tg, yg);
        t(j+1)= t(j) + Dt;
        y(:,j+1)= y(:,j) + (1-w)*k1 + w*k2;
    endfor


    figure(1)
    plot(t,y(1,:),'r')
    grid on


endfunction
