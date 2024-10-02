

function parcial

    fs=2500;
    Dt=1/fs;
    segundos=10;
    NDt=fs*segundos;
    w=1;

    f = @(t,y,b) [ 0*y(1) + 1*y(2) + 0*b ; -25*y(1) - 1.6*y(2) + 1*b ];
    t= zeros(1,NDt);
    y= zeros(2,NDt);
    g= load("p1_3k10_01.txt","-ascii");

    t(1)=0;
    y(:,1)=0;
    for j=1:NDt-1
        k1=Dt*f(t(j),y(:,j),g(j));
            tg= t(j) + Dt/(2*w);
            yg= y(:,j) + k1/(2*w);
        k2=Dt*f(tg,yg(:),g(j));
        y(:,j+1)= y(:,j) + (1-w)*k1 + w*k2;
        t(j+1)= t(j) + Dt;
    endfor

    figure(1)
    plot(t,y(1,:),'b')
    grid on

    % integral Igg
    t1=9.5*fs+1

    t_t1= t(t1)
    x1_t1= y(1,t1)
    x2_t2= y(2,t1)
    Igg=0;
    for j=1:NDt
        Igg = Igg + Dt * (g(j)^2 + g(j+1)^2)/2;
    endfor

    Igg
# t1=9.5 segundos
#x1_t1 = 0.8727
#x2_t1 = -2.1021
#Igg = 101.55
#I12 = 0.4349

endfunction
