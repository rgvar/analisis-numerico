function ejercicio_parcial_1

    fs=5000;
    Dt=1/fs;
    segundos=4;
    NDt=fs*segundos;
    w=0.5;

    f= @(t,y,b) [ 0*y(1) + 1*y(2) + 0*b ; -14400*y(1) - 2.4*y(2) + 1*b ];

    t=zeros(1,NDt);
    y=zeros(2,NDt);
    g=load("p1_3k9_01.txt","-ascii");

    t(1)=0;
    y(:,1)=0;

    for j=1:NDt-1
        k1= Dt * f(t(j), y(:,j), g(j));
            tg= t(j) + Dt/(2*w);
            yg= y(:,j) + k1/(2*w);
        k2= Dt * f(tg, yg, g(j));
        y(:,j+1)= y(:,j) + (1-w)*k1 + w*k2;
        t(j+1)= t(j) + Dt;
    endfor


    # función en un punto
    t1 = 1 * fs;
    fprintf("t = %.3f \n  x1(t1) = %.3f \n  x2(t1) = %.3f \n",t(t1),y(1,t1),y(2,t1));

    figure(1)
    plot(t,y(1,:),'r')
    grid on

    #integral trapecio

    # cálculo integral I
    FIt= @(y1,y2) Dt * (y1^2 + y2^2) / 2;
    # aplico derivada dx1(t)/dt (primera parte)
    f1 = @(y) 0*y(1) + 1*y(2);

    Igg=0;
    for j=1: NDt-1
        Igg= Igg + FIt(g(j),g(j+1));
    endfor

    I11=0;
    for j=1: NDt-1
        I11= I11 + FIt(y(1,j),y(1,j+1));
    endfor

    I22=0;
    for j=1: NDt-1
        I22= I22 + FIt(f1(y(:,j)),f1(y(:,j+1)));
    endfor

    fprintf("A22 = %f \n  Igg = %f \n  I11 = %f \n  I22 = %f",-1.2,Igg,I11,I22);

    # derivada en un punto

    # t=1s
    n=1*fs;
    df= -y(1,n-1)/(2*Dt) + y(1,n+1)/(2*Dt);
    df





endfunction
