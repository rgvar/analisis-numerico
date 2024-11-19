function global_3
    clc,clear
    graficar=1;

    ro=1.17;
    dt=0.02;
    N=2688;
    tf=(N-1)*dt;

    # 1) F DISCRETA

    g=load("registro_231123.txt","-ascii");

    t = zeros(1,N);
    for j=1:N-1
        t(1,j+1) = t(j) + dt;
    endfor

    if (graficar == 0)
        figure(1)
        plot(t,g,'b')
        grid on
    endif

    # 2) EULER EXPLÍCITO

    gu = ones(1,N);

    f = @ (y,g) [ -ro*y(1) + g ; y(1) + -ro*y(2) ; y(2) + -ro*y(3) ];

    y=zeros(3,N);
    k1=zeros(3,1);

    for j=1: N-1
        k1(:,1) = dt * f(y(:,j), g(j));
        y(:,j+1) = y(:,j) + k1(:,1);
    endfor

    if (graficar==0)
        figure(2)
        plot(t,y(3,:),'b')
        grid on
    endif

    # 3) F DISCRETA h(t)
    # h(t) = e^(-ro*t)

    for j=1: N
        h(j) = exp(-ro*t(j));
    endfor

    if (graficar==0)
        figure(3)
        plot(t,h,'b')
        grid on
    endif

    # 4) CONVOLUCIÓN : h3(t)
    # h3(t) = h3(t) o g(t)

    h3_ = dt * conv(h,h);
    h3 = dt * conv(h,h3_);

    if (graficar==0)
        figure(4)
        plot(t,h3(1:N),'b')
        grid on
    endif

    # 5) CONVOLUCIÓN : x3(t)
    x3 = dt * conv(h3,g);

    if (graficar==1)
        figure(5)
        plot(t,x3(1:N),'b')
        grid on
    endif

    # 6) TRAPECIOS

    x3_2 = y(3,:);
    x3_5 = x3;

    I1 = 0;
    Id = 0;
    for j=1:N-1
        I1 = I1 + dt * (x3_2(j)^2 + x3_2(j+1)^2) / 2;
        Id = Id + dt * ((x3_2(j) - x3_5(j))^2 + (x3_2(j+1) - x3_5(j+1))^2) / 2;
    endfor
    format long
    I1,Id






endfunction
