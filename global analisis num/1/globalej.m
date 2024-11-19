

function globalej
    clc,clear
    graficas=1;

    ro = 1.17;
    dt = 0.02;
    N = 2688;
    tf = (N - 1) * dt;

    f = @(y,g) [ -ro*y(1) + g ; y(1) + -ro*y(2) ; y(2) + -ro*y(3) ];

    # 1) EULER EXPLÍCITO

    t=zeros(1,N);
    y=zeros(3,N);
    k1=zeros(3,1);
    h3=zeros(1,N);

    g=ones(1,N);

    t(1)=0;
    y(:,1) = 0;
    h3(1)=0;

    for j=1 : N-1
        k1 = dt * f( y(:,j) , g(j) );
        t(1,j+1) = t(j) + dt;
        y(:,j+1) = y(:,j) + k1(:,1);
        h3(1,j+1) = y(2,j) + -ro*y(3,j);

    endfor

    if (graficas == 0 )
        figure(1)
        plot(t,y(3,:),'b')
        grid on
        title("x3(t)")

        figure(2)
        plot(t,h3(1,:),'b')
        grid on
        title("h3(t)")
    endif

    # 2) F DISCRETA
    gr = load("registro_231123.txt", "-ascii");

    if (graficas == 0)
        figure(3)
        plot(t,gr,'b')
        grid on
        title("gr(t)")
    endif

    # 3) TDF h3(t)

    H3 = fft(h3,N);

    for k=1:N
        H3_mod(k) = abs(H3(k));
    endfor

    if (graficas == 0)
        figure(4)
        stem(H3_mod(1:70),'bo-')
        grid on
        title('modulo tdf h3')
    endif


    # 4) TDF gr(t)

    Gr = fft(gr,N);

    for k=1:N
        Gr_mod(k) = abs(Gr(k));
    endfor

    if (graficas == 0)
        figure(5)
        stem(Gr_mod(1:N/2),'bo-')
        grid on
        title('modulo tdf gr')
    endif

    # 5) Ycf(k)

    for k=1:N
        Ycf(k) = Gr(k) * H3(k);
        Ycf_mod(k) = abs(Ycf(k));
    endfor

    if (graficas == 0)
        figure(6)
        stem(Ycf_mod(1:N/20),'bo-')
        grid on
        title('modulo Ycf')
    endif

    # 6) ycf(t)
    ycf = dt*ifft(Ycf);

    if (graficas ==1)
        figure(7)
        plot(t,ycf,'b')
        grid on
        title('modulo ycf')
    endif









endfunction
