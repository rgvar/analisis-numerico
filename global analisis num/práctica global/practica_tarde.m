function practica_tarde
    clc,clear

    fs = 16000;
    Dt = 1/fs;
    fprintf("Dt = %.8f \n",Dt);

    gr = load("audiovoz05.txt","-ascii");
    N = length(gr)
    t = (0:N-1)*Dt;

    Dw = (2*pi) / (N*Dt);
    fprintf("Dw = %.3f rad/seg \n",Dw);

    Gr = fft(gr);
    Gr_mod = abs(Gr);

    ## gráfico gr(t) ##
    figure(1)
    plot(t,gr,'b')
    grid on
    title("gr(t)")
    xlim([0 0.7])

    ## gráfico módulo Gr(k) ##
    figure(2)
    stem(Gr_mod(1:2500),'bo-')
    grid on
    title("módulo Gr(k)")

    ## FILTRO ##

    p = 200*Dw;

    for j=1:N-1
        h1(j) = exp(-p*t(j));
        h2(j) = t(j)*exp(-p*t(j));
    endfor

    h = Dt * conv(h1,h2);

    ## gráfico gr(t) ##
    figure(3)
    plot(t,h(1:N),'r')
    grid on
    title("h(t)")
    xlim([0 0.02])

    # encontrar H(k) tdf de h(t)
    H = fft(h,N);
    H_mod = abs(H);

    ## gráfico módulo Gr(k) ##
    figure(4)
    stem(H_mod(1:500),'bo-')
    grid on
    title("módulo H(k)")

    ## encontrar yf(t)

    yf = Dt * conv(gr, h);

    ## gráfico gr(t) ##
    figure(5)
    plot(t,yf(1:N),'r')
    grid on
    title("yf(t)")
    xlim([0 0.7])

    # encontrar YF(k) tdf de yf(t)
    YF = fft(yf,N);
    YF_mod = abs(YF);

    ## gráfico módulo YF(k) ##
    figure(6)
    stem(YF_mod(1:2500),'bo-')
    grid on
    title("módulo YF(k)")





endfunction
