function global_2
    clc,clear

    fs = 16000;
    Dt = 1/fs;
    fprintf("Dt = %.8f \n",Dt)

    gr = load("audiovoz05.txt", "-ascii");
    N = length(gr)

    Dw = (2*pi) / (N*Dt)

    t = (1:N)*Dt;

    figure(1)
    plot(t,gr,'b')
    grid on
    title("gr(t)")

    Gr = fft(gr);
    Gr_mod = abs(Gr);

    figure(2)
    stem(Gr_mod(1:2500), 'bo-')
    grid on
    title("Mod Gr(k)")

    p = 200*Dw;

    for j=1:N-1
        h1(j) = exp(-p*t(j));
        h2(j) = t(j)*exp(-p*t(j));
    endfor

    h = Dt * conv(h1,h2);

    figure(3)
    plot(t,h(1:N),'r')
    grid on
    title("h(t)")
    xlim([0 0.02])

    H = fft(h);
    H_mod = abs(H);

    figure(4)
    stem(H_mod(1:500), 'bo-')
    grid on
    title("Mod H(k)")

    # yf(t)
    yf = Dt * conv(h,gr);

    figure(5)
    plot(t,yf(1:N),'r')
    grid on
    title("yf(t)")

    # YF(t)
    YF = fft(yf,N);
    YF_mod = abs(YF);

    figure(6)
    stem(YF_mod(1:2500), 'bo-')
    grid on
    title("Mod YF(k)")







endfunction
