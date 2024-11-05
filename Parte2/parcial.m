


function parcial2
    clc,clear
    clf reset;

    g = load('registro-13-nov.txt');
    N = length(g);

    Dt = 0.16755;
    A1 = 1/4;

    % F discreta

    for j=1:N
        t(j) = (j-1)*Dt;
    endfor

    figure(1)
    plot(t,g,'b.')
    grid on
    title('g(t)');

    % Màximo g en t < t(N/3)

    t_max = 0;
    g_max = 0;
    for j=1:N/3
        g_abs(j) = abs(g(j));
        if (g_abs(j) > g_max)
            t_max = t(j);
            g_max = g_abs(j);
        endif
    endfor
    t_max,g_max

    % TDF y Mòdulo

    G_tdf = fft(g,N);
    for j=1: N
        G_tdf_mod(j) = abs(G_tdf(j));
    endfor
    G_max = max(G_tdf_mod)

    for k=1:10
        if G_tdf_mod(k) > (0.05 * G_max)
            G_tdf_kc = G_tdf_mod(k);
            kc = k;
        endif
    endfor
    kc, G_tdf_kc

    figure(2)
    plot(t,G_tdf_mod,'bo-')
    grid on
    title('mòdulo TDF G')

    % Aproximaciòn min2 en base trigonomètrica

    m = kc;
    w0 = (2*pi)/N;
    FI = zeros(N, m*2+1);
    FI(:,1) = 1;
    for j=1:N
        for k=1:m
            FI(j,2*k) = cos(k*w0*(j-1));
            FI(j,2*k+1) = sin(k*w0*(j-1));
        endfor
    endfor

    D = diag(diag(FI' * FI));
    b = FI' * g;

    a=zeros(2*m+1,1);
    for k=1:2*m+1
        a(k,1) = b(k) / D(k,k);
    endfor

    Pa = FI*a;

    tPa_max = 0;
    Pa_max = 0;
    for k=1:N/3
        if (Pa(k) > Pa_max)
            tPa_max = t(k);
            Pa_max = Pa(k);
        endif
    endfor
    tPa_max,Pa_max

    figure(3)
    plot(t,Pa,'r-')
    grid on
    title('Pa(t)')

    % Convoluciòn

    Tp = Dt * N;
    Dw = (2*pi) / Tp;
    p = kc * Dw;

    for i=1:N
        hn(i) = A1*exp((-p)*t(i));
    endfor

    yc = Dt * conv(hn,g);

    figure(4)
    plot(t,yc(1:N),'g')
    title('Convoluciòn')
    grid on

    tyc_max=0;
    yc_max=0;
    for i=1:N/3
        if (yc(i) > yc_max)
            tyc_max=t(i);
            yc_max=yc(i);
        endif
    endfor
    tyc_max,yc_max

    figure(10)
    plot(t,g,'b',t,Pa,'r-')
    grid on
    title('comparaciòn g(t), Pa(t), c(t)')

endfunction
