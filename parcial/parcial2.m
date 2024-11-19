

function parcial2
    clc,clear
    Dt = 0.13464;
    A1 = 1/6;

    g = load("parcial_2_3k10_T1.txt");
    N = length(g);

    % f discreta

    for j=1:N
        t(j) = (j-1)*Dt;
    endfor

    figure(1)
    plot(t,g,'r-')
    grid on
    title('g(t)')

    tg_max = 0;
    g_max = 0;
    for j=1:N/2
        if g(j) > g_max
            tg_max = t(j);
            g_max = g(j);
        endif
    endfor
    tg_max,g_max

    % tdf

    G_tdf = fft(g,N);

    for k=1:N
        G_tdf_mod(k) = abs(G_tdf(k));
    endfor

    G_tdf_MAX = max(G_tdf_mod)

    figure(2)
    stem(G_tdf_mod(1:150),'bo-')
    grid on
    title('modulo tdf')

    G_tdf_kc=0;
    kc=0;
    for k=1:10
        if G_tdf_mod(k) > (0.1 * G_tdf_MAX)
            G_tdf_kc = G_tdf_mod(k);
            kc = k;
        endif
    endfor
    G_tdf_kc,kc

    % min2
    m=3;
    Tp = Dt*N;
    Dw = (2*pi)/Tp;
    w0 = (2*pi)/N;

    FI = zeros(N,m*2+1);

    FI(:, 1) = 1;
    for k=1:m
        for j=1:N
            FI(j,2*k)=cos(k*w0*(j-1));
            FI(j,2*k+1)=sin(k*w0*(j-1));
        endfor
    endfor

    D = diag(diag(FI'*FI));
    b = FI' * g;

    a=zeros(m*2+1,1);
    for k = 1:(2 * m) + 1
        a(k, 1) = b(k) / D(k, k);
    end

    Pa = FI*a;

    tPa_max=0;
    Pa_max=0;
    for i=1:(N/2)
        if(Pa(i)>Pa_max)
            Pa_max = Pa(i);
            tPa_max = t(i);
        endif
    endfor
    tPa_max,Pa_max

    figure(3)
    plot(t,Pa,'r-')
    grid on
    title('aprox min2')

    % convolucion

    kc = 7;
    Tp = Dt*N;
    Dw = (2*pi)/Tp;
    p=kc*Dw;

    for j=1:N
        hn(i)=A1*exp((-p)*t(j));
    endfor
    yc=Dt*conv(hn,g);

    figure(4)
    plot(t,yc(1:N),'g')
    grid on
    title('conv')

    tyc_max=0;
    yc_max=0;
    for j=1:N/2
        if yc(j) > yc_max
            tyc_max = t(j);
            yc_max = yc(j);
        endif
    endfor
    tyc_max,yc_max

    % grafica en conjunto

    figure(5)
    plot(t,g,'r',t,yc(1:N),'g',t,Pa,'b')
    grid on
    title('comparaciòn')


endfunction
