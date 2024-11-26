function globalVargas
    clc, clear
    grafico=1;

    # Parámetros

    fs = 16000;
    Dt = 1/fs
    gr = load("audiovoz02.txt", "-ascii");
    N = length(gr)

    t = (1:N)*Dt;

    if (grafico==0)
        figure(1)
        plot(t,gr,'b')
        grid on
        title("gr(t)")
    endif

    # TDF de gr(t)

    Dw = (2*pi)/(N*Dt)

    Gr = fft(gr);
    Gr_mod = abs(Gr);

    [Adm, k_max] = max(Gr_mod);
    fm = (k_max-1)*Dw;
    Adm,fm

    if (grafico==0)
        figure(2)
        stem(Gr_mod(1:N/2), 'bo-')
        grid on
        title("Mod Gr(k)")
    endif

    ## FILTROS

    zita = 0.4;
    wn = fm;

    # filtro con entrada función pulso unitario

    x = zeros(3,N);
    dx_dt = zeros(3,1);
    k1 = zeros(3,1);

    % definición función delta de dirac
    g = zeros(1,N);
    g(1,1) = 1/Dt;

    for j=1:N-1
        dx_dt(1) = x(2,j) + g(1,j);
        dx_dt(2) = (-(wn^2))*x(1,j) + ((-2)*zita*wn)*x(2,j);
        dx_dt(3) = x(1,j) + (-wn)*x(3,j);
        k1(:) = Dt * dx_dt(:);
        x(:,j+1) = x(:,j) + k1(:);
    endfor

    h = x(3,:);

    if (grafico==0)
        figure(3)
        plot(t,h,'b')
        grid on
        title("h(t)")
    endif

    # TDF de h(t)

    H = fft(h);
    H_mod = abs(H);

    Ah = H_mod(k_max) % (fm/Dw)+1

    if (grafico==0)
        figure(4)
        stem(H_mod(1:N/2),'bo-')
        grid on
        title("Mod H(k)")
    endif

    # filtro con entrada función gr(t)
    x = zeros(3,N);
    dx_dt = zeros(3,1);
    k1 = zeros(3,1);

    for j=1:N-1
        dx_dt(1) = x(2,j) + gr(j);
        dx_dt(2) = (-(wn^2))*x(1,j) + ((-2)*zita*wn)*x(2,j);
        dx_dt(3) = x(1,j) + (-wn)*x(3,j);
        k1(:) = Dt * dx_dt(:);
        x(:,j+1) = x(:,j) + k1(:);
    endfor

    yf = x(3,:);

    if (grafico==0)
        figure(5)
        plot(t,yf,'b')
        grid on
        title("yf(t)")
    endif

    YF = fft(yf);
    YF_mod = abs(YF);

    Ay = YF_mod(k_max)

    if (grafico==1)
        figure(6)
        stem(YF_mod(1:N/2),'bo-')
        grid on
        title("Mod yf(k)")
    endif











endfunction
