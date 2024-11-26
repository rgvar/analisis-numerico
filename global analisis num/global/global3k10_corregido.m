function global3k10
    clc,clear

    fs=16000;
    gr = load("audiovoz02.txt", "-ascii");

    # Dt, N, Dw
    Dt = 1/fs
    N = length(gr)
    Dw = 2*pi / (N*Dt)

    t = (0:N-1)*Dt;

    # GRÁFICO gr(t)
    figure(1)
    plot(t,gr,'b')
    grid on
    title("función gr")

    # Gr(k) = TDF{gr(t)} y módulo
    Gr = fft(gr,N);
    Gr_mod = abs(Gr);

    figure(2)
    stem(Gr_mod,'bo-')
    grid on
    title("módulo Gr")

    [Adm, k_max] = max(Gr_mod(1:N/2))
    fm = (k_max - 1) * (fs / N)

    % Sistema EDO
    zita = 0.4;
    wn = fm;


    ## Encontrar h(t)

    % definición función pulso unitario
    g = zeros(1,N);
    g(1) = 1;

    # resolución EDO con filtro, entrada pulso unitario (g(k))

    y = zeros(3,N);
    k1 = zeros(3,1);
    for j=1 : N-1
        x_dt(1) = y(2,j) + g(1,j);
        x_dt(2) = (-(wn)^2) * y(1) + ((-2)*zita*wn)* y(2) ;
        x_dt(3) =  1 * y(1) + (-wn) * y(3);

        k1 = Dt * x_dt(:);

        y(:,j+1) = y(:,j) + k1(:,1);
    endfor

    h = y;

    figure(3)
    plot(t,y,'b')
    grid on
    title("h(t)")

    % H(k) = TDF{h(k)}
    H = fft(h, N);
    H_mod = abs(H);

    figure(4)
    stem(H_mod,'bo-')
    grid on
    title("módulo tdf H")

    Ah = H_mod(k_max-1)

    % Encontrar yf(t)

    # resolución EDO con filtro, entrada función gr(k)
    y2 = zeros(3,N);

    for j=1: N-1
        k1 = Dt * 1;
        y2(:,j+1) = y2(:,j) + k1(:,1);
    endfor

    yf = y2(3,:);

    figure(5)
    plot(t,yf,'b')
    grid on
    title("yf(t)")

    % YF(k) = TDF{yf(k)}

    YF = fft(yf,N);
    YF_mod = abs(YF);

    Ay = YF_mod(k_max-1)

    figure (6)
    stem(YF_mod,'bo-')
    grid on
    title("módulo tdf YF")


endfunction



