function global3k10
    clc,clear

    fs=16000;
    gr = load("audiovoz02.txt", "-ascii");
    N = length(gr)

    % calculo Dt
    Dt = 1/fs


    t = zeros(1,N);

    for j=1: N-1
        t(j+1) = t(j) + Dt;
    endfor

    figure(1)
    plot(t,gr,'b')
    grid on
    title("funciòn gr")

    Gr = fft(gr,N);

    for k=1:N
        Gr_mod(k) = abs(Gr(k));
    endfor

    figure(2)
    stem(Gr_mod,'bo-')
    grid on
    title("mòdulo tdf Gr")

    Dw = 2 * pi / (N * Dt)

    % Mòdulo Gr; Adm y fm
    [Adm, k_max] = max(Gr_mod(2:N/2))
    fm = (k_max - 1) * (fs / N)



    % EDO
    zita = 0.4;
    wn = fm;

    f = @(y,g) [ y(2) + g ;
                -(wn)^2 * y(1) + -2*zita*wn*y(2) ;
                y(1) + (-wn * y(3))];

    y = zeros(3,N);
    k1 = zeros(3,1);
    g = zeros(1,N);
    g(1,1) = 1;
    h= zeros(1,N);
    for j=1 : N-1
        k1 = Dt * f( y(:,j) , g(j) );
        y(:,j+1) = y(:,j) + k1(:,1);
    endfor

    h = y(3,:);

    figure(3)
    plot(t,h,'b')
    grid on
    title("h(t)")

    % TDF de h

    H = fft(h, N);
    H_mod = abs(H);

    Ah = max(H_mod)

    figure(4)
    stem(H_mod,'bo-')
    grid on
    title("modulo tdf H")

    % TDF YF

    y2 = zeros(3,N);
    for j=1: N-1
        k1= Dt * f (y(:,j) , gr(j) );
        y2(:,j+1) = y2(:,j) + k1(:,1);
    endfor

    yf = y2(1,:);

    YF = fft(yf,N);
    YF_mod = abs(YF);

    Ay = max(YF_mod)






endfunction
