function ejercicio_global
    clc,clear

    fs = 16000;
    Dt = 1/fs

    gr = load("audiovoz02.txt", "-ascii");
    N = length(gr)

    t = (0:N-1)*Dt;

    figure(1)
    plot(t,gr,'b')
    grid on
    title("gr(t)")

    w0 = 2*pi / N;
    Dw = w0 / Dt;

    Gr = fft(gr);
    Gr_mod = abs(Gr);

    figure(2)
    stem(Gr_mod(1:N/2), 'bo-')
    grid on
    title("Módulo de Gr(k)")

    [Adm, k] = max(Gr_mod(1:N/2));
    fm = (k-1) * Dw

    % filtro

    % EDO

    wn = fm;
    zita = 0.4;

    dx_dt = zeros(3, 1);
    y = zeros(3, N);

    g = zeros(1,N);
    g(1) = 1/Dt;

    for j=1:N-1
        dx_dt(1) = y(2,j) + g(j);
        dx_dt(2) = (-(wn)^2)*y(1,j) + (-2*zita*wn)*y(2,j);
        dx_dt(3) = y(1,j) + (-wn)*y(3,j);
        y(:,j+1) = y(:,j) + Dt * dx_dt(:);
    endfor

    h_t = y(3,:);

    figure(3)
    plot(t,h_t,'b')
    grid on
    title("h(t)")

    H_k = fft(h_t);
    H_k_mod = abs(H_k);
    Ah = H_k_mod(k)

    figure(4)
    stem(H_k_mod(1:N/2),'bo-')
    grid on
    title("Mod H(k)")

    # yf(t); Yf(t); Ay

    y_2 = zeros(3, N);
    dx_dt_2 = zeros(3,1);

    for j=1:N-1
        dx_dt_2(1) = y_2(2,j) + gr(j);
        dx_dt_2(2) = (-(wn)^2)*y_2(1,j) + (-2*zita*wn)*y_2(2,j);
        dx_dt_2(3) = y_2(1,j) + (-wn)*y_2(3,j);
        y_2(:,j+1) = y_2(:,j) + Dt * dx_dt_2(:);
    endfor

    yf = y_2(3,:);

    figure(5)
    plot(t,yf,'b')
    grid on
    title("yf(t)")

    Yf = fft(yf);
    Yf_mod = abs(Yf);
    Ay = Yf_mod(k)

    figure(6)
    stem(Yf_mod(1:N/2),"bo-")
    grid on
    title("Mod Yf(t)")

endfunction



