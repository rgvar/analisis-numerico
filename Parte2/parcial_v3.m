
function parcial3
    clc,clear

    Dt = 0.083776;
    g = load('registro_4_nov.txt');
    N = length(g)

    % f discreta
    t = zeros(N,1);
    for j=1:N
        t(j) = (j-1)*Dt;
    endfor

    norma=0;
    for j=1:N
        norma = norma + g(j)^2;
    endfor
    norma = norma^(1/2)

    figure(1)
    plot(t,g,'r-')
    grid on

    % tdf
    G_tdf = fft(g,N);

    for j=1:N
        mod_G(j) = abs(G_tdf(j));
    endfor

    Norma_2_G = 0;
    for k=1:N
        Norma_2_G = Norma_2_G + mod_G(k)^2;
    endfor
    Norma_2_G = Norma_2_G^(1/2);
    fprintf("Norma_2_G = %.4f \n", Norma_2_G);

    figure(2)
    stem(mod_G(1:110),'bo-')
    title('modulo de G')

    % convoluciòn
    Tp = Dt*N;
    Dw = 2*pi/Tp
    kc=8;
    A1 = 1/3;
    a = 0.3;
    w = kc*Dw;
    p = a * w;

    for j=1:N
        h(j) = A1 * (exp((-p)*t(j))/w) * sin(w*t(j));
    endfor

    yc = Dt * conv(h,g)

    figure(3)
    plot(t, yc(1:N), 'b-')
    grid on

    % grafico conjunto

    figure(40)
    plot(t,g,'r-',t,yc(1:N),'b-')
    grid ons



endfunction
