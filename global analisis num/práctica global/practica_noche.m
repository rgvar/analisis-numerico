function practica_noche
    clc,clear

    gr = load("audiovoz02.txt", "-ascii");
    N = length(gr)

    fs = 16000;
    Dt = 1/fs;
    fprintf("Dt = %.3f*10^(-4) \n",(Dt*10^4))

    t = (0:N-1)*Dt;

    ## gráfico gr(t) ###
    figure(1)
    plot(t,gr,'b')
    grid on
    title("gr(t)")
    ####

    dw = 2*pi/(N*Dt);
    fprintf("dw = %.3f \n",dw)

    ### Gr(k) TDF de gr(t)

    Gr = fft(gr,N);
    Gr_mod = abs(Gr);

    ## gráfico Módulo Gr(k) ###
    figure(2)
    stem(Gr_mod(1:N/2),'bo-')
    grid on
    title("Mod Gr(k)")
    ####

    ## Amplitud máxima Adm y fm
    [Adm, k] = max(Gr_mod);
    fm = (k-1)*dw;
    fprintf("Adm = %.2f \n",Adm);
    fprintf("fm = %.2f \n",fm);

    ###### FILTRO #######

    zita = 0.4;
    wn = fm;

    # impulso unitario
    g = zeros(1,N);
    g(1,1) = 1/Dt;

    # euler explícito para resolver la EDO y encontrar h(t) con entrada impulso unitario
    dx_dt = zeros(3,1);
    k1 = zeros(3,1);
    x = zeros(3,N);
    for j=1:N-1
        dx_dt(1) = x(2,j) + g(1,j);
        dx_dt(2) = (-(wn)^2)*x(1,j) + (-2*zita*wn)*x(2,j);
        dx_dt(3) = x(1,j) + (-wn)*x(3,j);
        k1(:) = Dt * dx_dt(:);

        x(:,j+1) = x(:,j) + k1(:);
    endfor

    h = x(3,:);

    ## gráfico h(t) ##
    figure(3)
    plot(t,h,'b')
    grid on
    title("h(t)")
    ####

    ### H(k) tdf de h(t)
    H = fft(h);
    H_mod = abs(H);

    ## gráfico módulo H(k) ##
    figure(4)
    stem(H_mod(1:N/2),'bo-')
    grid on
    title("Mód H(k)")
    ####

    ## admplitud Ah para frecuencia fm
    Ah = H_mod(k)
    fprintf("Ah = %.3f*10^(-6) \n",(Ah*10^6));

    # euler explícito para resolver la EDO y encontrar yf(t) con entrada función gr(t)
    dx_dt = zeros(3,1);
    k1 = zeros(3,1);
    x_2 = zeros(3,N);
    for j=1:N-1
        dx_dt(1) = x_2(2,j) + gr(j);
        dx_dt(2) = (-(wn)^2)*x_2(1,j) + (-2*zita*wn)*x_2(2,j);
        dx_dt(3) = x_2(1,j) + (-wn)*x_2(3,j);
        k1(:) = Dt * dx_dt(:);

        x_2(:,j+1) = x_2(:,j) + k1(:);
    endfor

    yf = x_2(3,:);

    ## gráfico yf(t) ##
    figure(5)
    plot(t,yf,'b')
    grid on
    title("yf(t)")
    ####

    ### YF(k) tdf de yf(t)
    YF = fft(yf);
    YF_mod = abs(YF);

    ## gráfico módulo YF(k) ##
    figure(6)
    stem(YF_mod(1:N/2),'bo-')
    grid on
    title("Mód YF(k)")
    ####

    ## admplitud Ay para frecuencia fm
    Ay = YF_mod(k)
    fprintf("Ay = %.2f*10^(-8) \n",(Ay*10^8));




endfunction



















