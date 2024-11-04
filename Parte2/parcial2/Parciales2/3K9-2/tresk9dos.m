


function tresk9dos
    clc,clear

    gn=load('parcial_2_3k9_2023.txt','-ascii');
    N= length(gn)
    Dt=0.02244
    w0=2*pi/N
    Tp=N*Dt
    Dw=w0/Dt

    %% FUNCIÓN DISCRETA

    for k=1:N
        tn(k)=Dt*(k-1);
    endfor

    figure(1)
    plot(tn,gn,'b')
    grid on;
    title('gn discreta');

    tg_max=0;
    g_max=0;
    for k=1:N/2
        gn_abs = abs(gn(k));
        if (gn_abs>g_max)
            tg_max=tn(k);
            g_max=gn_abs;
        endif
    endfor
    tg_max
    g_max

    %% TRANSFORMADA DE FOURIER

    G_tdf=fft(gn,N);
    for k=1:N
        G_mod(k)=abs(G_tdf(k));
    endfor
    G_tdf_max = max(G_mod)

    for k=1:10
        if G_mod(k) > (0.1*G_tdf_max)
            kc=k;
            G_tdf_kc=G_mod(k);
        endif
    endfor
    G_tdf_kc
    kc
    %EJERCICIO 3
    m=3;
    dim = 2*m+1;
    phi = zeros(N,dim);
    %%Frecuencias para cada Kw
    kw0(1) = 0;
    kw(1) = 0;
    %Calculo la matriz phi
    phi(:, 1) = 1;
    for k=1:m
    kw0(k + 1) = k * w0;
    kw(k + 1) = k * Dw;
    for j=1:N
    phi(j,2*k)=cos(k*w0*(j-1));
    phi(j,2*k+1)=sin(k*w0*(j-1));
    endfor
    endfor
    %Calculo los coeficientes alfa y b
    D = diag(diag(phi' * phi));
    b = phi' * gn;
    # vector a
    a=zeros(dim,1);
    for k = 1:(2 * m) + 1
    a(k, 1) = b(k) / D(k, k);
    end
    Pa=phi*a;
    tPa_max=0;
    Pa_max=0;
    for i=1:(N/2)
    if(Pa(i)>Pa_max)
    Pa_max = Pa(i);
    tPa_max = tn(i);
    endif
    endfor
    tPa_max
    Pa_max
    figure(3)
    plot(tn,Pa,'b')
    title('Aproximación con Pa')
    grid on
    title('figura 3')
    %EJERCICIO 4
    A1=1;
    p=kc*Dw;
    for k=1:N
    h(k)=A1*exp(-p*tn(k));
    endfor
    yc=Dt*conv(h,gn);
    yc_max=0;
    tyc_max=0;
    for k=1:N/2
    if (yc(k)>yc_max)
    yc_max=yc(k);
    tyc_max=tn(k);
    endif
    endfor
    yc_max
    tyc_max
    %Ejercicio 5
    figure(3)
    plot(tn,gn,'-r',tn,Pa,'b',tn,yc(1:N),'g')
    grid on;
    legend('g(t)','conv','Aprox Min')
endfunction
