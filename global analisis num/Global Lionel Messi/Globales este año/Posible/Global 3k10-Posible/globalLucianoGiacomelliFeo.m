function globalLucianoGiacomelliFeo
  g=load('audiovoz02.txt','-ascii');
  N=length(g)
  t0=0;
  fs=16000;
  Dt=1/fs
  w0 = 2*pi/N;
  Dw = 2*pi / (N*Dt)
  for k=1:N
    tr(k)=t0+(k-1)*Dt;
  endfor
  figure(1)
  plot(tr,g(1:N),'b');
  title('Grafica del txt')
  printf("El valor de Dt = %.3f\n", Dt*10^4);
  printf("El valor de Dw = %.3f\n", Dw);
  Gr=fft(g,N);
    for k=1:N
    Gmod(k)=abs(Gr(k));
  endfor
[Adm,max_index]=max(Gmod);
fm=(max_index-1)*Dw;
printf("El valor de Adm = %.2f\n", Adm);
printf("El valor de fm = %.2f Hz\n", fm);
  figure(2);
  stem(Gmod(1:N),'b');
  grid on;
  title('Grafica del modulo del registro')
  % Parámetros iniciales
zeta = 0.4;
wn = fm;
N = length(g);
Dt = 1 / fs;
t(1) = 0;
x1(1) = 0;
x2(1) = 0;
x3(1) = 0;
h_t = zeros(1, N);
escalon = zeros(N);

    %es raro, esta es la escalonada
    for i = 1:N

        if tr(i) < N / (fs * 2)
            escalon(i) = 0;
        else
            escalon(i) = 1;
        end

    end



% Método de Euler explícito
for j = 1:N
    t(j+1) = t(j) + Dt;
    x1(j+1) = x1(j) + Dt * (x2(j) + escalon(j));
    x2(j+1) = x2(j) + Dt * ((-(wn*wn))*x1(j) - 2*zeta*wn*x2(j));
    x3(j+1) = x3(j) + Dt * (1*x1(j) - wn*x3(j));
    h_t(j) = x3(j);  %Salida de la funcion
endfor
derivada = zeros(N, 1);

    for i = 1:N

        switch (i)
            case 1
                derivada(1) = (1 / (2 * Dt)) * (-3 * h_t(1) + 4 * h_t(2) - h_t(3));
            case N
                derivada(N) = (1 / (2 * Dt)) * (3 * h_t(N) - 4 * h_t(N - 1) + h_t(N - 2));
            otherwise
                derivada(i) = (1 / (2 * Dt)) * (-h_t(i - 1) + h_t(i + 1));

        end

    end
figure(3)
plot(tr,derivada,'b')
title('Grafico de h')
HK=fft(derivada,N);
    for k=1:N
    Hmod(k)=abs(HK(k));
  endfor
    % Graficar el módulo de H(k)
    figure(4);
    stem(Hmod(1:N),'b');
    title('Módulo de la Transformada de Fourier de h(t)');
    grid on;
    Ah = Hmod(max_index) %Necesito sacar este modulo con el indice que obtuve para sacar la amplitud maxima
  printf("El valor de Ah en f_m = %.3f\n", Ah * 10^6);
    yf = zeros(1, N);
for j = 1:N
    t(j+1) = t(j) + Dt;
    x1(j+1) = x1(j) + Dt * (x2(j) + g(j));
    x2(j+1) = x2(j) + Dt * ((-(wn*wn))*x1(j) - 2*zeta*wn*x2(j));
    x3(j+1) = x3(j) + Dt * (1*x1(j) - wn*x3(j));
    yf(j) = x3(j);  %Salida de la funcion
endfor
  % Graficar yf(t_n)
  figure(5);
  plot(tr, yf, 'b');
  title('Función Filtrada yf(t)');
  xlabel('Tiempo [s]');
  ylabel('Amplitud');
  grid on;
  % Transformada de Fourier de yf(t_n)
  YF = fft(yf);
  for k=1:N
  YFmod(k) = abs(YF(k));
  endfor
  % Graficar el módulo de YF(k)
  figure(6);
  stem(YFmod, 'b');
  title('Módulo de la Transformada de Fourier de yf(t)');
  grid on;
 %  Amplitud en f_m para YF(k)
  Ay = YFmod(max_index)
  printf("El valor de Ay en f_m = %.2f\n", Ay * 10^8);
  %Elegi el metodo de euler explicito para resolver las EDO porque es simple, eficiente y adecuado para sistemas discretos. Y no tiene tanto error del resultado original
  endfunction
