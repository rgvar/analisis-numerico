function globalTarde
  g=load('audiovoz05.txt','-ascii');
  N=length(g)
  t0=0;
  fs=16000;
  Dt=1/fs;
  w0 = 2*pi/N;
  Dw = (2*pi*fs)/N;
  p=200*Dw;
  for k=1:N
    tr(k)=t0+(k-1)*Dt;
  endfor
    for k=1:N
    h1(k) = exp(-p*tr(k));
    h2(k) = tr(k) * exp (-p*tr(k));
  endfor
  figure(1);
  plot(tr,g(1:N),'b');
  grid on;
  G=fft(g,N);
for k=1:N
  Gmod(k)=abs(G(k));
  end
  figure(2);
  stem(Gmod(1:2500),'b');
  grid on;
  title('Grafica del modulo del registro')
  h=Dt*conv(h1,h2);
  H=fft(h,N);
  for k=1:N
  Hmod(k)=abs(H(k));
end
figure(3)
plot(tr,h(1:N),'r');
grid on;
title('h(t) respuesta a impulso unitario del filtro')
  xlim([0 0.02]); % Limitar el eje x a 0-0.02 segundos
figure(4)
stem(Hmod(1:N),'b');
grid on;
title('Mod H')
  xlim([0 500]);
yf=Dt*conv(h,g)
YF=fft(yf);
for k=1:N
  YFmod(k)=abs(YF(k));
end
figure(5)
plot(tr,yf(1:N),'b');
grid on;
title('yf');


figure(6)
stem(YFmod(1:2500),'b');
grid on;
title('Modulo de YF');

  endfunction
