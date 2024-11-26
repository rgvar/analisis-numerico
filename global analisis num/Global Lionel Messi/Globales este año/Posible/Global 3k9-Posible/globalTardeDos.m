function globalTardeDos
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
  h=Dt*conv(h1,h2);
  figure(1);
  plot(tr,g(1:N),'b');
  grid on;

G = fft(g,N);
for k=1:N
  Gmod (k) = abs(G(k));
  end

H= fft(h,N);

for k=1:N
  Hmod(k) = abs(H(k));
  end


for k=1:N
  X(k) = G(k) * H(k);
endfor

x = Dt * ifft(X,N);

for k=1:N
  Xmod (k) = abs(X(k));
endfor

figure(1)
plot(tr,x(1:N),'b');
title('Misma que 5');
grid on;

figure(2);
stem(Xmod(1:2500),'b');
title('Misma figura que 6')
grid on;

figure(3);
stem(Gmod(1:2500),'b');
title('Grafica de modulo de G');
grid on;

figure(4);
stem(Hmod(1:2500),'b');
title('Grafica del modulo de H');
grid on;
  endfunction
