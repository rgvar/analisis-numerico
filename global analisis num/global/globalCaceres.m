function globalCaceres
    clc,clear
  g = load("audiovoz02.txt");
  N = length(g);
  fs = 16000;
  Dt = 1/fs;
  Dt
  t = (1:N-1)*Dt;
  w0 = 2*pi/N;
  Dw = w0/Dt;
  Dw
  w = (0:N/2 -1)*Dw;
  G = abs((fft(g,N)))(1:N/2);
  0.02
  (2*pi/12)/0.02


  [G_max, w_max] = max(G);
  G_max
  fm = w_max*Dw

  zita = 0.4;
  wn = fm;

  h1 = exp(-wn*t);
  h2 = (1/wn)*exp(-wn*zita*t).*sin(wn*t);
  h = Dt*conv(h1,h2);
  H = abs(fft(h,N));
  Ah = H(w_max)

  y = Dt*conv(h,g);
  Y = abs(fft(y, N));
  Ay = Y(w_max)

  figure(2)
  plot(t,h)

  h2 = [1,2,3,4,5,0,0,0,0,0,0,0];
  H2 = abs(fft(h2));
  H2

endfunction
