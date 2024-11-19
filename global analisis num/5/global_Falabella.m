function global_Falabella

  %1_LEER-----------------------------------------------------------------------

  g = load("registro_231121.txt");
  N = 2688; s = 1.17;  dt = 0.02; tf = (N-1)*dt;

  for i = 1:N
    t(i) = dt*(i-1);
  endfor

  figure(1)
  plot(t, g, 'b')
  grid on

  %2_RESOLVER-------------------------------------------------------------------

  x3_2(1,1) = 0;
  x3_2(2,1) = 0;
  x3_2(3,1) = 0;

  for i=1: N-1
    k1 = dt * f_pend(g(i), x3_2(:,i));
    x3_2(:,i+1) = x3_2(:,i) + k1;
  end

  figure(2)
  plot(t, x3_2(3,:) ,'b')
  grid on

  %3_OBTENER--------------------------------------------------------------------

  for i = 1:N
    h(i) = e^(-s*t(i));
  endfor

  figure(3)
  plot(t,h,'b')
  grid on

  %4_OBTENER--------------------------------------------------------------------

  for i = 1:N
    h2(i) = 0;
    for j = 1:i
      h2(i) = h2(i) + h(i+1-j)*h(j)*dt;
    endfor
  endfor

  for i = 1:N
    h3(i) = 0;
    for j = 1:i
      h3(i) = h3(i) + h2(i+1-j)*h(j)*dt;
    endfor
  endfor

  figure(4)
  plot(t,h3,'b')
  grid on

  %5_OBTENER--------------------------------------------------------------------

  for i = 1:N
    x3_5(i) = 0;
    for j = 1:i
      x3_5(i) = x3_5(i) + h3(i+1-j) * dt * g(j);
    endfor
  endfor

  figure(5)
  plot(t,x3_5,'b')
  grid on

  %6_COMPARAR-------------------------------------------------------------------
  for i=1:N
    y1(i) = (x3_2(i))^2;
  endfor

  for i=1:N
    y2(i)=(x3_2(i)-x3_5(i))^2;
  endfor

  %I1-------------------------------
  iTrap =0;
  for k = 1:N-1
    iTrap = iTrap +dt*(y1(k)+y1(k+1))/2;
  endfor

  i1(1) = -iTrap;
  for k = 2:N-1
    i1(k) = i1(k-1)+dt*(y1(k)+y1(k+1))/2;
  endfor

  I1=0;
  for i = 1:N-1
    I1 = I1 + i1(i);
  endfor

  %Id-------------------------------

  iTrap1 = 0;
  for k = 1:N-1
    iTrap1 = iTrap1 +dt*(y2(k)+y2(k+1))/2;
  endfor

  id(1) = -iTrap1;
  for k = 2:N-1
    id(k) = id(k-1)+dt*(y2(k)+y2(k+1))/2;
  endfor

  Id=0;
  for i = 1:N-1
    Id = Id + id(i);
  endfor

    disp(['I1 = ',num2str(I1,'%.2f')]);
	  disp(['Id = ',num2str(Id,'%.2f')]);

endfunction

function [fy] = f_pend(x,z)
  s = 1.17;
  fy(1,1) = -s * z(1) + x(1);
  fy(2,1) = 1*z(1) + (-s*z(2)) + 0;
  fy(3,1) = 1*z(2) + (-s*z(3)) + 0;
endfunction
