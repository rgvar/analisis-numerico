
function integral_primitiva(N,w)
    t0=0;
    tf=5/w;
    dt=(tf-t0) / N;
    dim = N+1;

    for k=1:dim
        t(k) = (k-1) * dt;
        y(k) = e^(-w*t(k));
    endfor


	figure(1)
    plot(t,y,'ob')
    grid on

	yp(1) = 0;
	for k=1:N
		yp(k+1) = yp(k) + dt * (y(k) + y(k+1)) / 2;
	endfor


	figure(2)
	plot(t,yp, 'ob')
	grid on

	format long
	t = t(dim)
	y = y(dim)
	yp = yp(dim)

endfunction
