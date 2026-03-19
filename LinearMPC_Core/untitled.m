Ac = [0 1;0 0];
Bc = [0;1];
Ac = blkdiag(Ac,Ac,Ac);
Bc = blkdiag(Bc,Bc,Bc);
Ts = 0.05;
ssd = c2d(ss(Ac,Bc,[],[]),Ts);
A = ssd.A;
B = ssd.B;
Kt = -place(A,B,[0.98,0.96,0.98,0.96,0.98,0.96]);
alpha = 3;
[V,beta,St]= MPCSim_tmnl(A,B,Kt,alpha);



