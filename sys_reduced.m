function dxdt = sys_reduced(t,x,Kp,phi1,KI,phi2,om0,om1,om2,uf,ug,uh,alpha)

mu0_tmp = @(s0,s2) s0./(1+s0).*s2./(Kp+s2);
mu1_tmp = @(s1,s2) phi1*s1./(1+s1).*1./(1+KI*s2);
mu2_tmp = @(s2) phi2*s2./(1+s2);

dxdt = [x(1)*(mu0_tmp(-x(1)+uf,-om2*x(1)+om1*x(2)-x(3)+uh)-alpha); x(2)*(mu1_tmp(om0*x(1)-x(2)+ug,-om2*x(1)+om1*x(2)-x(3)+uh)-alpha); x(3)*(mu2_tmp(-om2*x(1)+om1*x(2)-x(3)+uh)-alpha)];

end
