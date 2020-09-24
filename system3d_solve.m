clear;
% close all;
%parameters
Y0=0.019;
Y1=0.04;
Y2=0.06;
km0=29;
km1=26;
km2=35;
KS0=0.053;
KS1=0.302;
KS2=2.5e-5;
L0=1e-6;
KI2=3.5e-6;
om0=(KS0/KS1)*(224/208)*(1-Y0);
om1=(KS1/KS2)*(32/224)*(1-Y1);
om2=(16/208)*(KS0/KS2);
phi1=(km1/km0)*(Y1/Y0);
phi2=(km2/km0)*(Y2/Y0);
Kp=L0/KS2;
KI=KS2/KI2;

% om0=0.1854055;
% om1=1656.6857;
% om2=163.07692;
% phi1=1.8874773;
% phi2=3.8112523;
% Kp=0.04;
% KI=7.1428571;


mu0 = @(s0,s2) s0./(1+s0).*s2./(Kp+s2);
mu1 = @(s1,s2) phi1*s1./(1+s1).*1./(1+KI*s2);
mu2 = @(s2) phi2*s2./(1+s2);

% % Before Hopf (in uf)
% uf=0.538;
% ug=0;
% uh=0.1;
% alpha=0.01;
% 
% x0_int = 1 + uf + 1/(Kp*(phi2-alpha)+alpha-1);
% x1_int = om0*x0_int + ug + 1 + phi1/(alpha*(1+KI*alpha/(phi2-alpha))-phi1);
% x2_int = -om2*x0_int + om1*x1_int + uh - alpha/(phi2-alpha);
% 
% % Checking the value at the int equilibrium
% mu0(-x0_int+uf,-om2*x0_int+om1*x1_int-x2_int+uh);
% mu1(om0*x0_int-x1_int+ug,-om2*x0_int+om1*x1_int-x2_int+uh);
% mu2(-om2*x0_int+om1*x1_int-x2_int+uh);
% 
% options = odeset('RelTol',1e-8,'AbsTol',1e-8);
% ic_per = [0.33; 0.0528; 5];
% % 0.344611244,0.058465896,40.75884261,
% % ic_per=[0.3806711559962231; 0.06638815809539414; 48.00340445405122];
% ic_eq = [0.55;0.1;0.1];
% [t1,x_per] = ode15s(@(t,X) sys_reduced(t,X,Kp,phi1,KI,phi2,om0,om1,om2,uf,ug,uh,alpha),[0 30000],ic_per,options);
% [t2,x_eq] = ode15s(@(t,X) sys_reduced(t,X,Kp,phi1,KI,phi2,om0,om1,om2,uf,ug,uh,alpha),[0 500],ic_eq,options);
% [t3,x_per_orb] = ode15s(@(t,X) sys_reduced(t,X,Kp,phi1,KI,phi2,om0,om1,om2,uf,ug,uh,alpha),[0 10000],[x_per(end,1);x_per(end,2);x_per(end,3)],options);
% figure(1);
% plot3(x_per(:,1),x_per(:,2),x_per(:,3),'LineWidth',2);
% hold on;
% plot3(x_eq(:,1),x_eq(:,2),x_eq(:,3),'LineWidth',2);
% plot3(x_per_orb(:,1),x_per_orb(:,2),x_per_orb(:,3),'k','LineWidth',3);
% plot3(x0_int,x1_int,x2_int,'g*','LineWidth',3);
% scatter3(0.33,0.0528,5,'r','filled');
% scatter3(0.55,0.1,0.1,'r','filled');
% disp([x_per(end,1),x_per(end,2),x_per(end,3)]);
% xlabel('x0');
% ylabel('x1');
% zlabel('x2');
% legend('IC: [0.33, 0.0528, 5]', 'IC: [0.55, 0.1, 0.1]','Periodic orbit','Interior equilibrium','Location','Best');
% grid on;

% %After Hopf (in uf)
% uf=0.6;
% ug=0;
% uh=0.1;
% alpha=0.01;
% 
% x0_int = 1 + uf + 1/(Kp*(phi2-alpha)+alpha-1);
% x1_int = om0*x0_int + ug + 1 + phi1/(alpha*(1+KI*alpha/(phi2-alpha))-phi1);
% x2_int = -om2*x0_int + om1*x1_int + uh - alpha/(phi2-alpha);
% 
% % Checking the value at the int equilibrium
% mu0(-x0_int+uf,-om2*x0_int+om1*x1_int-x2_int+uh);
% mu1(om0*x0_int-x1_int+ug,-om2*x0_int+om1*x1_int-x2_int+uh);
% mu2(-om2*x0_int+om1*x1_int-x2_int+uh);
% 
% options = odeset('RelTol',1e-6,'AbsTol',1e-8);
% ic_per = [0.33; 0.0528; 5];
% % ic_eq = [0.1;0.01;0.3];
% ic_eq = [0.002106;0.008478;13.81];
% [t1,x_per] = ode15s(@(t,X) sys_reduced(t,X,Kp,phi1,KI,phi2,om0,om1,om2,uf,ug,uh,alpha),[0 4000],ic_per,options);
% [t2,x_eq] = ode15s(@(t,X) sys_reduced(t,X,Kp,phi1,KI,phi2,om0,om1,om2,uf,ug,uh,alpha),[0 5000],ic_eq,options);
% plot3(x_per(:,1),x_per(:,2),x_per(:,3),'LineWidth',2);
% hold on;
% plot3(x_eq(:,1),x_eq(:,2),x_eq(:,3),'LineWidth',2);
% plot3(x0_int,x1_int,x2_int,'g*','LineWidth',3);
% scatter3(0.33,0.0528,5,'r','filled');
% scatter3(0.002106,0.008478,13.81,'r','filled');
% % scatter3(x_per(end,1),x_per(end,2),x_per(end,3),'g','filled','LineWidth',10);
% disp([x_per(end,1),x_per(end,2),x_per(end,3)]);
% disp([x_eq(end,1),x_eq(end,2),x_eq(end,3)]);
% xlabel('x0');
% ylabel('x1');
% zlabel('x2');
% legend('IC: [0.33, 0.0528, 5]', 'IC: [0.002106, 0.008478, 13.81]','Interior equilibrium','Location','Best');
% grid on;
% hold off;

% --------------------------------------------------------------------------------------------------------------------------------

% % Before Hopf (in ug)
% uf=0.5;
% ug=0.0006;
% uh=0.1;
% alpha=0.01;
% 
% x0_int = 1 + uf + 1/(Kp*(phi2-alpha)+alpha-1);
% x1_int = om0*x0_int + ug + 1 + phi1/(alpha*(1+KI*alpha/(phi2-alpha))-phi1);
% x2_int = -om2*x0_int + om1*x1_int + uh - alpha/(phi2-alpha);
% 
% % Checking the value at the int equilibrium
% mu0(-x0_int+uf,-om2*x0_int+om1*x1_int-x2_int+uh);
% mu1(om0*x0_int-x1_int+ug,-om2*x0_int+om1*x1_int-x2_int+uh);
% mu2(-om2*x0_int+om1*x1_int-x2_int+uh);
% 
% options = odeset('RelTol',1e-4,'AbsTol',1e-8);
% ic_per = [0.5; 0.05; 1];
% ic_eq = [0.5;0.12;1];
% [t1,x_per] = ode15s(@(t,X) sys_reduced(t,X,Kp,phi1,KI,phi2,om0,om1,om2,uf,ug,uh,alpha),[0 50000],ic_per,options);
% [t2,x_eq] = ode15s(@(t,X) sys_reduced(t,X,Kp,phi1,KI,phi2,om0,om1,om2,uf,ug,uh,alpha),[0 5000],ic_eq,options);
% [t3,x_per_orb] = ode15s(@(t,X) sys_reduced(t,X,Kp,phi1,KI,phi2,om0,om1,om2,uf,ug,uh,alpha),[0 10000],[x_per(end,1);x_per(end,2);x_per(end,3)],options);
% plot3(x_per(:,1),x_per(:,2),x_per(:,3),'LineWidth',2);
% hold on;
% plot3(x_eq(:,1),x_eq(:,2),x_eq(:,3),'LineWidth',2);
% plot3(x_per_orb(:,1),x_per_orb(:,2),x_per_orb(:,3),'k','LineWidth',3);
% plot3(x0_int,x1_int,x2_int,'g*','LineWidth',3);
% scatter3(0.5,0.05,1,'r','filled');
% scatter3(0.5,0.12,1,'r','filled');
% disp([x_per(end,1),x_per(end,2),x_per(end,3)]);
% xlabel('x0');
% ylabel('x1');
% zlabel('x2');
% legend('IC: [0.5, 0.05, 1]', 'IC: [0.5, 0.12, 1]','Periodic orbit','Interior equilibrium','Location','Best');
% grid on;
% hold off;

% %After Hopf (in ug)
% uf=0.5;
% ug=0.0008;
% uh=0.1;
% alpha=0.01;
% 
% x0_int = 1 + uf + 1/(Kp*(phi2-alpha)+alpha-1);
% x1_int = om0*x0_int + ug + 1 + phi1/(alpha*(1+KI*alpha/(phi2-alpha))-phi1);
% x2_int = -om2*x0_int + om1*x1_int + uh - alpha/(phi2-alpha);
% 
% % Checking the value at the int equilibrium
% mu0(-x0_int+uf,-om2*x0_int+om1*x1_int-x2_int+uh);
% mu1(om0*x0_int-x1_int+ug,-om2*x0_int+om1*x1_int-x2_int+uh);
% mu2(-om2*x0_int+om1*x1_int-x2_int+uh);
% 
% options = odeset('RelTol',1e-4,'AbsTol',1e-8);
% ic_per = [0.5; 0.05; 1];
% ic_eq = [0.5;0.12;1];
% [t1,x_per] = ode15s(@(t,X) sys_reduced(t,X,Kp,phi1,KI,phi2,om0,om1,om2,uf,ug,uh,alpha),[0 40000],ic_per,options);
% [t2,x_eq] = ode15s(@(t,X) sys_reduced(t,X,Kp,phi1,KI,phi2,om0,om1,om2,uf,ug,uh,alpha),[0 5000],ic_eq,options);
% plot3(x_per(:,1),x_per(:,2),x_per(:,3),'LineWidth',2);
% hold on;
% plot3(x_eq(:,1),x_eq(:,2),x_eq(:,3),'LineWidth',2);
% plot3(x0_int,x1_int,x2_int,'g*','LineWidth',3);
% scatter3(0.5,0.05,1,'r','filled');
% scatter3(0.5,0.12,1,'r','filled');
% disp([x_per(end,1),x_per(end,2),x_per(end,3)]);
% disp([x_eq(end,1),x_eq(end,2),x_eq(end,3)]);
% xlabel('x0');
% ylabel('x1');
% zlabel('x2');
% legend('IC: [0.5, 0.05, 1]', 'IC: [0.5, 0.12, 1]','Interior equilibrium','Location','Best');
% grid on;
% hold off;

% --------------------------------------------------------------------------------------------------------------------------------

% % Before Hopf (in uh)
% uf=0.5;
% ug=0.0006;
% uh=0.05;
% alpha=0.01;
% 
% x0_int = 1 + uf + 1/(Kp*(phi2-alpha)+alpha-1);
% x1_int = om0*x0_int + ug + 1 + phi1/(alpha*(1+KI*alpha/(phi2-alpha))-phi1);
% x2_int = -om2*x0_int + om1*x1_int + uh - alpha/(phi2-alpha);
% 
% % Checking the value at the int equilibrium
% mu0(-x0_int+uf,-om2*x0_int+om1*x1_int-x2_int+uh);
% mu1(om0*x0_int-x1_int+ug,-om2*x0_int+om1*x1_int-x2_int+uh);
% mu2(-om2*x0_int+om1*x1_int-x2_int+uh);
% 
% options = odeset('RelTol',1e-4,'AbsTol',1e-8);
% ic_per = [0.5; 0.07; 1];
% ic_eq = [0.5;0.1;1];
% [t1,x_per] = ode15s(@(t,X) sys_reduced(t,X,Kp,phi1,KI,phi2,om0,om1,om2,uf,ug,uh,alpha),[0 50000],ic_per,options);
% [t2,x_eq] = ode15s(@(t,X) sys_reduced(t,X,Kp,phi1,KI,phi2,om0,om1,om2,uf,ug,uh,alpha),[0 5000],ic_eq,options);
% [t3,x_per_orb] = ode15s(@(t,X) sys_reduced(t,X,Kp,phi1,KI,phi2,om0,om1,om2,uf,ug,uh,alpha),[0 10000],[x_per(end,1);x_per(end,2);x_per(end,3)],options);
% plot3(x_per(:,1),x_per(:,2),x_per(:,3),'LineWidth',2);
% hold on;
% plot3(x_eq(:,1),x_eq(:,2),x_eq(:,3),'LineWidth',2);
% plot3(x_per_orb(:,1),x_per_orb(:,2),x_per_orb(:,3),'k','LineWidth',3);
% plot3(x0_int,x1_int,x2_int,'g*','LineWidth',3);
% scatter3(0.5,0.07,1,'r','filled');
% scatter3(0.5,0.1,1,'r','filled');
% disp([x_per(end,1),x_per(end,2),x_per(end,3)]);
% xlabel('x0');
% ylabel('x1');
% zlabel('x2');
% legend('IC: [0.5, 0.07, 1]', 'IC: [0.5, 0.1, 1]','Periodic orbit','Interior equilibrium','Location','Best');
% grid on;
% hold off;

%After Hopf (in uh)
uf=0.5;
ug=0.0006;
uh=0.3;
alpha=0.01;

x0_int = 1 + uf + 1/(Kp*(phi2-alpha)+alpha-1);
x1_int = om0*x0_int + ug + 1 + phi1/(alpha*(1+KI*alpha/(phi2-alpha))-phi1);
x2_int = -om2*x0_int + om1*x1_int + uh - alpha/(phi2-alpha);

% Checking the value at the int equilibrium
mu0(-x0_int+uf,-om2*x0_int+om1*x1_int-x2_int+uh);
mu1(om0*x0_int-x1_int+ug,-om2*x0_int+om1*x1_int-x2_int+uh);
mu2(-om2*x0_int+om1*x1_int-x2_int+uh);

options = odeset('RelTol',1e-4,'AbsTol',1e-8);
ic_per = [0.5; 0.07; 1];
ic_eq = [0.5;0.1;1];
[t1,x_per] = ode15s(@(t,X) sys_reduced(t,X,Kp,phi1,KI,phi2,om0,om1,om2,uf,ug,uh,alpha),[0 50000],ic_per,options);
[t2,x_eq] = ode15s(@(t,X) sys_reduced(t,X,Kp,phi1,KI,phi2,om0,om1,om2,uf,ug,uh,alpha),[0 5000],ic_eq,options);
plot3(x_per(:,1),x_per(:,2),x_per(:,3),'LineWidth',2);
hold on;
plot3(x_eq(:,1),x_eq(:,2),x_eq(:,3),'LineWidth',2);
plot3(x0_int,x1_int,x2_int,'g*','LineWidth',3);
scatter3(0.5,0.07,1,'r','filled');
scatter3(0.5,0.1,1,'r','filled');
disp([x_per(end,1),x_per(end,2),x_per(end,3)]);
disp([x_eq(end,1),x_eq(end,2),x_eq(end,3)]);
xlabel('x0');
ylabel('x1');
zlabel('x2');
legend('IC: [0.5, 0.07, 1]', 'IC: [0.5, 0.1, 1]','Interior equilibrium','Location','Best');
grid on;
hold off;

% --------------------------------------------------------------------------------------------------------------------------------

% uf=0.5;
% ug=0.0006;
% uh=0.3;
% alpha=0.02;
% 
% x0_int = 1 + uf + 1/(Kp*(phi2-alpha)+alpha-1);
% x1_int = om0*x0_int + ug + 1 + phi1/(alpha*(1+KI*alpha/(phi2-alpha))-phi1);
% x2_int = -om2*x0_int + om1*x1_int + uh - alpha/(phi2-alpha);
% 
% % Checking the value at the int equilibrium
% mu0(-x0_int+uf,-om2*x0_int+om1*x1_int-x2_int+uh);
% mu1(om0*x0_int-x1_int+ug,-om2*x0_int+om1*x1_int-x2_int+uh);
% mu2(-om2*x0_int+om1*x1_int-x2_int+uh);
% 
% options = odeset('RelTol',1e-4,'AbsTol',1e-8);
% % ic_per = [0.5; 0.07; 1];
% % ic_per = [0.2786;0.0352;13.1317];
% % ic_per = [0.296915; 0.03004; 0];
% ic_per = [0.4775; 0.0471; 0];
% ic_eq = [0.5;0.1;1];
% [t1,x_per] = ode15s(@(t,X) sys_reduced(t,X,Kp,phi1,KI,phi2,om0,om1,om2,uf,ug,uh,alpha),[0 50000],ic_per,options);
% [t2,x_eq] = ode15s(@(t,X) sys_reduced(t,X,Kp,phi1,KI,phi2,om0,om1,om2,uf,ug,uh,alpha),[0 5000],ic_eq,options);
% plot3(x_per(:,1),x_per(:,2),x_per(:,3),'LineWidth',2);
% hold on;
% plot3(x_eq(:,1),x_eq(:,2),x_eq(:,3),'LineWidth',2);
% plot3(x0_int,x1_int,x2_int,'g*','LineWidth',3);
% scatter3(0.5,0.07,1,'r','filled');
% scatter3(0.5,0.1,1,'r','filled');
% disp([x_per(end,1),x_per(end,2),x_per(end,3)]);
% disp([x_eq(end,1),x_eq(end,2),x_eq(end,3)]);
% xlabel('x0');
% ylabel('x1');
% zlabel('x2');
% legend('IC: [0.5, 0.07, 1]', 'IC: [0.5, 0.1, 1]','Interior equilibrium','Location','Best');
% grid on;
% hold off;

% Y0=0.019                                                                                                    
% Y1=0.04                                                                                                     
% Y2=0.06                                                                                                     
% km0=29                                                                                                      
% km1=26                                                                                                      
% km2=35                                                                                                      
% KS0=0.053                                                                                                   
% KS1=0.302                                                                                                   
% KS2=2.5e-5                                                                                                  
% L0=1e-6                                                                                                     
% KI2=3.5e-6                                                                                                  
% om0=(KS0/KS1)*(224/208)*(1-Y0)                                                                              
% om1=(KS1/KS2)*(32/224)*(1-Y1)                                                                               
% om2=(16/208)*(KS0/KS2)                                                                                      
% phi1=(km1/km0)*(Y1/Y0)                                                                                      
% phi2=(km2/km0)*(Y2/Y0)                                                                                      
% Kp=L0/KS2                                                                                                   
% KI=KS2/KI2                                                                                                  
% x0'=((s0/(1+s0))*(s2/(Kp+s2))-alph)*x0                                                                      
% x1'=(phi1*(s1/1+s1)*(1/(1+KI*s2))-alph)*x1                                                                  
% x2'=(phi2*s2/(1+s2)-alph)*x2                                                                                
% s0'=alph*(uf-s0)-(s0/(1+s0))*(s2/(Kp+s2))*x0                                                                
% s1'=alph*(ug-s1)+om0*(s0/(1+s0))*(s2/(Kp+s2))*x0-(phi1*(s1/(1+s1))*(1/(1+KI*s2)))*x1                        
% s2'=alph*(uh-s2)-om2*(s0/(1+s0))*(s2/(Kp+s2))*x0+om1*(phi1*(s1/(1+s1))*(1/(1+KI*s2)))*x1-(phi2*s2/(1+s2))*x2
% 
% s0=(-x0+uf)
% s1=(om0*x0-x1+ug)
% s2=(-om2*x0+om1*x1-x2+uh)

% x0=0.58953892
% x1=0.058655684
% x2=0.001
% s0=0.010461084
% s1=0.050648074
% s2=1.1338417