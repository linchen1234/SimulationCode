%欠驱动无人艇，Critic-actor架构，强化学习
%对比仿真实验---ANNC方法

clear
clc
close all
tic

global center11 center12 center_f1 
global center21 center22 center_f2
global Node1 Node2 Node3 Node4 Nodef1 Nodef2 
global A B w 

global e1_bar_0 e1_bar_inf e1_under_0 e1_under_inf kappa_1  
global e2_bar_0 e2_bar_inf e2_under_0 e2_under_inf kappa_2 
global center1 center2 Node11 Node12

%系统状态和设计参数赋初值
xs1=[0;30;0.2;0;0;0]; 
A=65; 
B=45;  
w=0.041;
x_d0=0;   
y_d0=0;

hat_tau_wr0=10; 
hat_tau_wu0=10;  
rho1=20; 
rho2=20;
rho=[rho1;rho2];

beta1=2.0;  
beta2=4.0;  
beta3=3.5; 
beta4=0.8;  

beta=[beta1;beta2;beta3;beta4];
alpha_r0=0;
alpha_u0=0;
alpha0=[alpha_r0,alpha_u0]; 
tau_rm=0;
tau_um=0;
tau0=[tau_rm,tau_um]; 
c=[0.01;0.01];
hbv=0;
gammaf1=0.01;  
gammaf2=0.01;
gammaf=[gammaf1;gammaf2];

sigmaf1=0.1; 
sigmaf2=0.1; 
sigmaf=[sigmaf1;sigmaf2];  
gammav=2;
sigmav=0.01;

Node1=72;
Node2=72;
Node3=72;
Node4=72;

width1=0.1; 
width2=0.1; 
width3=0.8;
width4=0.01;
width=[width1,width2,width3,width4];
widthf1=1.8;
widthf2=0.7; 
widthf=[widthf1,widthf2];

gamma_c1=1;          
gamma_a1=8;  
gamma_c2=2.7;  
gamma_a2=10;  
gamma_c3=3;   
gamma_a3=12;   
gamma_c4=4;   
gamma_a4=20;  
gamma_c=[gamma_c1,gamma_c2,gamma_c3,gamma_c4];
gamma_a=[gamma_a1,gamma_a2,gamma_a3,gamma_a4];

center1_1=zeros(1,Node1);
for i=-4:1:4
   for j=1:1:8
     center1_1(j+(i+4)*8)=i*1.8;      
   end    
end
center11=center1_1;
center1_2=zeros(1,Node2);
for i=-4:1:4
   for j=1:1:8
     center1_2(j+(i+4)*8)=i*1.8; 
   end    
end
center12=center1_2;
center2_1=zeros(1,Node3);
for i=-4:1:4
   for j=1:1:8
     center2_1(j+(i+4)*8)=i*1.8;    
   end    
end
center21=center2_1;
center2_2=zeros(1,Node4);
for i=-4:1:4
   for j=1:1:8
     center2_2(j+(i+4)*8)=i*1.8;    
   end    
end
center22=center2_2;

Wc10=0.02*ones(1,Node1);  
Wa10=0.01*ones(1,Node1);
Wc20=0.2*ones(1,Node2); 
Wa20=0.1*ones(1,Node2);
Wc30=0.1*ones(1,Node3); 
Wa30=0.1*ones(1,Node3);
Wc40=0.2*ones(1,Node4); 
Wa40=0.2*ones(1,Node4);

Nodef1=188;
Nodef2=163;
Wf10=zeros(1,Nodef1);
Wf20=zeros(1,Nodef2); 

center_f11=zeros(1,Nodef1);
center_f12=zeros(1,Nodef1);
center_f13=zeros(1,Nodef1);
for i=-5:1:5
   for j=1:1:8
      center_f11(j+(i+5)*8)=i*1.4;  
      center_f12(j+(i+5)*8)=i*1.8; 
      center_f13(j+(i+5)*8)=i*1.2; 
   end
end
center_f1=[center_f11;center_f12;center_f13];
center_f21=zeros(1,Nodef2);
center_f22=zeros(1,Nodef2);
center_f23=zeros(1,Nodef2);
for i=-5:1:5
   for j=1:1:8
       center_f21(j+(i+5)*8)=1.8*i;  
       center_f22(j+(i+5)*8)=1.8*i;
       center_f23(j+(i+5)*8)=1.3*i;
   end
end
center_f2=[center_f21;center_f22;center_f23];

tfinal=90;
x0=[xs1', hat_tau_wr0, hat_tau_wu0, alpha0, tau0, Wc10, Wa10, Wc20, Wa20, Wc30, Wa30, Wc40, Wa40, Wf10, Wf20,alpha0,hbv]';
length(x0)

[~,x]=ode45(@Underactuated_ship_RL_function,0:0.1:tfinal,x0,[],c,width,widthf,beta,gamma_c,gamma_a,gammaf,sigmaf,rho,gammav,sigmav);
length(x0)

Node11=88;
Node12=63;

center11_nn=zeros(1,Node11); 
for i=-5:1:5
    for j=1:1:8
        center11_nn(j+(i+5)*7)=0.02*i;
    end
end
center1=center11_nn;

center21_nn=zeros(1,Node12);
for i=-10:1:10
    for j=-1:1:1
            center21_nn(j+2+(i+10)*3)=0.5*i;
    end
end
center2=center21_nn;

width_nn=[0.01;0.01];
sigma_nn=[0.01;0.01];
gamma_nn=[0.01;0.01];

W1_0=zeros(Node11,1);
W2_0=zeros(Node12,1);

xi_3_nn=0;            
xi_4_nn=0;

k_nn=[0.1;0.1;   
    0.08;0.08; 
    0.01;0.01]; 
hat_d_wr=0;
hat_d_wu=0;
tau_u_nn=0;
tau_r_nn=0;

e1_bar_0=pi/2;
e1_bar_inf=0.06;
e1_under_0=pi/2;
e1_under_inf=0.02;
kappa_1=0.06;

psi_d0=atan2(xs1(2)-y_d0,xs1(1)-x_d0);
psi_e0=xs1(3)-psi_d0;
D_e1_bar_0=-kappa_1*(e1_bar_0 - e1_bar_inf);
D_e1_under_0=-kappa_1*(e1_under_0 - e1_under_inf);
gamma_10=e1_under_0/e1_bar_0;
D_gamma_10=(D_e1_under_0*e1_bar_0-e1_under_0*D_e1_bar_0)/(e1_bar_0^2);

chi_10=0.5*(1/(e1_under_0+psi_e0)+1/(e1_bar_0-psi_e0));
phi_10=0.5*D_e1_under_0/(e1_under_0+psi_e0)-0.5*D_e1_bar_0/(e1_bar_0-psi_e0)-0.5*D_gamma_10/gamma_10; 
z_10=0.5*log(1+psi_e0/e1_under_0)-0.5*log(1-psi_e0/e1_bar_0);

e2_bar_0=45;
e2_bar_inf=2.5;
e2_under_0=12;
e2_under_inf=0.2;
kappa_2=0.04;   

z_e0=sqrt((xs1(1)-x_d0)^2+(xs1(2)-y_d0)^2);
D_e2_bar_0=-kappa_2*(e2_bar_0 - e2_bar_inf);
D_e2_under_0=-kappa_2*(e2_under_0 - e2_under_inf);
gamma_20=e2_under_0/e2_bar_0;
D_gamma_20=(D_e2_under_0*e2_bar_0-e2_under_0*D_e2_bar_0)/(e2_bar_0^2);

chi_20=0.5*(1/(-e2_under_0+z_e0)+1/(e2_bar_0-z_e0));
phi_20=0.5*D_e2_under_0/(e2_under_0-z_e0)+0.5*D_e2_bar_0/(-e2_bar_0+z_e0)-0.5*D_gamma_20/gamma_20; % 非对称的情况下
z_20=0.5*log(-1+z_e0/e2_under_0)-0.5*log(1-z_e0/e2_bar_0);

D_x_d0=A*w;
D_y_d0=0;
D_psi_d0=(xs1(4)*sin(psi_e0)+xs1(5)*cos(psi_e0)-D_y_d0*cos(psi_d0)+D_x_d0*sin(psi_d0))/z_e0;

alpha_fr0=(-k_nn(1)*z_10-phi_10)/chi_10+D_psi_d0;
alpha_fu0= (-(k_nn(2)*z_20+phi_20)/chi_20+xs1(5)*sin(psi_e0)+D_x_d0*cos(psi_d0)+D_y_d0*sin(psi_d0))/cos(psi_e0);
alpha_f0=[alpha_fr0,alpha_fu0];

y0=[xs1' xi_3_nn xi_4_nn hat_d_wr hat_d_wu tau_r_nn tau_u_nn alpha_f0 W1_0' W2_0' alpha_f0]';
[t,y]=ode45(@Underactuated_ship_position_error_NN_OB_20180330_function,0:0.1:tfinal,y0,[],k_nn,c,width_nn,gamma_nn,sigma_nn);


e1_bar=(e1_bar_0 - e1_bar_inf)*exp(-kappa_1*t)+e1_bar_inf;
e1_under=(e1_under_0 - e1_under_inf)*exp(-kappa_1*t)+e1_under_inf;
e2_bar=(e2_bar_0 - e2_bar_inf)*exp(-kappa_2*t)+e2_bar_inf;
e2_under=(e2_under_0 - e2_under_inf)*exp(-kappa_2*t)+e2_under_inf;

x_d=A*sin(w*t);
y_d=B*(1-cos(w*t));
psi_d=atan2(x(:,2)-y_d,x(:,1)-x_d);  
psi_d_n=atan2(y(:,2)-y_d,y(:,1)-x_d);  

len=length(t);
psi_e=zeros(len,1);    
psi_e_nn=zeros(len,1);    
for i=1:len
      psi_e(i)=x(i,3)-psi_d(i);
       psi_e_nn(i)=y(i,3)-psi_d_n(i);
       
    while abs(psi_e(i))>pi   
        if x(i,3)>0
            x(i,3)=x(i,3)-2*pi;   
        else
            x(i,3)=x(i,3)+2*pi;
        end
        psi_e(i)=x(i,3)-psi_d(i);
    end
        while abs(psi_e_nn(i))>pi   
        if y(i,3)>0
            y(i,3)=y(i,3)-2*pi;  
        else
            y(i,3)=y(i,3)+2*pi;
        end
        psi_e_nn(i)=y(i,3)-psi_d_n(i);
    end
end

z_e=sqrt((x(:,2)-y_d).^2+(x(:,1)-x_d).^2);
z_e_nn=sqrt((y(:,2)-y_d).^2+(y(:,1)-x_d).^2);

figure(2);
plot(x(:,1),x(:,2),'g',y(:,1),y(:,2),'-b',x_d,y_d,'-r','LineWidth',1.5);
legend('RLOTC','ANNC','Reference trajectory');
hold off;
xlabel('x(m)');
ylabel('y(m)');
axis equal

figure(3);
plot(t,psi_e,'g',t,psi_e_nn,'-b','LineWidth',1.5);
legend('RLOTC','ANNC');
hold on
hold off
xlabel('Time(s)');
ylabel('Tracking error z_{1}');

figure(4);
plot(t,z_e,'g',t,z_e_nn,'-b','LineWidth',1.5);
legend('RLOTC','ANNC');
hold on
hold off
xlabel('Time(s)');
ylabel('Tracking error z_{3}');

L=length(t);
hat_tau_wr=zeros(L,1);
hat_tau_wu=zeros(L,1);
alpha_fr=zeros(L,1);
alpha_fu=zeros(L,1);
alpha_fr_nn=zeros(L,1);
alpha_fu_nn=zeros(L,1);
alpha_r=zeros(L,1);
alpha_u=zeros(L,1);
tau_r=zeros(L,1);
tau_u=zeros(L,1);
Wc1_norm=zeros(L,1);
Wa1_norm=zeros(L,1);
Wc2_norm=zeros(L,1);
Wa2_norm=zeros(L,1);
Wc3_norm=zeros(L,1);
Wa3_norm=zeros(L,1);
Wc4_norm=zeros(L,1);
Wa4_norm=zeros(L,1);
Wf1_norm=zeros(L,1);
Wf2_norm=zeros(L,1);
 
for i=2:L
    hat_tau_wr(i)=x(i,7);
    hat_tau_wu(i)=x(i,8);

    alpha_fr(i)=x(i,9);
    alpha_fu(i)=x(i,10);  
    
    alpha_fr_nn(i)=y(i,13);
    alpha_fu_nn(i)=y(i,14);  
      
        Wc1_norm(i)=norm(x(i,13:12+Node1));
        Wa1_norm(i)=norm(x(i,13+Node1:12+Node1+Node1));
    
        Wc2_norm(i)=norm(x(i,157:156+Node2));
        Wa2_norm(i)=norm(x(i,157+Node2:156+Node2+Node2));
    
        Wc3_norm(i)=norm(x(i,301:300+Node3));
        Wa3_norm(i)=norm(x(i,301+Node3:300+Node3+Node3));
    
        Wc4_norm(i)=norm(x(i,445:444+Node4));
        Wa4_norm(i)=norm(x(i,445+Node4:444+Node4+Node4));
    
        Wf1_norm(i)=norm(x(i,589:588+Nodef1));
        Wf2_norm(i)=norm(x(i,589+Nodef1:588+Nodef1+Nodef2));
    
end

r_e=x(6)-alpha_fr;
u_e=x(4)-alpha_fu;
r_e_nn=y(6)-alpha_fr_nn;
u_e_nn=y(4)-alpha_fu_nn;
alpha_r_nn=zeros(L,1);
alpha_u_nn=zeros(L,1);
tau_r_nn=zeros(L,1);
tau_u_nn=zeros(L,1);

for i=2:L
    alpha_r(i)=(x(i,840)-x(i-1,840))/(t(i)-t(i-1));
    alpha_u(i)=(x(i,841)-x(i-1,841))/(t(i)-t(i-1));
    tau_r(i)=(x(i,11)-x(i-1,11))/(t(i)-t(i-1));
    tau_u(i)=(x(i,12)-x(i-1,12))/(t(i)-t(i-1));
    alpha_r_nn(i)=(y(i,166)-y(i-1,166))/(t(i)-t(i-1));
    alpha_u_nn(i)=(y(i,167)-y(i-1,167))/(t(i)-t(i-1));
    tau_r_nn(i)=(y(i,11)-y(i-1,11))/(t(i)-t(i-1));
    tau_u_nn(i)=(y(i,12)-y(i-1,12))/(t(i)-t(i-1));
end

r11=zeros(L,1);
for i=2:L
   r11(i)=psi_e(i)^2+alpha_r(i)^2; 
end
r12=zeros(L,1);
for i=2:L
   r12(i)=r_e(i)^2+tau_r(i)^2; 
end
r21=zeros(L,1);
for i=2:L
   r21(i)=z_e(i)^2+alpha_u(i)^2; 
end
r22=zeros(L,1);
for i=2:L
   r22(i)=u_e(i)^2+tau_u(i)^2; 
end

r11_nn=zeros(L,1);
for i=2:L
   r11_nn(i)=psi_e_nn(i)^2+alpha_r_nn(i)^2; 
end
r12_nn=zeros(L,1);
for i=2:L
   r12_nn(i)=r_e_nn(i)^2+tau_r_nn(i)^2; 
end
r21_nn=zeros(L,1);
for i=2:L
   r21_nn(i)=z_e_nn(i)^2+alpha_u_nn(i)^2; 
end
r22_nn=zeros(L,1);
for i=2:L
   r22_nn(i)=u_e_nn(i)^2+tau_u_nn(i)^2; 
end

figure(5);
plot(t,tau_r,'g',t,tau_r_nn,'b','LineWidth',1.5);
xlabel('Time(s)');
ylabel('Control input \tau_{r}');
legend('RLOTC','ANNC');

figure(6);
plot(t,tau_u,'g',t,tau_u_nn,'b','LineWidth',1.5);
xlabel('Time(s)');
ylabel('Control input \tau_{u}');
legend('RLOTC','ANNC');

figure(7)
plot(t,Wc1_norm,'b',t,Wc2_norm,'y',t,Wc3_norm,'g',t,Wc4_norm,'m','LineWidth',1.5);
h=legend('1','2','3','4');
set(h,'Interpreter', 'latex','string',{'$||\hat{W}_{c1}||$','$||\hat{W}_{c2}||$','$||\hat{W}_{c3}||$','$||\hat{W}_{c4}||$'});
xlabel('Time(s)');
ylabel('Norms of critic NN weights');

figure(8)
plot(t,Wa1_norm,'b',t,Wa2_norm,'y',t,Wa3_norm,'g',t,Wa4_norm,'m','LineWidth',1.5);
h=legend('1','2','3','4');
set(h,'Interpreter', 'latex','string',{'$||\hat{W}_{a1}||$','$||\hat{W}_{a2}||$','$||\hat{W}_{a3}||$','$||\hat{W}_{a4}||$'});
xlabel('Time(s)');
ylabel('Norms of actor NN weights');

figure(9)
plot(t,Wf1_norm,'b',t,Wf2_norm,'-m','LineWidth',1.5);
h=legend('1','2','3');
set(h,'Interpreter', 'latex','string',{'$||\hat{W}_{f1}||$','$||\hat{W}_{f2}||$'});
xlabel('Time(s)');
ylabel('Norms of NN weight estimates');

figure(10);
plot(t,r11+r12+r21+r22,'g',t,r11_nn+r12_nn+r21_nn+r22_nn,'b','LineWidth',1.5);
xlabel('Time(s)');
ylabel('Total cost function r_{1}+r_{2}+r_{3}+r_{4}');
legend('RLOTC','ANNC');
