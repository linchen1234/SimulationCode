function xdot=Underactuated_ship_RL_function(t,X,c,width,widthf,beta,gamma_c,gamma_a,gammaf,sigmaf,rho,gammav,sigmav)
%RLOTC

global center11 center12 center_f1 
global center21 center22 center_f2
global Node1 Node2 Node3 Node4 Nodef1 Nodef2 
global A B w 

x=X(1);
y=X(2);
psi=X(3);
u=X(4);
v=X(5);
r=X(6);
hat_tau_wr=X(7);  
hat_tau_wu=X(8);
alpha_fr=X(9);
alpha_fu=X(10);
hb_v=X(742);

m11=120*10^3;
m22=177.9*10^3;
m33=636*10^5;

x_d=A*sin(w*t);
y_d=B*(1-cos(w*t));
psi_d=atan2(y-y_d,x-x_d);

D_x_d=A*w*cos(w*t);
D_y_d=B*w*sin(w*t);

d_u=215*10^2;
d_u2=0.2*d_u;
d_u3=0.1*d_u;

d_v=147*10^3;
d_v2=0.2*d_v;
d_v3=0.1*d_v;

d_r=802*10^4;
d_r2=0.2*d_r;
d_r3=0.1*d_r;

f_u=-d_u*u/m11-d_u2*abs(u)*u/m11-d_u3*u^3/m11;
f_v=-d_v*v/m22-d_v2*abs(v)*v/m22-d_v3*v^3/m22;
f_r=-d_r*r/m33-d_r2*abs(r)*r/m33-d_r3*r^3/m33;

tau_w=diag([10000,20000,30000])...
          *[1.5+sin(0.1*t)+cos(0.01*t); 
           1.5+cos(0.1*t)+sin(0.01*t);  
           1.5-sin(0.1*t)+cos(0.01*t)];

rho_1=rho(1);
rho_2=rho(2);
z_e=sqrt((x-x_d)^2+(y-y_d)^2);
z_3=z_e;

psi_e=psi-psi_d;
z_1=psi_e;

while abs(z_1)>pi
    if psi>0
        psi=psi-2*pi;
    else
        psi=psi+2*pi;
    end
    z_1=psi-psi_d;
end

D_psi_d=(u*sin(psi_e)+v*cos(psi_e)-D_y_d*cos(psi_d)+D_x_d*sin(psi_d))/z_e;

fz=v*sin(psi_e)   +D_x_d*cos(psi_d)+D_y_d*sin(psi_d);
r_e=r-alpha_fr;
u_e=u-alpha_fu;
c_1=c(1);
c_2=c(2);

beta1=beta(1);
beta2=beta(2);
beta3=beta(3);
beta4=beta(4);

gamma_c1=gamma_c(1);
gamma_c2=gamma_c(2);
gamma_c3=gamma_c(3);
gamma_c4=gamma_c(4);

gamma_a1=gamma_a(1);
gamma_a2=gamma_a(2);
gamma_a3=gamma_a(3);
gamma_a4=gamma_a(4);

gamma_f1=gammaf(1);
gamma_f2=gammaf(2);


sigma_f1=sigmaf(1);
sigma_f2=sigmaf(2);

gamma_v=gammav;
sigma_v=sigmav;
width11=width(1);
width12=width(2);
width21=width(3);
width22=width(4);

width_f1=widthf(1);
width_f2=widthf(2);

Wc1=X(13:12+Node1);  
Wa1=X(13+Node1:12+Node1+Node1);

dS1=zeros(Node1,1);
Z1=z_1; 
for i=1:Node1
    dS1(i)=exp( -(Z1-center11(i))^2/(width11^2) ); 
end

omega1=-dS1*(beta1*z_1+0.5*Wa1'*dS1+D_psi_d); 

dWc1=-gamma_c1/(1+norm(omega1)^2)*omega1*( omega1'*Wc1-(beta1^2-1)*z_1^2-2*beta1*z_1*D_psi_d+0.25*norm(Wa1'*dS1)^2 );
dWa1=0.5*dS1*z_1+0.25*gamma_c1/(1+norm(omega1)^2)*(dS1)*dS1'*Wa1*omega1'*Wc1-gamma_a1*(dS1)*dS1'*Wa1;

alpha_r=-beta1*z_1-0.5*Wa1'*dS1;
y_r=alpha_fr-alpha_r;
D_alpha_fr=-y_r/c_1;

Wf1=X(589:588+Nodef1);
Sf1=zeros(Nodef1,1);
Zf1=[u;v;r]; 
for i=1:Nodef1
    Sf1(i)=exp( -sum((Zf1-center_f1(:,i)).^2 )/(width_f1^2) );  
end
Wc2=X(13+Node1+Node1:12+Node1+Node1+Node2);
Wa2=X(13+Node1+Node1+Node2:12+Node1+Node1+Node2+Node2);

dS2=zeros(Node2,1);
Z2=r_e; 
for i=1:Node2
    dS2(i)=exp( -(Z2-center12(i))^2/(width12^2));  
end
omega2=-dS2*( beta2*r_e+0.5*Wa2'*dS2+D_alpha_fr);
dWf1=gamma_f1*( gamma_c2/(1+norm(omega2)^2)*Sf1*(Sf1')*Wf1*omega2'*Wc2-2*gamma_c2/(1+norm(omega2)^2)*Sf1*D_alpha_fr*omega2'*Wc2+Sf1*r_e-sigma_f1*Wf1 );
F1_nn=Wf1'*Sf1;
dWc2=-gamma_c2/(1+norm(omega2)^2)*omega2*( -(beta2^2-1)*r_e^2-2*beta2*r_e*D_alpha_fr+2*hat_tau_wr*(0.5*hat_tau_wr-D_alpha_fr)+0.25*norm(Wa2'*dS2)^2+norm(Wf1'*Sf1)^2+2*Wf1'*Sf1*(hat_tau_wr-D_alpha_fr)+omega2'*Wc2 );
dWa2=0.5*dS2*r_e+0.25*gamma_c2/(1+norm(omega2)^2)*(dS2)*dS2'*Wa2*omega2'*Wc2+2*gamma_c2/(1+norm(omega2)^2)*omega2*beta2*r_e*hat_tau_wr+2*gamma_c2/(1+norm(omega2)^2)*omega2*hat_tau_wr^2-gamma_a2*(dS2)*dS2'*Wa2;
dhat_tau_wr=-rho_1*hat_tau_wr+2*gamma_c2/(1+norm(omega2)^2)*(beta2*r_e+hat_tau_wr)*(Wc2'-Wa2')*omega2+2*gamma_c2/(1+norm(omega2)^2)*Wf1'*Sf1*Wc2'*omega2;

Wc3=X(301:300+Node3);
Wa3=X(301+Node3:300+Node3+Node3);
dS3=zeros(Node3,1);
Z3=z_3;
for i=1:Node3
   dS3(i)=exp( -(Z3-center21(i))^2/(width21^2) );
end

omega3=-dS3*( (beta3*z_3+0.5*Wa3'*dS3)*cos(psi_e)^2+ fz);
dWc3=-gamma_c3/(1+norm(omega3)^2)*omega3*( -(beta3^2*cos(psi_e)^2-1)*z_3^2+0.25*norm(Wa3'*dS3)^2*cos(psi_e)^2-2*beta3*z_3*fz+Wc3'*omega3 );
dWa3=0.5*dS3*z_3*cos(psi_e)^2+0.25*gamma_c3/(1+norm(omega3)^2)*(dS3)*dS3'*Wa3*omega3'*Wc3*cos(psi_e)^2-gamma_a3*(dS3)*dS3'*Wa3;
dhb_v=gamma_v*( abs(z_3*sin(psi_e))-sigma_v*hb_v );

alpha_u=(-beta3*z_3-0.5*Wa3'*dS3)*cos(psi_e);
y_u=alpha_fu-alpha_u;
D_alpha_fu=-y_u/c_2;
Wc4=X(301+Node3+Node3:300+Node3+Node3+Node4);
Wa4=X(301+Node3+Node3+Node4:300+Node3+Node3+Node4+Node4);

dS4=zeros(Node4,1);
Z4=u_e; 
for i=1:Node4
    dS4(i)=exp( -(Z4-center22(i))^2/(width22^2) ); %lie
end

Wf2=X(677:676+Nodef2);
Sf2=zeros(Nodef2,1);
Zf2=[u;v;r]; 
for i=1:Nodef2
    Sf2(i)=exp( -sum((Zf2-center_f2(:,i)).^2 )/(width_f2^2) ); %lie 
end
omega4=-dS4*( beta4*u_e+0.5*Wa4'*dS4+D_alpha_fu );
dWf2=gamma_f2*( gamma_c4/(1+norm(omega4)^2)*Sf2*(Sf2')*Wf2*(omega4')*Wc4-2*gamma_c4/(1+norm(omega4)^2)*Sf2*D_alpha_fu*(omega4')*Wc4+Sf2*u_e-sigma_f2*Wf2  );
F2_nn=Wf2'*Sf2;
dWc4=-gamma_c4/(1+norm(omega4)^2)*omega4*( -(beta4^2-1)*u_e^2-2*beta4*u_e*D_alpha_fu+2*hat_tau_wu*(0.5*hat_tau_wu-D_alpha_fu)+0.25*norm(Wa4'*dS4)^2+norm(Wf2'*Sf2)^2+2*Wf2'*Sf2*(hat_tau_wu-D_alpha_fu)+omega4'*Wc4 );
dWa4=0.5*dS4*u_e+0.25*gamma_c4/(1+norm(omega4)^2)*(dS4)*dS4'*Wa4*omega4'*Wc4+2*gamma_c4/(1+norm(omega4)^2)*omega4*beta4*u_e*hat_tau_wu+2*gamma_c4/(1+norm(omega4)^2)*omega4*hat_tau_wu^2-gamma_a4*(dS4)*dS4'*Wa4;
dhat_tau_wu=-rho_2*hat_tau_wu+2*gamma_c4/(1+norm(omega4)^2)*(beta4*u_e+hat_tau_wu)*(Wc4'-Wa4')*omega4+2*gamma_c4/(1+norm(omega4)^2)*Wf2'*Sf2*Wc4'*omega4;

tau_rm=-beta2*r_e-hat_tau_wr-F1_nn-0.5*Wa2'*dS2;
tau_um=-beta4*u_e-hat_tau_wu-F2_nn-0.5*Wa4'*dS4;
tau_r=m33*tau_rm;
tau_u=m11*tau_um;
dx=u*cos(psi)-v*sin(psi);
dy=u*sin(psi)+v*cos(psi);
dpsi=r;
du=m22*v*r/m11+f_u+tau_u/m11+tau_w(1)/m11;
dv=-m11*u*r/m22+f_v+tau_w(2)/m22; 
dr=(m11-m22)*u*v/m33+f_r+tau_r/m33+tau_w(3)/m33;

dalpha_fr=(alpha_r-alpha_fr)/c_1;
dalpha_fu=(alpha_u-alpha_fu)/c_2;

xdot=[dx;dy;dpsi;du;dv;dr;dhat_tau_wr;dhat_tau_wu;dalpha_fr;dalpha_fu;tau_r;tau_u;dWc1;dWa1;dWc2;dWa2;dWc3;dWa3;dWc4;dWa4;dWf1;dWf2;alpha_r;alpha_u;dhb_v];
t
end

