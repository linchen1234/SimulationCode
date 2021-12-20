
function xdot=Underactuated_ship_position_error_NN_OB_20180330_function(t,X,k,c,width,gamma,sigma)

%ANNC¿ØÖÆ

global e1_bar_0 e1_bar_inf e1_under_0 e1_under_inf kappa_1  
global e2_bar_0 e2_bar_inf e2_under_0 e2_under_inf kappa_2 
global A B w 
global center1 center2  
global   Node11 Node12

k_z1=k(1);
k_z2=k(2);
k_3=k(3);
k_4=k(4);
k_dr=k(5);
k_du=k(6);
c_1=c(1);
c_2=c(2);

width1=width(1);
width2=width(2);

gamma1=gamma(1);
gamma2=gamma(2);

sigma1=sigma(1);
sigma2=sigma(2);
x=X(1);
y=X(2);
psi=X(3);
u=X(4);
v=X(5);
r=X(6);
xi_3=X(7);
xi_4=X(8);
alpha_fr=X(13);
alpha_fu=X(14);

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
       
e1_bar=(e1_bar_0 - e1_bar_inf)*exp(-kappa_1*t)+e1_bar_inf;
e1_under=(e1_under_0 - e1_under_inf)*exp(-kappa_1*t)+e1_under_inf;
gamma_1=e1_under/e1_bar;

e2_bar=(e2_bar_0 - e2_bar_inf)*exp(-kappa_2*t)+e2_bar_inf;
e2_under=(e2_under_0 - e2_under_inf)*exp(-kappa_2*t)+e2_under_inf;
gamma_2=e2_under/e2_bar;

D_e1_bar=-kappa_1*(e1_bar_0 - e1_bar_inf)*exp(-kappa_1*t);
D_e1_under=-kappa_1*(e1_under_0 - e1_under_inf)*exp(-kappa_1*t);
D_gamma_1=(D_e1_under*e1_bar-e1_under*D_e1_bar)/(e1_bar^2);

D_e2_bar=-kappa_2*(e2_bar_0 - e2_bar_inf)*exp(-kappa_2*t);
D_e2_under=-kappa_2*(e2_under_0 - e2_under_inf)*exp(-kappa_2*t);
D_gamma_2=(D_e2_under*e2_bar-e2_under*D_e2_bar)/(e2_bar^2);

psi_e=psi-psi_d;
while abs(psi_e)>pi
    if psi>0
        psi=psi-2*pi;
    else
        psi=psi+2*pi;
    end
    psi_e=psi-psi_d;
end
z_1=0.5*log(1+psi_e/e1_under)-0.5*log(1-psi_e/e1_bar);

z_e=sqrt((x-x_d)^2+(y-y_d)^2);
z_2=0.5*log(-1+z_e/e2_under)-0.5*log(1-z_e/e2_bar);

chi_1=0.5*(1/(e1_under+psi_e)+1/(e1_bar-psi_e));
phi_1=0.5*D_e1_under/(e1_under+psi_e)-0.5*D_e1_bar/(e1_bar-psi_e)-0.5*D_gamma_1/gamma_1; 

chi_2=0.5*(1/(-e2_under+z_e)+1/(e2_bar-z_e));
phi_2=0.5*D_e2_under/(e2_under-z_e)+0.5*D_e2_bar/(-e2_bar+z_e)-0.5*D_gamma_2/gamma_2; 
D_psi_d=(u*sin(psi_e)+v*cos(psi_e)-D_y_d*cos(psi_d)+D_x_d*sin(psi_d))/z_e;

alpha_r=D_psi_d+(-k_z1*z_1-phi_1)/chi_1;
alpha_u=(-(k_z2*z_2+phi_2)/chi_2+v*sin(psi_e)+D_x_d*cos(psi_d)+D_y_d*sin(psi_d))/cos(psi_e);

r_e=r-alpha_fr;
u_e=u-alpha_fu;

y_r=alpha_fr-alpha_r;
y_u=alpha_fu-alpha_u;

hat_d_wr=xi_3+k_dr*r_e;
hat_d_wu=xi_4+k_du*u_e;

W1=X(15:14+Node11);
Z1=r;
S1=zeros(Node11,1);
for i=1:Node11
    S1(i,1)=exp(-sum((Z1-center1(:,i)).^2)/(width1^2));
end
dW1=gamma1*(r_e*S1-sigma1*abs(r_e)*W1);
F1_nn=W1'*S1;

W2=X(15+Node11:14+Node11+Node12);
Z2=u;
S2=zeros(Node12,1);
for i=1:Node12
    S2(i,1)=exp(-sum((Z2-center2(:,i)).^2)/(width2^2));
end
dW2=gamma2*(u_e*S2-sigma2*abs(u_e)*W2);
F2_nn=W2'*S2;

tau_r=m33*(-k_3*r_e-z_1*chi_1-(m11-m22)*u*v/m33-F1_nn-hat_d_wr-y_r/c_1);
tau_u=m11*(-k_4*u_e-z_2*chi_2*cos(psi_e)-m22*v*r/m11-F2_nn-hat_d_wu-y_u/c_2);

dx=u*cos(psi)-v*sin(psi);
dy=u*sin(psi)+v*cos(psi);
dpsi=r;
du=m22*v*r/m11+f_u+tau_u/m11+tau_w(1)/m11;
dv=-m11*u*r/m22+f_v+tau_w(2)/m22;
dr=(m11-m22)*u*v/m33+f_r+tau_r/m33+tau_w(3)/m33;

dxi_3=r_e-k_dr*(y_r/c_1+(m11-m22)*u*v/m33+F1_nn+tau_r/m33+hat_d_wr);
dxi_4=u_e-k_du*(y_u/c_2+m22*v*r/m11+F2_nn+tau_u/m11+hat_d_wu);

dalpha_fr=(alpha_r-alpha_fr)/c_1;
dalpha_fu=(alpha_u-alpha_fu)/c_2;

xdot=[dx dy dpsi du  dv dr dxi_3 dxi_4 hat_d_wr hat_d_wu tau_r tau_u dalpha_fr dalpha_fu dW1' dW2' alpha_r alpha_u]';

t
end

