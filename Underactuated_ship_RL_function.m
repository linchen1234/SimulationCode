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


