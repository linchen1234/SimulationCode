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
