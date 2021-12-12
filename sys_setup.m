function [A, B, C, D, G, B_tot, D_tot] = sys_setup()

n = sqrt(398600/6778^3);
A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     3*n^2 0 0 0 2*n 0;
     0 0 0 -2*n 0 0;
     0 0 -n^2 0 0 0];
B = [0 0 0;
     0 0 0;
     0 0 0;
     1 0 0;
     0 1 0;
     0 0 1];
C = [1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 1 0 0 0];
D = [0 0 0;
     0 0 0;
     0 0 0];
G = [0;0;0;0;1;0]; % Disturbance matrix (disturbance only affects in track dynamics)
B_tot = [B,G];
D_tot = [D,[0;0;0]];