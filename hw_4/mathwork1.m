close all; clear all; clc
syms J a0 v0 p0 pf T
J = 20*(pf-p0-v0*T-0.5*a0*T^2)^2/(T^2);
J_diff = diff(J,T);
sol = solve(J_diff == 0, T);
J_min = subs(J,T,sol);