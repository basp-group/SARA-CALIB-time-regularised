clc; clear all; close all;
format compact;

K1 = 21;
K2 = K1;
K = [K1, K2];
S = 4;

point = [4,4];

center = floor(K2/2)*K1 + floor(K1/2)+ 1;
[ll_center,I_center] = getSupportRect(center, S, K1, K2);

[ll_p,I_p] = getSupportRect((point(2)-1)*K1 + point(1), S, K1, K2);
