clc; clear all; close all;
format compact;

% fftshift 2D
K1 = 3;
K2 = 14;
a = reshape(1:K1*K2, [K1,K2]);
b1 = fftshift(a);
c1 = reshape(fftshiftId2D(K1,K2,a(:)), [K1, K2]);
isequal(b1,c1)

% ifftshift 2D
b2 = ifftshift(a);
c2 = reshape(ifftshiftId2D(K1,K2,a(:)), [K1, K2]);
isequal(b2,c2)