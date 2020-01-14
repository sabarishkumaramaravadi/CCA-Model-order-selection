clear all;
close all;
clc;
X = 0: 0.2: 15;
Y = chi2pdf (X, 4);
plot(X,Y);