close all
clear all
clc

X = randn(100,1);
x = detect_event(X,0.5);
subplot(211), plot(X), hold on, stem(x)

Y = randn(100,1);
y = detect_event(Y,0.5);
subplot(212), plot(Y), hold on, stem(y)

[N,S] = nec_suf(x,y);
disp([N S])