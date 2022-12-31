clear all;
clc;

syms s a b;

A = [a b;-b a];

S = [s 0; 0 s];

C = S - A

D = inv(C)

E = ilaplace(D)



