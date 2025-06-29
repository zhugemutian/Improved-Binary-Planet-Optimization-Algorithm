clc;
clear;

load('kp_sc_1000');

T = 200;
pop = 30;
UB = 100;
LB = -100;
Rmin = (UB-LB)*n/1000;
kkk = 30; %算法重复次数

%% omega

E1 = [P./W;1:n];
E2 = sortrows(E1',1,'descend')';

P = P(E2(2,:));
W = W(E2(2,:));
Tweight = 0;
%greedy algorithm
for i = 1:n
    if Tweight + W(i) <= C
        Tweight = Tweight+W(i);
    else
        b = i;
        break
    end
end
XT = zeros(1,n);
XT(1:b-1) = UB;
XT(b:n) = LB;
% residual capacity
r = C-Tweight;
% omega
Omega1 = [];
for i = 1:b-1
    Omega1 = [Omega1, floor(r*P(b)/(P(i)*W(b)-P(b)*W(i)))+1];
end
Omega2 = [inf];
for i = b+1:n
    Omega2 = [Omega2, floor(r*P(b)/(P(b)*W(i)-P(i)*W(b)))+1];
end
Omega = [Omega1,Omega2];

%convergence_all = [];
solution_all = [];
Time_all = [];
for order_transfer = 1:20
    [solution, Time] = BPOATF(order_transfer,n,P,W,C,kkk);
    solution_all = [solution_all; solution];
    Time_all = [Time_all; Time];
end
[solution, Time] = BPOATF2(n,P,W,C,b,Omega,kkk);
solution_all = [solution_all;solution];
Time_all = [Time_all;Time];