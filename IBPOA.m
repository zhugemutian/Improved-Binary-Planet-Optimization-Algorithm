clc;
clear;
%% parameter setting
load('kp_sc_1000');
%the Instances should include
% profits of all items $P = [p_1,p_2,\dots,p_n]$
% weights of all items $W = [w_1,w_2,\dots,w_n]$
% knapsack capacity $C$
% number of items $n$

T = 200;
pop = 30;
UB = 100;
LB = -100;
Rmin = (UB-LB)*n/1000;

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

%% initial solution
X = rand(pop,n)*(UB-LB)+LB;
X(1,:) = XT;
t = 1;
while t <= T
    
    for i = 1:pop
        for j = 1:b-1
            RAND(j) = 1/(1+exp(-X(i,j)))+1/Omega(j); %+(b-n/2)/n
            if rand < 1/(1+exp(-X(i,j)))+1/Omega(j); %+(b-n/2)/n
                Y(i,j) = 1;
            else
                Y(i,j) = 0;            
            end
        end
        for j = b:n
            RAND(j) = 1/(1+exp(-X(i,j)))-1/Omega(j); %+(b-n/2)/n
            if rand < 1/(1+exp(-X(i,j)))-1/Omega(j); %+(b-n/2)/n
                Y(i,j) = 1;
            else
                Y(i,j) = 0;            
            end
        end
    end

    %% GRO
    for i = 1:pop
        for j = n:-1:1
            if W*Y(i,:)'>C
                Y(i,j) = 0;
            else
                break
            end
        end
        for j = 1:n
            if W*Y(i,:)'+(1-Y(i,j))*W(j)<=C
                Y(i,j) = 1;
            end
        end
    end

%% moment of the plant
    FY = P*Y';
    [t1,t2] = max(FY);
    
    record(t) = t1;
    if t>1
        record(t) = max(t1,record(t-1));
    end
    
    Xsun = X(t2,:);
    Ysun = Y(t2,:);
    alpha = abs(min(FY)-P*Ysun');

    a = 2;
    for i = 1:pop
        mass(i) = a^(alpha/(alpha+min(FY)-FY(i)+1));
        Risun(i) = sum((Y(i,:)-Ysun).^2);
    end
    mass_sun = mass(t2);

    for i = 1:pop
        M(i) = mass(i)*mass_sun/Risun(i);
    end
    Mmax = max(M);
    beta = M./Mmax;

    %% search 
    for i = 1:pop
        %% global search
        if Risun(i) > Rmin
            Xt = X + beta(i)*rand*(Xsun-X(i,:));
        else
        %% local search
            c = 2 - t/T;
            mu = 0.5;
            sigma = 0.2;
            Xt = X+ c*rand*((mu+sigma*randn)*Xsun-X(i,:));
        end
    end
    X = Xt;
    t = t+1;
end
%% local search

plot(1:T,record)

