function [solution, Time] = BPOATF(transferFuncNum,n,P,W,C,kkk)

% parameters
MaxIter = 200;
pop = 30;
UB = 100;
LB = -100;
Rmin = (UB-LB)*n/1000;

% Transfer function selection
transferFuncs = {
    @(v) 1./(1+exp(-2*v));                  % S1
    @(v) 1./(1+exp(-v));                    % S2
    @(v) 1./(1+exp(-v/2));                  % S3
    @(v) 1./(1+exp(-v/3));                  % S4
    @(v) erf(sqrt(pi)/2*v);                 % V1
    @(v) abs(tanh(v));                      % V2
    @(v) v/sqrt(1+v^2);                     % V3
    @(v) abs(2/pi*atan(pi*2*v))             % V4
    @(v) abs(v)^1.5;                        % U1
    @(v) abs(v)^2;                          % U2
    @(v) abs(v)^3;                          % U3
    @(v) abs(v)^4;                          % U4
    @(v) sqrt(abs(v))/sqrt(UB);             % T1
    @(v) abs(v)/UB;                         % T2
    @(v) abs(v)^(1/3)/UB^(1/3);             % T3
    @(v) abs(v)^(1/4)/UB^(1/4);             % T4
    @(v) sqrt(1-2^v);                       % Z1
    @(v) sqrt(1-5^v);                       % Z2
    @(v) sqrt(1-8^v);                       % Z3
    @(v) sqrt(1-10^v);                      % Z4
};

if transferFuncNum < 1 || transferFuncNum > length(transferFuncs)
    error('Transfer function number must be between 1 and %d', length(transferFuncs));
end
    
transferFunc = transferFuncs{transferFuncNum};


% Display selected transfer function
%fprintf('Using transfer function #%d\n', transferFuncNum);

Ztt = 0;
ZZ = [];
ZZ_con_1 = [];

tic
% 算法重复次数 kk
for kk = 1:kkk

% Initialization
X = -5*ones(pop,n);

% main loop
for t = 1:MaxIter
    % 计算选择概率
    for i = 1:pop
        for j = 1:n
            SV(i,j) = transferFunc(X(i,j));
        end
    end
    % 离散化
    for i = 1:pop
        for j = 1:n
            if rand <= SV(i,j);
                Y(i,j) = 1;
            else
                Y(i,j) = 0;
            end
        end
    end
    
    for i = 1:pop
        if Y(i,:)*W' > C
            Y(i,:) = zeros(1,n);
            X(i,:) = -5*ones(1,n);
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
            c = 2 - t/MaxIter;
            mu = 0.5;
            sigma = 0.2;
            Xt = X+ c*rand*((mu+sigma*randn)*Xsun-X(i,:));
        end
    end
    X = Xt;
end
ZZ = [ZZ,max(record)];
%ZZ_con_1 = [ZZ_con_1;record];
end
Ztt = toc;
%convergence = sum(ZZ_con_1/kk);
solution = ZZ;
Time = Ztt;