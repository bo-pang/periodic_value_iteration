function Triple_VI_TAC_V2_beta(sf)
% This example comes from
% Control of General Dynamic Systems With Periodically Varying Parameters
% Via Liapunov-Floquet Transformation. Journal of Dynamic Systems, 
% Measurement, and Control 116.4 (1994): 650-658.
global A B Q R n m s0 nx ww N0 wp theta phi nP nK;

% Construt the LPT system
kbar = 1;
cbar = 0.5;
wp = 1;
T = 2*pi/wp;
PIk1 = 1;
PIk2 = 2; % 2 unstable, 0.7 stable
gamma = @(t) PIk1+PIk2*cos(wp*t);
alpha = 1.0;

Ao = @(t) [zeros(3),eye(3);
    kbar*(gamma(t)-3), kbar*(3-gamma(t)), -kbar, -cbar, 0, 0;
    kbar*(4-gamma(t)), 2*kbar*(gamma(t)-3), kbar*(3-gamma(t)), cbar, -cbar, 0;
    -kbar, kbar*(4-gamma(t)), -kbar*(3+gamma(t)*(alpha-2)), 0, cbar, -cbar;];
Bo = @(t) [zeros(3);
    1, -1, 0;
    -1, 2, -1;
    0, -1, 2;];
da = @(t) [1+sin(3*t),0,0;
    0,1+sin(3*t),0;
    0,0,1+sin(3*t);];
DA = @(t) [
    zeros(3),zeros(3);
    da(t),zeros(3);];
A = @(t) Ao(t)+DA(t);
B = @(t) Bo(t);
[n,m] = size(B(0));
Q = @(t) eye(n);
R = @(t) eye(m);

N0 = 6;
N = 2*N0+1; % number of terms in Fourier
nx = n*(n+1)/2;
nP = nx*N;
nK = m*n*N;


% collect data
theta = [];
phi = [];
i = 0;
RK = nP+nK;
M0 = 30;
h = 0.2;    % sampling interval
M = 20*M0+1;     % sampling number
x = zeros(n+RK,1);
x(1:n) = zeros(n,1);
K = @(t) zeros(m,n);
ww = -500+(500-(-500))*rand(500,m);
n_reset = 0;

size_data = 0;
i = 0;
while size_data<M
    ts = i*h;
    tf = (i+1)*h;
    [~,y] = ode45(@(t,y) simpleSYS(t,y,K),[ts,tf],x);
    xs = y(1,:);
    xf = y(end,:);
    phi = [phi;kronv(xf(1:n))'-kronv(xs(1:n))'];
    theta = [theta;xf(n+1:n+nP),xf(n+nP+1:end)];
    [size_data,~] = size(phi);
    if norm(x(1:n))>10
        x(1:n) = zeros(n,1);
        i = 0;
        ww = -500+(500-(-500))*rand(500,m);
        n_reset = n_reset + 1;
    else
        x(1:n) = xf(1:n);
    end
    x(n+1:end) = zeros(RK,1);
    i = i+1;
end
M
n_reset
RK
rk = rank(theta'*theta);
if RK~=rk
    disp('Rank defficient!');
end

% learning process
s0 = 0;   % starting time

% VI ADP
Xhat0 = zeros(nP+nK,1);
% ode_opt = odeset('MaxStep',0.3);
tspan = [sf:-0.1:s0];
[That,Xhat_v] = ode45(@(s,y) vi_adp(s,y),tspan,Xhat0);

lt = length(That);
Xhat_m = zeros(m*n,N,lt);
Khat = zeros(m,n,lt);
for i=1:lt
    Xhat_m(:,:,i) = reshape(Xhat_v(i,nP+1:end),[m*n,N]);
    Khat(:,:,i) = reshape(Xhat_m(:,:,i)*F(That(i)),[m,n]);  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try value iteration
P0 = zeros(n);
Pvi0 = sm2vec(P0);
[~,Pvi] = ode45(@(s,y) vi(s,y),That,Pvi0);
PVI = zeros(n,n,lt);
KVI = zeros(m,n,lt);
for j=1:lt
    PVI(:,:,j) = vec2sm(Pvi(j,:)',n);
    KVI(:,:,j) = R(That(j))\(B(That(j))'*PVI(:,:,j));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try VI to obatin optimal solution
P0 = zeros(n);
Pvi_uniform0 = sm2vec(P0);
Tuniform = [100:-0.1:40.1,That'];
[~,Pvi_uniform] = ode45(@(s,y) vi(s,y),Tuniform,Pvi_uniform0);
PVI_uniform = zeros(n,n,length(Tuniform));
KVI_uniform = zeros(m,n,length(Tuniform));
for j=1:length(Tuniform)
    PVI_uniform(:,:,j) = vec2sm(Pvi_uniform(j,:)',n);
    KVI_uniform(:,:,j) = R(Tuniform(j))\(B(Tuniform(j))'*PVI_uniform(:,:,j));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fit the final value matrix with fourier
KVI = flip(KVI,3);
KVI_uniform = flip(KVI_uniform,3);
Tuniform = flip(Tuniform);
Xhat_m = flip(Xhat_m,3);
Khat = flip(Khat,3);
That = flip(That,1);
% T0 = 20;
% for i=1:lt
%     if That(i)>T0
%         k0 = i;
%         break;
%     end
% end
k0 = floor(lt/3);
k0
WW = zeros(k0,m*n);
UU = zeros(k0,N);

for i=1:k0
    WW(i,:) = reshape(Khat(:,:,i),[1 m*n]);
    UU(i,:) = F(That(i));
end

XbarK = UU\WW;

Kfinal = @(t) R(t)\reshape(XbarK'*F(t),[m,n]);

figure(1);
yplot_VIadp_norm = zeros(lt,1);
KVI_norm = zeros(lt,1);
for i=1:lt
    yplot_VIadp_norm(i) = norm(Kfinal(That(i)),'fro');
    KVI_norm(i) = norm(KVI(:,:,i),'fro');
end
plot(That,yplot_VIadp_norm,'r-');hold on;
plot(That,KVI_norm,'b*','markersize',5);hold on;
Khat_norm = zeros(lt,1);
for i=1:lt
    Khat_norm(i) = norm(Khat(:,:,i),'fro');
end
plot(That,Khat_norm,'mo','markersize',5);hold on;
xlabel('Algorithmic time (s)');
ylabel('Control gains');
legend({'$\Vert\bar{K}(\cdot)\Vert$','$\Vert K(\cdot)\Vert$',...
    '$\Vert\hat{K}_k\Vert$'},'interpreter','latex');


save(['Khat_response',num2str(sf),'.mat'],'That','Khat_norm');
% save the Xbar in the disturbance-free case
%save('Dfree.mat','XbarK');
% save the Xbar in the disturbance-free case
save('DXbarX.mat','XbarK');
disp('Learning complete');

figure(2);
hat_true_diff = zeros(lt,1);
Xhat_norm =zeros(lt,1);
hat_optimal_diff = zeros(lt,1);
bar_optimal_diff = zeros(lt,1);
for i=1:lt
    hat_true_diff(i) = norm(Khat(:,:,i)-KVI(:,:,i));
    Xhat_norm(i) = norm(Xhat_m(:,:,i));
    hat_optimal_diff(i) = norm(KVI_uniform(:,:,i)-Khat(:,:,i));
    bar_optimal_diff(i) = norm(Kfinal(That(i))-KVI_uniform(:,:,i));
end
subplot(2,2,1);
plot(That,hat_true_diff,'bo','MarkerSize',2);
xlabel('Algorithmic Time');
ylabel('$\Vert\hat{K}_k-K(s_k,s_k)\Vert$','interpreter','latex');
title('(a)');
subplot(2,2,2);
plot(That,bar_optimal_diff);
xlabel('System Evolution Time');
ylabel('$\Vert\bar{K}(t)-K^*(t)\Vert$','interpreter','latex');
title('(b)');
subplot(2,2,4);
plot(That,Xhat_norm);hold on;
plot(That,norm(XbarK)*ones(size(That)),'--');
xlabel('Algorithmic Time');
legend({'$\Vert\hat{W}^{K}\Vert$','$\Vert\bar{W}^{K}\Vert$'},'interpreter','latex');
title('(d)');
subplot(2,2,3);
plot(That,hat_optimal_diff,'ro','MarkerSize',2);
xlabel('Algorithmic Time');
ylabel('$\Vert \hat{K}_k-K^*(s_k)\Vert$','interpreter','latex');
title('(c)');

max_diff = max(bar_optimal_diff)

% Compare response
x0 = [0.8;0.5;0.3;0;0;0];
ts = 0;
tf = 20;
%[tini,xini] = ode45(@(t,y) A(t)*y,[0,300],x0);
[tfinal,xfinal] = ode45(@(t,y) (A(t)-B(t)*Kfinal(t))*y,[ts,tf],x0);
ufinal = [];
for i=1:length(tfinal)
    ufinal = [ufinal -Kfinal(tfinal(i))*xfinal(i,:)'];
end
tf = 100;
figure(3);
[t,x] = ode45(@(t,y) A(t)*y,[ts,tf],x0);
plot(t,x);
xlabel('System Evolution Time (s)');
ylabel('States');
legend({'$\eta_1$','$\eta_2$','$\eta_3$','$\dot{\eta}_1$',...
    '$\dot{\eta}_2$','$\dot{\eta}_3$'},'interpreter','latex');
figure(4);
plot(tfinal,xfinal);hold on;
xlabel('System Evolution Time (s)');
ylabel('States');
legend({'$\eta_1$','$\eta_2$','$\eta_3$','$\dot{\eta}_1$',...
    '$\dot{\eta}_2$','$\dot{\eta}_3$'},'interpreter','latex');
axis([0 20 -0.4 0.9]);
end

function dx = vi(t,x)
global A B Q R n;
xm = vec2sm(x,n);
dxm = A(t)'*xm+xm*A(t)-xm*B(t)*(R(t)\B(t)')*xm+Q(t);
dxm = -dxm;
dx = sm2vec(dxm);
end

function dx = vi_adp(t,x)
global theta phi Q R nP nx N0 m n;

term1 = reshape(x(1:nP),[nx,2*N0+1])*F(t);
term2 = sm2vec(Q(t));
tmp = reshape(reshape(x(nP+1:end),[m*n,2*N0+1])*F(t),[m,n]);
term3 = sm2vec(tmp'*(R(t)\tmp));
%dx = (theta'*theta)\(theta'*phi)*(-term1-term2+term3);
dx = theta\(phi*(-term1-term2+term3));
end

function dx = simpleSYS(t,x,K)
    global A B ww;
    [n,m] = size(B(0));
    
    % exploration noise
    %e1 = 0.02*sin(100*t)+0.282*exp(-0.04*t)*sin(20*t);
    e = zeros(m,1);
    for i=1:m
        e(i) = 0.2*sum(sin(ww(:,i)*t));
    end
    
    u = -K(t)*x(1:n)+e;
    
    % derivatives
    dx = A(t)*x(1:n)+B(t)*u;
    tildex = kronv(x(1:n));
    dx = [dx;kron(F(t)',tildex')'];
    dx = [dx;kron(F(t)',kron(x(1:n)',2*u'))'];
end

function Fn = F(t,n)
    global N0 wp;
    if nargin<2
        n=N0;
    end
    Fn = 1;
    for i = 1:n
        Fn = [Fn;cos(wp*i*t);sin(wp*i*t)];
    end
end

function X = kronv(x)
len = length(x);
X = [];
for i=1:len
    for j=i:len
        if i==j
            X(end+1) = x(i)*x(j);
        else
            X(end+1) = sqrt(2)*x(i)*x(j);
        end
    end
end
X = X';
end

function x = sm2vec(X)
[n,~] = size(X);
N = n*(n+1)/2;
x = zeros(N,1);
k = 1;
for i=1:n
    for j=i:n
        if i==j
            x(k) = X(i,j);
        else
            x(k) = sqrt(2)*X(i,j);
        end
        k = k+1;
    end
end
end

function X = vec2sm(x,n)
X = zeros(n);
num = flip(1:n);
for i=1:n
    index = 0;
    for k=1:i-1
        index = index+num(k);
    end
    for j=0:n-i
        if j~=0
            X(i,i+j)=x(index+j+1)/sqrt(2);
            X(j+i,i)=X(i,j+i);
        else
            X(i,j+i)=x(index+j+1);
        end
    end
end
end
