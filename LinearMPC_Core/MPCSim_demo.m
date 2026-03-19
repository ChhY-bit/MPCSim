clear,clc
%% 仿真设置
dt = 1e-4;
T = 20;

A = [0 1;0 0];
B = [0 ; 1];

xr_fun = @(t) [0.5+sin(0.5*t);...
               0+0.5*cos(0.5*t)];
ur_fun = @(x,t) -0.25*sin(0.5*t);
x_ini = [0;0];
%% 初始化
prob = MPCSim_init(A,B,diag([1000,1]),0.1,20,0.05,[1.5;1.2],-[1.5;1.2],1,-1);
N = prob.Horizen;
Ts = prob.Ts;
tspan = 0:dt:T;
N_max = length(tspan);

tspan_r = 0:dt:(T+N*Ts);
xr = xr_fun(tspan_r);
ur = ur_fun(xr,tspan_r);
x = zeros(2,N_max);
x(:,1) = x_ini;
u = zeros(1,N_max);
%% 终端预计算
Kt = -place(prob.Ad,prob.Bd,[0.95,0.9]);
alpha = 10;
terminal_term = MPCSim_tmnl(prob,Kt,alpha);
%% 仿真计算
dynamic_fun = @(x,u,t)A*x + B*u;
for k = 1:N_max-1
    if mod(tspan(k),Ts) == 0    % 采样周期
    % =========================== 在此处编辑控制器 ========================
        x_measured = x(:,k);
        error = xr(:,k) - x_measured;
        % 注意：必须取预测时域长度的参考值送入MPC--------------------------
        ur_part = ur(:,k:Ts/dt:k+Ts/dt*(N-1));     % 取预测时域长度，用于约束
        xr_part = xr(:,k+Ts/dt:Ts/dt:k+Ts/dt*N);   % 取预测时域长度，用于约束
        % -----------------------------------------------------------------
        u(:,k) = ur(:,k) - MPCSim_solv(prob,error,xr_part,ur_part,terminal_term);
        % u(:,k) = ur(:,k) - MPCSim_solv(prob,error,xr_part,ur_part,[]);
        % u(:,k) = ur(:,k);
    % =====================================================================
    else                        % 零阶保持
        u(:,k) = u(:,k-1);
    end

    % 动态更新
    x(:,k+1) = update_rk4(dynamic_fun,x(:,k),u(:,k),[],dt);
end
%% 结果
figure(1)
subplot(2,1,1);
plot(tspan,x(1,:))
hold on
plot(tspan,xr(1,1:N_max),'r--');
subplot(2,1,2);
plot(tspan,x(2,:))
hold on
plot(tspan,xr(2,1:N_max),'r--');

figure(2)
plot(tspan(1:N_max-1),u(1:N_max-1))
hold on
plot(tspan(1:N_max-1),ur(1:N_max-1),'r--');

figure(3)
subplot(2,1,1);
plot(tspan,xr(1,1:N_max)-x(1,:))
hold on
subplot(2,1,2);
plot(tspan,xr(2,1:N_max)-x(2,:))
hold on