function [Fxz,Gxz] = cal_Fzz_mothed_3(mod, gobs, nodex, nodez)
% 初始化输出矩阵
Fxz = zeros(size(gobs, 1), size(mod, 1) * size(mod, 2));

% 物理常数
G = 6.67430e-11;  % 万有引力常数
scale_factor = 1e12; % 与C代码中的1e12对应

% 奇点处理阈值
eps_val = 1e-10;     % 奇点判断阈值
dz_eps = 1e-10;      % 水平边阈值

% 获取网格尺寸
[grid_rows, grid_cols] = size(mod);
notez = grid_rows + 1;  % 水平边行数 (NZ+1)
notex = grid_cols + 1;  % 垂直边列数 (NX+1)

% 主循环
tic;
for k = 1:size(gobs, 1) % 观测点循环
    
    x_obs = gobs(k, 1);
    z_obs = gobs(k, 2);
    
    % 预分配边贡献矩阵
    H_edge = zeros(notez, grid_cols);  % 水平边贡献 (NOTEZ × NX)
    V_edge = zeros(grid_rows, notex);  % 垂直边贡献 (NZ × NOTEX)
    
    % 1. 计算所有水平边贡献 (上下边)
    for j = 1:notez  % j: 1..NZ+1
        for i = 1:grid_cols  % i: 1..NX
            % 获取边的两个端点坐标
            x1 = nodex(j, i) - x_obs;
            z1 = nodez(j, i) - z_obs;
            x2 = nodex(j, i+1) - x_obs;
            z2 = nodez(j, i+1) - z_obs;
            
            % 奇点判断 - 跳过包含观测点的边
            if (abs(x1) < eps_val && abs(z1) < eps_val) || ...
               (abs(x2) < eps_val && abs(z2) < eps_val)
                H_edge(j, i) = 0;
                continue;
            end
            
            % 水平边判断 (dz ≈ 0)
            dz_edge = z2 - z1;
            if abs(dz_edge) < dz_eps
                % 计算水平边贡献
                H_edge(j, i) = 0;
            else
                % 一般情况公式 (保持与原始代码一致)
                c1 = atan2(z1, x1);
                c2 = atan2(z2, x2);
                if c1 < 0, c1 = c1 + 2*pi; end
                if c2 < 0, c2 = c2 + 2*pi; end
                
                % 处理角度差环绕问题
                delta_c = c2 - c1;
                if delta_c > pi
                    delta_c = delta_c - 2*pi;
                elseif delta_c < -pi
                    delta_c = delta_c + 2*pi;
                end
                
                % 计算对数项
                c3 = log((x2^2 + z2^2) / (x1^2 + z1^2));
                
                % 计算边向量
                dx = x2 - x1;
                dz = z2 - z1;
                denominator = dx^2 + dz^2;
                
                % 避免除以零
                if denominator < 1e-12
                    H_edge(j, i) = 0;
                    continue;
                end
                
                % 计算分子项
                numerator = dz;
                
                % 计算当前边的贡献
                factor = dz * delta_c + 0.5 * dx * c3;
                H_edge(j, i) = 2.0 * G * (numerator / denominator) * factor * scale_factor;
            end
        end
    end
    
    % 2. 计算所有垂直边贡献 (左右边)
    for j = 1:grid_rows  % j: 1..NZ
        for i = 1:notex  % i: 1..NX+1
            % 获取边的两个端点坐标
            x1 = nodex(j, i) - x_obs;
            z1 = nodez(j, i) - z_obs;
            x2 = nodex(j+1, i) - x_obs;
            z2 = nodez(j+1, i) - z_obs;
            
            % 奇点判断 - 跳过包含观测点的边
            if (abs(x1) < eps_val && abs(z1) < eps_val) || ...
               (abs(x2) < eps_val && abs(z2) < eps_val)
                V_edge(j, i) = 0;
                continue;
            end
            
            % 垂直边公式 (使用x1作为参考点)
            c1 = atan2(z1, x1);
            c2 = atan2(z2, x2);
            if c1 < 0, c1 = c1 + 2*pi; end
            if c2 < 0, c2 = c2 + 2*pi; end
            
            % 处理角度差环绕问题
            delta_c = c2 - c1;
            if delta_c > pi
                delta_c = delta_c - 2*pi;
            elseif delta_c < -pi
                delta_c = delta_c + 2*pi;
            end
            V_edge(j, i) = 2 * G * delta_c * scale_factor;
        end
    end
    
    % 3. 组合所有网格的贡献 (修正索引)
    n = 1; % 网格索引
    for i = 1:grid_cols % 网格列循环 (1..NX)
        for j = 1:grid_rows  % 网格行循环 (1..NZ)
            % 组合四条边的贡献:
            %   上边: H_edge(j, i)
            %   右边: V_edge(j, i+1)
            %   下边: H_edge(j+1, i) 取负
            %   左边: V_edge(j, i) 取负
            
            f = H_edge(j, i) ...      % 上边
                + V_edge(j, i+1) ...  % 右边
                - H_edge(j+1, i) ...  % 下边 (负贡献)
                - V_edge(j, i);       % 左边 (负贡献)
            
            Fxz(k, n) = f;
            n = n + 1;
        end
    end
end

% F*m=g
Gxz = Fxz * reshape(mod,[],1);

% 计算耗时
computation_time = toc;
fprintf('计算重力核函数矩阵Fzz完成\n');
fprintf('总耗时: %.6f 秒\n', computation_time);
end
