function [Fx,Gx] = cal_Fx_mothed_4(mod, gobs, nodex, nodez)
% 初始化输出矩阵
Fx = zeros(size(gobs, 1), size(mod, 1) * size(mod, 2));

% 物理常数
G = 6.67430e-11;        % 万有引力常数
scale_factor = 1e11;     % 缩放因子

% 阈值设置
eps_val = 1e-10;         % 奇点判断阈值
dz_eps = 1e-10;          % 水平边阈值

% 获取网格尺寸
[NZ, NX] = size(mod);
NOTEX = NX + 1;
NOTEZ = NZ + 1;
num_points = size(gobs,1);  % 测点数量

% 预分配水平边和垂直边存储数组
H_edge = zeros(NOTEZ, NX);  % 水平边贡献 (NOTEZ x NX)
V_edge = zeros(NZ, NOTEX);  % 垂直边贡献 (NZ x NOTEX)

% ====== 记录开始时间 ======
start_time = tic;
% ===========================

for jj = 1:num_points
    
    xx = gobs(jj, 1);
    zz = gobs(jj, 2);
    
    % 预分配节点计算值
    node_atan = zeros(NOTEZ, NOTEX); % 存储每个节点的atan2值
    node_ln = zeros(NOTEZ, NOTEX);   % 存储每个节点的log值
    dis_x = zeros(NOTEZ, NOTEX);     % 存储每个节点的x距离
    dis_z = zeros(NOTEZ, NOTEX);     % 存储每个节点的z距离
    
    % 预计算每个节点的atan2和log值
    for j = 1:NOTEZ
        for i = 1:NOTEX
            dis_x(j, i) = nodex(j, i) - xx;
            dis_z(j, i) = nodez(j, i) - zz;
            
            % 奇点判断
            if abs(dis_x(j, i)) < eps_val && abs(dis_z(j, i)) < eps_val
                node_atan(j, i) = 0.0;
                node_ln(j, i) = 0.0;
                continue;
            end
            
            % 计算atan2值（注意MATLAB的atan2参数顺序：y, x）
            node_atan(j, i) = atan2(dis_z(j, i), dis_x(j, i));
            
            % 计算ln值
            node_ln(j, i) = log(dis_x(j, i)^2 + dis_z(j, i)^2);
        end
    end
    
    % 计算所有上下边贡献 (j: 1..NOTEZ, i: 1..NX)
    for j = 1:NOTEZ
        for i = 1:NX
            x1 = dis_x(j, i);
            z1 = dis_z(j, i);
            x2 = dis_x(j, i+1);
            z2 = dis_z(j, i+1);
            
            % 奇点判断
            if (abs(x1) < eps_val && abs(z1) < eps_val) || ...
                    (abs(x2) < eps_val && abs(z2) < eps_val)
                H_edge(j, i) = 0.0;
                continue;
            end
            
            % 水平边判断
            dz_edge = z2 - z1;
            is_level = abs(dz_edge) < dz_eps;
            
            if is_level
                % 水平边公式
                % 计算对数项差值
                c3 = node_ln(j, i+1) - node_ln(j, i);
                
                % 水平边贡献
                H_edge(j, i) = -G * z1 * c3 * scale_factor;
            else
                % 一般情况公式
                c1 = node_atan(j, i);
                c2 = node_atan(j, i+1);
                
                % 调整角度到0~2π范围
                if c1 < 0, c1 = c1 + 2*pi; end
                if c2 < 0, c2 = c2 + 2*pi; end
                
                % 处理角度差环绕问题
                delta_c = c2 - c1;
                if delta_c > pi
                    delta_c = delta_c - 2*pi;
                elseif delta_c < -pi
                    delta_c = delta_c + 2*pi;
                end
                
                % 计算对数项差值
                c3 = node_ln(j, i+1) - node_ln(j, i);
                
                % 计算边向量
                dx = x2 - x1;
                dz = z2 - z1;
                denominator = dx^2 + dz^2;
                
                % 避免除以零
                if denominator < 1e-12
                    H_edge(j, i) = 0.0;
                    continue;
                end
                
                % 计算分子项
                numerator = x1 * z2 - x2 * z1;
                
                % 计算当前边的贡献
                factor = dz * delta_c + 0.5 * dx * c3;
                H_edge(j, i) = 2.0 * G * (numerator / denominator) * factor * scale_factor;
            end
        end
    end
    
    % 计算所有垂直边贡献 (j: 1..NZ, i: 1..NOTEX)
    for j = 1:NZ
        for i = 1:NOTEX
            x1 = dis_x(j, i);
            z1 = dis_z(j, i);
            x2 = dis_x(j+1, i);
            z2 = dis_z(j+1, i);
            
            % 奇点判断
            if (abs(x1) < eps_val && abs(z1) < eps_val) || ...
                    (abs(x2) < eps_val && abs(z2) < eps_val)
                V_edge(j, i) = 0.0;
                continue;
            end
            
            % 垂直边公式
            c1 = node_atan(j, i);
            c2 = node_atan(j+1, i);
            
            % 调整角度到0~2π范围
            if c1 < 0, c1 = c1 + 2*pi; end
            if c2 < 0, c2 = c2 + 2*pi; end
            
            % 处理角度差环绕问题
            delta_c = c2 - c1;
            if delta_c > pi
                delta_c = delta_c - 2*pi;
            elseif delta_c < -pi
                delta_c = delta_c + 2*pi;
            end
            
            % 计算竖直边贡献
            V_edge(j, i) = 2 * G * x1 * delta_c * scale_factor;
        end
    end
    
    % 遍历所有单元组合核函数（列优先顺序）
    n = 1; % 网格索引
    for i = 1:NX        % 列循环 (X方向)
        for j = 1:NZ    % 行循环 (Z方向)
            % 组合四条边的贡献：
            %   上边: H_edge(j, i)
            %   右边: V_edge(j, i+1)
            %   下边: H_edge(j+1, i) 取负
            %   左边: V_edge(j, i) 取负
            
            f = H_edge(j, i) ...      % 上边 (正贡献)
                + V_edge(j, i+1) ...    % 右边 (正贡献)
                - H_edge(j+1, i) ...    % 下边 (负贡献)
                - V_edge(j, i);         % 左边 (负贡献)
            
            Fx(jj, n) = f;
            n = n + 1;
        end
    end
end

% F*m=g
Gx = Fx * reshape(mod,[],1);

% ====== 输出结果 ======
computation_time = toc(start_time);
fprintf('计算重力核函数矩阵Fx完成\n');
fprintf('总耗时: %.6f 秒\n', computation_time);

end