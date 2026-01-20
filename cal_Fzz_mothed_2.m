function [Fzz,Gzz] = cal_Fzz_mothed_2(mod, gobs, nodex, nodez)
% 初始化输出矩阵
Fzz = zeros(size(gobs, 1), size(mod, 1) * size(mod, 2));

% 物理常数
G = 6.67430e-11;  % 万有引力常数
scale_factor = 1e12; % 与C代码中的1e12对应

% 奇点处理阈值
eps = 1e-10;        % 奇点判断阈值
dx_eps = 1e-10;     % 垂直边阈值
dz_eps = 1e-10;     % 水平边阈值

% 主循环
tic;
for k = 1:size(gobs, 1) % 观测点循环
    
    x_obs = gobs(k, 1);
    z_obs = gobs(k, 2);
    n = 1; % 网格索引初始化
    
    for j = 1:size(mod, 2)   % 网格列循环
        for i = 1:size(mod, 1)   % 网格行循环
            % 计算网格四个角点相对于观测点的位移
            x_dis = [
                nodex(i, j) - x_obs;
                nodex(i, j+1) - x_obs;
                nodex(i+1, j+1) - x_obs;
                nodex(i+1, j) - x_obs;
                nodex(i, j) - x_obs;  % 闭合多边形
            ];
            
            z_dis = [
                nodez(i, j) - z_obs;
                nodez(i, j+1) - z_obs;
                nodez(i+1, j+1) - z_obs;
                nodez(i+1, j) - z_obs;
                nodez(i, j) - z_obs;  % 闭合多边形
            ];
            
            f = 0; % 初始化当前网格的贡献值
            
            for edge = 1:4 % 单个网格4条边循环
                % 获取当前边的两个端点
                idx1 = edge;
                idx2 = edge + 1;
                
                xi = x_dis(idx1);
                zi = z_dis(idx1);
                xii = x_dis(idx2);
                zii = z_dis(idx2);
                
                % 1. 奇点判断 - 如果端点接近观测点，跳过该边
                if (abs(xi) < eps && abs(zi) < eps) || (abs(xii) < eps && abs(zii) < eps)
                    continue;
                end
                
                % 计算边向量
                dx = xii - xi;
                dz = zii - zi;
                
                % 2. 垂直边判断 (dx ≈ 0)
                if abs(dx) < dx_eps
                    % 垂直边公式
                    % 计算角度 (使用atan2处理象限)
                    c1 = atan2(zi, xi);
                    c2 = atan2(zii, xii);
                    
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
                    
                    % 计算当前边的贡献
                    factor = delta_c;
                    f_val = 2 * G * factor * scale_factor;
                    f = f + f_val;
                    continue; % 处理完后跳过后续计算
                end
                
                % 3. 水平边判断 (dz ≈ 0)
                if abs(dz) < dz_eps
                    % 水平边公式
                    f_val = 0;
                    f = f + f_val;
                    continue; % 处理完后跳过后续计算
                end
                
                % 4. 一般情况 (非垂直非水平边)
                % 计算角度 (使用atan2处理象限)
                c1 = atan2(zi, xi);
                c2 = atan2(zii, xii);
                
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
                
                % 计算对数项
                c3 = log((xii^2 + zii^2) / (xi^2 + zi^2));
                
                % 计算分母
                denominator = dx^2 + dz^2;
                
                % 避免除以零 (退化边处理)
                if denominator < 1e-12
                    continue; % 跳过退化边
                end
                
                % 计算分子项
                numerator = dz;
                
                % 计算当前边的贡献
                factor = dz * delta_c + 0.5 * dx * c3;
                f_val = 2 * G * (numerator / denominator) * factor * scale_factor;
                f = f + f_val;
            end
            
            % 存储当前网格对当前观测点的贡献
            Fzz(k, n) = f;
            n = n + 1;
        end
    end
end

% F*m=g
Gzz = Fzz * reshape(mod,[],1);

% 计算耗时
computation_time = toc;
fprintf('计算重力核函数矩阵Fzz完成\n');
fprintf('总耗时: %.6f 秒\n', computation_time);

end