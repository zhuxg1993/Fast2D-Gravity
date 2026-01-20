function [F_Uxz,Uxz] = cal_Fxz_rectangle_mothed_1(mod,gobs,nodex,nodez)

G = 6.67430e-11;
F_Uxz = zeros(size(gobs,1),size(mod,1)*size(mod,2));

tic;
for xi = 1:size(gobs,1)
    
    k = 1;
    x0 = gobs(xi,1); % 观测点坐标x0，z0
    z0 = gobs(xi,2);
    
    for j = 1:size(mod,2)
        for i = 1:size(mod,1)
            
            xx1 = nodex(i,j); % 网格的上下左右边界对应坐标
            xx2 = nodex(i,j+1);
            zz1 = nodez(i,j);
            zz2 = nodez(i+1,j);
            
            x1 = xx1 - x0; % 网格四个边界到观测点对应轴的坐标之差
            x2 = xx2 - x0;
            z1 = zz1 - z0;
            z2 = zz2 - z0;
            
            LN11 = log(x1^2 + z1^2);
            LN12 = log(x1^2 + z2^2);
            LN21 = log(x2^2 + z1^2);
            LN22 = log(x2^2 + z2^2);
            
            % 奇点判断
            if abs(z1) < eps
                if abs(x1) < eps
                    LN11 = 0;
                end
                if abs(x2) < eps
                    LN21 = 0;
                end
            end
            
            if abs(z2) < eps
                if abs(x1) < eps
                    LN12 = 0;
                end
                if abs(x2) < eps
                    LN22 = 0;
                end
            end
            
            K_Uxz_11 = -G.*LN11; % 计算Uxz核函数
            K_Uxz_12 = -G.*LN12;
            K_Uxz_21 = -G.*LN21;
            K_Uxz_22 = -G.*LN22;
            
            F_Uxz(xi,k) = (K_Uxz_22 - K_Uxz_12 - K_Uxz_21 + K_Uxz_11).*10^12; % F_Uxz为核函数
            
            k = k + 1;
            
        end
    end
%     fprintf('正在计算第 %d 个测点\r', xi);
end
toc;

Uxz = F_Uxz*reshape(mod,[],1); % F_Uxz为重力梯度张量核函数
