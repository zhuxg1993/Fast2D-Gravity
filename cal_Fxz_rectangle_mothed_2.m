function [F_Uxz,Uxz] = cal_Fxz_rectangle_mothed_2(mod,gobs,nodex,nodez)

G = 6.67430e-11;
F_Uxz = zeros(size(gobs,1),size(mod,1)*size(mod,2));

LN = zeros(size(mod,1)+1,size(mod,2)+1);
dx = zeros(size(mod,1)+1,size(mod,2)+1);
dz = zeros(size(mod,1)+1,size(mod,2)+1);

tic;
for xi = 1:size(gobs,1)
    
    k = 1;
    x0 = gobs(xi,1); % 观测点坐标x0，z0
    z0 = gobs(xi,2);
    
    for j = 1:size(mod,2)+1
        for i = 1:size(mod,1)+1
            
            dx(i,j) = nodex(i,j) - x0; % 网格四个边界到观测点对应轴的坐标之差
            dz(i,j) = nodez(i,j) - z0;
            
            LN(i,j) = log(dx(i,j)^2 + dz(i,j)^2);
            
            % 奇点判断
            if abs(dz(i,j)) < eps
                if abs(dx(i,j)) < eps
                    LN(i,j) = 0;
                end
            end
            
        end
    end
    
    for j = 1:size(mod,2)
        for i = 1:size(mod,1)
            
            LN11 = LN(i,j);
            LN12 = LN(i+1,j);
            LN21 = LN(i,j+1);
            LN22 = LN(i+1,j+1);
            
            K_Uxz_11 = -G.*LN11; % 计算Uxz核函数
            K_Uxz_12 = -G.*LN12;
            K_Uxz_21 = -G.*LN21;
            K_Uxz_22 = -G.*LN22;
            
            F_Uxz(xi,k) = (K_Uxz_22 - K_Uxz_12 - K_Uxz_21 + K_Uxz_11).*10^12; % F_Uxz为重力异常核函数
            
            k = k + 1;
            
        end
    end
%     fprintf('正在计算第 %d 个测点\r', xi);
end
toc;

Uxz = F_Uxz*reshape(mod,[],1); % F_Uxz为重力异常
