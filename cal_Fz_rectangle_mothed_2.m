function [F_Uz,Uz] = cal_Fz_rectangle_mothed_2(mod,gobs,nodex,nodez)

G = 6.67430e-11;
F_Uz = zeros(size(gobs,1),size(mod,1)*size(mod,2));

LN = zeros(size(mod,1)+1,size(mod,2)+1);
ATAN = zeros(size(mod,1)+1,size(mod,2)+1);
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
            ATAN(i,j) = atan2(dx(i,j),dz(i,j));
            
            % 奇点判断
            if abs(dz(i,j)) < eps
                ATAN(i,j) = 0;
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
            ATAN11 = ATAN(i,j);
            ATAN12 = ATAN(i+1,j);
            ATAN21 = ATAN(i,j+1);
            ATAN22 = ATAN(i+1,j+1);
            x1 = dx(i,j);
            x2 = dx(i,j+1);
            z1 = dz(i,j);
            z2 = dz(i+1,j);
            
            K_Uz_11 = G.*(x1.*LN11 + 2.*z1.*ATAN11); % 计算Uz核函数
            K_Uz_12 = G.*(x1.*LN12 + 2.*z2.*ATAN12);
            K_Uz_21 = G.*(x2.*LN21 + 2.*z1.*ATAN21);
            K_Uz_22 = G.*(x2.*LN22 + 2.*z2.*ATAN22);
            
            F_Uz(xi,k) = (K_Uz_22 - K_Uz_12 - K_Uz_21 + K_Uz_11).*10^11; % F_Uz为重力异常核函数
            
            k = k + 1;
            
        end
    end
%     fprintf('正在计算第 %d 个测点\r', xi);
end
toc;

Uz = F_Uz*reshape(mod,[],1); % F_Uz为重力异常
