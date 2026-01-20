function [F_Ux,Ux] = cal_Fx_rectangle_mothed_2(mod,gobs,nodex,nodez)

G = 6.67430e-11;
F_Ux = zeros(size(gobs,1),size(mod,1)*size(mod,2));

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
            ATAN(i,j) = atan2(dz(i,j),dx(i,j));
            
            % 奇点判断
            if abs(dx(i,j)) < eps
                ATAN(i,j) = 0;
                if abs(dz(i,j)) < eps
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
            
            K_Ux_11 = G.*(z1.*LN11 + 2.*x1.*ATAN11); % 计算Ux核函数
            K_Ux_12 = G.*(z2.*LN12 + 2.*x1.*ATAN12);
            K_Ux_21 = G.*(z1.*LN21 + 2.*x2.*ATAN21);
            K_Ux_22 = G.*(z2.*LN22 + 2.*x2.*ATAN22);
            
            F_Ux(xi,k) = (K_Ux_22 - K_Ux_12 - K_Ux_21 + K_Ux_11).*10^11; % F_Ux为重力异常核函数
            
            k = k + 1;
            
        end
    end
%     fprintf('正在计算第 %d 个测点\r', xi);
end
toc;

Ux = F_Ux*reshape(mod,[],1); % F_Ux为重力异常
