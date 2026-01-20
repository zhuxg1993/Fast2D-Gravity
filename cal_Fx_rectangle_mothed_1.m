function [F_Ux,Ux] = cal_Fx_rectangle_mothed_1(mod,gobs,nodex,nodez)

G = 6.67430e-11;
F_Ux = zeros(size(gobs,1),size(mod,1)*size(mod,2));

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
            
%             x1 = xx1 - x0 + 1e-6; % 网格四个边界到观测点对应轴的坐标之差
%             x2 = xx2 - x0 + 1e-6;
%             z1 = zz1 - z0 + 1e-6;
%             z2 = zz2 - z0 + 1e-6;
            
            LN11 = log(x1^2 + z1^2);
            LN12 = log(x1^2 + z2^2);
            LN21 = log(x2^2 + z1^2);
            LN22 = log(x2^2 + z2^2);
            ATAN11 = atan2(z1,x1);
            ATAN12 = atan2(z2,x1);
            ATAN21 = atan2(z1,x2);
            ATAN22 = atan2(z2,x2);
            
            % 奇点判断
            if abs(x1) < eps
                ATAN11 = 0;
                ATAN12 = 0;
                if abs(z1) < eps
                    LN11 = 0;
                end
                if abs(z2) < eps
                    LN12 = 0;
                end
            end
            
            if abs(x2) < eps
                ATAN21 = 0;
                ATAN22 = 0;
                if abs(z1) < eps
                    LN21 = 0;
                end
                if abs(z2) < eps
                    LN22 = 0;
                end
            end
            

            K_Ux_11 = G.*(z1.*LN11 + 2.*x1.*ATAN11); % 计算Ux核函数
            K_Ux_12 = G.*(z2.*LN12 + 2.*x1.*ATAN12);
            K_Ux_21 = G.*(z1.*LN21 + 2.*x2.*ATAN21);
            K_Ux_22 = G.*(z2.*LN22 + 2.*x2.*ATAN22);
            
            F_Ux(xi,k) = (K_Ux_22 - K_Ux_12 - K_Ux_21 + K_Ux_11).*10^11; % F_Ux为核函数
            
            k = k + 1;
            
        end
    end
%     fprintf('正在计算第 %d 个测点\r', xi);
end
toc;

Ux = F_Ux*reshape(mod,[],1); % F_Ux为重力水平分量核函数
