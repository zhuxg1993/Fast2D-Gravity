function [F_Uzz,Uzz] = cal_Fzz_rectangle_mothed_1(mod,gobs,nodex,nodez)

G = 6.67430e-11;
F_Uzz = zeros(size(gobs,1),size(mod,1)*size(mod,2));

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
            
            ATAN11 = atan2(z1,x1);
            ATAN12 = atan2(z2,x1);
            ATAN21 = atan2(z1,x2);
            ATAN22 = atan2(z2,x2);
            
            % 奇点判断
            if abs(x1) < eps
                ATAN11 = 0;
                ATAN12 = 0;
            end
            
            if abs(x2) < eps
                ATAN21 = 0;
                ATAN22 = 0;
            end
            
            K_Uzz_11 = 2*G.*ATAN11; % 计算Uzz核函数
            K_Uzz_12 = 2*G.*ATAN12;
            K_Uzz_21 = 2*G.*ATAN21;
            K_Uzz_22 = 2*G.*ATAN22;
            
            F_Uzz(xi,k) = (K_Uzz_22 - K_Uzz_12 - K_Uzz_21 + K_Uzz_11).*10^12; % F_Uzz为重力梯度张量核函数
            
            k = k + 1;
            
        end
    end
%     fprintf('正在计算第 %d 个测点\r', xi);
end
toc;

Uzz = F_Uzz*reshape(mod,[],1); % F_Uzz为重力梯度张量核函数
