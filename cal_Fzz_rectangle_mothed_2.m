function [F_Uzz,Uzz] = cal_Fzz_rectangle_mothed_2(mod,gobs,nodex,nodez)

G = 6.67430e-11;
F_Uzz = zeros(size(gobs,1),size(mod,1)*size(mod,2));

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
            
            ATAN(i,j) = atan2(dz(i,j),dx(i,j));
            
            % 奇点判断
            if abs(dx(i,j)) < eps
                ATAN(i,j) = 0;
            end
            
        end
    end
    
    for j = 1:size(mod,2)
        for i = 1:size(mod,1)
            
            ATAN11 = ATAN(i,j);
            ATAN12 = ATAN(i+1,j);
            ATAN21 = ATAN(i,j+1);
            ATAN22 = ATAN(i+1,j+1);
            
            K_Uzz_11 = 2*G.*ATAN11; % 计算Uzz核函数
            K_Uzz_12 = 2*G.*ATAN12;
            K_Uzz_21 = 2*G.*ATAN21;
            K_Uzz_22 = 2*G.*ATAN22;
            
            F_Uzz(xi,k) = (K_Uzz_22 - K_Uzz_12 - K_Uzz_21 + K_Uzz_11).*10^12; % F_Uzz为重力异常核函数
            
            k = k + 1;
            
        end
    end
%     fprintf('正在计算第 %d 个测点\r', xi);
end
toc;

Uzz = F_Uzz*reshape(mod,[],1); % F_Uzz为重力异常
