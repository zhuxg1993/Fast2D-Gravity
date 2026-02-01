function [Fx,Gx] = cal_Fx_mothed_1(mod, gobs, nodex, nodez)
% Initialize the output matrix
Fx = zeros(size(gobs, 1), size(mod, 1) * size(mod, 2));

% Physical constants
G = 6.67430e-11;
scale_factor = 1e11;

% Main iteration loop
tic;
for k = 1:size(gobs, 1) % Loop over observation points
    
    x_obs = gobs(k, 1);
    z_obs = gobs(k, 2);
    n = 1; % Initialize grid indices
    
    for j = 1:size(mod, 2)
        for i = 1:size(mod, 1)
            % Compute the distances from the four grid corner points to the observation point
            x_dis = [
                nodex(i, j) - x_obs;
                nodex(i, j+1) - x_obs;
                nodex(i+1, j+1) - x_obs;
                nodex(i+1, j) - x_obs;
                nodex(i, j) - x_obs;  % Closed polygon
            ];
            
            z_dis = [
                nodez(i, j) - z_obs;
                nodez(i, j+1) - z_obs;
                nodez(i+1, j+1) - z_obs;
                nodez(i+1, j) - z_obs;
                nodez(i, j) - z_obs;  % Closed polygon
            ];
            
            f = 0; % Initialize the contribution of the current grid cell
            
            for edge = 1:4 % Loop over the four edges of a single grid cell
                % Get the two endpoints of the current edge
                idx1 = edge;
                idx2 = edge + 1;
                
                xi = x_dis(idx1);
                zi = z_dis(idx1);
                xii = x_dis(idx2);
                zii = z_dis(idx2);
                
                % Singularity check - skip this edge if an endpoint is close to the observation point
                if (abs(xi) < eps && abs(zi) < eps) || (abs(xii) < eps && abs(zii) < eps)
                    continue;
                end
                
                % Compute the angle (use atan2 to handle quadrant ambiguity)
                c1 = atan2(zi, xi);
                c2 = atan2(zii, xii);
                
                % Adjust the angle to the range [0, 2¦Ð]
                if c1 < 0
                    c1 = c1 + 2*pi;
                end
                if c2 < 0
                    c2 = c2 + 2*pi;
                end
                
                % Handle angle difference wrapping
                delta_c = c2 - c1;
                if delta_c > pi
                    delta_c = delta_c - 2*pi;
                elseif delta_c < -pi
                    delta_c = delta_c + 2*pi;
                end
                
                % Compute the logarithmic term
                c3 = log((xii^2 + zii^2) / (xi^2 + zi^2));
                
                % Compute the edge vector
                dx = xii - xi;
                dz = zii - zi;
                denominator = dx^2 + dz^2;
                
                % Avoid division by zero (degenerate edge handling)
                if denominator < 1e-12
                    continue; % Skip degenerate edge
                end
                
                % Compute the numerator term
                numerator = xi * zii - xii * zi;
                
                % Compute the contribution of the current edge
                factor = dz * delta_c + 0.5 * dx * c3;
                f_val = 2 * G * (numerator / denominator) * factor * scale_factor;
                f = f + f_val;
            end
            
            % Store the contribution of the current grid cell to the current observation point
            Fx(k, n) = f;
            n = n + 1;
        end
    end
end

% F*m = g
Gx = Fx * reshape(mod,[],1);

% Compute elapsed time
computation_time = toc;
fprintf('Computation of gravity kernel matrix Fx completed, time is: %.6f seconds.\n', computation_time);

end
