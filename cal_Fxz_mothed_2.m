function [Fxz,Gxz] = cal_Fxz_mothed_2(mod, gobs, nodex, nodez)
% Initialize the output matrix
Fxz = zeros(size(gobs, 1), size(mod, 1) * size(mod, 2));

% Physical constants
G = 6.67430e-11;
scale_factor = 1e12;

% Singularity handling thresholds
eps = 1e-10;        % Singularity check threshold
dx_eps = 1e-10;     % Vertical edge threshold
dz_eps = 1e-10;     % Horizontal edge threshold

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
                
                % Compute edge vector
                dx = xii - xi;
                dz = zii - zi;
                
                % 2. Vertical edge check (dx ¡Ö 0)
                if abs(dx) < dx_eps
                    % Vertical edge formula
                    c3 = log((xii^2 + zii^2) / (xi^2 + zi^2));
                    f_val = -G * c3 * scale_factor;
                    f = f + f_val;
                    continue; % Skip remaining calculations after handling
                end
                
                % 3. Horizontal edge check (dz ¡Ö 0)
                if abs(dz) < dz_eps
                    % Horizontal edge formula
                    f_val = 0;
                    f = f + f_val;
                    continue; % Skip remaining calculations after handling
                end
                
                % 4. General case (non-vertical and non-horizontal edge)
                % Compute angle (use atan2 to handle quadrant ambiguity)
                c1 = atan2(zi, xi);
                c2 = atan2(zii, xii);
                
                % Adjust the angle to the range 0¨C2¦Ð
                if c1 < 0, c1 = c1 + 2*pi; end
                if c2 < 0, c2 = c2 + 2*pi; end
                
                % Handle angle difference wrapping
                delta_c = c2 - c1;
                if delta_c > pi
                    delta_c = delta_c - 2*pi;
                elseif delta_c < -pi
                    delta_c = delta_c + 2*pi;
                end
                
                % Compute the logarithmic term (ln)
                c3 = log((xii^2 + zii^2) / (xi^2 + zi^2));
                
                % Compute denominator
                denominator = dx^2 + dz^2;
                
                % Avoid division by zero (degenerate edge handling)
                if denominator < 1e-12
                    continue; % Skip degenerate edge
                end
                
                % Compute numerator term
                numerator = dz;
                
                % Compute the contribution of the current edge
                factor = dx * delta_c - 0.5 * dz * c3;
                f_val = 2 * G * (numerator / denominator) * factor * scale_factor;
                f = f + f_val;
            end
            
            % Store the contribution of the current grid cell for the current observation point
            Fxz(k, n) = f;
            n = n + 1;
        end
    end
end

% F*m = g
Gxz = Fxz * reshape(mod,[],1);

% Compute elapsed time
computation_time = toc;
fprintf('Computation of gravity kernel matrix Fxz completed, time is: %.6f seconds.\n', computation_time);

end
