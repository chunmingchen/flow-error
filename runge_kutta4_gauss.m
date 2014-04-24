%runge_kutta4(vecAll_fit, prex+1, pret+1, dt);
% stdAll: gaussian error standard deviation
function [y v] = runge_kutta4_gauss(vecAll, x, t0, dt, varAll)
    f1 = f(vecAll, x, t0);
    f2 = f(vecAll, x + dt * f1 / 2, t0 + dt / 2);
    f3 = f(vecAll, x + dt * f2 / 2, t0 + dt / 2);
    f4 = f(vecAll, x + dt * f3, t0 + dt);
    y = x + dt * (f1 + 2.0 * f2 + 2.0 * f3 + f4) / 6.0;
    
    v = [interp2(squeeze(varAll(1,:,:,:)), x(2)+1, x(1)+1) 
         interp2(squeeze(varAll(2,:,:,:)), x(2)+1, x(1)+1)] * abs(dt)  ; %%!! approximate
    %varAll(1:2, floor(x(1)), floor(x(2))) * abs(dt);
% %     
%      stderr_v = [interp2(squeeze(err_fitted(1,:,:,:)), prex(2)+1, prex(1)+1)*dt
%                  interp2(squeeze(err_fitted(2,:,:,:)), prex(2)+1, prex(1)+1)*dt];
end

function y = f(vecAll, x, t)
    y=zeros(2,1);
    y(1) = interp3(squeeze(vecAll(1,:,:,:,:)),   x(2)+1, x(1)+1, t+1); % interp3 accesses (y,x,z,...)
    y(2) = interp3(squeeze(vecAll(2,:,:,:,:)),   x(2)+1, x(1)+1, t+1);
    %interp_v(3) = interp3(squeeze(vecAll_fit(3,:,:,:,:)),   prex(2)+1, prex(1)+1, pret+1)*dt;
    if isnan(y(1))
        disp('Runge Kutta: nan')
    end
%     disp(sprintf('vec(%f %f %f)=(%f %f)', x(1), x(2), t, y(1), y(2)))
end

