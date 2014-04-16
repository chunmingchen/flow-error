%runge_kutta4(vecAll_fit, prex+1, pret+1, dt);
function y = runge_kutta4(vecAll, x, t0, dt)
    f1 = f(vecAll, x, t0);
    f2 = f(vecAll, x + dt * f1 / 2, t0 + dt / 2);
    f3 = f(vecAll, x + dt * f2 / 2, t0 + dt / 2);
    f4 = f(vecAll, x + dt * f3, t0 + dt);
    y = x + dt * (f1 + 2.0 * f2 + 2.0 * f3 + f4) / 6.0;
end

function y = f(vecAll, x, t)
    y=zeros(2,1);
    y(1) = interp3(squeeze(vecAll(1,:,:,:,:)),   x(2)+1, x(1)+1, t+1);
    y(2) = interp3(squeeze(vecAll(2,:,:,:,:)),   x(2)+1, x(1)+1, t+1);
    %interp_v(3) = interp3(squeeze(vecAll_fit(3,:,:,:,:)),   prex(2)+1, prex(1)+1, pret+1)*dt;
    if isnan(y(1))
        disp('nan')
    end
end

