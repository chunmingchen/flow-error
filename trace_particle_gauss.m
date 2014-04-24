% input:
%   vecAll
%   dt: stepsize 
%   x0: [2,1]
% output:
%   traced_x_ary: [4, n]
function [traced_x_ary, traced_x_var] = trace_particle_gauss(vecAll, x0, t0, t1, dt, varAll)
    t_ary = t0:dt:t1;

    %%% trace particle
    traced_x_ary = zeros(4, length(t_ary));
    traced_x_ary(:,1) = [x0; 0; t0];
    traced_x_var = zeros(2, length(t_ary));
    pre_x = x0;
    pre_t = t0;
    pre_v = [0;0]; % variance
%     pre_x
    for i=2:length(t_ary)
        t = t_ary(i);
        [x, v] = runge_kutta4_gauss(vecAll, pre_x, pre_t, t-pre_t, varAll);
        traced_x_ary(:, i) = [x;0;t];
        v = v + pre_v;
        traced_x_var(:, i) = v;
        pre_x = x; pre_t = t; pre_v = v;
%         x
    end
end