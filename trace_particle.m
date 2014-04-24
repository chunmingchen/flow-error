% input:
%   vecAll
%   dt: stepsize 
%   x0: [2,1]
% output:
%   traced_x_ary: [4, n]
function traced_x_ary = trace_particle(vecAll, x0, t0, t1, dt)
    t_ary = t0:dt:t1;

    %%% trace particle
    traced_x_ary = zeros(4, length(t_ary));
    traced_x_ary(:,1) = [x0; 0; t0];
    pre_x = x0;
    pre_t = t0;
%     pre_x
    for i=2:length(t_ary)
        t = t_ary(i);
        x = runge_kutta4(vecAll, pre_x, pre_t, t-pre_t);
        traced_x_ary(:, i) = [x;0;t];
        pre_x = x; pre_t = t;
%         x
    end
end