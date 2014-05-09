true_list_file = '/data/flow/isabel_all/2d/all24.list';
fitted_list_file = '/data/flow/isabel_all/2d/fitted_quad_endpoint24/all24.list';
stderr_file = '/data/flow/isabel_all/2d/fitted_quad_endpoint24/sampling24_00_stderr.raw';
global data_w data_d data_t data_h 

auto =0
% pathline_loader
% output: trace_list

% get fields
if 1
    [true_list, data_w, data_h, data_d, data_t, scaling] = load_list(true_list_file);
    fitted_list = load_list(fitted_list_file);
    
    vecAll_fit  = zeros(3, data_w, data_h, 1, data_t);
    vecAll      = zeros(3, data_w, data_h, 1, data_t);
    z=2
    for t=1:data_t
       vec_fit  = load_vec(fitted_list{t});
       vec      = load_vec(true_list{t});

       % trim
       vecAll_fit(:,:,:,:,t)    = vec_fit(:,:,:,z);
       vecAll(:,:,:,:,t)        = vec(:,:,:,z);
    end

    % get error field => err_fitted
    fin = fopen(stderr_file, 'rb');
    var_fitted = fread(fin, [3 data_w*data_h], 'float32').^2;
    var_fitted = reshape(var_fitted, [3 data_w data_h]);
    fclose(fin);
end

data_t = 25; % 9
DT = min(0.5, data_t/20);
    
traced_err_list = [];
filtered_err_list = [];
quad_fit_err_list = [];
linear_fit_err_list = [];
origin_list = [];

fig=figure;

for trace = trace_list
    true_raw_ary = trace{1};
    true_raw_ary = trace_list{27}; % 27

    if true_raw_ary(4,end) < data_t-1.1
        continue;
        disp('skipped')
    end
    
    origin_list(:,end+1) = true_raw_ary(:,1);
   
    % sampling timing
    t_ary = 0:DT:data_t-DT-1; 
    t_ary = t_ary';
    
    
    % recompute true trace
    if 1
        disp('recompute true trace')
        true_raw_ary1 = trace_particle(vecAll, true_raw_ary(1:2,1), t_ary(1), t_ary(end), 0.5);
        plot(true_raw_ary(1,:), true_raw_ary(2,:), true_raw_ary1(1,:), true_raw_ary1(2,:));
        legend('c++ computed', 'matlab computed')    
        true_raw_ary = true_raw_ary1;
        drawnow
    end
    % process trace_list - even trace:    
    uneven_t_ary = true_raw_ary(4,:)';
    true_x_ary = [interp1q(uneven_t_ary, true_raw_ary(1,:)', t_ary) ...
                  interp1q(uneven_t_ary, true_raw_ary(2,:)', t_ary)];
    true_x_ary = true_x_ary';
    
    
    % linear fit
    linterp_x_ary = ([interp1q([t_ary(1);t_ary(end)], [true_x_ary(1,1);true_x_ary(1,end)], t_ary) ...
                    interp1q([t_ary(1);t_ary(end)], [true_x_ary(2,1);true_x_ary(2,end)], t_ary)] )' ;
    
    % quad fit
    [errstd1, errsum1, yy1, a, b, c, rmserr1] = quadfit_endpoint(t_ary,true_x_ary(1,:)',t_ary);
    [errstd2, errsum2, yy2, a, b, c, rmserr2] = quadfit_endpoint(t_ary,true_x_ary(2,:)',t_ary);

    z_ary = [yy1 yy2];
    z_ary=z_ary';
    z_var = [rmserr1^2; rmserr2^2]* ones(1,length(t_ary));

    %%% trace particle
    disp('particle tracing')
    if 1
        traced_fitted_x_ary = trace_particle(vecAll_fit, true_raw_ary(1:2,1), t_ary(1), t_ary(end), 0.5);
        uneven_t_ary = traced_fitted_x_ary(4,:)';
        traced_fitted_x_ary = [interp1q(uneven_t_ary, traced_fitted_x_ary(1,:)', t_ary) ...
                                interp1q(uneven_t_ary, traced_fitted_x_ary(2,:)', t_ary)];
        traced_fitted_x_ary = traced_fitted_x_ary';
    else % save time
        traced_fitted_x_ary = true_x_ary;
    end
    
    disp('filtering')
%     traced_fitted_x_raw_ary=true_raw_ary'; %%[11.000000 11.000000 2.000000 0.000000; 15.784504 11.436985 2.000000 1.000000; 17.084600 11.512656 2.000000 1.250000; 18.414061 11.578578 2.000000 1.500000; 19.773125 11.641154 2.000000 1.750000; 21.158314 11.705994 2.000000 2.000000; 22.555178 11.764643 2.000000 2.250000; 23.968336 11.818457 2.000000 2.500000; 25.417698 11.870479 2.000000 2.750000; 26.905905 11.917877 2.000000 3.000000; 28.409637 11.972345 2.000000 3.250000; 29.909431 12.012804 2.000000 3.500000; 31.437849 12.008892 2.000000 3.750000; 33.005356 11.998375 2.000000 4.000000; 34.573181 11.995104 2.000000 4.250000; 36.156193 11.963609 2.000000 4.500000; 37.761539 11.896452 2.000000 4.750000; 39.372112 11.805254 2.000000 5.000000; 40.969849 11.706592 2.000000 5.250000; 42.560917 11.593554 2.000000 5.500000; 44.161194 11.454815 2.000000 5.750000; 45.768669 11.272562 2.000000 6.000000; 47.391045 11.029262 2.000000 6.250000; 49.026974 10.745399 2.000000 6.500000; 50.657955 10.442777 2.000000 6.750000; 52.266922 10.100547 2.000000 7.000000; 53.861134 9.671535 2.000000 7.250000; 55.444309 9.137972 2.000000 7.500000; 57.030910 8.526555 2.000000 7.750000; 58.645233 7.841821 2.000000 8.000000];
%     traced_fitted_x_ary = [interp1q(traced_fitted_x_raw_ary(:,4), traced_fitted_x_raw_ary(:,1), t_ary) ...
%                            interp1q(traced_fitted_x_raw_ary(:,4), traced_fitted_x_raw_ary(:,2), t_ary)];
%     traced_fitted_x_ary=traced_fitted_x_ary'; 

if 0
 
elseif 1
    if 0  % kalman filter fw bw
        KR_DIR = 1;
        trace_filtering
        fw_post_x_ary = post_x_ary;
        fw_post_x_var = post_x_var;

        KR_DIR = -1;
        trace_filtering
        bw_post_x_ary = post_x_ary;
        bw_post_x_var = post_x_var;
    else  % fw bw
        [fw_post_x_ary fw_post_x_var] = trace_particle_gauss(vecAll_fit, true_x_ary(1:2,1), t_ary(1), t_ary(end), DT, var_fitted);
        [bw_post_x_ary bw_post_x_var] = trace_particle_gauss(vecAll_fit, true_x_ary(1:2,end), t_ary(end), t_ary(1), -DT, var_fitted);
        fw_post_x_ary = fw_post_x_ary(1:2,:);
        bw_post_x_ary = bw_post_x_ary(1:2,length(bw_post_x_ary):-1:1);  % reverse order
        bw_post_x_var = bw_post_x_var(:, length(bw_post_x_ary):-1:1); % reverse order
        
        if 1
            n= length(bw_post_x_ary);
            z_ary = bw_post_x_ary(n:-1:1);
            z_var = bw_post_x_var(:,n:-1:1);
            KR_DIR = 1;
            trace_filtering
            kf_fw_post_x_ary = post_x_ary;
            kf_fw_post_x_var = post_x_var;
            
            n = length(fw_post_x_ary);
            z_ary = fw_post_x_ary(n:-1:1);
            z_var = fw_post_x_var(:,n:-1:1);
            KR_DIR = -1;
            trace_filtering
            bw_post_x_ary = post_x_ary;
            bw_post_x_var = post_x_var;
            
            fw_post_x_ary = kf_fw_post_x_ary;
            fw_post_x_var = kf_fw_post_x_var;
            
        end
    end
        
    % merge
    merged_sum_var = fw_post_x_var+ bw_post_x_var;
    merged_x_ary = (bw_post_x_var.*fw_post_x_ary + fw_post_x_var.*bw_post_x_ary) ./ merged_sum_var;
    merged_x_ary(:,1) = fw_post_x_ary(:,1);
    merged_x_ary(:,end) = bw_post_x_ary(:,end);
    merged_x_var= fw_post_x_var.*bw_post_x_var./merged_sum_var;
    merged_x_var(:,1) = [0;0];
    merged_x_var(:,end) = [0;0];
    
    % statistics
    meanerr_fitted     = mean(sqrt(sum( (true_x_ary(1:2,:)- z_ary).^2 )))
    meanerr_filtered   = mean(sqrt(sum( (true_x_ary(1:2,:)- merged_x_ary).^2)))
    meanerr_traced     = mean(sqrt(sum( (true_x_ary(1:2,:)- traced_fitted_x_ary(1:2,:) ).^2)))
    meanerr_linear_fit = mean(sqrt(sum( (true_x_ary(1:2,:) - linterp_x_ary).^2 )))
    
    
    plot(fw_post_x_ary(1,:), fw_post_x_ary(2,:), bw_post_x_ary(1,:), bw_post_x_ary(2,:), ... 
        true_x_ary(1,:), true_x_ary(2,:),z_ary(1,:), z_ary(2,:), merged_x_ary(1,:), merged_x_ary(2,:), traced_fitted_x_ary(1,:), traced_fitted_x_ary(2,:) ...
        );
    legend('fw', 'bw', ...
        'true xy', sprintf('quad fitted xy, mean err=%f', meanerr_fitted), ...
        sprintf('merged filtered xy, mean err=%f', meanerr_filtered), ...
        sprintf('traced from fitted field, mean err=%f', meanerr_traced) ...
        );
    xlabel('x')
    ylabel('y')
    if auto==0
        disp('press key..')
        pause
    end

    hold on
    for i=1:length(merged_x_var)
        tmp = diag(merged_x_var(:,i));
        tmpx = merged_x_ary(:,i);
        if det(tmp) > 0
            error_ellipse(tmp, tmpx, 'style', 'm');
        end
        tmp = diag(fw_post_x_var(:,i));
        tmpx = fw_post_x_ary(:,i);
        if det(tmp) > 0
            error_ellipse(tmp, tmpx, 'style', 'b');
        end
        tmp = diag(bw_post_x_var(:,i));
        tmpx = bw_post_x_ary(:,i);
        if det(tmp) > 0
            error_ellipse(tmp, tmpx, 'style', 'g');
        end
    end
    hold off
    drawnow
    if auto==0
        disp('press key...')
        pause
    end
    

    
    if 0
        z_ary = merged_x_ary;
        z_var = merged_x_var;
        
        KR_DIR = 1;
        
        trace_filtering
        
        
    end
    
    if auto==1
        traced_err_list(end+1) = meanerr_traced;
        filtered_err_list(end+1) = meanerr_filtered;
        quad_fit_err_list(end+1) = meanerr_fitted;
        linear_fit_err_list(end+1) = meanerr_linear_fit;
        if auto==1
            saveas(fig, sprintf('trace_x%d_y%d.jpg', true_x_ary(1,1), true_x_ary(2,1)));
            disp('saved')
        end
    end
else
    n = length(z_ary);
    s = repmat(struct('P',zeros(2)), 1, n );
    for iter=1:4
        KR_DIR = 1;
        
        trace_filtering

        diff = post_x_ary(:,end)-z_ary(:,end);
        disp(sprintf('plus %f  ', diff))
        if (norm(diff) < 0.001)
            disp('end-point match')
            break;
        end

        % reverse
        KR_DIR = -1;
        z_ary = post_x_ary;
        z_ary(:,end) = true_x_ary(:,end);
        z_var = post_x_var + 3*norm(diff) * [1;1] * (1:-1/(n-1):0); 
%         z_stderr = post_x_stderr + 3*repmat(abs(diff), 1, length(post_x_stderr));
%         z_stderr(:,end) = [0;0];
        
        trace_filtering

        diff = post_x_ary(:,1)-z_ary(:,1);
        disp(sprintf('plus %f', abs(diff)))
        if (norm(diff) < 0.01)
            disp('end-point match')
            break;
        end
            
        % reverse
        z_ary = post_x_ary;
        z_ary(:,1) = true_x_ary(:,1);
        z_var = post_x_var + 3* norm(diff) * [1;1] * (1:-1/(n-1):0); %repmat(abs(diff), 1, length(post_x_stderr));
%         z_stderr = post_x_stderr + 3*repmat(abs(diff), 1, length(post_x_stderr));
%         z_stderr(:,1) = [0;0];
        

%         trace_smoothing
        
    end
    disp('*******************8')
end
    
%     plot(true_raw_ary(1,:), true_raw_ary(2,:))
%     disp('press key...')
%     pause
end

