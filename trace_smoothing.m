% kalman filter testing

% input : 
%   sampling - total n=sampling+1
%   flow field  : V(k) = f(k+1) (f is one based)
%   flow field error
%   trace x_ary : vector
%   trace time t_ary : vector
%   trace error xerrstd: single value
% output :
%   updated trace
%   updated error : vector

if 0
    % get fields
    vecAll_fit = zeros(3, data_w, data_h, 1, data_t);
    z=2
    for t=1:data_t
       vec_fit = load_vec(fitted_list{t});

       % trim
       vecAll_fit(:,:,:,:,t) = vec_fit(:,:,:,z);
    end

    % get error field => err_fitted
    fin = fopen(stderr_file, 'rb');
    err_fitted = fread(fin, [3 data_w*data_h], 'float32');
    err_fitted = reshape(err_fitted, [3 data_w data_h]);
    fclose(fin);
end

%% start smoothing 
% input: post_x_ary, post_x_stderr, z_ary

if length(t_ary)~=length(z_ary)
    error ('time series length mismatch');
end
pret = length(t_ary)-1;
smooth_x_ary = zeros(2,length(t_ary));
smooth_x_ary(:,end) = z_ary(:,end);

smooth_x_stderr = zeros(2,2,length(t_ary));
smooth_x_stderr(:,:,end) = zeros(2);
post_x_stderr(:,end)=[0.0001 0.0001];
for i=length(t_ary)-1:-1:2
    t=t_ary(i);
    z=z_ary(:,i);
    dt = t-pret; 
    %  Vq = interp3(V,Xq,Yq,Zq) assumes a default grid of sample points. 
    %  The default grid points cover the region, X=1:n, Y=1:m, Z=1:p, where [m,n,p] = size(V). 
    % note : as z is useless here, the 3rd dimension is time
    
    A = eye(2);
    L = diag(post_x_stderr(:,i)) *A* inv(diag(post_x_stderr(:,i+1)))
    smooth_x_ary(:,i) = post_x_ary(:,i) + L*( smooth_x_ary(:,i+1) - post_x_ary(:,i+1))
    smooth_x_stderr(:,:,i) = diag(post_x_stderr(:,i)) + L * (smooth_x_stderr(:,:,i+1) - diag(post_x_stderr(:,i+1))) *L'
    
    pret = t;
end
smooth_x_ary(:,1) = z_ary(:,1);  % first element
smooth_x_stderr(:,1) = zeros(2,1);

% comput hausdorff
% haus_fitted     = HausdorffDist(true_x_ary(1:2,:), z_ary);
% haus_filtered   = HausdorffDist(true_x_ary(1:2,:), post_x_ary);
% haus_traced     = HausdorffDist(true_x_ary(1:2,:), traced_fitted_x_ary(1:2,:));
% meanerr_fitted     = mean(sqrt(sum( (true_x_ary(1:2,:)- z_ary).^2 )));
% meanerr_filtered   = mean(sqrt(sum( (true_x_ary(1:2,:)- post_x_ary).^2)));
% meanerr_traced     = mean(sqrt(sum( (true_x_ary(1:2,:)- traced_fitted_x_ary(1:2,:) ).^2)));
meanerr_smoothed     = mean(sqrt(sum( (true_x_ary(1:2,:)- smooth_x_ary(1:2,:) ).^2)));
% x

plot(true_x_ary(1,:), true_x_ary(2,:),z_ary(1,:), z_ary(2,:), post_x_ary(1,:), post_x_ary(2,:), traced_fitted_x_ary(1,:), traced_fitted_x_ary(2,:), smooth_x_ary(1,:), smooth_x_ary(2,:));
legend('true xy', sprintf('interpolated xy, mean err=%f', meanerr_fitted), ...
    sprintf('filtered xy, mean err=%f', meanerr_filtered), ...
    sprintf('traced from fitted field, mean err=%f', meanerr_traced), ...
    sprintf('smoothed, mean err=%f', meanerr_smoothed));
xlabel('x')
ylabel('y')
disp('press key..')
pause

hold on
for i=1:length(t_ary)-1
    si=s(i);
%     post_x_stderr
    if det(si.Q) > 0
        error_ellipse(si.Q*3, si.x, 'style', 'k');
    end
    if det(si.P) > 0
        error_ellipse(si.P*3, si.x, 'style', 'r');
    end
    if det(si.R) > 0
        error_ellipse(si.R*3, si.z, 'style', 'g');
    end
    quiver(s(i).x(1), s(i).x(2), s(i).u(1), s(i).u(2), 'color', 'k');
    
    if (det(smooth_x_stderr(:,:,i))>0)
        error_ellipse(smooth_x_stderr(:,:,i)*3, smooth_x_ary(:,i), 'style', 'c');
    end
end
hold off

    