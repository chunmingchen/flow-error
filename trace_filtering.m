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

%% start filtering

if length(t_ary)~=length(z_ary)
    error ('time series length mismatch');
end
pret = 0;
prex = z_ary(:,1);
clear s
s.P = zeros(2);
post_x_ary = zeros(2,length(t_ary));
post_x_ary(:,1) = prex;
for i=2:length(t_ary)-1
    t=t_ary(i);
    z=z_ary(:,i);
    dt = t-pret; 
    %  Vq = interp3(V,Xq,Yq,Zq) assumes a default grid of sample points. 
    %  The default grid points cover the region, X=1:n, Y=1:m, Z=1:p, where [m,n,p] = size(V). 
    % note : as z is useless here, the 3rd dimension is time
%     interp_v=zeros(2,1);
%     interp_v(1) = interp3(squeeze(vecAll_fit(1,:,:,:,:)),   prex(2)+1, prex(1)+1, pret+1)*dt;
%     interp_v(2) = interp3(squeeze(vecAll_fit(2,:,:,:,:)),   prex(2)+1, prex(1)+1, pret+1)*dt;
    %interp_v(3) = interp3(squeeze(vecAll_fit(3,:,:,:,:)),   prex(2)+1, prex(1)+1, pret+1)*dt;
    stderr_v = zeros(2,1);
    stderr_v(1) = interp2(squeeze(err_fitted(1,:,:,:)), prex(2)+1, prex(1)+1)*dt;
    stderr_v(2) = interp2(squeeze(err_fitted(2,:,:,:)), prex(2)+1, prex(1)+1)*dt;
    %stderr_v(3) = interp2(squeeze(err_fitted(3,:,:,:)), prex(2)+1, prex(1)+1)*dt;
        
    % build kalmanf input
    s(i-1).x = prex;
    s(i-1).A = eye(2);
%     if i==2
%         s(i-1).Q = zeros(2);
%     else
        s(i-1).Q = diag(stderr_v); % variance of fitted flow field
%     end
    s(i-1).H = eye(2);
    s(i-1).R = diag(z_stderr(:,i)); % variance of fitted trace
    s(i-1).B = eye(2);
    s(i-1).u = runge_kutta4(vecAll_fit, prex, pret, dt)-prex-1;  %interp_v;
    s(i-1).z = z_ary(:,i);
        
    s(i)=kalmanf(s(i-1));
    
    prex = s(i).x;
    post_x_ary(:,i) = s(i).x;
    pret = t;
end
post_x_ary(:,length(t_ary)) = z_ary(:,length(t_ary)); % last element

% comput hausdorff
% haus_fitted     = HausdorffDist(true_x_ary(1:2,:), z_ary);
% haus_filtered   = HausdorffDist(true_x_ary(1:2,:), post_x_ary);
% haus_traced     = HausdorffDist(true_x_ary(1:2,:), traced_fitted_x_ary(1:2,:));
meanerr_fitted     = mean(sqrt(sum( (true_x_ary(1:2,:)- z_ary).^2 )));
meanerr_filtered   = mean(sqrt(sum( (true_x_ary(1:2,:)- post_x_ary).^2)));
meanerr_traced     = mean(sqrt(sum( (true_x_ary(1:2,:)- traced_fitted_x_ary(1:2,:) ).^2)));
% x

plot(true_x_ary(1,:), true_x_ary(2,:),z_ary(1,:), z_ary(2,:), post_x_ary(1,:), post_x_ary(2,:), traced_fitted_x_ary(1,:), traced_fitted_x_ary(2,:));
legend('true xy', sprintf('interpolated xy, mean err=%f', meanerr_fitted), ...
    sprintf('filtered xy, mean err=%f', meanerr_filtered), ...
    sprintf('traced from fitted field, mean err=%f', meanerr_traced));
xlabel('x')
ylabel('y')
pause
hold on
post_x_stderr=zeros(size(post_x_ary));
for i=1:length(t_ary)-1
    si=s(i);
    post_x_stderr(1,i)=si.P(1,1);
    post_x_stderr(2,i)=si.P(2,2);
    post_x_stderr
    if det(si.P) > 0
        error_ellipse(si.P*3, si.x, 'style', 'r');
    end
    if det(si.R) > 0
        error_ellipse(si.R*3, si.z, 'style', 'g');
    end
    quiver(s(i).x(1), s(i).x(2), s(i).u(1), s(i).u(2), 'color', 'k');
end
hold off

    