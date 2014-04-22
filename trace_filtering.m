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
KR_DIR

if length(t_ary)~=length(z_ary)
    error ('time series length mismatch');
end
n = length(t_ary);
if KR_DIR~=-1
    pret = t_ary(1);
    prex = z_ary(:,1);
    i_list = 2:n;
    clear s
    s.P = zeros(2);
    post_x_ary = zeros(2,n);
    post_x_ary(:,1) = prex;
    post_x_var = zeros(2,n);
else
    pret = t_ary(end);
    prex = z_ary(:,end);    
    i_list = n-1: -1: 1;
    s = repmat(struct('P',zeros(2)), 1, n );
    post_x_ary = zeros(2,n);
    post_x_ary(:,end) = prex;
    post_x_var = zeros(2,n);
end
for i=i_list
    t=t_ary(i);
    z=z_ary(:,i);
    dt = t-pret; 
    %  Vq = interp3(V,Xq,Yq,Zq) assumes a default grid of sample points. 
    %  The default grid points cover the region, X=1:n, Y=1:m, Z=1:p, where [m,n,p] = size(V). 
    % note : as z is useless here, the 3rd dimension is time
    
    [rk4_y, rk4_var] = runge_kutta4_gauss(vecAll_fit, prex, pret, dt, var_fitted);  %interp_v;
%     if (KR_DIR==1 && i==2 || KR_DIR==-1 && i==n-1)
%         rk4_var=[0;0];
%     end
        
    % build kalmanf input
    s(i-KR_DIR).x = prex;
    s(i-KR_DIR).A = eye(2);
    s(i-KR_DIR).Q = diag(rk4_var); % variance of fitted flow field
    s(i-KR_DIR).H = eye(2);
    s(i-KR_DIR).R = diag(z_var(:,i)); % variance of fitted trace
    s(i-KR_DIR).B = eye(2);
    s(i-KR_DIR).u = rk4_y-prex;  %interp_v;
    s(i-KR_DIR).z = z_ary(:,i);
        
    s(i)=kalmanf(s(i-KR_DIR));
    
    prex = s(i).x;
    post_x_ary(:,i) = s(i).x;
    post_x_var(:,i) = [s(i).Q(1,1); s(i).Q(2,2)];
    pret = t;
end
if KR_DIR~=-1
%     post_x_ary(:,end) = z_ary(:,end); % last element
%     post_x_var(:,end) = zeros(2,1);
else
%     post_x_ary(:,1) = z_ary(:,1); % last element
%     post_x_var(:,1) = zeros(2,1);
%     s(1).Q = zeros(2);
end

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
disp('press key..')
pause

hold on
post_x_stderr=zeros(size(post_x_ary));
if KR_DIR~=-1
    i_list = 1:n-1;
else
    i_list = 2:n;
end
for i=i_list
    si=s(i);
%     post_x_stderr
    if det(si.Q) > 0
        error_ellipse(si.Q, si.x, 'style', 'k');
    end
    if det(si.P) > 0
        error_ellipse(si.P, si.x, 'style', 'r');
    end
    if det(si.R) > 0
        error_ellipse(si.R, si.z, 'style', 'g');
    end
    quiver(s(i).x(1), s(i).x(2), s(i).u(1), s(i).u(2), 'color', 'k');
end
hold off
disp('press key...')
pause

    