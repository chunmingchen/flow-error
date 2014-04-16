base_list = '/data/flow/isabel_all/all.list';
skipped_list = '/data/flow/isabel_all/all8.list';
fitted_list = '/data/flow/isabel_all/fitted_quad/all.list';
fitted_endpoint_list = '/data/flow/isabel_all/fitted_quad_endpoint/all.list';
global data_w data_d data_t data_h 

base_files = load_list(base_list);
skipped_files = load_list(skipped_list);
fitted_list = load_list(fitted_list);
fitted_endpoint_list = load_list(fitted_endpoint_list);

data_t = 9;

% fitted
err = zeros(3, data_w, data_h, 1, data_t);
for t=1:data_t
   vec = load_vec(base_files{t});
   vec_fitted = load_vec(fitted_list{t});
      
   % trim
   vec = vec(:,:,:,50);
   vec_fitted = vec_fitted(:,:,:,50);

   err(:,:,:,:,t) = vec-vec_fitted;
   
end
err_std = std(err, 1, 5);
err_std_overtime = std(reshape(err(1, :,:,:,:), data_w*data_h, data_t))

img = squeeze(err_std(1,:,:));
imagesc(img');
axis image;
axis xy;
colormap(jet);


% fitted_endpoint
err = zeros(3, data_w, data_h, 1, data_t);
for t=1:data_t
   vec = load_vec(base_files{t});
   vec_fitted_endpoint = load_vec(fitted_endpoint_list{t});
      
   % trim
   vec = vec(:,:,:,50);
   vec_fitted_endpoint = vec_fitted_endpoint(:,:,:,50);

   err(:,:,:,:,t) = vec-vec_fitted_endpoint;
   
end
err_std = std(err, 1, 5);
err_std_overtime = std(reshape(err(1, :,:,:,:), data_w*data_h, data_t))

img = squeeze(err_std(1,:,:));
imagesc(img');
axis image;
axis xy;
colormap(jet);


% gen interp
vec1=load_vec(skipped_files{1});
vec2=load_vec(skipped_files{2});
for t=1:data_t
    vec = load_vec(base_files{t});
    vec_interp = vec1*(data_t-t)/(data_t-1) + vec2*(t-1)/(data_t-1);
    
    % trim
    vec = vec(:,:,:,50);
    vec_interp = vec_interp(:,:,:,50);
    
    err(:,:,:,:,t) = vec-vec_fitted;
end
err_std = std(err, 1, 5);

figure
img = squeeze(err_std(1,:,:));
imagesc(img');
axis image;
axis xy;
colormap(jet);

% gen polyfit
slice_polyfit = 
for x=1:data_w
    for y=1:data_h
        
vec_polyfit = 
for t=1:data_t
    vec = load_vec(base_files{t});
    