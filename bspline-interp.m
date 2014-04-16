base_list = '/data/flow/isabel_all/all.list';
skipped_list = '/data/flow/isabel_all/all8.list';
fitted_list = '/data/flow/isabel_all/fitted_quad/all.list';
global data_w data_d data_t data_h 

base_files = load_list(base_list);
skipped_files = load_list(skipped_list);
fitted_list = load_list(fitted_list);

data_t = 9;

err_fit = zeros(3, data_w, data_h, 1, data_t);
z=50
for t=1:data_t
   vec = load_vec(base_files{t});
   vec_fitted = load_vec(fitted_list{t});
      
   % trim
   vec = vec(:,:,:,z);
   vec_fitted = vec_fitted(:,:,:,z);

   err_fit(:,:,:,:,t) = vec-vec_fitted;
end


  for dd=1:3
    errOverTimeFit = reshape(err_fit(dd, :,:,:,:), 1, data_t*data_w*data_h);
    figure
%
    len = length(errOverTimeFit);
    min_err = min(errOverTimeFit);
    max_err = max(errOverTimeFit);
    % just plot
    subplot(2,2,1)
    plot(errOverTimeFit);
    title('error over time')

    % kde
    [f,xi] = ksdensity(errOverTimeFit);

    subplot(2,2,2)                
    plot(xi,f);
    title('KDE')

    % linear interp
    d = (max_err-min_err)/50;
    xx=min_err-d: d :max_err+d;
    yy=zeros(1,length(xx));
    for t=1:len-1
        if mod(t,data_t)==0
            continue
        end
        small = min(errOverTimeFit(t+1), errOverTimeFit(t));
        large = max(errOverTimeFit(t+1), errOverTimeFit(t));
        m = (len-1)*(large-small);
        range = xx>=small & xx<=large;
        yy(find(range)) = yy(find(range)) + 1/m;

    end

    subplot(2,2,3)
    plot(xx,yy)
    title('linearly interpolated')

    % histogram
    subplot(2,2,4)
    hist(errOverTimeFit,40)
    title('histogram')
  end
  
return


%%%%%%%%%%%%%%%%%%%% show for every point
figure
for z=1:1
%     for y=1:data_h
%         for x=1:data_w
x=13
y=13
            
            for dd=1:3
                
                errOverTimeFit = squeeze(err_fit(dd,x,y,z,:));

                len = length(errOverTimeFit);
                min_err = min(errOverTimeFit);
                max_err = max(errOverTimeFit);
                % just plot
                subplot(2,3,1)
                plot(errOverTimeFit);
                title('error over time')
                
                % kde
                [f,xi] = ksdensity(errOverTimeFit);
                
                subplot(2,3,2)                
                plot(xi,f);
                title('KDE')
                
                % linear interp
                d = (max_err-min_err)/50;
                xx=min_err-d: d :max_err+d;
                yy=zeros(1,length(xx));
                for t=1:len-1
                    small = min(errOverTimeFit(t+1), errOverTimeFit(t));
                    large = max(errOverTimeFit(t+1), errOverTimeFit(t));
                    m = (len-1)*(large-small);
                    range = xx>=small & xx<=large;
                    yy(find(range)) = yy(find(range)) + 1/m;
                    
                end
                
                subplot(2,3,3)
                plot(xx,yy)
                title('linearly interpolated')
                
                % histogram
                subplot(2,3,4)
                hist(errOverTimeFit)
                title('histogram')
                
                %
                
                xx=1:len/100:len;
                yy=spline(1:len, errOverTimeFit, xx)
                subplot(2,3,5)
                plot(xx,yy,'g');
                
                %
                subplot(2,3,6)
                hist(yy,20)
                title('histogram of b-spline')
                
                pause
            end
%         end
%     end
end


% err_std = std(err, 1, 5);
% 
% img = squeeze(err_std(1,:,:));
% imagesc(img');
% axis image;
% axis xy;
% colormap(jet);
% 
% % gen interp
% vec1=load_vec(skipped_files{1});
% vec2=load_vec(skipped_files{2});
% for t=1:data_t
%     vec = load_vec(base_files{t});
%     vec_interp = vec1*(data_t-t)/(data_t-1) + vec2*(t-1)/(data_t-1);
%     
%     % trim
%     vec = vec(:,:,:,50);
%     vec_interp = vec_interp(:,:,:,50);
%     
%     err(:,:,:,:,t) = vec-vec_fitted;
% end
% err_std = std(err, 1, 5);
% 
% figure
% img = squeeze(err_std(1,:,:));
% imagesc(img');
% axis image;
% axis xy;
% colormap(jet);