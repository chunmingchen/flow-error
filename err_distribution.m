global data_w data_d data_t data_h 
dataid=2
if dataid==1
    base_list = '/data/flow/isabel_all/2d/all.list';
    % skipped_list = '/data/flow/isabel_all/all8.list';
    fitted_list = '/data/flow/isabel_all/2d/fitted_quad_endpoint8/all8.list';
    sampling = 8; % sampling+1
elseif dataid==2
    base_list = '/data/flow/isabel_all/2d/all.list';
    fitted_list = '/data/flow/isabel_all/2d/fitted_quad_endpoint24/all.list';
    sampling = 24; % sampling+1
end

base_files = load_list(base_list);
% skipped_files = load_list(skipped_list);
fitted_list = load_list(fitted_list);

data_t = sampling+1;

vecAll = zeros(3, data_w, data_h, 1, data_t);
vecAll_fit = zeros(3, data_w, data_h, 1, data_t);
z=3
for t=1:data_t
   vec = load_vec(base_files{t});
   vec_fitted = load_vec(fitted_list{t});
      
   % trim
   vecAll(:,:,:,:,t) = vec(:,:,:,z);
   vecAll_fit(:,:,:,:,t) = vec_fitted(:,:,:,z);

end
err_fit(:,:,:,:,t) = vecAll-vecAll_fit;


%  for dd=1:3
%     errOverTimeFit = reshape(err_fit(dd, :,:,:,:), 1, data_t*data_w*data_h);
%     figure
% %
%     len = length(errOverTimeFit);
%     min_err = min(errOverTimeFit);
%     max_err = max(errOverTimeFit);
%     % just plot
%     subplot(2,2,1)
%     plot(errOverTimeFit);
%     title('error over time')
% 
%     % kde
%     [f,xi] = ksdensity(errOverTimeFit);
% 
%     subplot(2,2,2)                
%     plot(xi,f);
%     title('KDE')
% 
%     % linear interp
%     d = (max_err-min_err)/50;
%     xx=min_err-d: d :max_err+d;
%     yy=zeros(1,length(xx));
%     for t=1:len-1
%         if mod(t,data_t)==0
%             continue
%         end
%         small = min(errOverTimeFit(t+1), errOverTimeFit(t));
%         large = max(errOverTimeFit(t+1), errOverTimeFit(t));
%         m = (len-1)*(large-small);
%         range = xx>=small & xx<=large;
%         yy(find(range)) = yy(find(range)) + 1/m;
% 
%     end
% 
%     subplot(2,2,3)
%     plot(xx,yy)
%     title('linearly interpolated')
% 
%     % histogram
%     subplot(2,2,4)
%     hist(errOverTimeFit,40)
%     title('histogram')
%   end
%   
% return


%%%%%%%%%%%%%%%%%%%% show for every point

% hist for every block
new_samples = 2400;
figure
bw = 10;
x0=277
y0=188
for dd=1:1  %dim
    for yy=y0:bw:data_h-bw
        for xx=x0:bw:data_w-bw
    
    refx = xx+floor(bw/2);
    refy = yy+floor(bw/2);
    refz = 1;
    refVecOverTime = squeeze(vecAll(dd,refx,refy,refz,:));
    refVecOverTimeFit = squeeze(vecAll_fit(dd,refx,refy,refz,:));
    tt = 1:(data_t-1)/(new_samples-1):data_t; tt=tt';
    t0 = 1:data_t; t0=t0';
    % lenear fit
    refVecOverTimeInterp = interp1q(t0, refVecOverTime, tt);

    % poly fit
    p = polyfit(t0, refVecOverTimeFit, 2);
    refVecOverTimeInterpFit = polyval(p, tt);
    % error between poly fit and linear fit
    refErrOverTimeFit = abs(refVecOverTimeInterp - refVecOverTimeInterpFit);
    refErrOverTimeFitPosNeg = refVecOverTimeInterp - refVecOverTimeInterpFit;
    
    [val, mappingPos] = sort(refErrOverTimeFitPosNeg);
%     [val, mappingPos] = sort(refVecOverTimeInterpFit);
            
    refrms = sqrt(mean(refErrOverTimeFitPosNeg.*refErrOverTimeFitPosNeg));

    subplot(5,1,1)
    hold off
    subplot(5,1,2)
    hold off
    subplot(5,1,3)
    hold off
    for z=1:1
         for y=yy:yy+bw-1
             for x=xx:xx+bw-1
                vecOverTime = squeeze(vecAll(dd,x,y,z,:));
                vecOverTimeFit = squeeze(vecAll_fit(dd,x,y,z,:));
                tt = 1:(data_t-1)/(new_samples-1):data_t; tt=tt';
                t0 = 1:data_t; t0=t0';
                
                % lenear fit
                vecOverTimeInterp = interp1q(t0, vecOverTime, tt);
                
                % poly fit
                p = polyfit(t0, vecOverTimeFit, 2);
                vecOverTimeInterpFit = polyval(p, tt);
                
                % error between poly fit and linear fit
                errOverTimeFit = abs(vecOverTimeInterp - vecOverTimeInterpFit);
                errOverTimeFitPosNeg = vecOverTimeInterp - vecOverTimeInterpFit;
                
                rms = sqrt(mean(errOverTimeFitPosNeg.*errOverTimeFitPosNeg));
                ratio = rms / refrms;
                
                % residual
                subplot(5,1,1);
                plot(tt,  errOverTimeFitPosNeg);
                hold on  
                
                subplot (5,1,2);
                [f,xi] = ksdensity(errOverTimeFitPosNeg);
                plot(xi,f,'.')
                hold on
%                 plot(tt,errOverTimeFitPosNeg, '.r');

                subplot (5,1,3);
                X1 = errOverTimeFitPosNeg;  Y1=refErrOverTimeFitPosNeg;
                X2 = vecOverTimeInterp; Y2 = refVecOverTimeInterp;
                plot(X1, Y1, '.', X2, Y2, '.');
%                 plot(errOverTimeFitPosNeg, refVecOverTimeInterpFit, '.');
                covariance1 = sum( (X1-mean(X1)).*(Y1-mean(Y1)) ) / (length(Y1)-1);
                correlation1 = covariance1 / std(Y1) / std(X1);
                covariance2 = sum( (X2-mean(X2)).*(Y2-mean(Y2)) ) / (length(Y2)-1);
                correlation2 = covariance2 / std(Y2) / std(X2);
                title (sprintf('correlation-err=%f correlation-val=%f', correlation1, correlation2))
                
%                 
%                 % sort
%                 [val mappingPos2] = sort(errOverTimeFitPosNeg);
% %                 [val mappingPos2] = sort(errOverTimeFitPosNeg);
%                 orderedErrOverTimeFitPosNeg = zeros(size(errOverTimeFitPosNeg));
%                 orderedErrOverTimeFitPosNeg(mappingPos) = val;
%                 subplot (5,1,4);
%                 
%                 
%                 plot(tt, errOverTimeFitPosNeg, tt, orderedErrOverTimeFitPosNeg, tt, refErrOverTimeFitPosNeg, tt, refErrOverTimeFitPosNeg*ratio);
%                 legend('truth', 'ordered sample', 'ref', 'rms ratio')
%                 title (sprintf('mean error: ordered: %f, ref: %f, rms ratio: %f', ...
%                     mean(abs(orderedErrOverTimeFitPosNeg-errOverTimeFitPosNeg)), ...
%                     mean(abs(refErrOverTimeFitPosNeg-errOverTimeFitPosNeg)), ...
%                     mean(abs(refErrOverTimeFitPosNeg*ratio - errOverTimeFitPosNeg)) ))
%                 
%                 % sort from gaussian
%                 subplot(5,1,5)
%                 orderedErrOverTimeFitPosNegList = zeros(length(mappingPos), 20);
%                 for i=1:20
%                     val = sort( rms * randn(new_samples, 1) );
% 
%                     orderedErrOverTimeFitPosNegList(mappingPos, i) = val;
%                     
%                     orderedErrOverTimeFitPosNeg1 = orderedErrOverTimeFitPosNeg;
%                     LWIDTH=5;
%                     for j=1+LWIDTH:length(orderedErrOverTimeFitPosNeg)-LWIDTH
%                         orderedErrOverTimeFitPosNeg(j)=sum(orderedErrOverTimeFitPosNeg1(j-LWIDTH:j+LWIDTH))/(LWIDTH*2+1);
%                     end
%                                     
%                 end     
%                 orderedErrOverTimeFitPosNegGauss = mean(orderedErrOverTimeFitPosNegList, 2);
%                 
%                 % sort from histogram
%                 minval = min (errOverTimeFitPosNeg);
%                 maxval = max (errOverTimeFitPosNeg);
%                 bins = minval: (maxval-minval)/20: maxval;
%                 hist = histc(errOverTimeFitPosNeg, bins);
%                 orderedErrOverTimeFitPosNegList = zeros(length(mappingPos), 20);
%                 for i=1:20
%                     val = sort( sample_histogram(new_samples, bins', hist) );
%                     orderedErrOverTimeFitPosNegList(mappingPos, i) = val;
% %                     plot(tt,errOverTimeFitPosNeg, tt,orderedErrOverTimeFitPosNegList(:,i));
% %                     hold on
%                 end
%                 orderedErrOverTimeFitPosNegHist = mean(orderedErrOverTimeFitPosNegList, 2);
%                 
%                 plot(tt, errOverTimeFitPosNeg, tt, orderedErrOverTimeFitPosNegGauss, tt, orderedErrOverTimeFitPosNegHist);
%                 legend('truth', 'gaussian', 'histogram')
%                 title (sprintf('mean error: gaussian: %f, histogram: %f', mean(abs(orderedErrOverTimeFitPosNegGauss - errOverTimeFitPosNeg)), ...
%                     mean(abs(orderedErrOverTimeFitPosNegHist - errOverTimeFitPosNeg))));
%                 hold on
%                 hold off
%                 
%                 
                pause
             end
         end
    end
%     title(sprintf('dimension %d', dd))
    hold off
    
    pause
        end
    end
end

new_samples = 25;
figure
rmse_list=[];
count = 1;
for z=1:1
     for y=11:20:data_h
         for x=11:20:data_w
% x=277;
% y=188;
% x=11;
% y=11;
            for kind=1:2
                
                for dd=1:2
                    tt = 0:(data_t-1)/(new_samples-1):data_t-1; tt=tt';
                    t0 = 0:data_t-1; t0=t0';
                    if kind==1
                        vecOverTime = squeeze(vecAll(dd,x,y,z,:));
                        vecOverTimeFit = squeeze(vecAll_fit(dd,x,y,z,:));
                    else % radius
                        vecOverTimeCartesian = squeeze(vecAll(:,x,y,z,:));
                        magOverTime = sqrt(sum(vecOverTimeCartesian.*vecOverTimeCartesian));
                        thetaOverTime = atan(vecOverTimeCartesian(2,:)./ vecOverTimeCartesian(1,:));
                        [errstd1, errsum1, magOverTimeFit, a, b, c] = quadfit_endpoint(t0,magOverTime',tt); 
%                         magOverTimeFit = magOverTime';
                        magOverTimeFit=magOverTimeFit';
                        [errstd1, errsum1, thetaOverTimeFit, a, b, c] = quadfit_endpoint(t0,thetaOverTime',tt);
%                         thetaOverTimeFit=thetaOverTime';
                        thetaOverTimeFit=thetaOverTimeFit';
                        vecOverTimeFitCartesian(1,:) = magOverTimeFit.*cos(thetaOverTimeFit);
                        vecOverTimeFitCartesian(2,:) = magOverTimeFit.*sin(thetaOverTimeFit);
                        
                        
                        if dd==1
                            vecOverTime = magOverTime';
                            vecOverTimeFit = magOverTimeFit';
                        else
                            vecOverTime = thetaOverTime';
                            vecOverTimeFit = thetaOverTimeFit';
                        end
                        
                    end

                    % lenear fit
                    vecOverTimeInterp = interp1q(t0, vecOverTime, tt);

                    % poly fit
                    p = polyfit(t0, vecOverTimeFit, 2);
                    vecOverTimeInterpFit = polyval(p, tt);

                    % error between poly fit and linear fit
                    errOverTimeInterpFit = abs(vecOverTimeInterp - vecOverTimeInterpFit);
                    errOverTimeInterpFitPosNeg = vecOverTimeInterp - vecOverTimeInterpFit;

    %                 % Cubic B-spline fit
    %                 vecOverTimeSpline = spline(t0, vecOverTime, tt)
    %                 
                    if kind==1
                        subplot(3,2,dd)
                        plot(t0, vecOverTime, 'o-', tt, vecOverTimeInterpFit, '--')
                        title(sprintf('vec x/y over time. mse=%g', mean(errOverTimeInterpFit.^2)))
                    else
                        subplot(3,2,2+dd)
                        
                        plot(t0, vecOverTime, 'o-', tt, vecOverTimeInterpFit, '--')
                        title(sprintf('mag/theta over time. mse=%g', mean(errOverTimeInterpFit.^2)))
                        
                        subplot(3,2,4+dd)                        
                        plot(t0, vecOverTimeCartesian(dd,:), 'o-', tt, vecOverTimeFitCartesian(dd,:), '--')
                        title(sprintf('converted vec x/y over time, mse=%g', mean((vecOverTimeCartesian(dd,:) - vecOverTimeFitCartesian(dd,:)).^2)))
                    end

                    legend('Original', 'Quadratic')
                    
                end
            end
            pause
         end
     end
end

new_samples = 24;
figure
rmse_list=[];
count = 1;
for z=1:1
     for y=11:20:data_h
         for x=11:20:data_w
x=277;
y=188;
% x=11;
% y=11;
            
            for dd=1:2
                vecOverTime = squeeze(vecAll(dd,x,y,z,:));
                vecOverTimeFit = squeeze(vecAll_fit(dd,x,y,z,:));
                tt = 1:(data_t-1)/(new_samples-1):data_t; tt=tt';
                t0 = 1:data_t; t0=t0';
                
                % lenear fit
                vecOverTimeLinear = interp1q(t0, vecOverTime, tt);
                
                % poly fit
                p = polyfit(t0, vecOverTimeFit, 2);
                vecOverTimeInterpFit = polyval(p, tt);
                
                % error between poly fit and linear fit
                errOverTimeFit = abs(vecOverTimeLinear - vecOverTimeInterpFit);
                errOverTimeFitPosNeg = vecOverTimeLinear - vecOverTimeInterpFit;
                
%                 % Cubic B-spline fit
%                 vecOverTimeSpline = spline(t0, vecOverTime, tt)
%                 
                subplot(3,3,1)
                plot(t0, vecOverTime, 'o-', tt, vecOverTimeInterpFit, '--')
                
                legend('Original', 'Quadratic')
                title('vec over time')
                

                len = length(errOverTimeFitPosNeg);
                min_err = min(errOverTimeFitPosNeg);
                max_err = max(errOverTimeFitPosNeg);
                % just plot
                subplot(3,3,2)
                plot(errOverTimeFitPosNeg, '.');
                title('error over time')
                
                % abs
                subplot(3,3,3)
                hist(errOverTimeFit)
                hold on
                
                [f,xi] = ksdensity(errOverTimeFit);
                plot(xi,f, 'r');
                hold off
                title('residual')
                
                % residual
                subplot(3,3,4)                
                hist(errOverTimeFitPosNeg)
                hold on
                
                [f,xi] = ksdensity(errOverTimeFitPosNeg);
                plot(xi,f, 'r');
                hold off
                title('KDE')

                % historgram inference
                subplot(3,3,5)
                minvec = min(vecOverTime);
                maxvec = max(vecOverTime);
                xi =  minvec:(maxvec-minvec)/256:maxvec;
                bincount = histc(vecOverTime, xi);
                
                minerr = min(errOverTimeFitPosNeg);
                maxerr = max(errOverTimeFitPosNeg);
                err_xi =  minerr:(maxerr-minerr)/256:maxerr;
                
                estErrDist = MCErrHist(xi, bincount, p, err_xi);
                
                errDist = histc(errOverTimeFitPosNeg, err_xi);
                errDist = errDist / sum(errDist);
                plot(err_xi, errDist, err_xi, estErrDist)
                
                %
                Y = vecOverTimeLinear;
                fX = vecOverTimeInterpFit;
                covarianceY_fX = sum( (Y-mean(Y)).*(fX-mean(fX)) ) / (length(Y)-1)
                correlationY_fX = covarianceY_fX / std(Y) / std(fX)
                
                scatter(Y,fX)
                title (sprintf('cov = %f, cor = %f', covarianceY_fX, correlationY_fX))
                xlabel('Y')
                ylabel('f(X)')

                 
                % Use RMS as std
                rms = sqrt(mean(errOverTimeFitPosNeg.^2));
                max_err = max(abs(errOverTimeFitPosNeg));
                xi =  -max_err:max_err/40:max_err;
                
                subplot(3,3,6)
                pdf = normpdf(xi, 0, rms);
%                 hist(errOverTimeFitPosNeg, 30);
                bincount = histc(errOverTimeFitPosNeg, xi);
%                 hold on
                plot(xi, bincount, '.', xi, pdf, '-r');
                
                % vs real mean and var
                hold on
                realmean = mean(errOverTimeFitPosNeg);
                realstd = std(errOverTimeFitPosNeg, 1);
                pdf = normpdf(xi, realmean, realstd);
                plot(xi, pdf, '--c');
                
                % kde
                [f,xi] = ksdensity(errOverTimeFitPosNeg);                
                plot(xi,f);
                
                hold off
                title (sprintf('mean: %f std:%f rms:%f', realmean, realstd, rms))
                
                % lenear fit
                vecOverTimeInterp = interp1q(t0, vecOverTime, tt);
                
                % poly fit
                p = polyfit(t0, vecOverTimeFit, 2);
                vecOverTimeInterpFit = polyval(p, tt);
                
                % error between poly fit and linear fit
                errOverTimeFit = abs(vecOverTimeInterp - vecOverTimeInterpFit);
                errOverTimeFitPosNeg = vecOverTimeInterp - vecOverTimeInterpFit;
                % goodness of fit
                n = length(errOverTimeFit);
                count1rms = sum(errOverTimeFit<=rms) / n
                count2rms = sum(errOverTimeFit<=rms*2) / n
                count3rms = sum(errOverTimeFit<=rms*3) / n
                
                
                % rmse 
                subplot(3,3,7);
                vecOverTimeLinear2pt = interp1q([t0(1); t0(end)], [vecOverTime(1); vecOverTime(end)], (t0(1):t0(end))');
                err = vecOverTime - vecOverTimeLinear2pt;
                rmse1 = sqrt(mean(err.*err));
                disp(sprintf('Order %d, err mean=%f, mean square=%f', 1, mean(err), rmse1))
                rmse_list(count, 1) = 1;
                err = vecOverTime - vecOverTimeFit;
                rmse = sqrt(mean(err.*err));
                disp(sprintf('Order %d, err mean=%f, mean square=%f', 2, mean(err), rmse))
                rmse_list(count, 2) = rmse/rmse1;
                
                
                % 2nd, 3rd, 4th, 5th order
                for i=1:4
                    [p errEst] = polyfit(t0, vecOverTime, i);
                    [f delta] = polyval(p,t0, errEst);
                    plot(t0,vecOverTime,'.',...
                        t0,f,'g-',...
                        t0,f+2*delta,'r:',...
                        t0,f-2*delta,'r:');
                    err = vecOverTime - f;
                    mu = mean(err);
                    rmse = sqrt(mean(err.*err));
                    disp(sprintf('Order %d, err mean=%f, mean square=%f', i, mu, rmse))
                    if i>2
                        rmse_list(count, i)=rmse/rmse1;
                    end
                end               
                count = count +1;
                
                % error bound
                subplot(3,3,8)
                rms = sqrt(mean(errOverTimeFitPosNeg.^2));
                max_err = max(abs(errOverTimeFitPosNeg));
                [f,xi] = ksdensity(errOverTimeFitPosNeg);   
                
                % Chebyshev Bounds - Pr[|X-E[X]| >= lumbda] <= var(X)/lumbda^2
                s = std(errOverTimeFitPosNeg);
                m = mean(errOverTimeFitPosNeg);
                lumbda = abs(xi-m);
                Chebyshev = s*s./lumbda./lumbda;
                Chebyshev(find(Chebyshev>1))=1;
                
                % Chernoff Bound
                % Pr[X>(1+d)mu] <= exp(-d^2 mu/3)
                % Pr[X<(1-d)mu] <= exp(-d^2 mu/2)
                posmean = mean(errOverTimeFit);
                idx = find(abs(xi)>=posmean);
                d = abs(xi(idx))/posmean-1;
                Chernoff = xi;
                Chernoff(idx) = exp(-d.*d.*m/3);
                idx = find(abs(xi)<posmean);
                d = 1-abs(xi(idx))/posmean;
                Chernoff(idx) = exp(-d.*d.*m/2);
                
                % rms
                Gaussian = normpdf(xi, 0, rms);
                                
                plot(xi,f, '-' , xi, Chebyshev, '.', xi, Chernoff, '.', xi, Gaussian, '.');             
                legend('kde', 'Chebyshev', 'Chernoff', 'Gaussian')
                title(sprintf('rms=%g, 95%%: %g', rms, rms*2));
                
                pause
            end
         end
     end
end

 
%                 % exponential
%                 xi = 0:max_err/20:max_err;
%                 bincount = histc(errOverTimeFit, xi);
%                 f = fit(xi',bincount,'exp1')
%                 
%                 subplot(2,3,5)
%                 plot(f,xi,bincount)
%                 title('Exponential Distribution')
                
%                 
%                 % linear interp
%                 d = (max_err-min_err)/50;
%                 xx=min_err-d: d :max_err+d;
%                 yy=zeros(1,length(xx));
%                 for t=1:len-1
%                     small = min(errOverTimeFit(t+1), errOverTimeFit(t));
%                     large = max(errOverTimeFit(t+1), errOverTimeFit(t));
%                     m = (len-1)*(large-small);
%                     range = xx>=small & xx<=large;
%                     yy(find(range)) = yy(find(range)) + 1/m;
%                     
%                 end
%                 
%                 subplot(2,3,6)
%                 plot(xx,yy)
%                 title('linearly interpolated')
                
                
                %
                
%                 xx=1:len/100:len;
%                 yy=spline(1:len, errOverTimeFit, xx)
%                 subplot(2,3,5)
%                 plot(xx,yy,'g');
%                 
%                 %
%                 subplot(2,3,6)
%                 hist(yy,20)
%                 title('histogram of b-spline')




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



