
clear all 

% settings
n_trial = 100;
mu = 0;
sigma = 1;
gridsize = 10;
range_mu = [-3, 3];
range_sigma = [0.01, 3];

%
p_lapse = 0.01;

% start trials
start_trial = 20;

%
figure('Color','w');

%
grp = 20;
posterior = NaN(grp,grp);
xc_p = linspace(range_sigma(1),range_sigma(2),grp);
yc_p = linspace(range_mu(1),range_mu(2),grp);

% 
x = NaN(n_trial,1);
r = NaN(n_trial,1);
time_taken = NaN(n_trial,1);
x(1) = mean(range_mu) + 0.01;

for t = 1:n_trial
    
    if t<= start_trial
        tic;
        x(t) = -2 + (4).*rand;
        lambda_hat = 0.05;
        time_taken(t) = toc;
    else
        tic;
        [mu_hat, sigma_hat, lambda_hat] = fit_p_r(x(1:t-1), r(1:t-1));
        
        % x(t) = compute_sweetpoint(mu_hat, sigma_hat, max(lambda_hat, 0.01),mod(t,2));
        
        % x(t) = compute_sweetpoint(mu_hat, sigma_hat, lambda_hat);
        
        if t <= start_trial*2
            x(t) = compute_sweetpoint(mu_hat, sigma_hat, 0.05); % ,mod(t,2)
        else
            x(t) = compute_sweetpoint(mu_hat, sigma_hat, lambda_hat);
        end
        time_taken(t) = toc;
    end
        
        
    % generate response
    r(t) = x(t)+sigma*randn >= mu;
    
    % add lapses
    if binornd(1,p_lapse)
        r(t) = abs(r(t)-1);
    end
    
    % likelihood for sigma & mean
    posterior = NaN(grp,grp);
    for x_i = 1:grp
        for y_i = 1:grp
             posterior(y_i,x_i) = exp(L_r(x(1:t), r(1:t), yc_p(y_i), xc_p(x_i), lambda_hat));
        end
    end
    
    % plot
    subplot(1,2,1)
    contourf(xc_p, yc_p, posterior);
    xlabel('sigma'); ylabel('mu')
    line([sigma,sigma],range_mu,'Color','k')
    line(range_sigma,[mu mu],'Color','k')
    title('Likelihood')
    
    subplot(1,2,2)
    plot(x,'o','Color','k'); hold on
    t_i = 1:t;
    plot(t_i(r(1:t)==1), x(r(1:t)==1),'o','Color','k', 'MarkerFaceColor', 'k'); hold off
    ylim([-3.5 3.5])
    xlim([-2 n_trial+3])
    xlabel('trial')
    ylabel('stimulus')
    title('Stimuli placements & responses')
    drawnow
    
end

% fpo=fopen('sim.log','w');
% for n=1:n_trial
%     fprintf(fpo,'%.4f\t%.4f\n',x(n),r(n));
% end
% fclose(fpo);
% system('/usr/local/bin/R CMD BATCH fit_2asym.R')
% fileID = fopen('par_hat.log','r');
% par_hat = fscanf(fileID,'%f');
% fclose(fileID);
% mu_hat = par_hat(1,1);
% sigma_hat = par_hat(2,1);

% final estimate
[mu_hat, sigma_hat, lambda_hat, L] = fit_p_r(x, r);

sprintf('estimates: mu=%.2f  sigma=%.2f lambda=%.2f \ncomputation time (mean) =%.2f sec.',[mu_hat, sigma_hat, lambda_hat, mean(time_taken)])

