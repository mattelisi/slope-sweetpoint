%
% This script performa a simple test of the procedure
% and illustrate how it can be used.
% It also plots the likelihood of different values of mu and sigma, and
% show how it evolves during testing.
%

clear all 

% settings
n_trial = 100;
plot_l = 1;

% true parameters
mu = 0;
sigma = 1;
p_lapse = 0.01; % AKA lambda

% start trials:
% this indicates for how man trials at the beginning the stimulus should be 
% selected at random before switching to the sweets points.
start_trial = 20;

% this is just for plotting
if plot_l
    grp = 20;
    likelihood2 = NaN(grp,grp);
    range_mu = [-3, 3];
    range_sigma = [0.01, 3];
    xc_p = linspace(range_sigma(1),range_sigma(2),grp);
    yc_p = linspace(range_mu(1),range_mu(2),grp);
    figure('Color','w', 'Position',[0 50 500 220]); %[left bottom width height]
end

% 
x = NaN(n_trial,1);
r = NaN(n_trial,1);
time_taken = NaN(n_trial,1);

for t = 1:n_trial
    
    if t<= start_trial
        tic;
        x(t) = -2 + (4).*rand;
        lambda_hat = 0.05;
        time_taken(t) = toc;
    else
        tic;
        [mu_hat, sigma_hat, lambda_hat] = fit_p_r(x(1:t-1), r(1:t-1));
        
        
        %%% Simple: %%%
        %x(t) = compute_sweetpoint(mu_hat, sigma_hat, lambda_hat);
        
        %%% This assumes a given minimum probability of lapses at the beginning %%%
        % of the session (it can be useful becasue it shifts the sweetpoints 
        % toward the threshold): 
        if t <= start_trial*2
            x(t) = compute_sweetpoint(mu_hat, sigma_hat, max(lambda_hat, 0.05)); % ,mod(t,2)
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
    
    % plot
    if plot_l
        
        % compute likelihood for sigma & mean (only for plotting)
        likelihood2 = NaN(grp,grp);
        for x_i = 1:grp
            for y_i = 1:grp
                likelihood2(y_i,x_i) = exp(L_r(x(1:t), r(1:t), yc_p(y_i), xc_p(x_i), lambda_hat));
            end
        end
        
        subplot(1,2,1)
        likelihood2 = likelihood2/sum(likelihood2(:));
        contourf(xc_p, yc_p, likelihood2);
        xlabel('sigma'); ylabel('mu')
        line([sigma,sigma],range_mu,'Color','k')
        line(range_sigma,[mu mu],'Color','k')
        title('Likelihood')
        
        subplot(1,2,2)
        plot(x,'o','Color','k'); hold on
        t_i = 1:t;
        % use filled dots for "+" responses
        plot(t_i(r(1:t)==1), x(r(1:t)==1),'o','Color','k', 'MarkerFaceColor', 'k'); hold off
        ylim([-3.5 3.5])
        xlim([-2 n_trial+3])
        xlabel('trial')
        ylabel('stimulus')
        title('Stimuli placements & responses')
        drawnow;
    end
    
end

% gives the final estimate
[mu_hat, sigma_hat, lambda_hat, L] = fit_p_r(x, r);

fprintf('\n\nEstimates: mu=%.2f  sigma=%.2f lambda=%.2f \ncomputation time (mean) =%.2f sec.',[mu_hat, sigma_hat, lambda_hat, mean(time_taken)]);
fprintf('\nThe true values were: mu=%.2f  sigma=%.2f lambda=%.2f.\n',[mu, sigma, p_lapse]);

