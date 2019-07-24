function [g, alpha, index, order] = KOR89(y,n,omega_d)
% This version does NOT include sinusoid 
%% data to be simulated

% n = linspace(0,150,151);
% omega   = [0.60  0.10  0.45 0.85  0.08];
% b       = [0.44  0.95  0.61 0.00 -0.59];
% c       = [1.55 -0.60 -0.88 0.95  0.10];
% y       = zeros(1,length(n));
% 
% for i = 1:5
%     y = b(i)*cos(omega(i)*n) + c(i)*sin(omega(i)*n)+y;
% end

figure;
plot(y,'.-b','linewidth',0.5);


prompt = 'What is the order of simulation?\n order = ';
order = input(prompt);
% order = 25;
tic;
%% searching
% omega_d = sigma;
% omega_d = linspace(0.01,1.00,100);

p = NaN(length(omega_d)+1,length(n));
p(1,:) = ones(1,length(n));

w = NaN(length(omega_d)+1,length(n));
w(1,:) = ones(1,length(n));

for i = 1:length(omega_d)
    p(i+1,:) = cos(omega_d(i)*n);
end

%%
alpha = NaN(order+1,order);

g = NaN(1,order+1);
g(1) = sum(y.*w(1,:))/sum(w(1,:).^2);

index = NaN(1,order);

%%

for m = 1:order
    temp = zeros(1,length(omega_d));
    for i = 1:length(omega_d)
        
        if (~isnan(omega_d(i))) % this is to skip the selected frequency.
        
        w(m+1,:) = p(i+1,:);
        
        for r = 1:m
        alpha(m+1,r) = mean(p(i+1,:).*w(r,:))/mean(w(r,:).^2);
        w(m+1,:) = w(m+1,:) - alpha(m+1,r)*w(r,:);
        end
        
        g(m+1) = mean(y.*w(m+1,:))/mean(w(m+1,:).^2);
%         if (g(m+1)<=0)
%             g(m+1) = 0;
%         end
        
        temp(i) = g(m+1)^2*mean(w(m+1,:).^2);
        else
        end
    end
    
    
    
    [~, I] = max(temp);
    omega_d(I) = NaN;
    index(1,m) = I;
    
    
    for i = I:I
        
        w(m+1,:) = p(i+1,:);
        
        for r = 1:(m)
        alpha(m+1,r) = mean(p(i+1,:).*w(r,:))/mean(w(r,:).^2);
        w(m+1,:) = w(m+1,:) - alpha(m+1,r)*w(r,:);
        end
        
        g(m+1) = mean(y.*w(m+1,:))/mean(w(m+1,:).^2);
%         if (g(m+1)>=1)
%             g(m+1) = 1;
%         end
    end
end



%% plot retrieved data to check
% hold on;
% yyy = 0;
% 
% for i = 1:(order+1)
% yyy = yyy+g(i)*w(i,:);
% end
% plot(n,yyy,'.-r');



end