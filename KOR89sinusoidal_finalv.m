% Annals of Biomedical Engineering, Vol. 17, pp. 219-231, 1989 Applications
% of Fast Orthogonal Search
%% initialization
clc; clear all; close all;

%% data to be simulated

n = linspace(0,150,151);
omega   = [0.60  0.10  0.45 0.85  0.08];
b       = [0.44  0.95  0.61 0.00 -0.59];
c       = [1.55 -0.60 -0.88 0.95  0.10];
y       = zeros(1,length(n));

for i = 1:5
    y = b(i)*cos(omega(i)*n) + c(i)*sin(omega(i)*n)+y;
end
figure;hold on;
for i = 1
    plot(b(i)*cos(omega(i)*n) + c(i)*sin(omega(i)*n),'-c','linewidth',2);
end
for i = 2
    plot(b(i)*cos(omega(i)*n) + c(i)*sin(omega(i)*n),'-k','linewidth',1.5);
end
for i = 3
    plot(b(i)*cos(omega(i)*n) + c(i)*sin(omega(i)*n),'--b','linewidth',2);
end
for i = 4
    plot(b(i)*cos(omega(i)*n) + c(i)*sin(omega(i)*n),'-.r','linewidth',2);
end
for i = 5
    plot(b(i)*cos(omega(i)*n) + c(i)*sin(omega(i)*n),'.-m','linewidth',2);
end

figure;
plot(n,y,'o','LineWidth',1);

%% searching
omega_d = linspace(0.01,1.00,100);

p = NaN(2*length(omega_d)+1,length(n)); % +1 is for the DC term
p(1,:) = ones(1,length(n)); % DC term is set to 1

w = NaN(2*length(omega_d)+1,length(n)); % orthogonalized basis
w(1,:) = ones(1,length(n));

for i = 1:length(omega_d)
    p((2*i),:) = cos(omega_d(i)*n);
    p((2*i+1),:) = sin(omega_d(i)*n);
end

%%
alpha = NaN(2*5+1,2*5);

g = NaN(1,2*5+1);
g(1) = sum(y.*w(1,:))/sum(w(1,:).^2);

index = NaN(1,5);

%%

for m = 1:5
    temp = zeros(1,length(omega_d));
    for i = 1:length(omega_d)
        
        if (~isnan(omega_d(i))) % this is to skip the selected frequency.
        
        w(2*m,:) = p((2*i),:);
        
        for r = 1:(2*m-1)
        alpha(2*m,r) = mean(p(2*i,:).*w(r,:))/mean(w(r,:).^2);
        w(2*m,:) = w(2*m,:) - alpha(2*m,r)*w(r,:);
        end
        
        g(2*m) = mean(y.*w(2*m,:))/mean(w(2*m,:).^2);
        
        w((2*m+1),:) = p((2*i+1),:);
        
        for r = 1:(2*m)
        alpha(2*m+1,r) = sum(p(2*i+1,:).*w(r,:))/sum(w(r,:).^2);
        w(2*m+1,:) = w(2*m+1,:) - alpha(2*m+1,r)*w(r,:);
        end
        
        g(2*m+1) = sum(y.*w(2*m+1,:))/sum(w(2*m+1,:).^2);
        
        temp(i) = g(2*m)^2*mean(w(2*m,:).^2) + g(2*m+1)^2*mean(w(2*m+1,:).^2);
        else
        end
    end
    
    [C, I] = max(temp);
    
    index(1,m) = I;
    omega_d(I) = NaN;
    
    for i = I:I
        
        w(2*m,:) = p((2*i),:);
        
        for r = 1:(2*m-1)
        alpha(2*m,r) = mean(p(2*i,:).*w(r,:))/mean(w(r,:).^2);
        w(2*m,:) = w(2*m,:) - alpha(2*m,r)*w(r,:);
        end
        
        g(2*m) = mean(y.*w(2*m,:))/mean(w(2*m,:).^2);
        
        w((2*m+1),:) = p((2*i+1),:);
        
        for r = 1:(2*m)
        alpha(2*m+1,r) = sum(p(2*i+1,:).*w(r,:))/sum(w(r,:).^2);
        w(2*m+1,:) = w(2*m+1,:) - alpha(2*m+1,r)*w(r,:);
        end
        
        g(2*m+1) = sum(y.*w(2*m+1,:))/sum(w(2*m+1,:).^2);
        
    end
end

a = zeros(1,length(g));

a(11) = g(11);

store_v = zeros(11,11);
store_v (11,11) = 1;
% according to Equations (15),(16),(17), the coefficients a are calculated
% backwards from coefficients g.
% store_v is for v functions in Equation (15)

for m = (length(g)-1):-1:1
    v = zeros(1,length(g));
        v(m) = 1;
        a(m) = g(m)*v(m);
    for i = (m+1):length(g)
        for r = m:(i-1);
            v(i) = v(i) - alpha(i,r)*v(r);
        end
        a(m) = a(m)+ g(i)*v(i);
    end
    store_v(m,:) = v;
end

%% plot retrieved data to check
hold on;
yyy = 0;

for i = 1:(2*5+1)
yyy = yyy+g(i)*w(i,:);
end
plot(n,yyy,'-m','linewidth',1);
xlabel('Samples');
ylabel('Amplitude');
