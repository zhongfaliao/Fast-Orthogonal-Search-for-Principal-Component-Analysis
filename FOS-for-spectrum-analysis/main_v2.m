% main function
clc; clear all; close all;
% the functions theory_sim_sigma.m and KOR89.m are called

%%
[ F, p_in, resol, sigma_bar, delta_L, neff] = theory_sim_sigma_v2;
% theory_sim_sigma.m returns F(x), p_in(i), resol, sigma(i) x.  The
% corresponding variables are y(n), b(i), omega(i), and n. Candidate
% frequencies are those in omega(i).

%%
y = F;
b = p_in*resol;
n = 1:1:250;

% omega_d = linspace(min(sigma),max(sigma),281);

omega_d = 2*pi*sigma_bar*neff*delta_L;
% omega_d = 2*pi*linspace(min(sigma_bar),max(sigma_bar),401)*neff*delta_L;
omega_d = (round(omega_d/0.000001))/(1/0.000001);
%%

[g, alpha,index, order] = KOR89(y,n,omega_d);
% using the other function KOR89.m, g(i) are calculated.  the
% coefficients g(i) and alpha_cos are returned.
toc;

%% plot spectrum
a = zeros(1,length(g));

a(length(g)) = g(length(g));

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

end


% figure;
spectro = zeros(1, length(omega_d));

for i = 1:(order);
%     if (i~=20)
    spectro(index(i)) = a(i+1);
%     end
end

% plot(omega_d,spectrum,'-');

%%
figure;

plot(spectro,'.-m');
hold on;
plot(p_in/max(p_in)*max(spectro),'linewidth',2);


