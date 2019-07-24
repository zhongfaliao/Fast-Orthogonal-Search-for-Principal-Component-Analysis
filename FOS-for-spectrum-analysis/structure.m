function [A_1, A_2, B,L1,delta_Li,L2] = structure(number,gamma_s,kappa_s,gamma_c,kappa_c,L_upper,delta_L,alpha)

% Redundant
% S_s = gamma_s * [sqrt(1-kappa_s) -1i*sqrt(kappa_s);-1i*sqrt(kappa_s) sqrt(1-kappa_s)];
% S_c = gamma_c * [sqrt(1-kappa_c) -1i*sqrt(kappa_c);-1i*sqrt(kappa_c) sqrt(1-kappa_c)];

% upper arm length kept constant.
L1      = ones(1,number)*L_upper;

% unbalance between two arms of i-th MZI
% delta_Li = (0:1:(number-1))*delta_L;
delta_Li = (1:1:number)*delta_L;
% lower arm length gradually increase
L2 = L1 + delta_Li;

%% Find A A B coefficients

gamma_di = exp(-alpha*L2);

% redundant
% Sdi = zeros(2,2,number);
% for i=1:number
%     Sdi(:,:,i) = gamma_di(i) * exp(-1i*beta*L2(i)) * [exp(-alpha*delta_Li(i))*exp(-1i*beta*delta_Li(i)) 0; 0 1];
% end
%%
% Equation (8) :
% 
% $$A_{1,i}=2\gamma_s^2\gamma_{d,i}^2\gamma_c^2[\kappa_s\kappa_c+(1-\kappa_s)(1-\kappa_c)e^{-2\alpha\Delta L_i}]$$
% 
%%
% Equation (9) :
% 
% $$A_{2,i}=2\gamma_s^2\gamma_{d,i}^2\gamma_c^2[\kappa_s(1-\kappa_c)+\kappa_c(1-\kappa_s)e^{-2\alpha\Delta L_i}]$$
% 
%%
% Equation (10) :
% 
% $$B_i=4\gamma_s^2\gamma_{d,i}^2\gamma_c^2[\kappa_s(1-\kappa_c)\kappa_c(1-\kappa_s)]^{\frac{1}{2}}e^{-\alpha\Delta L_i}$$
% 
A_1 = 2*gamma_s^2.*gamma_di.^2*gamma_c.* ( kappa_s*kappa_c + (1-kappa_s)*(1-kappa_c).*exp(-2*alpha*delta_Li));
A_2 = 2*gamma_s^2.*gamma_di.^2*gamma_c.* ( kappa_s*(1-kappa_c) + (1-kappa_s)*kappa_c.*exp(-2*alpha*delta_Li));
B   = 4*gamma_s^2.*gamma_di.^2*gamma_c.* sqrt( kappa_s*kappa_c*(1-kappa_s)*(1-kappa_c) ).*exp(-alpha*delta_Li);
% above three coefficients depend on alpha, which is a constant, and
% delta_Li, which increase from i=1 to i=number. In addition,
% they are wavelength-independent.
end