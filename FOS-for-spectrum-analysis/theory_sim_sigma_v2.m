function [ F, p_in, resol, sigma_bar, delta_L, neff] = theory_sim_sigma_v2
    
    clc;clear all;close all;
    
    %% Wavelength and Wavenumber
    lambda0 = 840.7e-9;
    sigma0  = floor(1/lambda0);
    delta_min =  20000;
    delta_max = 160000;

    lambda_min = 1/(1/lambda0-delta_min);
    lambda_max = 1/(1/lambda0-delta_max);

    neff = 3.49769;

    sigma_min = floor(1/lambda_max);
    sigma_max = floor(1/lambda_min);

    resol = 100;
    % units m-1, therefore 1 cm-1
    sigma = (sigma_max):-resol:(sigma_min);
    d_sigma = ((sigma0-sigma_max):resol:(sigma0-sigma_min))/100;
    % d_sigma is the Raman shift which is from 200 to 1600.
    
    sigma_bar = sigma-sigma_min; % note this is in [1/m].

    %% Generate and test spectrum
%     lambda = 1./(sigma);

%     resol_lambda = ((lambda_max+lambda_min)/2)^2*resol;

    %%
    resol=resol/100; % to generate spectrum, this is in cm-1.
    spectro = generate_spectrum_v2(d_sigma,resol);
    figure;
    plot(spectro,'-k','linewidth',1);

    %% Define the MZI array
    number = 250;
    x = 1:1:number; % index for plotting only
    gamma_s = 1.0;
    kappa_s = 0.5;
    gamma_c = 1.0;
    kappa_c = 0.5;
    alpha   = 1e-9;
    L_upper = 5000e-9;
    delta_L = 3/(neff*sigma_min);
    %%

    [A_1, A_2, B,L1,delta_Li,L2] = ...
       structure(number,gamma_s,kappa_s,gamma_c,kappa_c,L_upper,delta_L,alpha);

    %% Input and Output Spectrum
%     delta_Li = linspace(0,2.7771e-07,250);
    p_in = spectro;

    P_in = (sum(p_in)*resol); 
    %    p_in is in [W/cm-1], therefore resol must be in [cm-1].
    % or,p_in is in [#/cm-1], therefore resol must be in [cm-1].

    Pout_1 = zeros(1,number);

    for i=1:number
        Pout_1(i) = (1/2)*A_1(i)*P_in - (1/2) * B(i) * sum(resol.*p_in.*cos(2*pi*sigma_bar*neff*delta_Li(i)));
%         Pout_1(i) = (1/2)*A_1(i)*P_in - (1/2) * B(i) * sum(resol.*p_in.*cos(2*pi*sigma*neff*delta_Li(i)));
        % resol is in [cm-1], sigma is the *wavenumber*.
    end

    
%     2.04216415e-6
    %% Modified spatial interferogram.

    F = - 1./B.*(2*Pout_1 - A_1*P_in);
    figure;
    plot(F,'.-r');

end