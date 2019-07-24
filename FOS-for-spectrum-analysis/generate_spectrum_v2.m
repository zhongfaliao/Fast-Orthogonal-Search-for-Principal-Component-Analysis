
function spectro = generate_spectrum_v2(d_sigma,resol)
spectro = zeros(1,length(d_sigma));
% spectro = spectro + gaussian(d_sigma(155),0.4,45*resol,d_sigma);
% spectro = spectro + gaussian(d_sigma(155),0.02,15*resol,d_sigma);
spectro(240) = 1;
spectro(230) = 0.5;
spectro(235) = 1.5;
spectro(260) = 1.5;
spectro(210) = 0.2;

% spectro = spectro + gaussian(d_sigma(235),0.95,33*resol,d_sigma).^2;
% spectro = spectro + gaussian(d_sigma(320),0.97,40*resol,d_sigma).^2;
end