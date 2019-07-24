function f = gaussian(x0,mag,var,lambda)
f = mag * exp( -(lambda-x0).^2/(2*var^2) );
end
