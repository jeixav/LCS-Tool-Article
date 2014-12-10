function derivative_ = derivative(t,x,$\epsilon$,A,$\omega$)
...
a = $\epsilon$*sin($\omega$*t);
b = 1 - 2*$\epsilon$*sin($\omega$*t);
f = a*x(idx1).^2 + b*x(idx1);

derivative_(idx1) = -$\pi$*A*sin($\pi$*f).*cos($\pi$*x(idx2));
derivative_(idx2) = $\pi$*A*cos($\pi$*f).*sin($\pi$*x(idx2)).*(2*a*x(idx1) + b);
