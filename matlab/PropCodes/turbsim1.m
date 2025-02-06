 function simturb=turbsim1(rootSDF1)
%        simturb=turbsim1(rootSDF1)
%   
% Generate realization of field with SDF rootSDF.^2
%
len=length(rootSDF1);


% Rodrigo: When i tried to implement xi as hermitian, the beggining of the
% simulations presented low values of scintillation. Is it really necessary
% to ensure that xi is Hermitian?
%
% half_real_part_xi_hermitian = randn(1,len/2);
% half_imag_part_xi_hermitian = randn(1,len/2);
%
% xi_hermitian = [half_real_part_xi_hermitian,flip(half_real_part_xi_hermitian)] + ...
%     1j*[half_imag_part_xi_hermitian,-flip(half_imag_part_xi_hermitian)];
%
% simturb=real(fftshift(fft(fftshift(rootSDF1.*xi_hermitian))));

xi=(randn(1,len)+1i*randn(1,len));
xi_real = randn(1,len);
simturb=real(fftshift(fft(fftshift(rootSDF1.*xi_real))));

% Rodrigo: Diagonising the power leaking to the imaginary part of the
% realization.
figure;
hold on;
plot(real(fftshift((fft(fftshift(rootSDF1.*xi))))));
plot(imag(fftshift((fft(fftshift(rootSDF1.*xi))))));
hold off;

return