function signal = denoise(signal, tolerance)

signal(abs(imag(signal)) < tolerance) = real(signal(abs(imag(signal)) < tolerance));
%denoise real part:
signal(abs(real(signal)) < tolerance) = imag(signal(abs(real(signal)) < tolerance));
%denoise 
signal(abs(signal) < tolerance) = 0;
