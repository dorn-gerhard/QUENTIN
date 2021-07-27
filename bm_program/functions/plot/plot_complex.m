function im = plot_complex(Values, counter, colorbar1, colorbar2)



N = 128;
colorbar1 = plot_colorbar(N ,[],[],1);
colorbar2 = plot_colorbar2(N ,[],[],1);

x_max = max(abs(real(Values(:))));
y_max = max(abs(imag(Values(:))));



%scale real and imaginary part between 0 and N
D_x = ceil((real(Values) +10^-10+ x_max)/(2*(x_max+10^-10)) * N);
D_y = ceil((imag(Values) +10^-10+ y_max)/(2*(y_max+10^-10)) * N);

%[min(D_x(:)), max(D_x(:)), min(D_y(:)), max(D_y(:))]


%image = zeros(size(Values,1), size(Values,2), 3);

im = zeros(size(Values,1), size(Values,2),3);

for k = 1:3
    im(:,:,k) = reshape(1-colorbar1(D_x,k),size(Values,1),size(Values,2))/2 + reshape(1-colorbar2(D_y,k),size(Values,1),size(Values,2))/2;
end

xy_max = max(x_max, y_max);
 

%scale real and imaginary part between 0 and N
D_x = ceil((real(Values) +10^-10+ xy_max)/(2*(xy_max+10^-10)) * N);
D_y = ceil((imag(Values) +10^-10+ xy_max)/(2*(xy_max+10^-10)) * N);

%[min(D_x(:)), max(D_x(:)), min(D_y(:)), max(D_y(:))]


%image = zeros(size(Values,1), size(Values,2), 3);

im_equal = zeros(size(Values,1), size(Values,2),3);

for k = 1:3
    im_equal(:,:,k) = reshape(1-colorbar1(D_x,k),size(Values,1),size(Values,2))/2 + reshape(1-colorbar2(D_y,k),size(Values,1),size(Values,2))/2;
end


ht = guidata(gcf);

ht.C_Data_equal{counter} = 1-im_equal;
ht.C_Data{counter} = 1-im;
ht.x_max = x_max;
ht.y_max = y_max;
guidata(gcf,ht);


% plot legend
subplot(1,5,5)
x = linspace(-1,1,64);
[xx,yy] = meshgrid(x);
Values2 = xx + 1i*(yy);
x_max2 = max(abs(real(Values2(:))));
y_max2 = max(abs(imag(Values2(:))));

%scale real and imaginary part between 0 and N
D_x2 = ceil((real(Values2) +10^-10+ x_max2)/(2*(x_max2+10^-10)) * N);
D_y2 = ceil((imag(Values2) +10^-10+ y_max2)/(2*(y_max2+10^-10)) * N);

im2 = zeros(size(Values2,1), size(Values2,2),3);

for k = 1:3
    im2(:,:,k) = reshape(1-colorbar1(D_x2,k),size(Values2,1),size(Values2,2))/2 + reshape(1-colorbar2(D_y2,k),size(Values2,1),size(Values2,2))/2;
end
image(1-im2)
axis square
xlabel('real part')
ylabel('imag part')
set(gca,'xtick', [0.5, 32,64.5])
set(gca,'ydir', 'normal')
set(gca,'xticklabel', {num2str(-x_max,3), '0', num2str(x_max,3)})

set(gca,'ytick', [0.5, 32,64.5])
set(gca,'yticklabel', {num2str(-y_max,3), '0', num2str(y_max,3)})
set(gcf, 'color', 'w')
    
%plot image
subplot(1,5,[1,4])
image(1-im)
    