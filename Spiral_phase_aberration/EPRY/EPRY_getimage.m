function EPRY_getimage()
%{
    FPM simulation 
    generate low resolution image
%}


%% read object imaging
pix = 1540;
I = imread('resource\\intensity.png');
I = double(I(:,:,1));I = I-min(I(:));I = 1-I/max(max(I));I = imresize(I,[pix,pix]); %intensity


O = sqrt(I).*exp(1i*0); % complex amplitude of object

%% initializing enviroment
[lambda,n_LED,~,CTF_object,...
           pix_CCD,plane_wave,df]=ini_enviroment();

%% Simulation of the process of snapshot
F = fftshift(fft2(O));
I_camera = zeros(pix_CCD,pix_CCD,n_LED^2);

% figure();
for con = 1:n_LED^2
    fxc = round((pix+1)/2+(plane_wave(1,con)/lambda)/df);
    fyc = round((pix+1)/2+(plane_wave(2,con)/lambda)/df);
    
    fxl=round(fxc-(pix_CCD-1)/2);fxh=round(fxc+(pix_CCD-1)/2);
    fyl=round(fyc-(pix_CCD-1)/2);fyh=round(fyc+(pix_CCD-1)/2);
    
    F_sub = F(fyl:fyh,fxl:fxh) .* CTF_object; 
    I_camera(:,:,con)=abs(ifft2(ifftshift(F_sub))).^2;
end

I_camera = I_camera-min(I_camera(:));
I_camera = I_camera/max(I_camera(:));
figure();
imshow((I_camera(:,:,1)),[])
save Lseq I_camera


end