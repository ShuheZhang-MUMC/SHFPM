function EPRY_reconstruction()
% close all
%% initializing enviroment
[lambda,n_LED,CTF_object,~,...
           pix_CCD,plane_wave,df]=ini_enviroment();

load('Lseq.mat');

PIX = 1540;
S = imresize(sqrt(I_camera(:,:,1)),[PIX,PIX]);
%imshow(S,[])
objectFT = fftshift(fft2(S));
pupil=1;

xx=linspace(-1,1,PIX);
[x,y]=meshgrid(xx);
mask = (abs(x)<0.98).*(abs(y)<0.98);


for loop = 1:20
    loop
    figure(112);
    imshow(abs(ifft2(ifftshift(objectFT))).^2.*mask,[])
    for con = 1:n_LED^2
        fxc = round((PIX+1)/2+plane_wave(1,con)/lambda/df);
        fyc = round((PIX+1)/2+plane_wave(2,con)/lambda/df);
        fxl=round(fxc-(pix_CCD-1)/2);fxh=round(fxc+(pix_CCD-1)/2);
        fyl=round(fyc-(pix_CCD-1)/2);fyh=round(fyc+(pix_CCD-1)/2);
        
        lowResFT1 = objectFT(fyl:fyh,fxl:fxh).*CTF_object.*pupil;
        
        lowResIM = ifft2(ifftshift(lowResFT1));
        
        cc = sum(sum(abs(lowResIM).^2))/sum(sum(I_camera(:,:,con))); % intensity correction factor
        lowResIM = sqrt(cc*I_camera(:,:,con)).*exp(1i.*angle(lowResIM));
        lowResFT2 = fftshift(fft2(lowResIM)); 

        temp_objectFT = objectFT(fyl:fyh,fxl:fxh);
        CTF_system = CTF_object.*pupil;
        
%%       PIE-type factor   
        objectFT(fyl:fyh,fxl:fxh)=objectFT(fyl:fyh,fxl:fxh)...
               + abs(CTF_system).*conj(CTF_system)./(max(max(abs(CTF_system))).*(abs(CTF_system).^2+eps)).*(lowResFT2-lowResFT1);
           
        pupil = pupil ...
               + abs(temp_objectFT).*conj(temp_objectFT)./(max(max(abs(temp_objectFT))).*(abs(temp_objectFT).^2+eps)).*(lowResFT2-lowResFT1);      

%%       ePIE-type factor   
%         objectFT(fyl:fyh,fxl:fxh)=objectFT(fyl:fyh,fxl:fxh)...
%                + conj(CTF_system)./max(max(abs(CTF_system).^2)).*(lowResFT2-lowResFT1);
%            
%         pupil = pupil ...
%                + conj(temp_objectFT)./max(max(abs(temp_objectFT).^2)).*(lowResFT2-lowResFT1);      

    I_camera(:,:,con) = I_camera(:,:,con)*cc;
    end
end

[m_x,m_y]= meshgrid(linspace(-1,1,PIX));
mask = (abs(m_x)<0.98) .* (abs(m_y)<0.98);

str = date;
I = abs(ifft2(ifftshift(objectFT))).^2 .* mask;
I = I - min(I(:));
I = I / max(I(:));
imwrite(I,['output//I_',str,'.png'])

P1 = angle(ifft2(ifftshift(objectFT))) .* mask;
P1 = P1 - min(P1(:));
P1 = P1 / max(P1(:));
imwrite(P1,['output//P_',str,'.png'])


P = mod(angle(CTF_object.*pupil),2*pi);
P = P - min(P(:));
P = P / max(P(:));
imwrite(P,['output//Pupil_',str,'.png'])


figure();
subplot(221);imshow(I,[]);title('recovered intensity')
subplot(222);imshow(P1,[]);title('recovered phase')
subplot(223);imshow(abs(CTF_object.*pupil),[]);title('recovered pupil modulus')
subplot(224);imshow(mod(angle(CTF_object.*pupil),2*pi),[]);title('recovered pupil phase')

end