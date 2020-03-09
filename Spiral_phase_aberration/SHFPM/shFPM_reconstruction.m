function shFPM_reconstruction()

%% initializing enviroment
[lambda,n_LED,CTF_object,~,N_lens,~,~,...
                     pix_CCD,sub_pix,plane_wave,df]=ini_enviroment();

LowSeq_data=load('resource//Lseq.mat');
I_camera = LowSeq_data.I_camera;

PIX = 1540;
S = imresize(sqrt(I_camera(sub_pix*(N_lens-1)/2+1:sub_pix*(N_lens+1)/2,...
                           sub_pix*(N_lens-1)/2+1:sub_pix*(N_lens+1)/2,1)),[PIX,PIX]);

objectFT = fftshift(fft2(S)); % initial guess of object spectrum
pupil=ones(pix_CCD);          % initial guess of pupil function

for loop = 1:60
    figure(112);
    loop
    imshow(angle(CTF_object.*pupil),[])
    for con = 1:n_LED^2
        fxc = round((PIX+1)/2+plane_wave(1,con)/lambda/df);
        fyc = round((PIX+1)/2+plane_wave(2,con)/lambda/df);
        fxl=round(fxc-(pix_CCD-1)/2);fxh=round(fxc+(pix_CCD-1)/2);
        fyl=round(fyc-(pix_CCD-1)/2);fyh=round(fyc+(pix_CCD-1)/2);
        CTF_system = CTF_object .* pupil;
        F_sub0 = objectFT(fyl:fyh,fxl:fxh);

        %--------------sub-image constrain begin---------------------------%   
        for x_con = 1:N_lens
            for y_con = 1:N_lens
                sub_image = I_camera(sub_pix*(y_con-1)+1:sub_pix*y_con,sub_pix*(x_con-1)+1:sub_pix*x_con,con); % get the m-th subimage in n-th large image
                if sum(sum(sub_image))~=0
                    sub_pupil    = CTF_system (sub_pix*(y_con-1)+1:sub_pix*y_con,sub_pix*(x_con-1)+1:sub_pix*x_con);
                    sub_spectrum = F_sub0     (sub_pix*(y_con-1)+1:sub_pix*y_con,sub_pix*(x_con-1)+1:sub_pix*x_con);
                    lowResFT_old = sub_spectrum.*sub_pupil;
                    lowResIM_old = ifft2(ifftshift(lowResFT_old));  
                    icf = sum(sum(abs(lowResIM_old).^2))/sum(sum(sub_image)); %intensity correct factor
                    lowResFT_new = fftshift(fft2(sqrt(icf*sub_image).*exp(1i*angle(lowResIM_old))));
     
                    sub_spectrum_old  = sub_spectrum;
                    sub_spectrum = sub_spectrum_old + ...
                            abs(sub_pupil).*conj(sub_pupil)./(max(max(abs(sub_pupil)))*(abs(sub_pupil).^2+eps)).*(lowResFT_new-lowResFT_old);

                    sub_pupil = sub_pupil +...
                            abs(sub_spectrum_old).*conj(sub_spectrum_old)./(max(max(abs(sub_spectrum_old)))*(abs(sub_spectrum_old).^2+eps)).*(lowResFT_new-lowResFT_old);
                               
                    F_sub0  (sub_pix*(y_con-1)+1:sub_pix*y_con,sub_pix*(x_con-1)+1:sub_pix*x_con)     = sub_spectrum;
                    pupil   (sub_pix*(y_con-1)+1:sub_pix*y_con,sub_pix*(x_con-1)+1:sub_pix*x_con)     = sub_pupil;
                    I_camera(sub_pix*(y_con-1)+1:sub_pix*y_con,sub_pix*(x_con-1)+1:sub_pix*x_con,con) = sub_image * icf;
                end
            end
        end
        %--------------sub-image constrain finished------------------------%  
        

        %--------------large image constrain begin-------------------------
        objectFT(fyl:fyh,fxl:fxh) = F_sub0;
        CTF_system = CTF_object.*pupil;
        lowResFT1 = F_sub0.*CTF_system;  
        
        lowResIM = arrayfft2(lowResFT1,N_lens,N_lens,sub_pix,'ifft');         
        cc = sum(sum(abs(lowResIM).^2))/sum(sum(I_camera(:,:,con)));       %intensity correct factor
        lowResIM = sqrt(cc*I_camera(:,:,con)).*exp(1i.*angle(lowResIM));
        lowResFT2 = arrayfft2(lowResIM,N_lens,N_lens,sub_pix,'fft'); 

        
        objectFT(fyl:fyh,fxl:fxh)=objectFT(fyl:fyh,fxl:fxh)...
               + conj(CTF_system)./(max(max(abs(CTF_system)))*(abs(CTF_system)+eps)).*(lowResFT2-lowResFT1);
           
        pupil = pupil ...
               + conj(F_sub0)./(max(max(abs(F_sub0)))*(abs(F_sub0)+eps)).*(lowResFT2-lowResFT1);  

        I_camera(:,:,con) = cc*I_camera(:,:,con);
        %----------large image constrain finished--------------------------
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

