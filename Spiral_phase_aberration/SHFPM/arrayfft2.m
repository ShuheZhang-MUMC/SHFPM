function out = arrayfft2(in,N_lensx,N_lensy,sub_pix,type)
% array fft2, lens array

out=zeros(size(in));
if strcmp(type,'fft')
    for x_con = 1:N_lensx
        for y_con = 1:N_lensy
            temp_in = in(sub_pix*(y_con-1)+1:sub_pix*y_con,sub_pix*(x_con-1)+1:sub_pix*x_con);       
            out(sub_pix*(y_con-1)+1:sub_pix*y_con,...
                sub_pix*(x_con-1)+1:sub_pix*x_con)= fftshift(fft2(temp_in));
        end
    end
elseif strcmp(type,'ifft')
    for x_con = 1:N_lensx
        for y_con = 1:N_lensy
            temp_in = in(sub_pix*(y_con-1)+1:sub_pix*y_con,sub_pix*(x_con-1)+1:sub_pix*x_con);       
            out(sub_pix*(y_con-1)+1:sub_pix*y_con,...
                sub_pix*(x_con-1)+1:sub_pix*x_con)= ifft2(ifftshift(temp_in));
        end
    end
end
    

end