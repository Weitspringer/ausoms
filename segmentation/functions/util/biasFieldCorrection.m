function [correctedStack] = biasFieldCorrection(stack)
%BIASFIELDCORRECTION 2D bias field correction on a whole stack
%   MICO, [0, 256], double output
%   Chunming Li (2021). MRI segmentation and bias field correction  
%   (https://www.mathworks.com/matlabcentral/fileexchange/
%   59752-mri-segmentation-and-bias-field-correction), 
%   MATLAB Central File Exchange. Retrieved April 15, 2021.
%   Modified: L. Springer on April 15, 2021

iterNum = 5;
N_region=3;  q=1;
correctedStack = stack;
for i = 1:size(stack, 3)
    Img = double(im2uint8(correctedStack(:,:,i)));
    A = 255;
    Img_original = Img;
    [nrow,ncol] = size(Img);n = nrow*ncol;
    %load ROI
    ROI = (Img>10); ROI = double(ROI);
    
    Bas=getBasisOrder3(nrow,ncol);
    N_bas=size(Bas,3);
    for ii=1:N_bas
        ImgG{ii} = Img.*Bas(:,:,ii).*ROI;
        for jj=ii:N_bas
            GGT{ii,jj} = Bas(:,:,ii).*Bas(:,:,jj).*ROI;
            GGT{jj,ii} = GGT{ii,jj} ;
        end
    end
    
    energy_MICO = zeros(3,iterNum);
    
    b=ones(size(Img));
    for ini_num = 1:1
        C=rand(3,1);
        C=C*A;
        M=rand(nrow,ncol,3);
        a=sum(M,3);
        for k = 1 : N_region
            M(:,:,k)=M(:,:,k)./a;
        end
        
        [e_max,N_max] = max(M,[], 3);
        for kk=1:size(M,3)
            M(:,:,kk) = (N_max == kk);
        end
        
        M_old = M; chg=10000;
        energy_MICO(ini_num,1) = get_energy(Img,b,C,M,ROI,q);
        
        for n = 2:iterNum
            
            [M, b, C]=  MICO(Img,q,ROI,M,C,b,Bas,GGT,ImgG,1, 1);
            energy_MICO(ini_num,n) = get_energy(Img,b,C,M,ROI,q);
            
            if(mod(n,1) == 0)
                PC=zeros(size(Img));
                for k = 1 : N_region
                    PC=PC+C(k)*M(:,:,k);
                end
                img_bc = Img./b;  % bias field corrected image
            end
        end
    end
    
    [M,C]=sortMemC(M,C);
    seg=zeros(size(Img));
    for k = 1 : N_region
        seg=seg+k*M(:,:,k);   % label the k-th region 
    end
    I = double(uint8(img_bc .* ROI));
    I = I ./ 255;
    correctedStack(:,:,i) = I;
end
end

