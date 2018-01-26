function [Fresnel_subaper_image] = Fresnel_ModifiedSH(SH, ML, WF)

wvl = SH.LAMBDA;
delta1 = SH.PixelSize;%the sample sapce of the aperture
delta2 = SH.PixelSize/5;%the sample space of the obsersve space
k = 2*pi / SH.LAMBDA;
Dz = ML.focal-ML.defocus;
Fresnel_subaper_image=cell(1,length(ML.coor));
% Dz = ((ML.focal)*((ML.focal)+ (ML.defocus)))/(ML.defocus);
%-----------------
%caculate the fresnel propageiton of each one of the subaperture
WF_aperture = SH.pupil;
lenslet_aperture = ML.AmplitudeMask;
WF_lenslet_aperture = WF_aperture .* lenslet_aperture;
[x,y] = meshgrid((-(ML.spacing)/2:(ML.spacing)/2-1)*(SH.PixelSize));
[theta,r] = cart2pol(x,y);
phase_sub_lens = -r.^2/(2*ML.focal);
Ein_sub_aper=cell(1,length(ML.coor));%Preallocation
for lens_idx = 1 : length(ML.coor)
    WF_sub_aper{lens_idx}=WF(ML.coor(lens_idx,1):ML.coor(lens_idx,2), ML.coor(lens_idx,3):ML.coor(lens_idx,4));
    Phase_sub_aper = exp(1i * k * WF_sub_aper{lens_idx}).*exp(1i * k * phase_sub_lens);
    transmit_sub_aper = WF_lenslet_aperture(ML.coor(lens_idx,1):ML.coor(lens_idx,2), ML.coor(lens_idx,3):ML.coor(lens_idx,4));
    Ein(ML.coor(lens_idx,1):ML.coor(lens_idx,2), ML.coor(lens_idx,3):ML.coor(lens_idx,4)) = transmit_sub_aper .* Phase_sub_aper;
    Ein_sub_aper{lens_idx} = Ein(ML.coor(lens_idx,1):ML.coor(lens_idx,2), ML.coor(lens_idx,3):ML.coor(lens_idx,4));
    [x2,y2,Uout] ...
         = ang_spec_prop(Ein_sub_aper{lens_idx}, wvl, delta1, delta2, Dz);
    outputImage = (abs(Uout).^2);
    outputImage = imresize(outputImage,delta2/delta1,'bilinear'); %Scaling the pixel size of the propagated field to the correct one.
    remainingPixels = (size(outputImage,1) -  size(Phase_sub_aper,1))/2;
    if remainingPixels < 0 %padding
        if mod(remainingPixels,2) ~= 0
            outputImage = imtranslate(outputImage,[0.5,0.5]);%%????????
            remainingPixels = remainingPixels - 0.5;                
        end
        padsize = abs([remainingPixels remainingPixels]);
        outputImage = padarray(outputImage, padsize ,0,'both');
        clipImg = outputImage(1:end-1,1:end-1);

%         figure,imagesc(outputImage);

    elseif remainingPixels > 0 %crop
        if mod(remainingPixels,2) ~= 0
            outputImage = imtranslate(outputImage,[0.5,0.5]);
            remainingPixels = remainingPixels - 0.5;
        end
        outputImage = imcrop(outputImage,[remainingPixels+1,remainingPixels+1,size(S,1)-1,size(S,1)-1]);
    end
    Fresnel_subaper_image{lens_idx} = clipImg; 
end


%---------------------------
%rescale the x2-y2 to the size of CCD

%     outputImage = imresize(outputImage,(L2/length(propagated))/pixelSize,'bilinear'); %Scaling the pixel size of the propagated field to the correct one.
%     remainingPixels = (size(outputImage,1) -  size(S,1))/2;
% 
%     if remainingPixels < 0 %padding
% %         if mod(remainingPixels,2) ~= 0
% %             outputImage = imtranslate(outputImage,[0.5,0.5]);%%????????
% %             remainingPixels = remainingPixels - 0.5;                
% %         end
%         padsize = abs([remainingPixels remainingPixels]);
%         outputImage = padarray(outputImage, padsize ,0,'both');
% %     elseif remainingPixels > 0 %crop
% %         if mod(remainingPixels,2) ~= 0
% %             outputImage = imtranslate(outputImage,[0.5,0.5]);
% %             remainingPixels = remainingPixels - 0.5;
% %         end
% %         outputImage = imcrop(outputImage,[remainingPixels+1,remainingPixels+1,size(S,1)-1,size(S,1)-1]);
%     end
end


