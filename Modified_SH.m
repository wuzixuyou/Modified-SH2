clear, clc, close all;
%  definition of all of the parameter
modes=36;ZerValues(1:15)=rand(15,1); ZerValues(16:modes)=0.1*rand(length(16:modes),1);
resolution=2048; PixelSize=10e-6; LAMBDA=0.78; Lenses=7; focal=10e-3; Prop=1; factor=1e-8; 
bits=16; MLGeometry=1 ;paint=1;MLCentroidMethod=1; MLSharedArea=1; MLAberration=1; 
MLfieldDistortion=0;MLvignetting=1; ObjectMagnitude=0; PupilDiameter=4.2; Exposuretime=1e-2; 
BandWidth=88; QE=0.90; WellCapacity=18e3; PhotonNoise=1; ReadNoise=1; DarkCurrent=0.01; ...
ReadOutNoise=8;
SH_Ap = 'allPupil';
% SH_Ap = 'Circle';

%Create the configuration of the Shack-Hartmann and microlenses array in a struct.
[SH,ML]=createSH(modes,resolution,PixelSize,LAMBDA,Lenses,focal,Prop,factor,bits,MLGeometry,...
    MLCentroidMethod,MLSharedArea,MLAberration,MLfieldDistortion,MLvignetting,WellCapacity,...
    DarkCurrent,ReadOutNoise,paint,SH_Ap);

%Calculate number of photons in reaching detector
nFot = SH.resolution/(10^(ObjectMagnitude/2.5)); %number of photons from the object (fot/s/cm2/Angstrom)
CollectorAreaCm2=pi*((100*PupilDiameter/2)^2); %CollectorArea in cm^2
nFot = round(nFot*Exposuretime*CollectorAreaCm2*BandWidth); 
nFot=nFot*QE;%total photons reaching the detector.

%{
*************************************************************************
********************Calculation of Zernike matrix************************
*************************************************************************
First, search if a Zernike Matrix with the actual number of modes already exist.
%}
N = SH.resolution;
R = SH.radius;
delta = 2*R/N;

[x, y] = meshgrid((-N/2:N/2-1) * delta);
[theta, r] = cart2pol(x,y);

W = 5 * zernike(6, r/R, theta) .* SH.pupil;
% imagesc(W);

%{
*************************************************************************
**** Introducing the incoming wavefront in microns (NOLL NOTATION**********
Zernike polyomials and atmospheric turbulence J Op Soc Am. Vol 66, No 3 ,
************************************March 1976****************************
It is suppose to set the 3 first value to zero (piston and tip/tilt)
%}

WF=W.*1e-6;%conversion to meters
 
%Painting
[SH.pupil4paint]=PaintWavefront(SH,WF);

%{
*************************************************************************
***************Definition of the microlenses array***********************
*************************************************************************
%}
[ML]=MicroLenses(SH,ML,WF);

%{
By defect, Microlenses are perfect shperical lenses according to the
ideal Fraunhofer propagation. If some aberrations has to be added to
microlenses, choose the zernike to be implemented and the amount of microns
in the following form:
ML.Aberration=    %This select the Zernike polynomial to be added in Noll
notation (i.e. ML.Aberration=11)
ML.AberrationAmount=  %Amount of zernike polynomial in meters
(i.e. ML.AberrationAmount= 1e-6)
%}

WFSubpupil=(WF+ML.AmplitudeMask);%Wavefront in each microlent (phase of 
%microlenses itself is not yet taking into account

Lenses=cell(1,length(ML.coor));%Preallocation
for i=1:length(ML.coor)
    Lenses{i}=WFSubpupil(ML.coor(i,1):ML.coor(i,2), ML.coor(i,3):ML.coor(i,4));
end


%{
*************************************************************************
*************** Point Spread Function Calculation ***********************
*************************************************************************
%}
PSF=cell(1,length(Lenses));
PSFmaxLocal=zeros(1,length(Lenses));
if ML.Prop==0 %CASE OF NO PROPAGATION BETWEEN MICROLENSES AND CCD. Edge effects will exist, it should be implemented a pad with zeros to avoid it
    fprintf('Propagation between microlenses and CCD is not considered. Microlenses are considered as pure amplitude objects\n')
    for i=1:length(Lenses)
        PF=ML.Pupil.*exp(-1i*SH.k.*Lenses{i});%Pupil function
        PSF{i}=abs(ifftshift(ifft2(fftshift(PF)))).^2;%PSF
        PSFmaxLocal(i)=max(max(PSF{i}));
    end
else
    L=length(Lenses{1})*SH.PixelSize; %Size of the region of each microlenses
    for j=1:length(Lenses)
        S = ML.Pupil.*exp(-1i*SH.k.*Lenses{j}); %complex phase screen
        
        if isfield(ML,'Aberration') == 1 %If aberration o each microlens itself is considered
            aberration = zernike(ML.AberrationZernike,size(S,1));
            S = S.*exp(-1i*SH.k*aberration*ML.AberationValue*SH.LAMBDA);
        end
        
        PSF{j} = lensletSimulation(S,L,SH.LAMBDA,ML.focal,SH.PixelSize);
        
%         if isfield(ML,'fieldDistortion') == 1
%             PSF{j} = applySeidelDistortion(ML.focal,SH.PixelSize,0.1*SH.LAMBDA,size(S,1),PSF{j});
%         end
        
        if isfield(ML,'vignetting') == 1
            PSF{j} = applyVignetting(ML.focal,SH.PixelSize,size(S,1),1e-3,PSF{j});
        end
        
        PSFmaxLocal(j)=max(max(PSF{j}));
    end
end
%{
***************************************************************************
****** Quantization of the signal, introduction of photonic noise**********
****************** & introduction of read noise****************************
***************************************************************************
Photon noise, also known as Poisson noise, is a basic form of uncertainty 
associated with the measurement of light, inherent to the quantized nature 
of light and the independence of photon detections. Its expected magnitude 
constitutes the dominant source of image noise except in low-light conditions.
Individual photon detections can be treated as independent events that 
follow a random temporal distribution. As a result, photon counting is a 
classic Poisson process.Photon noise is signal  dependent, and its standard 
deviation grows with the square root of the signal. Contrary to popular 
belief, shot noise experienced by the detector IS related to the QE of the 
detector! Back-illuminated sensors with higher QE yields a better 
Signal/Shot Noise ratio. There is a simple intuitive explanation for this ?
 shot noise must be calculated from the signal represented by the number 
of photoelectrons in the sensor (electrons generated from photons falling 
on the sensor), NOT JUST from the number of incoming photons. Therefore, 
if an average of 100 photons hit a pixel, but the sensor has a QE of 50% at
the wavelength of these photons, then an average of 50 photoelectrons will 
be created ?the shot noise and Signal/Shot Noise must be calculated from 
this value.
%}
if SH.bits==0
    fprintf('Calculations performed with "double" precision of MATLAB. \n');
    for j=1:length(PSF)
        PSF{j}=PSF{j}/max(max(PSF{j}));
    end
else
    PSFmax=max(max(PSFmaxLocal));
    fprintf('Calculations performed for %2.0f bits. \n',SH.bits);
    for i=1:length(PSF)
        PSF{i}=round((PSF{i}/PSFmax)*((2^SH.bits)-1));
    end
end

% paintingBigPSF.mat included the code that create the PSF "seen" by each 
% microlens, even when there exist energy from surrounding microlenses. It
% also includes the introduction of photonic and read noise.
if ML.radiusPixels*2==length(PSF{1})
    [PSF,idealPSF]=paintingPSF(PSF,ML,SH,nFot,0,1,PhotonNoise,ReadNoise,Exposuretime);
    drawnow();
elseif ML.radiusPixels*2~=length(PSF{1})
    [PSF,idealPSF]=paintingBigPSF(PSF,ML,SH,nFot,0,1,PhotonNoise,ReadNoise,Exposuretime);
    drawnow();
end



%{
*************************************************************************
*************** Defocus image Calculation ***********************
*************************************************************************
%}