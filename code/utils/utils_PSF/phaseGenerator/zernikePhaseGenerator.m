%Generate random or specified zernike phase
%%ELiiiiiii, 20210624
function pupilPlaneMatrix = zernikePhaseGenerator(xysize, kMaxPixels, folderName_zernikePhaseSave, polyMode, polyCoef, stretch)
% Inputs:
%     xysize kMaxPixels
%     folderName_zernikePhaseSave
%     polyMode: the mode of zernike aberation that will be generated
%     polyCoef: the coefficient of zernike aberation that will be generated.
%     stretch: phaseMatrix = phaseMatrix*stretch for large phases

%% initialization
pupilDiameter_pixels = ceil(2*kMaxPixels);
xn=linspace(-1,1,pupilDiameter_pixels);   %normalized x-coordinates 
yn=linspace(-1,1,pupilDiameter_pixels);   %normalized y-coordinates 

%% Wave Aberration definition 
d=2;                %normalized pupil diameter 

if length(polyMode) ~= length(polyCoef)
    error('The input polyMode and polyCoef should be vectors of same size!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Compute Wave Aberration function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
phaseMatrix=zeros(length(xn),length(yn)); 
for modeCount=1:length(polyMode)
    j=polyMode(modeCount);
    n=ceil((-3+sqrt(9+8*j))/2);   %highest power or order of the radial polynomial term 
    m=2*j-n*(n+2);                %azimuthal frequency of the sinusoidal component 
    phaseMatrix=phaseMatrix+polyCoef(modeCount)*zernike(n,m,xn,yn,d); 
end 


phaseMatrix=phaseMatrix * stretch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Plot&Save Wave Aberration function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
gcf = figure('color',[1 1 1]); 
% subplot(2,1,1) 
% W=rot90(W); 
% W=flipud(W); 
imagesc(xn,yn,phaseMatrix); 
colormap gray 
axis xy 
axis square 
%set(gca, 'TickDir', 'out') 
title(['Wave Aberration Function'],'FontSize', 10); 
xlabel('Normalized x pupil coordinate'); 
ylabel('Normalized y pupil coordinate'); 
%colormap gray 
colormap jet 
colorbar
saveas(gcf,[folderName_zernikePhaseSave,'/generatedZernikePhase.jpg']);

% figure('color',[1 1 1]); 
% meshc(xn,yn,phaseMatrix) 
% %view(-37.5,45)  
% title(['Wave Aberration Function'],'FontSize', 10); 
% xlabel('Normalized x pupil coordinate'); 
% ylabel('Normalized y pupil coordinate'); 
% zlabel('RMS Wavefront Error ( \mum)'); 
% colormap('default') 

save([folderName_zernikePhaseSave,'/phaseMatrix.mat'],'phaseMatrix','-v7.3');
pupilPlaneMatrix = expandAndCentralizeAMatrix(phaseMatrix,[xysize,xysize],0);
save([folderName_zernikePhaseSave,'/pupilPlaneMatrix.mat'],'pupilPlaneMatrix','-v7.3');

% close all;