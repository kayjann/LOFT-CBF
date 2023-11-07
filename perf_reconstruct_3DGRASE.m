% function [gcbf,cbfdat] = perf_resconstruct_3D (Filename, FieldStrength, ASLType, FirstimageType, SubtractionOrder, SubtractionType, ThreshFlag, threshold, CBFFlag, MeanFlag, Delaytime, Labeltime, Slicetime, T1b, PASLMo, optionPCA, ASLscaling, M0scaling, diffs)
%  
% This MATLAB function is to reconstruct the raw perfusion images from EPI images by the subtraction between labelled images 
% and control images. Quantified CBF images can also be reconstructed by select the option. It is based on SPM2 and 
% MATLAB 6. It is also comparable with SPM5 and Matlab7.
%    
% The MATLAB code of this function will be found in: http://cfn.upenn.edu/perfusion/software.htm
%    
% All the images are 3D SPM ANALYZE formatted (.img and .hdr). All the results are also saved in SPM ANALYZE format; 
% The labelled and control images should be the data after motion correction.
%    
% The method used here are based on the "simple subtraction", "surround subtraction" and "sinc subtraction" approaches described in
% Aguirre GK et al (2002) Experimental design and the relative sensitivity of perfusion and BOLD fMRI, NeuroImage. 15:488-500. 
%    
% BOLD data (or whatever the underlying pulse sequence that was used) are generated in addition to the perfusion data
%    
% for CASL and pCASL,
% CBF data are calculated according to the formula from
% Wang J, Alsop DC, et al. (2003) Arterial transit time imaging with flow encoding arterial spin tagging (FEAST).
% Magn Reson Med. 50:599-607. Page600, formula [1]
% CBF_CASL (ml/100g/min) = 60*100*deltaM*¦Ë*R/(2*alp*Mo*(exp(-w*R)-exp(-(t+w)*R))
% where deltaM = raw ASL signal (Control-Label)
%        ¦Ë = blood/tissue water partition coefficient, R =longitudinal relaxation rate of blood,
%       alp = tagging efficiency, Mo =  equilibrium magnetization of brain, 
%       w = post-labeling delay, t = duration of the labeling pulse,  
% and we use the assumed parameters for calculation as ¦Ë=0.9g/ml, 
% for 3T, alp=0.68, T1b=1650ms, R=1/T1b=0.606sec-1. 
% for 1.5T, alp=0.71, T1b=1200ms, R=1/T1b=0.83sec-1.                                                      
%
% for PASL,
% CBF data are calculated according to the formula from
% Wang J, Aguirre GK, et al. (2003) Arterial Spin Labeling Perfusion fMRI With Very Low Task Frequency
% Magn Reson Med. 49:796-802. Page798, formula [1]
% CBF_PASL (ml/100g/min) = 60*100*deltaM*¦Ë/(2*alp*Mo*t*exp(-(t+w)*R))
% where deltaM = raw ASL signal (Label-control)
%        ¦Ë = blood/tissue water partition coefficient, R =longitudinal relaxation rate of blood,
%       alp = tagging efficiency, Mo =  equilibrium magnetization of brain, 
%       w = Post Inf Sat delay, t = Post IR delay (TI1)  
% and we use the assumed parameters for calculation as ¦Ë=0.9g/ml, 
% for 3T, alp=0.68, T1b=1650ms, R=1/T1b=0.606sec-1. 
% for 1.5T, alp=0.95, T1b=1200ms, R=1/T1b=0.83sec-1.                                                      
%
%  Inputs:
%    Firstimage - integer variable indicating the type of first image 
%    - 0:control; 1:labeled 
%   Select raw images (*.img, images in a order of control1.img, label1.img, control2.img, label2.img,....;
%   or images in a order of label1.img, control1.img, label2.img, control2.img, .... )
%    
%    SubtractionType - integer variable indicating which subtraction method will be used 
%    -0: simple subtraction; 1: surround subtraction;2: sinc subtractioin.
%    for CASL, suppose Perfusion = Control - Label;
%    if the raw images is: (L1, C1, L2, C2...)
%     the simple subtraction is: (C1-L1, C2-L2...)
%     the surround subtraction is: (C1-(L1+L2)/2, C2-(L2+L3)/2,...)
%     the sinc subtraction is: (C1-L3/2, C2-L5/2...)
%
%    for PASL, suppose Perfusion = Label - Control;
%    if the raw images is: (L1, C1, L2, C2...)
%     the simple subtraction is: (L1-C1, L2-C2...)
%     the surround subtraction is: ((L1+L2)/2-C1, (L2+L3)/2-C2,...)
%     the sinc subtraction is: (L3/2-C1, L5/2-C2...)
%    
%  Outputs:
%    BOld Images: Bold_*.img,Bold_*.hdr;  Mean_Bold.img, Mean_Bold.hdr; 
%    Perfusion Images: Perf_*.img, Perf_*.hdr; Mean_Perf.img, Mean_Perf.hdr;
%    CBF Images: CBF_*.img, CBF_*.hdr; Mean_CBF.img, Mean_CBF.hdr;
%    
%  By H.Y. Rao & J.J. Wang, @CFN, UPenn Med. 07/2004.
%  Updated for SPM5 comparable 12/2009
%  add ASLscaliong by Danny 08/2017


function [gcbf,cbfdat, mCBF, Mask] = perf_resconstruct_3DGRASE(Filename, FieldStrength, ASLType, FirstimageType, SubtractionOrder, SubtractionType, ThreshFlag, threshold, CBFFlag, MeanFlag, Delaytime, Labeltime, Slicetime, T1b, PASLMo, optionPCA, ASLscaling, M0scaling, diffs);


% close all;
spm('ver',[],1);

% try Filename;
%    ;
% catch
%   if strcmp(spm('ver',[],1),'SPM8')
%      Filename_lab=spm_select(Inf,'any','Select Labeling imgs', [],pwd,'.*img');
%      Filename_con=spm_select(Inf,'any','Select Control imgs', [],pwd,'.*img');
%   else
%     Filename_lab = spm_get(Inf,'*.img','Select Labeling imgs');
%     Filename_con = spm_get(Inf,'*.img','Select Control imgs');
%   end;
%   
%   if length(Filename_lab) ~= length(Filename_con), fprintf('the labeling image number is equal to the control'); return;end;
%   
%   for i=1:size(Filename_lab,1)
%       Filename(2*i-1,:) = Filename_lab(i,:);
%       Filename(2*i,:) = Filename_con(i,:);
%   end
% 
%   if isempty(Filename), fprintf('No images selected!\n');return;end;
% 
%paranum = 1; pos=1;
% end;
% CASLmask = spm_select(1,'image','Select mask image'); 

% FieldStrength = spm_input('Scanner Strength: 1:3T; 2:1.5T', '+1', 'e', 1);
%  paranum = paranum + 1;
%  
% ASLType = spm_input('ASLType: 2:pCASL; 1:CASL; 0:PASL', '+1', 'e', 2);
%   paranum = paranum + 1;

% if ASLType == 0 || ASLType == 2 || ASLType == 3;
%   if strcmp(spm('ver',[],1),'SPM8')
%     PASLMo =spm_select(1,'any','Select Mo imgs', [],pwd,'.*img');
%   else
%      PASLMo = spm_get(1,'*.img','Select PASL Mo image'); 
%      paranum = paranum + 1;
%    end;
   if isempty(PASLMo) fprintf('No PASL Mo images selected!\n');return;end;
 
% end;

 %FirstimageType = spm_input('Select 1st Image type? 0:control; 1:labeled', '+1', 'e', 1);
 %paranum = paranum + 1;
%  FirstimageType=1;

%SubtractionOrder = spm_input('Select SubtractionOrder', '+1', 'm',  ['*Even-Odd(Img2-Img1)|Odd-Even(Img1-Img2)'], [0 1], 0);
%   paranum = paranum + 1;
%SubtractionOrder=1;

%SubtractionType = spm_input('Selct SubtractionType', '+1', 'm',  ['*Simple |Surround|Sinc'], [0 1 2], 0);
% paranum = paranum + 1;
%  SubtractionType=0;

if SubtractionType==2, 			
  %Timeshift = spm_input('Time shift of sinc interpolation', '+1', 'e', 0.5);
  Timeshift=0.5
%  paranum = paranum + 1;
end;

 %CBFFlag = spm_input('Produce quanperf_resconstructtified CBF images? 0:no; 1:yes', '+1', 'e', 1);
 %paranum = paranum + 1;
  %CBFFlag=1;


 %ThreshFlag = spm_input('Threshold EPI images? 0:no; 1:yes', '+1', 'e', 1);
 %paranum = paranum + 1;
  %ThreshFlag=1;

% if ThreshFlag==1,
%   threshold =  spm_input('Input EPI Threshold value', '+1', 'e', 0.8);
%   paranum = paranum + 1;
% end;
% absthreshold=200;

 %MeanFlag = spm_input('Produce mean images? 0:no; 1:yes', '+1', 'e', 1);
 %paranum = paranum + 1;
 % MeanFlag=1;
  
% if CBFFlag==1,
  if ASLType==2 || ASLType==3 %pCASL
%    %Labeltime = spm_input('Enter Label time:sec', '+1', 'e', 1.2);
%    %Delaytime = spm_input('Enter Delay time:sec', '+1', 'e', 1.0);
%    %Slicetime = spm_input('Enter Slice acquisition time:msec', '+1', 'e', 0);
    alp = 0.80;   %pCasl tagging efficiency with correction of background suppression
    alp = 0.73      %new tagging efficiency for Siemens with background suppression
%    if FieldStrength == 1, R = 0.606; else R = 0.83; end;  %longitudinal relaxation rate of blood
%     paranum = paranum + 5;
   elseif   ASLType ==0 %CASL
%    %Labeltime = spm_input('Enter Label time:sec', '+1', 'e', 1.6);
%    %Delaytime = spm_input('Enter Delay time:sec', '+1', 'e', 1.2);
%    %Slicetime = spm_input('Enter slice acquisition time:msec', '+1', 'e', 0);
%    if FieldStrength == 1, alp = 0.68; else alp = 0.71; end;   %Casl tagging efficiency
%    if FieldStrength == 1, R = 0.606; else R = 0.83; end;  %longitudinal relaxation rate of blood
%     paranum = paranum + 5;
   else  %PASL
%     %Labeltime = spm_input('Enter Post IR Delay time:sec', '+1', 'e', 0.7); % TI1
%     %Delaytime = spm_input('Enter Post Inf Sat Delay time:sec', '+1', 'e', 1.2);
%     %Slicetime = spm_input('Enter slice acquisition time:msec', '+1', 'e', 42);
     alp = 0.95;   %PASL tagging efficiency
%     if FieldStrength == 1, R = 0.606; else R = 0.83; end;  %longitudinal relaxation rate of blood
%     paranum = paranum + 5;
%    end;
  end;
 
 %T1b = spm_input('Enter blood T1 : msec', '+1', 'e', 1650);  %you can input the updated blood T1
   R = 1000/T1b;

%  if ASLType ==2 %pCASL
%   alp =  spm_input('Enter label efficiency', '+1', 'e', 0.80);
% end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the main program
% [Finter,Fgraph,CmdLine] = spm('FnUIsetup','Perf Reconstruct',0);
% spm('FigName','Perf Reconstruct: working',Finter,CmdLine);
% spm('Pointer','Watch')

tmp=Filename(1,:);
prefix=[tmp(1:length(tmp)-3) 'txt']
clear tmp;

% Map images
V=spm_vol(deblank(Filename));

if ASLType==0 || ASLType == 2 || ASLType==3
  VMo = spm_vol(deblank(PASLMo)); 
  PASLModat = zeros([VMo.dim(1:2) 1]);
end;
 
 if length(V)==0, fprintf('no raw img files was selected'); return; end;
 if rem(length(V),2)==1, warning('the number of raw img files is not even, last img is ignored'); end;
 perfnum=fix(length(V)/2);
 
% Create output images...
VO = V(1:perfnum);
VB = V(1:perfnum);
VCBF=V(1:perfnum);
VMP = V(1);
VMCBF = V(1);
VMC = V(1);

  for k=1:length(V),
        [pth,nm,xt] = fileparts(deblank(V(k).fname)); %vr
        if SubtractionType==0, 
         VO(k).fname = fullfile(pth,['Perf_0' nm xt]);
         if CBFFlag==1, VCBF(k).fname = fullfile(pth,['CBF_0_' nm xt ]);end;
        end;
        if SubtractionType==1, 
         VO(k).fname = fullfile(pth,['Perf_1' nm xt]);
         if CBFFlag==1, VCBF(k).fname = fullfile(pth,['CBF_1_' nm xt ]);end;
        end;
        if SubtractionType==2, 
         VO(k).fname = fullfile(pth,['Perf_2' nm xt]);
         if CBFFlag==1, VCBF(k).fname = fullfile(pth,['CBF_2_' nm xt ]);end;
        end;
        VB(k).fname    = fullfile(pth,['Bold_' nm xt ]);
  end;

  for k=1:perfnum,
           VO(k)  = spm_create_vol(VO(k));
           VB(k)  = spm_create_vol(VB(k));
           VCBF(k)  = spm_create_vol(VCBF(k));
          %if strcmp(spm('ver',[],1),'SPM8')   
            VO(k).dt=[16,0]; VB(k).dt=[16,0]; VCBF(k).dt =[16,0];  %'float' type
         %else
           %VO(k).dim(4) = 16; VB(k).dim(4) = 16; VCBF(k).dim(4) = 16; %'float' type
          %end;
  end;
  
cdat = zeros([VO(1).dim(1:3) perfnum]);
ldat = zeros([VO(1).dim(1:3) perfnum]);
pdat = zeros([VO(1).dim(1:3) perfnum]);
bdat = zeros([VB(1).dim(1:3) perfnum]);

linear_cdat=zeros([VB(1).dim(1:3) 2*perfnum]);
linear_ldat=zeros([VB(1).dim(1:3) 2*perfnum]);
sinc_ldat=zeros([VB(1).dim(1:3) 2*perfnum]);
sinc_cdat=zeros([VB(1).dim(1:3) 2*perfnum]);


%-Start progress plot
%-----------------------------------------------------------------------
%spm_progress_bar('Init',perfnum,'Perf Reconstruct','Images completed');

% read raw data
temp=load_untouch_nii(Filename(1,:));
dat = spm_read_vols(V);
dat=dat/temp.hdr.dime.scl_slope;

threshvalue = zeros(1,length(V));

% threshold the EPI images 
  %Mask=zeros(V(1).dim(1:3));
  
% read the Mo data
%  if ASLType==0 || ASLType==2 || ASLType==3;
    tempMo=load_untouch_nii(PASLMo);

    PASLModat = spm_read_vols(VMo); 
    PASLModat=PASLModat/tempMo.hdr.dime.scl_slope;
    
    PASLModat=PASLModat./M0scaling;
    %Mask = Mask.*(PASLModat>threshold*mean(mean(mean(PASLModat))));
    %Mask = Mask.*(PASLModat>20);
%  end;
  
  if ThreshFlag ==1,
      Mask=zeros(size(PASLModat));
      Mask(find(PASLModat>(threshold*(max(PASLModat(:))))))=1;
      %Mask=repmat(Mask,[1 1 1 length(V)]);
%    for k=1:length(V),
%      Mask = Mask.*(dat(:,:,:,k)>threshold*mean(mean(mean(dat(:,:,:,k)))));
%      threshvalue(1, k) = max(20, threshold*mean(mean(mean(dat(:,:,:,k)))));
% %     Mask = Mask.*(dat(:,:,:,k)>absthreshold);
%    end;
  end;


 for k=1:length(V),
%    datamk= spm_read_vols(V(k));
%   datamk = datamk.*Mask;
    dat(:,:,:,k) = (dat(:,:,:,k).*Mask)./ASLscaling;   %Danny add ASLscaling 
 end;

% define the control and label images...
 for k=1:length(V),
  if SubtractionOrder==0, 
      if rem(k,2)== 1, ldat(:,:,:,(k+1)/2) = dat(:,:,:,k); end;
      if rem(k,2)== 0, cdat(:,:,:,k/2) = dat(:,:,:,k); end;
  end;
  if SubtractionOrder==1, 
      if rem(k,2)== 1, cdat(:,:,:,(k+1)/2) = dat(:,:,:,k); end;
      if rem(k,2)== 0, ldat(:,:,:,k/2) = dat(:,:,:,k); end;
  end;
 end;
 
 
 % obtained BOLD data
 for k=1:perfnum,
  bdat(:,:,:,k) = (dat(:,:,:,2*k-1) + dat(:,:,:,2*k))/2;
 end;
 
if optionPCA==0
 % do the simple subtraction...
if SubtractionType==0,
  for k=1:perfnum,
    pdat(:,:,:,k) = cdat(:,:,:,k) - ldat(:,:,:,k);
  end;
 spm_progress_bar('Set',k);
end;
 
  % do the linear interpolation...
  if SubtractionType==1,
     pnum=1:perfnum;
     lnum=1:0.5:perfnum;
     for x=1:V(1).dim(1),
      for y=1:V(1).dim(2),
       for z=1:V(1).dim(3),
        cdata = zeros(1,perfnum);
        ldata = zeros(1,perfnum);
        linear_cdata = zeros(1,length(V));
        linear_ldata = zeros(1,length(V));
         for k=1:perfnum, 
          cdata(k) = cdat(x,y,z,k);
          ldata(k) = ldat(x,y,z,k);
         end;
         linear_cdata=interp1(pnum,cdata,lnum);
         linear_ldata=interp1(pnum,ldata,lnum);
         for k=1:2*perfnum-1, 
          linear_cdat(x,y,z,k)= linear_cdata(k);
          linear_ldat(x,y,z,k)= linear_ldata(k);
         end;
        end; 
       end; 
      end; 

   
     % do the surround subtraction....
     if FirstimageType ==1; 
          pdat(:,:,:,1) = cdat(:,:,:,1) - ldat(:,:,:,1);
          spm_progress_bar('Set',1);
        for k=2:perfnum, 
          pdat(:,:,:,k) = linear_cdat(:,:,:,2*(k-1)) - ldat(:,:,:,k);
          spm_progress_bar('Set',k);
        end;
     end;
     if FirstimageType ==0; 
          pdat(:,:,:,1) = cdat(:,:,:,1) - ldat(:,:,:,1);
          spm_progress_bar('Set',1);
       for k=2:perfnum, 
          pdat(:,:,:,k) = cdat(:,:,:,k) - linear_ldat(:,:,:,2*(k-1));
          spm_progress_bar('Set',k);
        end;
     end;
end;


 % do the sinc interpolation...
  if SubtractionType==2,
     for x=1:V(1).dim(1),
       for y=1:V(1).dim(2),
         for z=1:V(1).dim(3),
           cdata = zeros(1,perfnum);
           ldata = zeros(1,perfnum);
           sinc_cdata = zeros(1,length(V));
           sinc_ldata = zeros(1,length(V));
           for k=1:perfnum, 
             cdata(k) = cdat(x,y,z,k);
             ldata(k) = ldat(x,y,z,k);
           end;
           sincnum = fix(perfnum/Timeshift);
           sinc_cdata=interpft(cdata,sincnum);
           sinc_ldata=interpft(ldata,sincnum);
           for k=1:2*perfnum, 
            sinc_cdat(x,y,z,k)= sinc_cdata(k);
            sinc_ldat(x,y,z,k)= sinc_ldata(k);
           end;
         end;  
       end;
     end;
 
      % do the sinc subtraction....
         if FirstimageType ==1; 
          pdat(:,:,:,1) = cdat(:,:,:,1) - ldat(:,:,:,1);
             for k=2:perfnum, 
               pdat(:,:,:,k) = sinc_cdat(:,:,:,2*(k-1)) - ldat(:,:,:,k);
               spm_progress_bar('Set',k);
             end;
          end;
         if FirstimageType ==0; 
           pdat(:,:,:,1) = cdat(:,:,:,1) - ldat(:,:,:,1);
             for k=2:perfnum, 
               pdat(:,:,:,k) = cdat(:,:,:,k) - sinc_ldat(:,:,:,2*(k-1));
               spm_progress_bar('Set',k);
           end;
         end;
  end;      
 elseif optionPCA==1
     diffs=(diffs.*permute(repmat(Mask,[1,1,1,size(diffs,1)]),[4 1 2 3]))./ASLscaling;
     pdat=permute(diffs,[2,3,4,1]);
 end

 % Write Bold and perfusion image...
   for k=1:perfnum,
      VO(k) = spm_write_vol(VO(k),pdat(:,:,:,k));
      VB(k) = spm_write_vol(VB(k),bdat(:,:,:,k));
   end;


 % calculated the mean image...
  if MeanFlag ==1,
    Mean_dat=zeros([V(1).dim(1:3)]);
    VMP.fname = fullfile(pth,['AMean_Perf' nm(1:5) xt]);
    VMP = spm_create_vol(VMP);
    %if strcmp(spm('ver',[],1),'SPM8')   
       VMP.dt=[16,0];   %'float' type
       VMP.pinfo(1)=1;
     % else
      %VMP.dim(4) = 16; %'float' type
    %end;
    
    for x=1:V(1).dim(1),
       for y=1:V(1).dim(2),
         for z=1:V(1).dim(3),
          Mean_dat(x,y,z) = mean(pdat(x,y,z,:));
         end;
       end;
    end;
    
    % Write mean perfusion image...
        VMP = spm_write_vol(VMP,Mean_dat);
  end;


  % calculated the CBF image...
  if CBFFlag ==1,
     
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    spm_progress_bar('Clear')
% 
%   [Finter,Fgraph,CmdLine] = spm('FnUIsetup','Perf Reconstruct',0);
%   spm('FigName','CBF Reconstruct: working',Finter,CmdLine);
%   spm('Pointer','Watch')
%   %-----------------------------------------------------------------------
%   spm_progress_bar('Init',perfnum,'CBF Reconstruct','Images completed');

     cbfdat = zeros([VO(1).dim(1:3), perfnum]);
     cmean_dat = zeros([VO(1).dim(1:3)]);

     for x=1:V(1).dim(1),
       for y=1:V(1).dim(2),
         for z=1:V(1).dim(3),
          cmean_dat(x,y,z) = mean(cdat(x,y,z,:));
         end;
       end;
     end;
   
     % Write mean BOLD/Control image...
   if MeanFlag ==1,
     VMC.fname = fullfile(pth,['AMean_BOLD' nm(1:5) xt]);
     VMC = spm_create_vol(VMC);
   %  if strcmp(spm('ver',[],1),'SPM8')   
            VMC.dt=[16,0];   %'float' type
            VMC.pinfo(1)=1;
    %  else
     %   VMC.dim(4) = 16; %'float' type
      % end;
       
     VMC = spm_write_vol(VMC,cmean_dat);
   end;
   
     for k=1:perfnum, 
       for x=1:V(1).dim(1),
       for y=1:V(1).dim(2),
       for z=1:V(1).dim(3),
        Dtime = Delaytime + Slicetime*z/1000;

        if ASLType ==2 || ASLType==3  % pCASL
         if cmean_dat(x,y,z)<mean(threshvalue)
          cbfdat(x,y,z,k)=0;
         else
          %cbfdat(x,y,z,k) = 2700*pdat(x,y,z,k)*R/alp/((exp(-Dtime*R)-exp(-(Dtime+Labeltime)*R))*cmean_dat(x,y,z));
          cbfdat(x,y,z,k) = 2700*pdat(x,y,z,k)*R/alp/((exp(-Dtime*R)-exp(-(Dtime+Labeltime)*R))*PASLModat(x,y,z));
         end;
        end;

        if ASLType ==0   % CASL
         if cmean_dat(x,y,z)<mean(threshvalue)
          cbfdat(x,y,z,k)=0;
         else
          cbfdat(x,y,z,k) = 2700*pdat(x,y,z,k)*R/alp/((exp(-Dtime*R)-exp(-(Dtime+Labeltime)*R))*cmean_dat(x,y,z));
         end;
        end;

        if ASLType ==1  %PASL
         if (PASLModat(x,y,z)<mean(threshvalue) | cmean_dat(x,y,z)<mean(threshvalue))
          cbfdat(x,y,z,k)=0;
         else
          cbfdat(x,y,z,k) = 2700*pdat(x,y,z,k)/Labeltime/alp/(exp(-(Dtime+Labeltime)*R)*PASLModat(x,y,z));
         end;
        end;
        
       end;
       end;
       end;
  
      % Write CBF images...
      VCBF(k) = spm_write_vol(VCBF(k),cbfdat(:,:,:,k));
      spm_progress_bar('Set',k);

     end;

  
    if MeanFlag ==1,
     Mean_cbfdat=zeros([VO(1).dim(1:3)]);
     VMCBF.fname = fullfile(pth,['AMean_CBF' nm(1:5) xt]);
     VMCBF = spm_create_vol(VMCBF);
     % if strcmp(spm('ver',[],1),'SPM8')   
        VMCBF.dt=[16,0];   %'float' type
        VMCBF.pinfo(1)=1;
      %else
       % VMCBF.dim(4) = 16; %'float' type
       %end;

    voxelnum=0;
    zeronum=0;
    globalCBF=0;
    meancontrol=0;
    
     for x=1:V(1).dim(1),
       for y=1:V(1).dim(2),
         for z=1:V(1).dim(3),
            Mean_cbfdat(x,y,z) = mean(cbfdat(x,y,z,:));
            if Mean_cbfdat(x,y,z) ==0,  
              zeronum = zeronum+1;
            else
              voxelnum = voxelnum+1;
              globalCBF = globalCBF+Mean_cbfdat(x,y,z);
              meancontrol = meancontrol+cmean_dat(x,y,z); 
            end;
         end;
       end;
     end;
 
   globalCBF = globalCBF/voxelnum;
   meancontrol = meancontrol/voxelnum;

     % Write mean CBf image...
     VMCBF = spm_write_vol(VMCBF, Mean_cbfdat);
    end;
    
end;

 mCBF=Mean_cbfdat(:); MASKm=Mask(:);
 mCBF=mean(mCBF(find(MASKm))); 

 gcbf = spm_global(VMCBF);
 gbold = spm_global(VMC);

save globalCBF globalCBF %EDIT(H)
 
  fprintf('\n\t Perfusion images written to: ''%s'' to \n',VO(1).fname)
  fprintf('\t  ''%s''. \n',VO(perfnum).fname)
  if MeanFlag ==1, 
    fprintf('\t Mean_perf image written to: ''%s''.\n\n',VMP.fname)
  end;  

  fprintf('\t BOLD images written to: ''%s''  to \n',VB(1).fname)
  fprintf('\t ''%s'' . \n',VB(perfnum).fname)
  if MeanFlag ==1, 
   fprintf('\t Mean_BOLD image written to: ''%s''.\n\n',VMC.fname)
  end;  

  if CBFFlag ==1, 
     fprintf('\t Quantified CBF images written to: ''%s'' to \n',VCBF(1).fname)
     fprintf('\t ''%s'' .\n',VCBF(fix(length(VCBF)/2)).fname)
     fprintf('\t Mean Quantified CBF image written to: ''%s'' \n\n',VMCBF.fname)

     fprintf('\t the spm global mean BOLD control signal is:')
     fprintf('\t %6.2f \n',gbold)
     fprintf('\t the spm global mean CBF signal is:')
     fprintf('\t %6.3f \n',gcbf)
     fprintf('\t ml/100g/min \n\n')
     
     fprintf('\t the calculated voxel number is:')
     fprintf('\t %8.1f\n',voxelnum)
     fprintf('\t the zero number is:')
     fprintf('\t %8.1f\n',zeronum)
     fprintf('\t the global mean BOLD control signal is:')
     fprintf('\t %6.2f \n',meancontrol)
     fprintf('\t the global mean CBF signal is:')
     fprintf('\t %6.3f',mCBF)
     fprintf('\t ml/100g/min \n\n')
  end;  

   fid=fopen(prefix,'w');
     fprintf(fid,' \n BIOCLINICA ASL BIOGEN trial \n \n');
     %fprintf(fid,' \n K. Jann @ University of Southern California, Los Angeles (version 1.1 Dez 2016) \n \n');
     fprintf(fid,' \n File created:    %s \n \n', datestr(clock));
     fprintf(fid,' \n TXT File stored under :    %s \n \n',prefix);

     fprintf(fid,'\t Quantified CBF images written to: ''%s'' \n',VCBF(1).fname);
     fprintf(fid,'\t Mean Quantified CBF image written to: ''%s'' \n\n',VMCBF.fname);

     fprintf(fid,'\t the spm global mean BOLD control signal is:');
     fprintf(fid,'\t %6.2f \n',gbold);
     
     fprintf(fid,'\t the calculated voxel number is:');
     fprintf(fid,'\t\t %8.1f\n',voxelnum);
     fprintf(fid,'\t the zero number is:');
     fprintf(fid,'\t\t\t %8.1f\n',zeronum);
     fprintf(fid,'\t the global mean BOLD control signal is:');
     fprintf(fid,'\t\t %6.2f \n\n',meancontrol);
     
     fprintf(fid,'\t [ml/100g/min] \n\n');
     fprintf(fid,'\t ---------------------------------------------------------------------\n');
     if ASLType ==0
     fprintf(fid,'\t -------continuous ASL [CASL]-----------------------------\n');
     end;
     if ASLType ==1
     fprintf(fid,'\t -------pulsed ASL [PASL]---------------------------------\n');
     end;
     if ASLType ==2
     fprintf(fid,'\t -------pseudoCASL [pCASL]--------------------------------\n');
     end;
     if ASLType ==3
     fprintf(fid,'\t -------3D pseudoCASL [pCASL]--------------------------------\n');
     end;
     fprintf(fid,'\t ---------------------------------------------------------------------\n\n');
     fprintf(fid,'\t the  global mean CBF signal is     :');
     fprintf(fid,'\t %6.3f \t [ml/100g/min] \n\n',mCBF);
     fprintf(fid,'\t ---------------------------------------------------------------------\n');
     
     
     
     
     fclose(fid);
  
 %spm_progress_bar('Clear')
  
 fprintf('......computing done.\n\n')

%  spm('Pointer');
%  spm('FigName','CBFReconstruct done',Finter,CmdLine);

% 
