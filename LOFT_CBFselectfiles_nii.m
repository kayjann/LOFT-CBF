dirpath=uigetdir('.', 'select subject folder')
cd(dirpath)

folders1.name='ASL'

%select ASL images

cd(folders1.name);
        result=dir('CBF*.nii');
        data = result; %regexp(result(1:end-1),'\n', 'split')';
        numel(data)
        u=1;
        for i=1:numel(data)
            data1{u}=[folders1.name '/' data(i).name];
            u=u+1;
        end
        cd(dirpath)
        clear result

%select M0 image
folders2.name='M0'

cd(folders2.name);
        result=dir('CBF*.nii');
        data = result; %regexp(result(1:end-1),'\n', 'split')';
        numel(data)
        for i=1%:numel(data)
            Mzero{i}=[folders2.name '/' data(i).name];
        end
        cd(dirpath)
        
        
%    params=ls('ASL4D_parms.mat');
%    x=open(deblank(params));
%    SS=x.parms.MRScaleSlope;
   
%LOFT_CBFquantify(images, Mzero, FieldStrength, ASLType, FirstimageType, SubtractionOrder, SubtractionType, PLD, LabelTime, Slicetime, threshold, optionPCA, ASLscaling, M0scaling)

   %change MoCSF value in perf_reconstruct_2D.m for different populations
   %siemens
   LOFT_CBFquantify(data1, Mzero, 1, 3, 1, 0, 0, 2, 1.5, 0, 0.1, 1,1 ,1 ) %change slicetiming  for different sequences
   %philips
   %LOFT_CBFquantify(data1, Mzero, 1, 3, 1, 1, 0, 2, 1.5, 0, 0.1, 1,1, 1) %change slicetiming  for different sequences

   
