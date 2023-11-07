% scale Mzero image due Grappa acceleration. 


function adjustMzero(Filename,PLD,slicetime)
% slicetime and PLD in milliseconds!!!!!
addpath('/Users/loft12/DATA_KJ/Utilities/SPM12/')


slicetime;
TI2=PLD;
T1=1600; %T1 gray matter at 3T 


VMo=spm_vol(deblank(Filename));
PASLModat = spm_read_vols(VMo); 
 

for slice=1:VMo.dim(3)
    Moadj(:,:,slice)=PASLModat(:,:,slice)/(1-exp(-((TI2+slicetime*slice)/T1))); 
end


[pth,nm,xt] = fileparts(deblank(VMo.fname));
VMo.fname = fullfile(pth,['adjusted_' nm xt]);    
VMo= spm_create_vol(VMo);
VMo= spm_write_vol(VMo,Moadj);

end
