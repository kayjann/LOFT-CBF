% CBF quantification toolbox. 
addpath('/ifs/loni/faculty/kjann/Utilities/NIFTI//')
addpath('/ifs/loni/faculty/kjann/Utilities/LOFT_CBF//')



%% select and prepare data

%workingdir= uigetdir('.' , 'Select Data Directory')
%cd(workingdir)
 
    
    clear;clc;close all;
    
    %expand 4D nii file into separate volumes
    [ASLfile, subjectfolder] = uigetfile('*.nii', 'Select 4D ASL image');
    cd(subjectfolder);
    display(['working on: ASLfile']) 
    display('expanding ASL file')
    
    %if file is  gz-zipped
    %unzip nii.gz file
       % zipfile=ASLfile;
       % ASLfile=ASLfile(1:end-3);
       % gunzip(zipfile);
        expand_nii_scan(ASLfile);
    
   
    %sort M0 and ASL data into subfolders
    rmdir('ASL','s')
    rmdir('M0','s')
    
    mkdir('ASL');
    mkdir('M0');
    ASLdata=dir([ASLfile(1:end-4) '_0*.nii'])
    for dataASL=1:2
        movefile(ASLdata(dataASL).name, './M0/')
    end
    for dataASL=3:numel(ASLdata)
        movefile(ASLdata(dataASL).name, './ASL/')
    end
    display('selecting data to quantify CBF')
            %select ASL images

            cd('ASL');
                result=dir([ASLfile(1:4), '*.nii']);
                data = result; %regexp(result(1:end-1),'\n', 'split')';
                numel(data);
                u=1;
                for i=1:numel(data)
                    data1{u}=['ASL/' data(i).name];
                    u=u+1;
                end
                cd([subjectfolder]);
                clear result
            %select M0 image
            cd('M0');
                result=dir([ASLfile(1:4), '*.nii']);
                data = result; %regexp(result(1:end-1),'\n', 'split')';
                numel(data);
                for i=1:numel(data)
                    Mzero{i}=['M0/' data(i).name];
                end
                cd([subjectfolder]);
                

    %% calculate CBF
    display('calculating CBF map') 
    
    % !!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!
    % please consult LOFT_CBFquantify stript or type 'help LOFT_CBFquantify' to select appropriate parameters for quantification
    LOFT_CBFquantify(data1, Mzero{1}, 1, 3, 1, 1, 0, 2, 1.5, 0, 0.1, 1, 1, 1)  %change slicetiming for different sequences
    % LOFT_CBFquantify(images, Mzero, FieldStrength, ASLType, FirstimageType, SubtractionOrder, SubtractionType, PLD, LabelTime, Slicetime, threshold, optionPCA, ASLscaling, M0scaling)
   
    cd(subjectfolder)
    display('DONE')
    

   


  
