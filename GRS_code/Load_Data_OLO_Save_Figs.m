clear all
close all

%% flags for what to run
whois='IF';
% subject and hemisphere list
hemi_list={'lh'}; h=1; % eventually do both hemispheres
sub_anat='mmSR20151113'; % anatomy data

nStim=72;
session_dir={['Videostim.sm0.', hemi_list{h}], ['Audiostim.sm0.', hemi_list{h}]};
% the first file is used as the sort index

%% set up directories
if 0
    code_dir='/project_space/Finelab/GRS_sightrecovery';
    anatomies_dir='/project_space/Finelab/GRS_sightrecovery/GRS_data/anatomies';
elseif strcmp(whois, 'IF')
    code_dir='C:\Users\Ione Fine\Documents\Work\Science\Projects\Tristram Savage\GRS_sightrecovery\GRS_code';
    anatomies_dir='C:\Users\Ione Fine\Documents\Work\Science\Projects\Tristram Savage\GRS_sightrecovery\GRS_data\anatomies';
    functional_dir='C:\Users\Ione Fine\Documents\Work\Science\Projects\Tristram Savage\GRS_sightrecovery\GRS_data\functional\mmSRGRS02\bold'
%    addpath(genpath('C:\Users\Ione Fine\Documents\Work\Science\Programming Utilities'))
end
addpath(genpath(code_dir))

%% vert2vox (from vertex space to full voxel space)
% pull in the label file and use it to subsample vertices based on them
% belonging to a single voxel
cd(anatomies_dir)
cd([sub_anat,filesep, 'label'])
disp('loading the label file')
[vertex,x, y, z]=readCortexLabels([hemi_list{h}, '.TS.VTC_first.label']);
% the size of vertex, x, y, z is smaller than max(vertex) because some
% vertices aren't labelled.
vertXYZ=[x y z]; % x yz position of each vertex
[vox, unique_vox_ind,restore_vert]=unique(round(vertXYZ), 'rows');
voxXYZ=round(vertXYZ(unique_vox_ind, :)); % x yz position of each voxel
clear jnk1 jnk2 jnk3 x y z

%% load in beta weights for first file in the list (from vertex space to voxel space)
cd([functional_dir, filesep, session_dir{1}]);
if strcmp(whois, 'IF')
    tmp=load('beta');
    beta(1).vox=squeeze(tmp.tmp.vol(1,unique_vox_ind, 1,1:nStim)); % number of unique vox by number of regressors
else
    tmp=MRIread('beta.nii.gz');
    beta(1).vox=squeeze(tmp.vol(1,unique_vox_ind, 1,1:nStim)); 
end

%% recursive reduction (from full to reduced voxel space)
thr=1;
thr_step=.05;
maxsize=2000; % how big a matrix do I want to use for OLO, set to Inf if you want the biggest matrix possible
disp(['matrix size is ', num2str(size(beta(1).vox,1))]);
keepIterating=1;
family=1:length(beta(1).vox);
while keepIterating
    if length(unique(family))<maxsize
        keepIterating=0;
    else% while ccmatrix still going to be too big
        thr=thr-thr_step; % gradually reduce the threshold for merging
        disp(['reducing interation. Correlation threshold = ',num2str(thr)]);
        beta(1).ccmap=corrcoef(beta(1).vox'); % cross correlation map
        family=zeros(size(beta(1).vox,1), 1); % entries are filled in as they join a family
        for t=1:length(beta(1).ccmap) % for each voxel
            if var(beta(1).vox(t, :))==0 % put duff voxels in a bad family
                family(t)=-1000;
            elseif family(t)==0
                same=find(beta(1).ccmap(t,:)>thr);
                family(same)=t;
            end
        end
        disp(['beta now of length ' num2str(length(unique(family)))])
    end
end

uvals=setdiff(unique(family), -1000);
for i=1:length(uvals)
    re_exp(i).ind=find(family==uvals(i));
    beta(1).mn(i, :)=nanmean(squeeze(beta(1).vox(re_exp(i).ind, :)), 1); % create a canonical beta weight profile for each voxel family
end

%% OLO (reduced voxel space)
cc_mat=corrcoef(beta(1).mn');
d = pdist(cc_mat);
tree = linkage(d);
leafOrder=optimalleaforder(tree, d);


%% Re-expansion of vertices (from reduced to full voxel space)
re_exp_index=[];
for r=1:length(leafOrder)
    re_exp_index=[re_exp_index; re_exp(leafOrder(r)).ind];
end

%% image the first OLO matrix (full voxel space)
subplot(1,length(session_dir),1)
imagesc(beta(1).ccmap(re_exp_index, re_exp_index));
title(session_dir{1})

%% load and image comparison beta files (full voxel space)
for n=2:length(session_dir)
    cd([functional_dir, filesep,session_dir{n}]);
    if strcmp(whois, 'IF')
        tmp=load('beta');
        beta(n).vox=squeeze(tmp.tmp.vol(1,unique_vox_ind, 1,1:nStim));
    else
        tmp=MRIread('beta.nii.gz');
        beta(n).vox=squeeze(tmp.vol(1,unique_vox_ind, 1,1:nStim));
    end
    beta(n).ccmap=corrcoef(beta(n).vox');
    subplot(1,length(session_dir),n)
    imagesc(beta(n).ccmap(re_exp_index, re_exp_index));
    title(session_dir{n})
end

%% overlap the OLO matrices (full voxel space)
figure(2)
allmat=.5*ones(length(re_exp_index), length(re_exp_index), 3);
for n=1:min([length(session_dir), 3]) % overlap the first three beta weight sets on the list
    allmat(:, :,n)=beta(n).ccmap(re_exp_index, re_exp_index);
end
imagesc(allmat);axis equal; axis square; drawnow; 

%% manually select clusters (full voxel space)
nc=input('How many clusters do you want to select? ...');
for n=1:nc
    [x,y] = ginput(2);
    cluster(n).vox=re_exp_index(round(x(1):x(2)));
    cluster(n).filename=input('filename for cluster ... ', 's');
end


%% convert to vertices

for n=1:nc
    cluster(n).vertices=[];
    cluster(n).vertXYZ=[];
    for v=1:length(cluster(n).vox)
        xyz=voxXYZ(v,:);
        ind=find(round(vertXYZ(:, 1))==xyz(1) & round(vertXYZ(:, 2))==xyz(2) & round(vertXYZ(:, 3))==xyz(3));
        cluster(n).vertices=cat(1,cluster(n).vertices, ind);
        cluster(n).vertXYZ=cat(1,  cluster(n).vertXYZ,  vertXYZ(ind,:));
      
    end
     cd(anatomies_dir)
        cd([sub_anat,filesep, 'label'])
        writeCortexLabels([hemi_list{h}, '.',cluster(n).filename, '.label'], cluster(n)); % write to a label file
end


