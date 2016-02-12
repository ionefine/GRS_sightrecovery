clear all
close all

whois='IF';

if 0
    code_dir='/project_space/Finelab/GRS_sightrecovery';
    anatomies_dir='/project_space/Finelab/GRS_sightrecovery/GRS_data/anatomies';
elseif strcmp(whois, 'IF')
    code_dir='C:\Users\Ione Fine\Documents\Work\Science\Projects\Tristram Savage\GRS_sightrecovery\GRS_code';
    anatomies_dir='C:\Users\Ione Fine\Documents\Work\Science\Projects\Tristram Savage\GRS_sightrecovery\GRS_data\anatomies';
    functional_dir='C:\Users\Ione Fine\Documents\Work\Science\Projects\Tristram Savage\GRS_sightrecovery\GRS_data\functional\mmSRGRS02\bold'
    addpath(genpath('C:\Users\Ione Fine\Documents\Work\Science\Programming Utilities'))
end
addpath(genpath(code_dir))

% subject and hemisphere list
hemi_list={'lh', 'rh'}; % eventually do both hemispheres
sub_anat='mmSR20151113'; % anatomy data

h=1;
session_dir={['Videostim.sm0.', hemi_list{h}], ['Audiostim.sm0.', hemi_list{h}]}; 
% the files we are going to look at using that sortfile index, using the
% first as the sorter

% pull in the label file and use it to subsample vertices based on them
% belonging to a single voxel
cd(anatomies_dir)
cd([sub_anat,filesep, 'label'])
disp('loading the label file')
[vertex,jnk1,x,jnk2, y,jnk3, z]=readCortexLabels([hemi_list{h}, '.TS.VTC_first.label']);
[vox, unique_vox_ind,restore_vert]=unique(round([x y z]), 'rows');

clear jnk1 jnk2 jnk3 x y z
% eureka! the size of vertex, x, y, z is smaller than the max vertex. That
% is because some vertices aren't labelled. But we probably don't care
% about those, so wtf

% configuration details
nStim=72;

% load in beta weights
cd(functional_dir);
cd(session_dir{1});
if strcmp(whois, 'IF')
    tmp=load('beta');
    beta_vox=squeeze(tmp.tmp.vol(1,unique_vox_ind, 1,1:nStim)); % number of unique vox by number of regressors
else
tmp=MRIread('beta.nii.gz');
beta_vox=squeeze(tmp.vol(1,unique_vox_ind, 1,1:nStim)); % number of unique vox by number of regressors
end

% so the default max array size for regression is limited, if beta is
% too big we are going to recursively reduce the number of voxels by until it's small
% enough to do corrcoef on.
% chunksize=10000; % the size of matrix at which Matlab starts to wedge 10,000-20,000 is good
% nchunks=ceil(size(beta_vox, 1)./chunksize);
%newbetasize=size(beta,1); % how far have we shrunk it?
thr=1;
thr_step=.05;
maxsize=2000; % how big a matrix do I want to use for OLO, set to Inf if you want the biggest matrix possible, remember olo will be slow
disp(['matrix size is ', num2str(size(beta_vox,1))]);
keepIterating=1;
family=1:length(beta_vox);
while keepIterating
    if length(unique(family))<maxsize % still works if maxsize inf because doesn't execute if try fails
        keepIterating=0;
    else% while ccmatrix still going to be too big
        thr=thr-thr_step; % gradually reduce the threshold for merging
        disp(['reducing interation. Correlation threshold = ',num2str(thr)])     
        beta_ccmap=corrcoef(beta_vox'); %cross correlation map
        family=zeros(size(beta_vox,1), 1); %Remake family matrix; entries are zeroed as their value is collapsed into other entries.
        for t=1:length(beta_ccmap) % for each voxel
            if var(beta_vox(t, :))==0
                family(t)=-1000;
            elseif family(t)==0  %family is vector os length beta+vox, if a entry is 0 it Nans the whole row of the
                same=find(beta_ccmap(t,:)>thr);
                family(t)=t;
                family(same)=t;
            end
        end
        disp(['beta now of length ' num2str(length(unique(family)))])
    end % did reduction work
end % did initial check show that reduction required

uvals=setdiff(unique(family), -1000);
for i=1:length(uvals)
    re_exp(i).ind=find(family==uvals(i));
    mnbeta(i, :)=nanmean(squeeze(beta_vox(re_exp(i).ind, :)), 1); % create a canonical beta weight profile for each voxel
end

%%%Olo stuff
cc_mat=corrcoef(mnbeta');
d = pdist(cc_mat);
tree = linkage(d);
leafOrder=optimalleaforder(tree, d);
%%%

%%%Re-expansion
beta_vox_reexp=[];
re_exp_index=[];

for r=1:length(leafOrder)
    %     beta_vox_reexp=[beta_vox_reexp;beta_vox(beta_coords3{leafOrder(r)},:)];
    %re_exp_index=[Re_exp_index; beta_coords3{leafOrder(r)}];
    re_exp_index=[re_exp_index; re_exp(leafOrder(r)).ind];
end
[x,y] = ginput(2);
clusterx=

beta_vox_reexp=beta_vox(re_exp_index,:);
sorted_full_ccmap=(corrcoef(beta_vox_reexp'));
return
% pull out a cluster
[x,y] = ginput(2)
beta_voxNZ=beta_vox;
[z, y]=find(beta_voxNZ(:,1)==0);
beta_voxNZ(z,:)=[]; %removes voxels with values = 0xels with values = 0
beta_ccmap=corrcoef(beta_voxNZ');




if SortCon==1
    
    audio_dir=(['Audiostim.sm0.', hemi_list{h}]);
    cd(functional_dir);
    cd(audio_dir);
    aud_struct=MRIread('beta.nii.gz');
    beta_vox_aud=squeeze(aud_struct.vol(1,unique_vox_ind, 1,1:nStim)); % number of unique vox by number of regressors
    beta_vox_aud_sorted=(beta_vox_aud(Re_exp_index,:));
    beta_vox_aud(z,:)=[];
    cc_map_aud=corrcoef(beta_vox_aud');
    cc_map_aud_sorted=corrcoef(beta_vox_aud_sorted');
    
    %         fig=figure(1);
    %         subplot(2,2,1)
    %         imagesc(beta_ccmap)
    %         subplot(2,2,2)
    %         imagesc(sorted_full_ccmap)
    %         subplot(2,2,3)
    %         imagesc(cc_map_aud)
    %         subplot(2,2,4)
    %         imagesc(cc_map_aud_sorted)
    %         saveas(fig,['AudStim_by_VidStim' hemi_list{h} '.png'])
    
elseif SortCon==2
    
    cd(functional_dir);
    video_dir=(['Videostim.sm0.', hemi_list{h}]);
    
    cd(functional_dir);
    cd(video_dir);
    vid_struct=MRIread('beta.nii.gz');
    beta_vox_vid=squeeze(vid_struct.vol(1,unique_vox_ind, 1,1:nStim)); % number of unique vox by number of regressors
    beta_vox_vid_sorted=(beta_vox_vid(Re_exp_index,:));
    beta_vox_vid(z,:)=[];
    cc_map_vid=corrcoef(beta_vox_vid');
    cc_map_vid_sorted=corrcoef(beta_vox_vid_sorted');
    %
    %         fig=figure(1)
    %         subplot(2,2,1)
    %         imagesc(beta_ccmap)
    %         subplot(2,2,2)
    %         imagesc(sorted_full_ccmap)
    %         subplot(2,2,3)
    %         imagesc(cc_map_vid)
    %         subplot(2,2,4)
    %         imagesc(cc_map_vid_sorted)
    %         saveas(fig,['VidStim_by_AudStim' hemi_list{h} '.png'])
    
end

%%%
%%%
%
% end


%%%imstuff done


%
% figure(1)
% subplot(2,2,1)
% imagesc(sorted_full_ccmap);
% subplot(2,2,2)
% imagesc(sorted_full_ccmap_test)
%
%
% % figure(1)
% % subplot(2,2,1)
% % imagesc(cc_mat2);
% % subplot(2,2,2)
% % imagesc(cc_mat2(leafOrder, leafOrder))
%
% subplot(2,2,3)
% imagesc(beta_ccmap)
% subplot(2,2,4)
% imagesc(sorted_full_ccmap)

%%%imstuff done




%%%%%%%%%%%%Below here are bits of test code

% Vert sum checker {for checkin' contents of large voxel matricies}

% zed=[];
% for r=1:length(beta_coords2)
%
% zed=zed+sum(beta_coords2{r});
% end
% zed






%%%duplicate finder
%
% for d=1:8000
%     zed2=[];
%     for z=1:length(beta_coords)
%         if max(beta_coords{z}(:)==d)==1
%             zed2=[zed2; z];
%         end
%     end
%
%     if length(zed2)>1
%         d
%     end
%
% end

%





