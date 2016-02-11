% clear all

code_dir='/project_space/Finelab/GRS_sightrecovery';
    anatomies_dir='/project_space/Finelab/GRS_sightrecovery/GRS_data/anatomies';
    addpath(genpath(code_dir))
    SortCon=1;
    
    h=1;
% for h=1:2
   clearvars -except h SortCon anatomies_dir code_dir savetest
    close all
    
    % general directory level above subjects
 
    
    
    % subject and hemisphere list
    hemi_list={'lh', 'rh'}; % eventually do both hemispheres
    sub_anat='mmSR20151113'; % anatomy data
    
    % pull in the label file and use it to subsample vertices based on them
    % belonging to a single voxel
    cd(anatomies_dir)
    cd([sub_anat,filesep, 'label'])
    disp('loading the label file')
    [vertex,jnk1,x,jnk2, y,jnk3, z]=readCortexLabels([hemi_list{h}, '.TS.VTC_first.label']);
    % eureka! the size of vertex, x, y, z is smaller than the max vertex. That
    % is because some vertices aren't labelled. But we probably don't care
    % about those, so wtf
    
    [vox, unique_vox_ind,restore_vert]=unique(round([x y z]), 'rows');
    % provides index into unique vertices, vox(restore_vert) takes you back to
    % vertex space
    
    % configuration details
    nStim=72;
    
    
    %Set sorting subject and targets
    % 1= V1AA,, 1=AVV,moew come later
    if SortCon==1
        functional_dir='/project_space/Finelab/GRS_sightrecovery/GRS_data/functional/mmSRGRS02/bold';
        session_dir={['Videostim.sm0.', hemi_list{h}]};
        
    elseif SortCon==2
        functional_dir='/project_space/Finelab/GRS_sightrecovery/GRS_data/functional/mmSRGRS02/bold';
        session_dir={['Audiostim.sm0.', hemi_list{h}]};
        
        % elseif SortCon==3
        %     functional_dir='/project_space/Finelab/GRS_sightrecovery/GRS_data/functional/mmSRGRS01/bold';
        %     session_dir={['Videostim.sm0.', hemi_list{1}]};
        
        % elseif SortCon==4
        %     functional_dir='/project_space/Finelab/GRS_sightrecovery/GRS_data/functional/mmSRGRS01/bold';
        %     session_dir={['Videostim.sm0.', hemi_list{1}]};
        
    else error('invalid SortCon selection')
    end
    
    
    
    % load in beta weights
    cd(functional_dir);
    cd(session_dir{1});
    tmp=MRIread('beta.nii.gz');
    
    beta_vox=squeeze(tmp.vol(1,unique_vox_ind, 1,1:nStim)); % number of unique vox by number of regressors
    
    
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
    beta_vox2=beta_vox;
    beta_coords=cell(size(beta_vox,1), 1);
    while keepIterating
        if size(beta_vox2, 1)<maxsize % still works if maxsize inf because doesn't execute if try fails
            keepIterating=0;
        else% while ccmatrix still going to be too big
            thr=thr-thr_step; % gradually reduce the threshold for merging
            disp(['reducing interation. Correlation threshold = ',num2str(thr)])
            
            beta_coords=cell(size(beta_vox,1), 1); %cell arrray which stores indecies of voxels which are highly correlated and compressed
            beta_ccmap=corrcoef(beta_vox'); %cross correlation map
            keep=ones(size(beta_vox,1), 1); %Remake keep matrix; entries are zeroed as their value is collapsed into other entries.
            
            for t=1:length(beta_ccmap) % for each voxel
                if keep(t) %keep is vector os length beta+vox, if a entry is 0 it Nans the whole row of the
                    
                    same=find(beta_ccmap(t,:)>thr);
                    
                    % only keep one voxel
                    beta_coords{t}=nonzeros(keep(same).*same'); %This line stores makes
                    %               sure each voxel number is storred in only an singlge cell-- no dobule-dipping.
                    %problem here here is that there is some bias in the
                    %compression based on voxel number.  For instance, if voxel
                    %2 were 99 % corrilated with voxel 4, but 85% correlated
                    %with voxel 1, on the last iteration voxel 4 would be
                    %"grouped" with voxel 1 instead of voxel 4 because it came
                    %first.
                    
                    
                    
                    keep(same)=0;  keep(t)=1;
                    
                end
                %                 keep(t)=0;
            end
            beta_vox2=beta_vox(logical(keep), :);
            disp(['beta now of length ' num2str(size(beta_vox2, 1))])
        end % did reduction work
        
    end % did initial check show that reduction required
    
    beta_coords2=beta_coords(logical(keep)); %kill empty cells {[]}.
    
    
    
    % %%%%% remove nan from cc_mat / zero values from
    % cc_mat=corrcoef(beta_vox2');
    % keep=ones(length(beta_vox2), 1);
    %
    %
    % for i=1:size(cc_mat, 1)
    %     if sum(isnan(cc_mat(i,:)))==length(beta_vox2)
    %         keep(i)=0;
    %
    %     end
    % end
    %
    % beta_vox3=beta_vox2(logical(keep),:);  %removes voxels with values = 0
    % beta_coords3=beta_coords2(logical(keep)); %Removes cells which are [0x1]
    
    
    
    %%%%% remove nan from cc_mat / zero values from alternate / faster
    beta_coords3=beta_coords2;
    beta_vox3=beta_vox2;
    [x, y]=find(beta_vox2(:,1)==0);
    beta_vox3(x,:)=[]; %removes voxels with values = 0
    beta_coords3(x,:)=[];%Removes cells which are [0x1]
    %%%%%%%
    
    
    
    
    
    
    
    %%%Olo stuff
    cc_mat2=corrcoef(beta_vox3');
    d = pdist(cc_mat2);
    tree = linkage(d);
    leafOrder=optimalleaforder(tree, d);
    %%%
    
    
    
    %%%Re-expansion
    beta_vox_reexp=[];
    Re_exp_index=[];
    
    for r=1:length(leafOrder)
        
        %     beta_vox_reexp=[beta_vox_reexp;beta_vox(beta_coords3{leafOrder(r)},:)];
        Re_exp_index=[Re_exp_index; beta_coords3{leafOrder(r)}];
        
        
    end
    
    beta_vox_reexp=beta_vox(Re_exp_index,:);
    sorted_full_ccmap=(corrcoef(beta_vox_reexp'));
   
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





