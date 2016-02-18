clear all
close all

%% flags for what to run
whois='IF';
% subject and hemisphere list
hemi_list={'rh'}; h=1; % eventually do both hemispheres
sub_anat='mmSR20151113'; % anatomy data

nStim=72;
session_dir={['Audiostim.sm0.', hemi_list{h}]};
roifilename={'.V1.label', '.V2.label',  '.MT.label', '.BA6.label', '.BA45.label'}%'.TS.VTC_first.label',
% the first file is used as the sort index
writeroifilename={'.V1_CA.label', '.V2_CA.label', '.TS.VTC_first_CA.label', '.MT_CA.label', '.BA6_CA.label',  '.BA45_CA.label'};
%% set up directories
if 0
    code_dir='/project_space/Finelab/GRS_sightrecovery';
    anatomies_dir='/project_space/Finelab/GRS_sightrecovery/GRS_data/anatomies';
elseif strcmp(whois, 'IF')
    code_dir='C:\Users\Ione Fine\Documents\Work\Science\Projects\Tristram Savage\GRS_sightrecovery\GRS_code';
    anatomies_dir='C:\Users\Ione Fine\Documents\Work\Science\Projects\Tristram Savage\GRS_sightrecovery\GRS_data\anatomies';
    functional_dir='C:\Users\Ione Fine\Documents\Work\Science\Projects\Tristram Savage\GRS_sightrecovery\GRS_data\functional\mmSRGRS02\bold';
    audio_dir='C:\Users\Ione Fine\Documents\Work\Science\Projects\Tristram Savage\GRS_sightrecovery\GRS_stimfiles';
end
addpath(genpath(code_dir))

%% read in and analyse audio files
cd(audio_dir)
audioFilenames=Stim_List_As_Presented_Nov_2016;
freq_bins=round(exp(linspace(log(50), log(12000), 16)));
freq_bins=[0,freq_bins, inf];
% find the power in each frequency bin
for a=1:length(audioFilenames)
    [y,Fs] = audioread(audioFilenames{a}) ; y=y(:, 1);
    maxt = length(y)/Fs;   
    t = linspace(0,maxt,length(y));
    Y = complex2real(fft(y),t);
    for f=1:length(freq_bins)-1
        ind=find(Y.freq>freq_bins(f) & Y.freq<freq_bins(f+1));
        A(a, f)=mean(Y.amp(ind));
    end
end

A=A(:, 1:end-1); %A is nStim  x nfrequency bands
for a=1:size(A, 2) % for each frequency band
    A(:, a)=scaleif(A(:, a)); % make all stimuli equally loud
end
  
% % the bootstrapped audio stimuli, scrambling over frequencies, this works
% too well
% for f=1:size(A, 2)
%     Ab(:, f)=A(randperm(size(A,1)), f); % for each freq bin, shuffle across stimuli
% end
% the bootstrapped audio stimuli, scrambling over stimuli


% for f=1:size(A, 2)
%     Ab(:, )=A(s, randperm(size(A,2))); % for each freq bin, shuffle across stimuli
% end
Ab=A(randperm(size(A,1)), :);
A=A(1:nStim,:);Ab=Ab(1:nStim,:);

%% load in the beta weights
cd([functional_dir, filesep, session_dir{1}]);
if strcmp(whois, 'IF')
    tmp=load('beta');
    beta=squeeze(tmp.tmp.vol(1,:, 1,1:nStim)); % number of unique vox by number of regressors
else
    tmp=MRIread('beta.nii.gz');
    beta=squeeze(tmp.vol(1,:, 1,1:nStim)); %   nvoxels x nStim
end

for r=1:length(roifilename)
    % load in the ROIs
    disp(['loading ', roifilename{r}])
    cd(anatomies_dir)
    cd([sub_anat,filesep, 'label'])
    [roi(r).vertex,x, y, z]=readCortexLabels([hemi_list{h}, roifilename{r}]);
    clear vox voxb
    
    for v=1:length(roi(r).vertex)
        if mod(v, 1000)==0
            disp(['processing voxel ', num2str(v), ' out of ', num2str(length(roi(r).vertex))])
        end
        if var(beta(roi(r).vertex(v)+1,:))>0
            Ab=A(randperm(size(A,1)), :);
        %    vox(v).amp=mean(beta(roi(r).vertex(v)+1,:));
            B=scaleif(beta(roi(r).vertex(v)+1,:));
            W=pinv(A)*B';  Wb=pinv(Ab)*B'; % find the weights for each frequency bin
            P=A*W; Pb=Ab*Wb;% predicted betas
            err=sum(sqrt((P'-B).^2)); errb=sum(sqrt((Pb'-B).^2));
            vox(v).B=B; vox(v).W=W; vox(v).P=P; vox(v).err=err;
            voxb(v).B=B; voxb(v).W=Wb; voxb(v).P=Pb; voxb(v).err=errb;
        
        else
            vox(v).B=NaN; vox(v).W=NaN; vox(v).P=NaN; vox(v).err=NaN;
            voxb(v).B=NaN; voxb(v).W=NaN; voxb(v).P=NaN; voxb(v).err=NaN;
        end   
    end
    subplot(2, ceil(length(roifilename)/2), r)
    
    statLims=prctile([voxb(:).err], [1, 99]); % do any voxels have fit values significantly higher or lower than expected?
    n=hist([vox(:).err],5:30);
    disp(sum(n))
    nb=hist([voxb(:).err],5:30);
    plot(5:30, n, 'r-', 'LineWidth', 2);hold on;
    plot(5:30, nb, 'k-', 'LineWidth', 1);
    plot([statLims(1) statLims(1)], [0 max(n)*1.1])
    set(gca, 'YLim', [0 max(n)*1.1])
%    text(15, max(n)*.7, num2str(r2p(mean([vox(:).amp]), 2)))
    title(roifilename{r})
    drawnow
    ind=find([vox(:).err]<statLims(1));
    W=[vox(ind).W]';
   return
%     
%     imagesc(W)
%     T = clusterdata(W,10.5); 
%     
%     for p=1:10
%         plot(freq_bins(2:end-1),mean(W(find(T==p), :), 1)); hold on
%     end
%     return
    
end


%     %% convert to vertices
%     cd(anatomies_dir)
%     cd([sub_anat,filesep, 'label'])
%     writeCortexLabels([hemi_list{h}, '.',writeroifilename, '.label'], cluster(n)); % write to a label file

