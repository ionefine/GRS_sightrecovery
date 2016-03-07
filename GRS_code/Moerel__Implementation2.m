
clear all
% close all

%% TODO
% right hemisphere ventral

% Marty's, Mello
% implement a pwelch
% read the Moerel paper
% run the whole brain, and save all the interesting voxels
%implement NSL bandpass filters ala Morele {http://www.isr.umd.edu/Labs/NSL/Software.htm;} in place of bins
    %moerel 2012 :128 --> 40 bins uniformly distributed on a logarithmic frequency axis; center frequency ranges from 186 to 6817 Hz


%% flags for what to run
whois='IF';
% subject and hemisphere list
hemi_list={'lh'}; h=1; % eventually do both hemispheres
sub_anat='mmSR20151113'; % anatomy data
sub_fun='msSrGRS01'
nStim=72;
session_dir={['Audiostim.sm0.', hemi_list{h}]};
roifilename={ 'lh.Aud.partial.label', 'lh.aud2.too_big.label', 'lh.MT.label', 'rh.MT.label', 'rh.TS.VTC_first.label' , 'lh.TS.VTC_first.label' };%'.V2.label', ,'.TS.VTC_first.label','.V1.label',, '.BA45.label'
% the first file is used as the sort index
writeroifilename={'.V1_CA.label', '.V2_CA.label', '.TS.VTC_first_CA.label', '.MT_CA.label', '.BA6_CA.label',  '.BA45_CA.label'};
%% set up directories

% %%%%Normal dirs
% % if 0
% %     code_dir='/project_space/Finelab/GRS_sightrecovery';
% %     anatomies_dir='/project_space/Finelab/GRS_sightrecovery/GRS_data/anatomies';
% % elseif strcmp(whois, 'IF')
% %     code_dir='C:\Users\Ione Fine\Documents\Work\Science\Projects\Tristram Savage\GRS_sightrecovery\GRS_code';
% %     anatomies_dir='C:\Users\Ione Fine\Documents\Work\Science\Projects\Tristram Savage\GRS_sightrecovery\GRS_data\anatomies';
% %     functional_dir='C:\Users\Ione Fine\Documents\Work\Science\Projects\Tristram Savage\GRS_sightrecovery\GRS_data\functional\mmSRGRS02\bold';
% %     audio_dir='C:\Users\Ione Fine\Documents\Work\Science\Projects\Tristram Savage\GRS_sightrecovery\GRS_stimfiles';
% % end

% Home deskrop dirs
%%Normal dirs

    code_dir='Z:\Dropbox\Git\GRS_sightrecovery\GRS_code';
    anatomies_dir='Z:\Dropbox\Git\GRS_sightrecovery\GRS_data\anatomies';
    functional_dir='Z:\Dropbox\Git\GRS_sightrecovery\GRS_data\functional';
    audio_dir='Z:\Dropbox\Git\GRS_sightrecovery\GRS_stimfiles';



addpath(genpath(code_dir))

%% read in and analyse audio files




% 

cd(audio_dir)
audioFilenames=Stim_List_As_Presented_Nov_2016;

para=[8.0000    8.0000   -2.0000   -0.4150]; %taken as default values from nsl toolbox
%may need to be changed is we customize much
rFs=12000;


for a=1:length(audioFilenames)  
    [y,Fs] = audioread(audioFilenames{a});% ; y=y(:, 1);
    yRsmp=resample(y,rFs,Fs);
    ySpect3=wav2aud(yRsmp',para,'fun');
    A(a,:)=mean(ySpect3);
    
    
end



% 
% 
% cd(audio_dir)
% audioFilenames=Stim_List_As_Presented_Nov_2016;
% freq_bins=round(exp(linspace(log(500), log(5000), 8)));
% freq_bins=[0,freq_bins, inf];

% 
% % find the power in each frequency bin
% for a=1:length(audioFilenames)  
%     [y,Fs] = audioread(audioFilenames{a}) ; y=y(:, 1);
%     maxt = length(y)/Fs;
%     t = linspace(0,maxt,length(y));
%     Y = complex2real2(fft(y),t);
%     for f=1:length(freq_bins)-1
%         ind=find(Y.freq>freq_bins(f) & Y.freq<freq_bins(f+1));
%         A(a, f)=mean(Y.amp(ind));
%     end
% end
% 
% 
% 



A=A(1:nStim,:);
figure
imagesc(A)
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
    [roi(r).vertex,x, y, z]=readCortexLabels([roifilename{r}]);
    clear vox voxb
    for rep=1:2
        AM=A;
        for v=1:length(roi(r).vertex)
            if var(beta(roi(r).vertex(v)+1,:))>0
                B=scaleif(beta(roi(r).vertex(v)+1,:), -1, 1);
                if rep>1
                    AM=AM(randperm(size(AM,1)), :);
                end   
                W=pinv(AM)*B';
                P=AM*W; % predicted betas
                err=sum(sqrt((P'-B).^2));
                vox(rep,v).B=B; vox(rep,v).W=W; vox(rep,v).P=P; vox(rep,v).err=err;
            else
                vox(rep,v).B=NaN; vox(rep,v).W=NaN; vox(rep,v).P=NaN; vox(rep,v).err=NaN;    
            end
        end
    end
    subplot(2, ceil(length(roifilename)/2), r)
    
    statLims=prctile([vox(2:end,:).err], [1, 5]); % do any voxels have fit values significantly higher or lower than expected?
    [n, bins]=hist([vox(1,:).err],5:40);n=n./sum(n);
    [nb, bins]=hist([vox(2:end,:).err],5:40);nb=nb./sum(nb);
    plot(bins, n, 'r-', 'LineWidth', 2);hold on;
    plot(bins, nb, 'k-', 'LineWidth', 1);
    plot([statLims(1) statLims(1)], [0 max(n)*1.1])
    set(gca, 'YLim', [0 max(n)*1.1])
    title(roifilename{r})
    drawnow
%     ind=find([vox(1,:).err]<statLims(1));
%     W=[vox(ind).W]';
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
