dirs = dir('spikes*.csv');
for datano = 1:length(dirs)
filename = dirs(datano).name;
T = readtable(filename);
disp(['............................................'])
disp(['Analyzing spike data obtained from CALIMA']);
tic
frames = height(T);
rois = width(T)-1;
A = table2array(T);
A(:,1) = [];
disp(['for ',sprintf('%d',rois),' ROIs with a diameter of 50 pixels:'])
    for i = 1:rois
        roi(i) = nnz(A(:,i)==1);
    end
wavesperroi = mean(roi);
freqpersec = wavesperroi/frames;
freqpermin = freqpersec*60;
wavesover = mean(roi(roi>0));
freqpersecover = wavesover/frames;
freqperminover = freqpersecover*60;
toc

disp(['Saving data...'])
waves = compose('%05.3f', wavesperroi);
fsec = compose('%05.3f', freqpersec);
fmin = compose('%05.3f', freqpermin);
waveso0 = compose('%05.3f', wavesover);
fseco0 = compose('%05.3f', freqpersecover);
fmino0 = compose('%05.3f', freqperminover);
results = [waves fsec fmin; waveso0 fseco0 fmino0];
re_mat = array2table(results);
re_mat.Properties.VariableNames(1:3) = {'Average No. waves per ROI','Frequency per second','Frequency per minute'};
name = ['results_' filename];
writetable(re_mat,name)
disp(['Data from ', sprintf('%s',filename), ' succesfully saved as ', sprintf('%s',name),'!'])
disp([' '])
disp([' '])
end

prompt = "> Would you like to summarise your data? y/N: ";
sum = input(prompt,"s");
if sum == 'y' || sum == 'Y'
    disp(['> Do you want to summarise for samples or ROIs?']);
    prompt = "Select s for samples or r for ROIs: ";
    sr = input(prompt, "s");
    if sr == 's' || sr == 'S'
        meanSampleAnalysis();
    elseif sr == 'r' || sr == 'R'
        meanROIAnalysis();
    else
        disp(['Invalid choice. Data will not be summarised.']);
    end
else
    disp(['Data not summarised.']);
    disp(['Data will be summarised manually.']);
end

prompt = "> Do you want to clean your workspace? y/N: ";
txt = input(prompt,"s");
if txt == 'y' || txt == 'Y'
        disp(['Goodbye!'])
        clear
else
    disp(['Goodbye!'])
end

function meanROIAnalysis()
warning('off','all')

tic 
xpro_all = [];
xpro_over = [];
tnf_all = [];
tnf_over = [];
ctrl_all = [];
ctrl_over = [];
Xpro = "_Xpro_";
TNF = "_TNF_";
CTRL = ["_Xpro_","_TNF_"];
resmat = dir('results_spikes*.csv');

% FOR TNF
for i = 1:length(resmat)
sample = readtable(resmat(i).name);
ss = table2array(sample);
freqall = ss(1,3);
freqover = ss(2,3);
    if contains(resmat(i).name,TNF)
        tnf_all(i,:) = freqall;
        tnf_over(i,:) = freqover;
    else
        tnf_all(i,:) = NaN;
        tnf_over(i,:) = NaN;
    end
end

% FOR XPRO
for i = 1:length(resmat)
sample = readtable(resmat(i).name);
ss = table2array(sample);
freqall = ss(1,3);
freqover = ss(2,3);
    if contains(resmat(i).name,Xpro)
        xpro_all(i,:) = freqall;
        xpro_over(i,:) = freqover;
    else
        xpro_all(i,:) = NaN;
        xpro_over(i,:) = NaN;    
    end
end

% FOR CONTROLS
for i = 1:length(resmat)
sample = readtable(resmat(i).name);
ss = table2array(sample);
freqall = ss(1,3);
freqover = ss(2,3);
TF = contains(resmat(i).name,CTRL);
    if TF == 0        
        ctrl_all(i,:) = freqall;
        ctrl_over(i,:) = freqover;    
    else
        ctrl_all(i,:) = NaN;
        ctrl_over(i,:) = NaN;    
    end
end

mean_control = mean(ctrl_all(ctrl_all>0));
mean_control_over0 = mean(ctrl_over(ctrl_over>0));
mean_tnf = mean(tnf_all(tnf_all>0));
mean_tnf_over0 = mean(tnf_over(tnf_over>0));
mean_xpro = mean(xpro_all(xpro_all>0));
mean_xpro_over0 = mean(xpro_over(xpro_over>0));
analysed_data = [mean_control mean_tnf mean_xpro; mean_control_over0 mean_tnf_over0 mean_xpro_over0];
final_data = array2table(analysed_data);
final_data.Properties.VariableNames(1:3) = {'Control','TNF 100ng/ml','Xpro 100ng/ml'};
final_data.Properties.RowNames(1:2) = {'All ROIs', 'Active ROIs'};
n = ['ROI_summarised_data.csv'];
writetable(final_data,n);
disp(['Data were successfully summarised!']);
toc
disp([' ']);
disp(['Your ROI_summarised_data.csv file has been saved on the current directory.'])
disp([' ']);
disp([' ']);
final_data

disp([' ']);
disp([' ']);
disp(['Calculating statistics for active ROIs...'])
Stats = [ctrl_over tnf_over xpro_over];
disp([' ']);
disp([' ']);
h = kstest(Stats);
if h == 0
    disp(['Running ONE-WAY ANOVA:']);
    [p,tbl,stats] = anova1(Stats)
    disp([' ']);
    disp([' ']);
    disp(['Running multiple t-tests:']);
    disp(['']);
    disp(['Comparing controls with TNF treatments...'])
    [h,p,ci,stats] = ttest2(Stats(:,1),Stats(:,2))
    disp(['Comparing controls with Xpro treatments...'])
    [h,p,ci,stats] = ttest2(Stats(:,1),Stats(:,2))
else
    disp(['Running Kruskal-Wallis ANOVA:']);
    [p,tbl,stats] = kruskalwallis(Stats)
    disp([' ']);
    disp([' ']);
    disp(['Running multiple t-tests:']);
    disp(['']);
    ctrl_over = ctrl_over(~isnan(ctrl_over));
    tnf_over = tnf_over(~isnan(tnf_over));
    xpro_over = xpro_over(~isnan(xpro_over));
    disp(['Comparing controls with TNF treatments...'])
    STATS = mwwtest(transpose(ctrl_over),transpose(tnf_over))
    disp(['Comparing controls with Xpro treatments...'])
    STATS = mwwtest(transpose(ctrl_over),transpose(xpro_over))
end 

prompt = "> Would you like to repeat the calculations for all ROIs? y/N: ";
answer = input(prompt,"s");
if answer == 'y' || answer == 'Y'
    disp([' ']);
    disp([' ']);
    disp(['Calculating statistics for all ROIs...'])
    Stats = [ctrl_all tnf_all xpro_all];
    disp([' ']);
    disp([' ']);
    h = kstest(Stats);
    if h == 0
        disp(['Running ONE-WAY ANOVA:']);
        [p,tbl,stats] = anova1(Stats)
        disp([' ']);
        disp([' ']);
        disp(['Running multiple t-tests:']);
        disp(['']);
        disp(['Comparing controls with TNF treatments...'])
        [h,p,ci,stats] = ttest2(Stats(:,1),Stats(:,2))
        disp(['Comparing controls with Xpro treatments...'])
        [h,p,ci,stats] = ttest2(Stats(:,1),Stats(:,2))
    else
        disp(['Running Kruskal-Wallis ANOVA:']);
        [p,tbl,stats] = kruskalwallis(Stats)
        disp([' ']);
        disp([' ']);
        disp(['Running multiple t-tests:']);
        disp(['']);
        ctrl_over = ctrl_over(~isnan(ctrl_over));
        tnf_over = tnf_over(~isnan(tnf_over));
        xpro_over = xpro_over(~isnan(xpro_over));
        disp(['Comparing controls with TNF treatments...'])
        STATS = mwwtest(transpose(ctrl_over),transpose(tnf_over))
        disp(['Comparing controls with Xpro treatments...'])
        STATS = mwwtest(transpose(ctrl_over),transpose(xpro_over))
    end 
else
    ctrl_all = ctrl_all(~isnan(ctrl_all));
    tnf_all = tnf_all(~isnan(tnf_all));
    xpro_all = xpro_all(~isnan(xpro_all));
end

prompt = "> Do you want to save the extracted data? y/N: ";
answer = input(prompt,'s');
if answer == 'y' || answer == 'Y'
    save xpro_activeROIs xpro_over
    save xpro_allROIs xpro_all
    save control_activeROIs ctrl_over
    save control_allROIs ctrl_all
    save tnf_activeROIs tnf_over
    save tnf_allROIs tnf_all
    disp([' ']);
    disp(['Group data saved on the current directory!']);
    disp([' ']);
else 
    disp(['Group data discarded.']);
    disp([' ']);
end

end

function meanSampleAnalysis()
warning('off','all')
disp(['Summarising per retina sample...']);
samples = dir('results_spikes_*.csv');
wt = [];
d1 = [];
rif = [];
thy = [];
wttnf = [];
d1tnf = [];
riftnf = [];
thytnf = [];
wtxpro = [];
d1xpro = [];
rifxpro = [];
thyxpro = [];

disp(['> Do you want to analyse for all or only active ROIs?']);
prompt = "Select L for all or A for only active ROIs: ";
response = input(prompt,'s');
if response == 'l' || response == 'L'
    var = 1;
elseif response == 'a' || response == 'A'
    var = 2;
else 
    disp(['Input not found. Data will be analysed for only active ROIs.']);
    disp([' ']);
    var = 2;
end

tic 
for nom = 1:length(samples)
    T = readtable(samples(nom).name);
    T = table2array(T);
    for nos = 1:length(samples)
        
        tfwt = "WT_" + [sprintf('%02d' ,nos)];
        tfd1 = "D1tdT_" + [sprintf('%02d' ,nos)];
        tfrif = "RIF_" + [sprintf('%02d' ,nos)];
        tfthy = "Thy_" + [sprintf('%02d' ,nos)];
        tfwttnf = "WT_TNF_" + [sprintf('%02d' ,nos)];
        tfd1tnf = "D1tdT_TNF_" + [sprintf('%02d' ,nos)];
        tfriftnf = "RIF_TNF_" + [sprintf('%02d' ,nos)];
        tfthytnf = "Thy_TNF_" + [sprintf('%02d' ,nos)];
        tfwtxpro = "WT_Xpro_" + [sprintf('%02d' ,nos)];
        tfd1xpro = "D1tdT_Xpro_" + [sprintf('%02d' ,nos)];
        tfrifxpro = "RIF_Xpro_" + [sprintf('%02d' ,nos)];
        tfthyxpro = "Thy_Xpro_" + [sprintf('%02d' ,nos)];
        
        if contains(samples(nom).name,tfwt) == 1
            wt(nom,nos) = T(var,3);
        elseif contains(samples(nom).name,tfd1) == 1
            d1(nom,nos) = T(var,3);
        elseif contains(samples(nom).name,tfrif) == 1
            rif(nom,nos) = T(var,3);
        elseif contains(samples(nom).name,tfthy) == 1
            thy(nom,nos) = T(var,3);
        elseif contains(samples(nom).name,tfwttnf) == 1
            wttnf(nom,nos) = T(var,3);
        elseif contains(samples(nom).name,tfd1tnf) == 1
            d1tnf(nom,nos) = T(var,3);
        elseif contains(samples(nom).name,tfriftnf) == 1
            riftnf(nom,nos) = T(var,3);
        elseif contains(samples(nom).name,tfthytnf) == 1
            thytnf(nom,nos) = T(var,3);
        elseif contains(samples(nom).name,tfwtxpro) == 1
            wtxpro(nom,nos) = T(var,3);
        elseif contains(samples(nom).name,tfd1xpro) == 1
            d1xpro(nom,nos) = T(var,3);
        elseif contains(samples(nom).name,tfrifxpro) == 1
            rifxpro(nom,nos) = T(var,3);
        elseif contains(samples(nom).name,tfthyxpro) == 1
            thyxpro(nom,nos) = T(var,3);
        end
    end
end

disp(['Data successfully summarised!']);
disp([' ']);
disp(['Calculating statistics...'])
disp([' ']);
disp(['Calculating average frequency (waves/min) of all the controls:'])

%For controls
WTs = [];
D1s = [];
RIFs = [];
THYs = [];

for s = 1:width(wt)
    meanwt(:,s) = mean(nonzeros(wt(:,s)));
    WTs = meanwt';
end
for s = 1:width(d1)
    meand1(:,s) = mean(nonzeros(d1(:,s)));
    D1s = meand1';
end
for s = 1:width(rif)
    meanrif(:,s) = mean(nonzeros(rif(:,s)));
    RIFs = meanrif';
end
for s = 1:width(thy)
    meanthy(:,s) = mean(nonzeros(thy(:,s)));
    THYs = meanthy';
end
controls = {WTs D1s RIFs THYs};
marray = max(cellfun(@numel, controls));
WT = nan(marray,1);
D = nan(marray,1);
RIF = nan(marray,1);
THY = nan(marray,1);
WT(1:length(WTs),1) = WTs;
D(1:length(D1s),1) = D1s;
RIF(1:length(RIFs),1) = RIFs;
THY(1:length(THYs),1) = THYs;
controls = cat(1,WT,D,RIF,THY);
controls = controls(~isnan(controls))
clear marray meanwt meand1 meanrif meanthy WT WTs D1s D RIFs RIF THYs THY

disp([' ']);
disp(['Calculating average frequency (waves/min) of all the TNF-treated retinas:'])
%For TNF treatments
WTs = [];
D1s = [];
RIFs = [];
THYs = [];

for s = 1:width(wttnf)
    meanwt(:,s) = mean(nonzeros(wttnf(:,s)));
    WTs = meanwt';
end
for s = 1:width(d1tnf)
    meand1(:,s) = mean(nonzeros(d1tnf(:,s)));
    D1s = meand1';
end
for s = 1:width(riftnf)
    meanrif(:,s) = mean(nonzeros(riftnf(:,s)));
    RIFs = meanrif';
end
for s = 1:width(thytnf)
    meanthy(:,s) = mean(nonzeros(thytnf(:,s)));
    THYs = meanthy';
end
TNFs = {WTs D1s RIFs THYs};
marray = max(cellfun(@numel, TNFs));
WT = nan(marray,1);
D = nan(marray,1);
RIF = nan(marray,1);
THY = nan(marray,1);
WT(1:length(WTs),1) = WTs;
D(1:length(D1s),1) = D1s;
RIF(1:length(RIFs),1) = RIFs;
THY(1:length(THYs),1) = THYs;
TNFs = cat(1,WT,D,RIF,THY);
TNFs = TNFs(~isnan(TNFs))
clear marray meanwt meand1 meanrif meanthy WT WTs D1s D RIFs RIF THYs THY

disp([' ']);
disp(['Calculating average frequency (waves/min) of all the TNF-treated retinas:'])
%For Xpro treatments
WTs = [];
D1s = [];
RIFs = [];
THYs = [];

for s = 1:width(wtxpro)
    meanwt(:,s) = mean(nonzeros(wtxpro(:,s)));
    WTs = meanwt';
end
for s = 1:width(d1xpro)  
    meand1(:,s) = mean(nonzeros(d1xpro(:,s)));
    D1s = meand1';
end
for s = 1:width(rifxpro)   
    meanrif(:,s) = mean(nonzeros(rifxpro(:,s)));
    RIFs = meanrif';
end
for s = 1:width(thyxpro)
    meanthy(:,s) = mean(nonzeros(thyxpro(:,s)));
    THYs = meanthy';
end
Xpro = {WTs D1s RIFs THYs};
marray = max(cellfun(@numel, Xpro));
WT = nan(marray,1);
D = nan(marray,1);
RIF = nan(marray,1);
THY = nan(marray,1);
WT(1:length(WTs),1) = WTs;
D(1:length(D1s),1) = D1s;
RIF(1:length(RIFs),1) = RIFs;
THY(1:length(THYs),1) = THYs;
Xpro = cat(1,WT,D,RIF,THY);
Xpro = Xpro(~isnan(Xpro))
clear marray meanwt meand1 meanrif meanthy WT WTs D1s D RIFs RIF THYs THY

data = {controls TNFs Xpro};
marray = max(cellfun(@numel, data));
CTRL = nan(marray,1);
TNF = nan(marray,1);
XPRO = nan(marray,1);
CTRL(1:length(controls),1) = controls;
TNF(1:length(TNFs),1) = TNFs;
XPRO(1:length(Xpro),1) = Xpro;
Stats = [CTRL TNF XPRO];

mean_control = mean(controls);
mean_tnf = mean(TNFs);
mean_xpro = mean(Xpro);
analysed_data = [mean_control mean_tnf mean_xpro];
final_data = array2table(analysed_data);
final_data.Properties.VariableNames(1:3) = {'Control','TNF 100ng/ml','Xpro 100ng/ml'};
n = ['sample_summarised_data.csv'];
writetable(final_data,n);
disp(['Data were successfully summarised!']);
toc
disp([' ']);
disp(['Your sample_summarised_data.csv file has been saved on the current directory.'])
disp([' ']);
disp([' ']);
final_data

h = kstest(Stats);
if h == 0
    disp(['Running ONE-WAY ANOVA:']);
    [p,tbl,stats] = anova1(Stats)
    disp([' ']);
    disp([' ']);
    disp(['Running multiple t-tests:']);
    disp(['']);
    disp(['Comparing controls with TNF treatments...'])
    [h,p,ci,stats] = ttest2(Stats(:,1),Stats(:,2))
    disp(['Comparing controls with Xpro treatments...'])
    [h,p,ci,stats] = ttest2(Stats(:,1),Stats(:,2))
else
    disp(['Running Kruskal-Wallis ANOVA:']);
    [p,tbl,stats] = kruskalwallis(Stats)
    disp([' ']);
    disp([' ']);
    disp(['Running multiple t-tests:']);
    disp(['']);
    disp(['Comparing controls with TNF treatments...'])
    STATS = mwwtest(transpose(controls),transpose(TNFs))
    disp(['Comparing controls with Xpro treatments...'])
    STATS = mwwtest(transpose(controls),transpose(Xpro))
end 

prompt = "> Do you want to save the grouped data? y/N: ";
answer = input(prompt,'s');
if answer == 'y' || answer == 'Y'
    save controls controls;
    save TNFs TNFs;
    save Xpro Xpro;
    disp([' ']);
    disp(['Group data saved on the current directory!']);
    disp([' ']);
else 
    disp(['Group data discarded.']);
    disp([' ']);
end

end