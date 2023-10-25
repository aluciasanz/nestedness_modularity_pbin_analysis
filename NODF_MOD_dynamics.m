% jan 2023 - ALSANZ Global nestednes and modularity and over time
clear 
clc
%% add path set working directory
PATH='./BiMat-master'; %write the path for BiMat-master
mkdir(strcat(PATH,'/myfolder'))   
GEN_PATH = genpath(PATH);
addpath(GEN_PATH);
cd(strcat(PATH,'/myfolder'));
PATH=cd;
%% Import DATA
% matrix_raw=readtable(strcat(PATH,'\Data\22-06-09-matrix-21d_bis.csv')); % PBIN located in Data 
matrix_raw=readtable(strcat(PATH,'/Data/22-06-09-matrix-21d_bis.csv')); % PBIN located in Data
hostID = matrix_raw.Properties.VariableNames(2:end); % Hosts are columns
phageID = matrix_raw.Phage; % Phage are rows
matrix=table2array(matrix_raw(:,2:end)); % convert table to matrix

%% Create the bipartite object
bp = Bipartite(matrix'); % hosts are rows, phage are columns
bp.row_labels =  hostID;
bp.col_labels = phageID;

%% divide by day
days_host=zeros(1,length(bp.row_labels));

for i=1:length(bp.row_labels)
    id=char(bp.row_labels(i));
    id2=strsplit(id,'_');
    id3=char(id2(1));
    id4=strsplit(id3,'T');
    day=str2num(char(id4(end)));
    days_host(i)=day;
end

days_phage=zeros(1,length(bp.col_labels));

for i=1:length(bp.col_labels)
    id=char(bp.col_labels(i));
    id2=strsplit(id,'_');
    id3=char(id2(1));
    id4=strsplit(id3,'T');
    day=str2num(char(id4(end)));
    days_phage(i)=day;
end
%% calculate modularity and nestedness by same-time window

day_vec=unique(days_phage); % number of days
mod=[]; %vector to save modularity over time
nes=[]; %vector to save nestedness over time

for i=1:length(day_vec)
    idx_h=logical(days_host==day_vec(i));
    idx_p=logical(days_phage==day_vec(i));
    
    %Delete all zeros rows and columns
    zero_p=sum(bp.matrix(idx_h,idx_p))==0;
    zero_h=sum(bp.matrix(idx_h,idx_p),2)==0;
    idx_h1=idx_h;
    idx_h1(find(idx_h==1))=idx_h(find(idx_h==1))-zero_h';
    idx_p1=idx_p;
    idx_p1(find(idx_p==1))=idx_p(find(idx_p==1))-zero_p;
    
    idx_h=idx_h1;
    idx_p=idx_p1;
    submatrix=bp.matrix(idx_h,idx_p);
    
    
    sub_bp=Bipartite(submatrix);
    sub_bp.row_labels =  bp.row_labels(idx_h);
    sub_bp.col_labels =  bp.col_labels(idx_p);
    sub_bp.community = LPBrim(sub_bp.matrix); %Uses LPBrim algorithm
    sub_bp.community.Detect();   
    
    sub_bp.nestedness.Detect();
    ntc=Nestedness.NTC(sub_bp.matrix); % Using a different algorithm to calculate nestedness
    mod(i)=sub_bp.community.Qb;
    nes(i)=sub_bp.nestedness.N;
    
    %%statistics
%     sub_bp.statistics.DoNulls(1000, @NullModels.AVERAGE); % Create 1000 random matrices using the average
%     sub_bp.statistics.TestCommunityStructure();
%     sub_bp.statistics.TestNestedness();
%     sub_bp.statistics.Print();
%     [a,b]=size(submatrix);
%     fprintf('Degrees of freedom: %d\n',a*b-1);
%     fprintf('Q_b p-value=%e\n', erfc(sub_bp.statistics.Qb_values.zscore));
%     fprintf('Q_r p-value=%e\n', erfc(sub_bp.statistics.Qr_values.zscore));
%     fprintf('NODF p-value=%e\n', erfc(sub_bp.statistics.N_values.zscore));
    
    % make a figure
    figure(i);
    sub_bp.plotter.font_size = 8; 
    sub_bp.plotter.cell_color=[1,1,0];
    sub_bp.plotter.back_color = [0.65,0.65,0.65];
    subplot(1,2,1)
    
    % Matlab command to change the figure window;
    set(gcf,'Position',[0+50 72 932 922]);
    sub_bp.plotter.use_isocline = true; %The NTC isocline will be plotted.
    sub_bp.plotter.isocline_color = 'red'; %Decide the color of the isocline.
    sub_bp.plotter.PlotNestedMatrix();
    title(['$NODF = $',num2str(sub_bp.nestedness.N)],'interpreter','latex','fontsize',23);
    subplot(1,2,2)
    sub_bp.plotter.use_isocline = true; %Although true is the default value
    sub_bp.plotter.PlotModularMatrix();
    title(['$Q_b =',num2str(sub_bp.community.Qb), ';Q_r =',num2str(sub_bp.community.Qr),'\\N = $', num2str(sub_bp.community.N)],'interpreter','latex','fontsize',23);
  
end
%% save calculated modularity and nestedeness values over time into a csv
cd Data
filename='nodf_mod_dynamics.csv';
fp = fopen(filename);
fprintf(fp,'time(days),lpbrim,nodf\n');
for i=1:length(mod)
    fprintf(fp,strcat([char(times(i)),',',num2str(mod(i)),',',num2str(nes(i)),'\n']));
end    
fclose(fp);
cd ..
%% Calculate GLOBAL modularity/nestedness
% Calculating Modularity
bp.community = LeadingEigenvector(bp.matrix);
% The next flag is exclusive of Newman Algorithm and what it does is to
% performn a final tuning after each sub−division (see Newman 2006).
bp.community.DoKernighanLinTunning = true; % Default value
bp.community.Detect();
fprintf('The modularity value Qb is %f\n', bp.community.Qb);
fprintf('The fraction inside modules Qr is %f\n',bp.community.Qr);


%% Calculating nestedness
bp.nestedness.Detect();
fprintf('The Nestedness value is %f\n', bp.nestedness.N);
bp.printer.PrintStructureValues();

%% PLotting raw matrix
figure(1);
% Matlab command to change the figure window;
set(gcf,'Position',[0 72 1751 922]);
bp.plotter.font_size = 8; %Change the font size of the rows and labels
% Use different color for each kind of interaction
bp.plotter.use_type_interaction = true; %
bp.plotter.color_interactions(1,:) = [1 1 1]; %Red color for clear lysis
bp.plotter.color_interactions(2,:) = [0 0 1]; %Blue color for turbid spots
bp.plotter.back_color = 'blue';
% After changing all the format we finally can call the plotting function.
bp.plotter.PlotMatrix();
title('original PBIN')
%% NESTEDNESS matrix
figure(2);
% Matlab command to change the figure window;
set(gcf,'Position',[0+50 72 932 922]);
bp.plotter.use_isocline = true; %The NTC isocline will be plotted.
bp.plotter.isocline_color = 'red'; %Decide the color of the isocline.
bp.plotter.PlotNestedMatrix();
title(['Nestedness, $NODF = $',num2str(bp.nestedness.N)],...
'interpreter','latex','fontsize',23);

%% MODULAR matrix 
bp.plotter.back_color = 'blue';
figure(3);
bp.plotter.font_size=10;
bp.community = LPBrim(bp.matrix); %Uses LPBrim algorithm
bp.plotter.use_isocline = true; %Although true is the default value
bp.plotter.PlotModularMatrix();
title(['$Q = $',num2str(bp.community.Qb),' $c = $', num2str(bp.community.N)],'interpreter','latex','fontsize',23);
%% Save data
mod_host=bp.community.row_modules;
mod_phage=bp.community.col_modules;
idx_host=bp.community.index_rows;
idx_phage=bp.community.index_cols;
mod_h=mod_host(idx_host);
mod_p=mod_phage(idx_phage);

IDX_host=[idx_host(mod_h==1);idx_host(mod_h==2);idx_host(mod_h==3);idx_host(mod_h==4)];
IDX_phage=[flipud(idx_phage(mod_p==1));flipud(idx_phage(mod_p==2));flipud(idx_phage(mod_p==3))];


%% save tables
cd 'Data'
writematrix(bp.matrix(IDX_host,IDX_phage),'22-06-09_d21_LPbrim.csv','Delimiter',',','QuoteStrings',true);
writecell(bp.row_labels(IDX_host),'hostID_d21_LPbrim.csv','Delimiter',',','QuoteStrings',true);
writecell(bp.col_labels(IDX_phage),'phageID_d21_LPbrim.csv','Delimiter',',','QuoteStrings',true);
writematrix(mod_h,'modules_host.csv','Delimiter',',','QuoteStrings',true);
writematrix(flipud(mod_p),'modules_phage.csv','Delimiter',',','QuoteStrings',true);
cd ..


