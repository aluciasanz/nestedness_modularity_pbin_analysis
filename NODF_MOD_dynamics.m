% jan 2023 - ALSANZ nestednes and modularity over time
clear 
clc
%% add path set working directory
PATH='~/Documents/ML_Meyer/JBorin/BiMat-master'; %write the path for BiMat-master
mkdir(strcat(PATH,'/myfolder'))   
GEN_PATH = genpath(PATH);
addpath(GEN_PATH);
cd(strcat(PATH,'/myfolder'))

%% Import DATA
matrix_raw=readtable(strcat(cd,'/Data/22-06-09-matrix-21d_bis.csv')); % PBIN located in Data
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

%% calculate modularity/nestedness in a 2 or 3 day overlapping window
day_vec=unique(days_phage);

mod=zeros(1,length(day_vec)-2);
nes=zeros(1,length(day_vec)-2);

for j=1:(length(day_vec)-2)
    i=day_vec(j);
    i2=day_vec(j+1);
%     idx_h=logical((days_host==i) + (days_host==i2));
%     idx_p=logical((days_phage==i) +(days_phage==i2));
    
    i3=day_vec(j+2);
    idx_h=logical((days_host==i) + (days_host==i2) +(days_host==i3));
    idx_p=logical((days_phage==i) +(days_phage==i2) +(days_phage==i3));
    
    submatrix=bp.matrix(idx_h,idx_p);
    
    sub_bp=Bipartite(submatrix);
    sub_bp.row_labels =  bp.row_labels(idx_h);
    sub_bp.col_labels =  bp.col_labels(idx_p);
    sub_bp.community = LPBrim(sub_bp.matrix); %Uses LPBrim algorithm
    sub_bp.community.Detect();
    
    mod(j)=sub_bp.community.Qb;    
    
    sub_bp.nestedness.Detect();
    nes(j)=sub_bp.nestedness.N;   
    ntc=Nestedness.NTC(sub_bp.matrix); % Using a different algorithm to calculate nestedness
    
    
%     figure(j);
%     % Matlab command to change the figure window;
%     set(gcf,'Position',[0 72 1751 922]);
%     sub_bp.plotter.font_size = 8; %Change the font size of the rows and labels
%     % Use different color for each kind of interaction
%     sub_bp.plotter.use_type_interaction = true; %
%     sub_bp.plotter.color_interactions(1,:) = [1 1 1]; %Red color for clear lysis
%     sub_bp.plotter.color_interactions(2,:) = [0 0 1]; %Blue color for turbid spots
%     sub_bp.plotter.back_color = 'blue';
%     % After changing all the format we finally can call the plotting function.
%     sub_bp.plotter.PlotMatrix();
    
    
%     figure(j+10);
%     sub_bp.plotter.PlotModularMatrix();
%     title(['$Q = $',num2str(sub_bp.community.Qb),' $c = $', num2str(sub_bp.community.N)],'interpreter','latex','fontsize',23);
%     
    %save modular matrix
%     writematrix(sub_bp.community.matrix,strcat(['modularity_T',num2str(i),'-',num2str(i3),'.csv']),'Delimiter',',','QuoteStrings',true);
%     writecell(sub_bp.row_labels(sub_bp.community.index_rows),strcat(['modularity_T',num2str(i),'-',num2str(i3),'id_host.csv']),'Delimiter',',','QuoteStrings',true);
%     writecell(sub_bp.col_labels(sub_bp.community.index_cols),strcat(['modularity_T',num2str(i),'-',num2str(i3),'id_phage.csv']),'Delimiter',',','QuoteStrings',true);  
    
%     figure(j+20);
%         % Matlab command to change the figure window;
%     set(gcf,'Position',[0+50 72 932 922]);
%     sub_bp.plotter.use_isocline = true; %The NTC isocline will be plotted.
%     sub_bp.plotter.isocline_color = 'red'; %Decide the color of the isocline.
%     sub_bp.plotter.PlotNestedMatrix();
%     title(['$NODF = $',num2str(sub_bp.nestedness.N),' $NTC = $',num2str(ntc.N),' $ connectance = $', num2str(ntc.connectance)],'interpreter','latex','fontsize',23);
%     
    %save NODF matrix
%     writematrix(sub_bp.nestedness.matrix,strcat(['NODF_T',num2str(i),'-',num2str(i3),'.csv']),'Delimiter',',','QuoteStrings',true);
%     writematrix(sub_bp.row_labels(sub_bp.community.index_rows),strcat(['modularity_T',num2str(i),'-',num2str(i3),'idx_host.csv']),'Delimiter',',','QuoteStrings',true);
%     writematrix(sub_bp.col_labels(sub_bp.community.index_cols),strcat(['modularity_T',num2str(i),'-',num2str(i3),'idx_phage.csv']),'Delimiter',',','QuoteStrings',true);  
%     
    

end
%% 3 day window
figure(101);plot(mod,'ko-','LineWidth',2); title('Modularity');xlabel('Day')
xticks([1:5]);
times={'3-9','6-12','9-15','12-18','15-21'};
xticklabels(times)
hold on;
plot(nes,'ro-','LineWidth',2); title('NODF');xlabel('Day')
xticks([1:5]);
xticklabels(times)
legend('LPbrim','NODF');
title('3 day window analysis')
%% save calculated modularity and nestedeness values over time into a csv
filename='nodf_mod_dynamics.csv';
fp = fopen(filename);
fprintf(fp,'time(days),lpbrim,nodf\n');
for i=1:length(mod)
    fprintf(fp,strcat([char(times(i)),',',num2str(mod(i)),',',num2str(nes(i)),'\n']));
end    
fclose(fp);
%% 2 Day window
figure(101);plot(mod,'ko-','LineWidth',2); title('Modularity');xlabel('Day')
xticks([1:6])
xticklabels({'3-6','6-9','9-12','12-15','15-18','18-21'})
hold on;
plot(nes,'ro-','LineWidth',2); title('NODF');xlabel('Day')
xticks([1:6])
xticklabels({'3-6','6-9','9-12','12-15','15-18','18-21'})
