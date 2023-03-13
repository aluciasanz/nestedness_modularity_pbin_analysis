% A Lucia-Sanz PBIN nestedness modularity analysis JUNE,30,2022

clear
clc
%% add path
PATH='~/Documents/ML_Meyer/JBorin/BiMat-master'; %write the path for BiMat-master
% mkdir(strcat(PATH,'/myfolder'))   
% GEN_PATH = genpath(PATH);
% addpath(GEN_PATH);
cd(strcat(PATH,'/myfolder'));

%% Import DATA
matrix_raw=readtable(strcat(cd,'/Data/22-06-09-matrix-21d_bis.csv')); % PBIN located in Data
hostID = matrix_raw.Properties.VariableNames(2:end); % Hosts are columns
phageID = matrix_raw.Phage; % Phage are rows
matrix=table2array(matrix_raw(:,2:end)); % convert table to matrix

%% Create the bipartite object
bp = Bipartite(matrix'); % hosts are rows, phage are columns
bp.row_labels =  hostID;
bp.col_labels = phageID;

%% Calculating Modularity
bp.community = LeadingEigenvector(bp.matrix);
% The next flag is exclusive of Newman Algorithm and what it does is to
% performn a final tuning after each sub−division (see Newman 2006).
bp.community.DoKernighanLinTunning = true; % Default value
bp.community.Detect();
fprintf('The modularity value Qb is %f\n', bp.community.Qb);
fprintf('The fraction inside modules Qr is %f\n',bp.community.Qr);

%%JAN 2023 until day 21
% Modularity:
% 	Used algorithm:             	  LeadingEigenvector
% 	N (Number of modules):      	                   5
% 	Qb (Standard metric):       	              0.2078
% 	Qr (Ratio of int/ext inter):	             -0.0552
% The modularity value Qb is 0.207758
% The fraction inside modules Qr is -0.055168
%% Calculating nestedness
bp.nestedness.Detect();
fprintf('The Nestedness value is %f\n', bp.nestedness.N);
bp.printer.PrintStructureValues();

%%JAN 2023 Until day 21
% 
% Nestedness NODF:
% 	NODF (Nestedness value):    	              0.6372
% 	NODF (Rows value):          	              0.5313
% 	NODF (Columns value):       	              0.9391
% The Nestedness value is 0.637200
% Modularity:
% 	Used algorithm:             	  LeadingEigenvector
% 	N (Number of modules):      	                   5
% 	Qb (Standard metric):       	              0.2078
% 	Qr (Ratio of int/ext inter):	             -0.0552
% Nestedness NODF:
% 	NODF (Nestedness value):    	              0.6372
% 	NODF (Rows value):          	              0.5313
% 	NODF (Columns value):       	              0.9391
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
