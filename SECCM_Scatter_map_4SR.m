clear
close all
clc

%% Experiment Details

droplet_size=200;
num_landings_y = 1; %Size of grid in y axis
num_landings_x = 11; %Size of grid in x axis
landing_separation_y=10; %Separation between landing points in y
landing_separation_x=10; %Separation between landing points in x
Lx=num_landings_x*landing_separation_x; %Length of x axis for map
Ly=num_landings_y*landing_separation_y; %Length of y axis for map

%% create x,y coordinates

%B = repmat(A,n) returns an array containing n copies of A in the row and column dimensions
x=repmat(0:landing_separation_x:landing_separation_x*(num_landings_x-1),num_landings_y*4,1);
y=repmat(0:landing_separation_y:landing_separation_y*(num_landings_y*4-1),num_landings_x,1);
x=x';
x=x(:);
y=y(:);

%% Extract data
% 200 mV/s
file1=('PureAl_phosphate_4scanrates_5s_CA_2_200mVs');
S1=load(file1);

Ecorr1=zeros(num_landings_x*num_landings_y,1);
OCP1=zeros(num_landings_x*num_landings_y,1);

for k=1:3*num_landings_x*num_landings_y

    format_spec1='Trace_1_%d_%d_%d';
    
    [m1,n1]=size(S1.(sprintf(format_spec1,k,1,1))); 
    %using the formatting operators specified by formatSpec and returns the resulting text in str.
   
    currents1=S1.(sprintf(format_spec1,k,1,1))(:,2);
    potentials1=S1.(sprintf(format_spec1,k,1,2))(:,2);
     
    if mod(k,3)==0 %PDP data
        yc1=log10(abs(currents1));
        [icorr1,I1]=min(yc1);
        Ecorr1(k/3)=S1.(sprintf(format_spec1,k,1,2))(I1,2);
    end
    if mod(k,3)==1 %OCP data
        OCP1((k+2)/3)=S1.(sprintf(format_spec1,k,1,2))(end,2);
    end
end

% 100 mV/s
file2=('PureAl_phosphate_4scanrates_5s_CA_2_100mVs');
S2=load(file2);

Ecorr2=zeros(num_landings_x*num_landings_y,1);
OCP2=zeros(num_landings_x*num_landings_y,1);

for k=1:3*num_landings_x*num_landings_y

    format_spec2='Trace_2_%d_%d_%d';
    
    [m2,n2]=size(S2.(sprintf(format_spec2,k,1,1))); 
   
    currents2=S2.(sprintf(format_spec2,k,1,1))(:,2);
    potentials2=S2.(sprintf(format_spec2,k,1,2))(:,2);
     
    if mod(k,3)==0 
        yc2=log10(abs(currents2));
        [icorr2,I2]=min(yc2);
        Ecorr2(k/3)=S2.(sprintf(format_spec2,k,1,2))(I2,2);
    end
    if mod(k,3)==1 
        OCP2((k+2)/3)=S2.(sprintf(format_spec2,k,1,2))(end,2);
    end
end

% 50 mV/s
file3=('PureAl_phosphate_4scanrates_5s_CA_2_50mVs');
S3=load(file3);

Ecorr3=zeros(num_landings_x*num_landings_y,1);
OCP3=zeros(num_landings_x*num_landings_y,1);

for k=1:3*num_landings_x*num_landings_y

    format_spec3='Trace_3_%d_%d_%d';
    
    [m3,n3]=size(S3.(sprintf(format_spec3,k,1,1))); 
   
    currents3=S3.(sprintf(format_spec3,k,1,1))(:,2);
    potentials3=S3.(sprintf(format_spec3,k,1,2))(:,2);
     
    if mod(k,3)==0
        yc3=log10(abs(currents3));
        [icorr3,I3]=min(yc3);
        Ecorr3(k/3)=S3.(sprintf(format_spec3,k,1,2))(I3,2);
    end
    if mod(k,3)==1
        OCP3((k+2)/3)=S3.(sprintf(format_spec3,k,1,2))(end,2);
    end
end

% 25 mV/s
file4=('PureAl_phosphate_4scanrates_5s_CA_2_25mVs');
S4=load(file4);

Ecorr4=zeros(num_landings_x*num_landings_y,1);
OCP4=zeros(num_landings_x*num_landings_y,1);

for k=1:3*num_landings_x*num_landings_y

    format_spec4='Trace_4_%d_%d_%d';
    
    [m4,n4]=size(S4.(sprintf(format_spec4,k,1,1))); 

    currents4=S4.(sprintf(format_spec4,k,1,1))(:,2);
    potentials4=S4.(sprintf(format_spec4,k,1,2))(:,2);
     
    if mod(k,3)==0
        yc4=log10(abs(currents4));
        [icorr4,I4]=min(yc4);
        Ecorr4(k/3)=S4.(sprintf(format_spec4,k,1,2))(I4,2);
    end
    if mod(k,3)==1
        OCP4((k+2)/3)=S4.(sprintf(format_spec4,k,1,2))(end,2);
    end
end

% Combine 4 rows
Ecorr=[Ecorr4;Ecorr3;Ecorr2;Ecorr1];
OCP=[OCP4;OCP3;OCP2;OCP1];

stdev_OCP = std(OCP)
average_OCP = mean(OCP)
CI_OCP = tinv([0.025  0.975],43)* std(OCP)
%% 
% OCP map
f1=figure;
scatter(x,y,droplet_size,OCP,'filled');

box on
axis([-landing_separation_y Lx -landing_separation_x Ly*4])
caxis([-1,-0.2])
pbaspect([num_landings_x num_landings_y*4 1])

h=colorbar('eastoutside');
set(gca,'FontSize',15,'LineWidth',1.5,'FontWeight','bold')
ylabel(h,'E_c_o_r_r(OCP) vs. SCE (V)','FontSize',15);
xlabel('X (\mum)','FontSize',16)
ylabel('Y (\mum)','FontSize',16)

%%
% Ecorr map
f2=figure;
scatter(x,y,droplet_size,Ecorr,'filled');

box on
axis([-landing_separation_y Lx -landing_separation_x Ly*4])
caxis([-1,-0.2])
pbaspect([num_landings_x num_landings_y*4 1])
set(gca,'layer','top')

h=colorbar('eastoutside');
set(gca,'FontSize',15,'LineWidth',1.5,'FontWeight','bold');
ylabel(h,'E_c_o_r_r(PDP) vs. SCE (V)','FontSize',15);
xlabel('X (\mum)','FontSize',16)
ylabel('Y (\mum)','FontSize',16)

%%  Sample OCP curves
figure
for k1=[4]
    for k=k1*3-2
        
    format_spec='Trace_1_%d_%d_%d';
   
    potentials=S1.(sprintf(format_spec,k,1,2))(:,2);
    times=S1.(sprintf(format_spec,k,1,2))(:,1);
       
    
    
    x=times(:);
    y=potentials(:);
    plot(x,y,'linewidth',2);  
      
    hold on   
  
    end
end

for k1=[7 11]
    for k=k1*3-2
        
    format_spec='Trace_3_%d_%d_%d';
   
    potentials=S3.(sprintf(format_spec,k,1,2))(:,2);
    times=S3.(sprintf(format_spec,k,1,2))(:,1);
       
    
    
    x=times(:);
    y=potentials(:);
    plot(x,y,'linewidth',2);  
      
    hold on   
  
    end
end

for k1=[2]
    for k=k1*3-2
        
    format_spec='Trace_4_%d_%d_%d';
   
    potentials=S4.(sprintf(format_spec,k,1,2))(:,2);
    times=S4.(sprintf(format_spec,k,1,2))(:,1);
       
    
    
    x=times(:);
    y=potentials(:);
    plot(x,y,'linewidth',2);  
      
    hold on   
  
    end
end

xlabel('t (s)')
ylabel('E vs. SCE (V)')
ylim([-1.8 -0.6])
set(gca,'FontSize',18,'linewidth',2) 
legend('1','2','3','4');
legend location southeast;

