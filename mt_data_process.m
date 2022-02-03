

clc;
close all;
clear all;
% load reference site in,远参考数据输入
[FileName,PathName]=uigetfile('CNB2000-20220120.mat');%('*.hor','Input the file you want rea
file=strcat(PathName,FileName);  
load(file);
KAKZ=BZ1;

H_NR=H_N
% Rn=length(H_NR);


[FileName,PathName]=uigetfile('*.mat');%('*.hor','Input the file you want rea
file=strcat(PathName,FileName);  
load(file);
%这里需要根据实际数据进行设置
a=6371;
cof=pi/180;
%90-x,x为地磁纬度
theta=(90+42.04)*cof;
%采样间隔为一小时
sample_frequency=1/(60*60);

% al=length(BZ1);bl=length(H_NR);
% al=[al bl;];
% m1=min(al);
% BZ1=BZ1(1:m1,1);
KAKZ=BZ1;
H_N=H_NR;
n=672;
% Dn=length(H_NR);

%远参考判断数据长度相等 
%if (Dn~=Rn)
%     error('data length is not equal')
% end


%根据数据长度进行切割
n_segment=38;
n_each_segment=5000;

n_decimation=8;

obj=mt_data_process_class1(sample_frequency,n_segment,n_each_segment,n_decimation);%%
[e_segment,b_segment,b_segmentR]=obj.data_partition(KAKZ',H_N,H_NR);


% obj.plot_data(e_plot,b_plot,1);

[e_segment,b_segment,b_segmentR]=obj.data_detrend(e_segment,b_segment,2,b_segmentR);

% plot the results for check

[e_segment,b_segment,b_segmentR]=obj.taper_window(e_segment,b_segment,b_segmentR);
% plot the results for check

% e_plot=reshape(e_segment(:,1,:),2,obj.n_each_segment);
% b_plot=reshape(b_segment(:,1,:),2,obj.n_each_segment);
% obj.plot_data(e_plot,b_plot,3);
% start to fourier transform.
[e_fourier,b_fourier,frequency_out,b_fourierR] ...
                =obj.data_fourier(e_segment,b_segment,b_segmentR);
            
% decide frequency band for output data
  band_point=obj.setBandLims(obj.n_transform,obj.n_decimation);
  n_band=length(band_point);
% stacking the data for diffrent frequency and time segment.
frequency=zeros(n_band,1);
impedance2=zeros(n_band,1);
c_error2=zeros(n_band,1);
bf=cell(1,n_band);
ef1=cell(1,n_band);
Cf=cell(1,n_band);

% KHfour=reshape(e_fourier,[n_decimation,16]);
% KZfour=reshape(b_fourier,[n_decimation,16]);
% for i=1:352
% C(i)=-a*tan(theta)/2*KZfour(i)*KHfour(i)'/(KHfour(i)*KHfour(i)');
% end

flag_method='least_square';    
% flag_method='remote_reference'flag_method='least_square'flag_method='robust'; 
 switch flag_method
     case 'least_square'
       for i=1:n_band
          flag_stackband=i;
          [e_stack,b_stack,datapoint_stack,frequency_band]= obj.frequency_time_stack ...
                (e_fourier,b_fourier,band_point,flag_stackband,n_band);
          frequency(i,1)=(frequency_band(1,2)-frequency_band(1,1))/2+frequency_band(1,1);
          % calculate the impedance 
          
          obj1=impedance_solution_class;
          [c,c_error]=obj1.least_square_method(e_stack,b_stack, ...
                datapoint_stack,theta,a);
          % save the magnetic field and impedance
          
          bf{i}=b_stack;
          ef1{i}=e_stack;
          Cf{i}=c;
          % output only impedance in xy yx direction and error.
         % impedance_arr(i,:,:)=[impedance(1,1),impedance(1,2);impedance(2,1),impedance(2,2)];
          %impedance2(i,:)=[impedance(1,2),impedance(2,1)];
         c_error2(i)=c_error;
          coh2(i)=abs((e_stack*b_stack')*(e_stack*b_stack'))*inv(b_stack*b_stack')*inv(e_stack*e_stack');
%          coh2(i)=coh;
       end 
        
 
  % calculate the apparent resistivity and phase
     case 'robust'
       for i=1:n_band
          flag_stackband=i;
          [e_stack,b_stack,datapoint_stack,frequency_band]= obj.frequency_time_stack ...
                (e_fourier,b_fourier,band_point,flag_stackband,n_band);
          frequency(i,1)=(frequency_band(1,2)-frequency_band(1,1))/2+frequency_band(1,1);
          % calculate the impedance 
          
          obj1=impedance_solution_class;
          [c,c_error]=obj1.robust_method(e_stack,b_stack, ...
                datapoint_stack,theta,a);
         %改// 
          bf{i}=b_stack;
          ef1{i}=e_stack;
          Cf{i}=c;
          c_error2(i)=c_error;
          coh2(i)=abs((e_stack*b_stack')*(e_stack*b_stack'))*inv(b_stack*b_stack')*inv(e_stack*e_stack');
         
         %//改
% % %           impedance2(i,:)=[impedance(1,2),impedance(2,1)];
% % %           impedance_arr(i,:,:)=[impedance(1,1),impedance(1,2);impedance(2,1),impedance(2,2)];
% % %           impedance_error2(i,:)=[impedance_error(1,2),impedance_error(2,1)];
% % %           disp(['Iterations :',num2str(i)]);
          
       end
     case 'remote_reference'
        for i=1:n_band
          flag_stackband=i;
          [e_stack,b_stack,datapoint_stack,frequency_band,b1_stack]= obj.frequency_time_stack ...
                (e_fourier,b_fourier,band_point,flag_stackband,n_band,b_fourierR);
          frequency(i,1)=(frequency_band(1,2)-frequency_band(1,1))/2+frequency_band(1,1);
          % calculate the impedance 
          
          obj1=impedance_solution_class;
          [c,c_error]=obj1.remote_reference(e_stack,b_stack, ...
                b1_stack,datapoint_stack,theta,a);
           Cf{i}=c;
           c_error2(i)=c_error;
           coh2(i)=abs((e_stack*b_stack')*(e_stack*b_stack'))*inv(b_stack*b_stack')*inv(e_stack*e_stack');
          disp(['Band :',num2str(i)]);
         
        end
       
         
     case 'robustreference'
        for i=1:n_band
          flag_stackband=i;
          [e_stack,b_stack,datapoint_stack,frequency_band,b1_stack]= obj.frequency_time_stack ...
                (e_fourier,b_fourier,band_point,flag_stackband,n_band,b1_fourier);
          frequency(i,1)=(frequency_band(1,2)-frequency_band(1,1))/2+frequency_band(1,1);
          % calculate the impedance 
          
          obj1=impedance_solution_class;
          [impedance,impedance_error]=obj1.robust_reference(e_stack,b_stack, ...
                b1_stack,datapoint_stack);
          impedance2(i,:)=impedance;
          impedance_arr(i,:)=impedance;
          impedance_error2(i,:)=[impedance_error];
          coh2(i)=abs((e_stack*b_stack')*(e_stack*b_stack'))*inv(b_stack*b_stack')*inv(e_stack*e_stack');
          disp(['Band :',num2str(i)]);
          
       end
 end
%  output the impedance and magnetic field  for test 


%   save( 'C-CNB-2000-2022.mat','frequency','Cf','c_error2','coh2');

%%
clf;
subplot(4,1,[2:3]);
T1=1./frequency/24/60/60;
axis([3 120 -1200 1200])
hold on
ylabel('C-respond')
xlabel('period/d')
% set(gca,'fontsize',10,'fontweight','bold','FontName','Times New Roman','xscal','log','YLim',[0,1]);% a1=cell2mat(Cf);
a1=cell2mat(Cf);
y11=real(a1);
y12=imag(a1); 
e4=sqrt(c_error2/2);

errorbar(T1,y11,e4,'--bd','LineWidth',1)
errorbar(T1,y12,e4,'--bd','LineWidth',1)
set(gca,'fontsize',10,'fontweight','bold','FontName','Times New Roman','xscal','log','YLim',[-1000,1500]);

subplot(4,1,1);
axis([3 120 0 1])
hold on
plot(T1,coh2);