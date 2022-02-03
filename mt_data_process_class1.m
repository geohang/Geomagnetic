classdef mt_data_process_class1<handle
    properties
        sample_frequency;
        n_segment;
        n_each_segment;
    
        n_transform;
        n_decimation;
        % total data points for each components
        
    end
    methods
        function obj=mt_data_process_class1(sample_frequency,n_segment,...
                n_each_segment,n_decimation)
            if nargin==4
                obj.sample_frequency=sample_frequency;
                obj.n_segment=n_segment;
                obj.n_each_segment=n_each_segment;
                obj.n_decimation=n_decimation;
                n_data_transform=2^nextpow2(obj.n_each_segment);
                obj.n_transform=n_data_transform;
            end
        end
        function [e_segment,b_segment,b1_segment]=  ...
                data_partition(obj,e_data,b_data,b1_data)
           
%            this part is used to partition the data into data segement
%            for later use, keep in mind
%            the data should be n_segement multiple the whole data 
%            otherwise special care needed.
           if nargin==3
              e_segment=zeros(obj.n_segment,obj.n_each_segment);
              b_segment=zeros(obj.n_segment,obj.n_each_segment);
              for i=1:obj.n_segment
                e_segment(i,:)=e_data(1+(i-1)*obj.n_each_segment...
                   :obj.n_each_segment*i);
% %                 e_segment(2,i,:)=e_data(2,1+(i-1)*obj.n_each_segment...
%                    :obj.n_each_segment*i);
                b_segment(i,:)=b_data(1+(i-1)*obj.n_each_segment...
                   :obj.n_each_segment*i);
% %                 b_segment(2,i,:)=b_data(2,1+(i-1)*obj.n_each_segment...
%                    :obj.n_each_segment*i);
              end
           elseif nargin==4
              e_segment=zeros(obj.n_segment,obj.n_each_segment);
              b_segment=zeros(obj.n_segment,obj.n_each_segment);
              b1_segment=zeros(obj.n_segment,obj.n_each_segment);
              for i=1:obj.n_segment
                e_segment(i,:)=e_data(1+(i-1)*obj.n_each_segment...
                   :obj.n_each_segment*i);
%                 e_segment(2,i,:)=e_data(2,1+(i-1)*obj.n_each_segment...
%                    :obj.n_each_segment*i);
                b_segment(i,:)=b_data(1+(i-1)*obj.n_each_segment...
                   :obj.n_each_segment*i);
%                 b_segment(2,i,:)=b_data(2,1+(i-1)*obj.n_each_segment...
%                    :obj.n_each_segment*i);
                b1_segment(i,:)=b1_data(1+(i-1)*obj.n_each_segment...
                   :obj.n_each_segment*i);
%                 b1_segment(2,i,:)=b1_data(2,1+(i-1)*obj.n_each_segment...
%                    :obj.n_each_segment*i);
              end
           end
           
        end
        
        function [e_segment,b_segment,b1_segment]= ...
                data_detrend(obj,e_segment1,b_segment1, ...
                flag_detrend, b1_segment1)
          % The function will be used to remove
            % the mean and linear trend
            % first demean
            if nargin==4
              e_segment=zeros(obj.n_segment,obj.n_each_segment);
              b_segment=zeros(obj.n_segment,obj.n_each_segment);
              if (flag_detrend==1)  
                 for i=1:1
                   e_data_mean=mean(e_segment1(:,:),2);
                   b_data_mean=mean(b_segment1(:,:),2);
                   for j=1:obj.n_segment
                     e_data_temp(j,:)=(e_segment1(j,:)-e_data_mean);
                     b_data_temp(j,:)=(b_segment1(j,:)-b_data_mean);
                   end
                   time=((1:1:obj.n_each_segment)*(1/obj.sample_frequency))';% 20 is sample time in s
                   for j=1:obj.n_segment
                     p=polyfit(time,e_data_temp(j,:)',1);
                     p1=polyfit(time,b_data_temp(j,:)',1);
                     e_relativedata=p(1)*time+p(2);
                     b_relativedata=p1(1)*time+p1(2);
                     e_data(j,:)=e_data_temp(j,:)-e_relativedata';
                     b_data(j,:)=b_data_temp(j,:)-b_relativedata';
                   end
                   e_segment(:,:)=e_data(:,:);
                   b_segment(:,:)=b_data(:,:);
                 end
              elseif(flag_detrend==2)
                   for j=1:obj.n_segment
                      e_segment(j,:)=detrend(e_segment1(j,:));
                      b_segment(j,:)=detrend(b_segment1(j,:));
                   end
              end
            elseif nargin==5
              e_segment=zeros(obj.n_segment,obj.n_each_segment);
              b_segment=zeros(obj.n_segment,obj.n_each_segment);
              b1_segment=zeros(obj.n_segment,obj.n_each_segment);
            
                for j=1:obj.n_segment
                   e_segment(j,:)=detrend(e_segment1(j,:));
                   b_segment(j,:)=detrend(b_segment1(j,:));
                   b1_segment(j,:)=detrend(b1_segment1(j,:));
                end
             
            end
            % demean and detrend is finished
            end
          
           
        
        function [e_segment,b_segment,b1_segment]= ...
                taper_window(obj,e_segment1,b_segment1,b1_segment1)
         % start to apply taper windows for the data
            % to remove the possible allies
            e_segment=zeros(obj.n_segment,obj.n_each_segment);
            b_segment=zeros(obj.n_segment,obj.n_each_segment);
            b1_segment=zeros(obj.n_segment,obj.n_each_segment);
            w_temp=tukeywin(obj.n_each_segment,.2);
            if nargin==3
              for i=1:1
                for j=1:obj.n_segment
                  e_segment(j,:)=reshape(e_segment1(j,:),1, ...
                     obj.n_each_segment).*w_temp';
                  b_segment(j,:)=reshape(b_segment1(j,:),1, ...
                     obj.n_each_segment).*w_temp';
                end
              end
            elseif nargin==4
              for i=1:1
                for j=1:obj.n_segment
                  e_segment(j,:)=reshape(e_segment1(j,:),1, ...
                     obj.n_each_segment).*w_temp';
                  b_segment(j,:)=reshape(b_segment1(j,:),1, ...
                     obj.n_each_segment).*w_temp';
                   b1_segment(j,:)=reshape(b1_segment1(j,:),1, ...
                     obj.n_each_segment).*w_temp';
                end
              end 
            end
            % window taper has been finished
              
        end
        function  [e_fourier,b_fourier,frequency,b1_fourier] ...
                =data_fourier(obj,e_segment1,b_segment1,b1_segment1)
           % this function is used to calculate the Fourier transform
           % using the data after preprocessed 
           % and will output the coefficients e_fourier,b_fourier
           % is used to output the fourier transformation coefficients 
           % for both electrical field and magnetic field
           % frequency is used to store the frequencies after transfomation
           %e_spectral,b_spectral are used to output the spectrum of both
           % electric field and magnetic field in magnitude rather than
           % in power spectrum.

         
           e_fourier1=complex(zeros(obj.n_segment,obj.n_transform));
           b_fourier1=complex(zeros(obj.n_segment,obj.n_transform));
           b1_fourier1=complex(zeros(obj.n_segment,obj.n_transform));
           e_fourier=complex(zeros(obj.n_segment,obj.n_transform/2));
           b_fourier=complex(zeros(obj.n_segment,obj.n_transform/2));
           b1_fourier=complex(zeros(obj.n_segment,obj.n_transform/2));
           if nargin==3
            % pad zeros to pow of 2
              %for i=1:1
                 for j=1:obj.n_segment
                   e_fourier1(j,:)=fft(e_segment1(j,:), ...
                     obj.n_transform)/obj.n_each_segment;
                   b_fourier1(j,:)=fft(b_segment1(j,:), ...
                     obj.n_transform)/obj.n_each_segment;
                   % half the frequencies is abundant and we only 
                   % keep the first half of the coefficients into 
                   % e_fourier and b_fourier
                   e_fourier(j,:)=e_fourier1(j,2:obj.n_transform/2+1);
                   b_fourier(j,:)=b_fourier1(j,2:obj.n_transform/2+1);
                 end
              %end
           elseif nargin==4
%               for i=1:1
                for j=1:obj.n_segment
                  e_fourier1(j,:)=fft(e_segment1(j,:), ...
                     obj.n_transform)/obj.n_each_segment;
                  b_fourier1(j,:)=fft(b_segment1(j,:), ...
                     obj.n_transform)/obj.n_each_segment;
                  b1_fourier1(j,:)=fft(b1_segment1(j,:), ...
                     obj.n_transform)/obj.n_each_segment;
                  % half the frequencies is abundant and we only 
                  % keep the first half of the coefficients into 
                  % e_fourier and b_fourier
                  e_fourier(j,:)=e_fourier1(j,2:obj.n_transform/2+1);
                  b_fourier(j,:)=b_fourier1(j,2:obj.n_transform/2+1);
                  b1_fourier(j,:)=b1_fourier1(j,2:obj.n_transform/2+1);
                end
%               end
           end
            % calculate the powerspectral for each component 
            % for output
           frequency1=obj.sample_frequency/2*linspace(0,1,obj.n_transform/2+1);
           % Fourier transform process is finished.         
           frequency=frequency1(2:end);
        end
        function [e_stack,b_stack,datapoint_stack,frequency_band,b1_stack]...
                = frequency_time_stack ...
                (obj,e_fourier,b_fourier,band_point,flag_stackband,n_band,b1_fourier)
            % this function is used to stack coefficents in to one
            % particular band (specified by flag_stackband)
            % and is used for impedance estimation, stored in one matrix.
            %lag_stackband is the band used to stack,n_decimation
            % is used to 
            % output e_stack,b_stack are the stacked e and b spectral
            % stored in one row and frequency_band keep the up and low
            % frequency for particular band.
           
            if (flag_stackband<=n_band)
               
                band_value=band_point(flag_stackband,:);
                band_interval=band_value(2)-band_value(1);
                n_total=band_interval+1;
                e_stack=complex(zeros(1,obj.n_segment*n_total));
                b_stack=complex(zeros(1,obj.n_segment*n_total));
                b1_stack=complex(zeros(1,obj.n_segment*n_total));
                if nargin==6
                 % for i=1:2
                    for j=1:obj.n_segment
                      e_stack((band_interval+1)*(j-1)+1:...
                          j*(band_interval+1))=e_fourier(j,band_value(1):...
                          band_value(2));
                      b_stack((band_interval+1)*(j-1)+1:...
                          j*(band_interval+1))=b_fourier(j,band_value(1):...
                          band_value(2));
                    end
                 % end
                elseif nargin==7
%                   for i=1:2
                    for j=1:obj.n_segment
                      e_stack((band_interval+1)*(j-1)+1:...
                          j*(band_interval+1))=e_fourier(j,band_value(1):...
                          band_value(2));
                      b_stack((band_interval+1)*(j-1)+1:...
                          j*(band_interval+1))=b_fourier(j,band_value(1):...
                          band_value(2));
                      b1_stack((band_interval+1)*(j-1)+1:...
                          j*(band_interval+1))=b1_fourier(j,band_value(1):...
                          band_value(2));
                    end
%                   end
                end
                datapoint_stack=(band_interval+1)*obj.n_segment;
                sample_time=1/obj.sample_frequency;
                frequency_band=band_value/(obj.n_transform*sample_time);
            else
                error('stack number is great than your set band');
            end
            
        end
        
          function [bl,n_band] = setBandLims(obj,N,ndec)
        %    define frequency limits for bands that are roughly evenly spaced in
        %    log frequency, ndec bands per decade.   Second argument is optional,
        %    defaults to 8
          if nargin ==1 
           ndec = 8;
          end
          l10Lims = log10(1/N):1/ndec:log10(1/2);
          lims = round(10.^(l10Lims)*N);
          lims = [lims(diff(lims)>0) lims(end)];
          lims = lims(lims>1);
          bl = [lims(1:end-1)' lims(2:end)'-1]+1;
          n_band=length(bl);
        end
        
        function plot_data(obj,e_data,b_data,figure_n)
            % this function is used to plot the data
            % the input data e_data has dimension of 
            % (2,n_each_segment),b_data has
            % dimension of (2,n_each_segment)
            figure(figure_n);
            start1=max(abs(e_data(1,:)))+max(abs(e_data(1,:)))/10;
            start2=2*start1+max(abs(e_data(2,:)))+max(abs(e_data(2,:)))/10;
            start3=start2+max(abs(e_data(2,:)))+max(abs(e_data(2,:)))/10+...
               max(abs(b_data(1,:)))+max(abs(b_data(1,:)))/10;
           start4=start3+max(abs(b_data(1,:)))+max(abs(b_data(1,:)))/10 ...
               +max(abs(b_data(2,:)))+max(abs(b_data(2,:)))/10;
           plot(b_data(2,:)+start4);
           hold on;  
           plot(b_data(1,:)+start3);      
           plot(e_data(2,:)+start2);         
           plot(e_data(1,:)+start1);
           xlabel('data points');
           n_point=size(e_data,2);
           text(-n_point/8,start1,'Ex');
           text(-n_point/8,start2,'Ey');
           text(-n_point/8,start3,'Hx');
           text(-n_point/8,start4,'Hy');
           
        end
        end  
end
