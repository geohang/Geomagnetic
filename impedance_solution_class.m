classdef impedance_solution_class<handle
    properties 
     flag_method;   
        
    end
    methods
    
        function [c,c_error]=least_square_method(obj,e_stack,b_stack, ...
                datapoint_stack,theta,a)
            %This function is used to calculate the impedance 
            % Using least square method.
            % e_stack and b_stack is 2*datapoint_stack
            vector_notation=1;
           % vector_notation =1 is using  the 2*N vector for both E and B
           % vector_notation=2 is using the N*2 vector for B and N*1 vector
           % for E which is regression
            c=complex(zeros(1,1));
            if vector_notation==1
             
               c=-a*tan(theta)/2*(e_stack*b_stack')*inv(b_stack*b_stack');
               
               % residual estimate
               %C(i)=-a*tan(theta)/2*KZfour(i)*KHfour(i)'/(KHfour(i)*KHfour(i)');
               data_residual=(-a*tan(theta)/2.*e_stack-c.*b_stack)/(a.*tan(theta)/2);
               variance_data=sum(abs(data_residual).^2,2)/...
                (datapoint_stack-2);
               variance_impedance=inv(b_stack*b_stack');
               variance_impedance=(-a*tan(theta)/2)^2*variance_impedance*variance_data;
%                variance_impedance(2,2).*variance_data(1,1); ...
%                variance_impedance(1,1).*variance_data(2,1),...
%                variance_impedance(2,2).*variance_data(2,1)];
            elseif vector_notation==2
                % estimate the impedance 
                impedance_temp1=inv(b_stack*b_stack')*b_stack*e_stack(1,:)';
                impedance_temp2=inv(b_stack*b_stack')*b_stack*e_stack(2,:)';
                impedance=[impedance_temp1';impedance_temp2'];
                % calculate the residual 
                data_residual=e_stack-impedance*b_stack;
                variance_data=sum(abs(data_residual).^2,2)/...
                (datapoint_stack-2);
                % caculate the variance of impedance
                variance_temp1=variance_data(1,1)*inv(b_stack*b_stack');
                variance_temp2=variance_data(2,1)*inv(b_stack*b_stack');
                variance_impedance=[variance_temp1(1,1),variance_temp1(2,2); ...
                   variance_temp2(1,1),variance_temp2(2,2) ];
            end
%           coh=abs((e_stack*b_stack')*(e_stack*b_stack'))*inv(e_stack*e_stack')*inv(b_stack*b_stack');
%           c_error=abs((e_stack*b_stack')*(e_stack*b_stack'))*inv(e_stack*e_stack')*inv(b_stack*b_stack');
            c_error= variance_impedance;


        end
        
        function [impedance,impedance_error] =remote_reference(obj,e_stack,b_stack, ...
                b_reference_stack,datapoint_stack,theta,a)
        
            vector_notation=1;
            impedance=complex(zeros(1,1));
            if vector_notation==1
               % calcualte the impedance
                impedance=-a*tan(theta)/2*(e_stack* b_reference_stack')*inv(b_stack* b_reference_stack');
               
%                impedance=e_stack*b_reference_stack'* ...
%                    inv(b_stack*b_reference_stack');
               % residual estimate
               %data_residual=e_stack-impedance*b_stack;
               data_residual=(-a*tan(theta)/2*e_stack-impedance*b_stack)/(a*tan(theta)/2);
               variance_data=sum(abs(data_residual).^2,2)/...
                (datapoint_stack-2);
              % calculate the variance of the impedance 
               variance_impedance=(-a*tan(theta)/2)^2*inv(b_stack*b_reference_stack' ...
                    *inv( b_reference_stack* b_reference_stack')* ...
                     b_reference_stack*b_stack')*variance_data;
%                variance_impedance=[variance_impedance(1,1).*variance_data(1,1), ...
%                    variance_impedance(2,2).*variance_data(1,1); ...
%                    variance_impedance(1,1).*variance_data(2,1),...
%                    variance_impedance(2,2).*variance_data(2,1)];
            elseif vector_notation==2
                % calculate the impedance using remote reference method
                impedance_temp1=inv(b_reference_stack*b_stack')*...
                    b_reference_stack*e_stack(1,:)';
                impedance_temp2=inv(b_reference_stack*b_stack')* ...
                    b_reference_stack*e_stack(2,:)';
                impedance=[impedance_temp1';impedance_temp2'];
                % data residual estimation
                data_residual=e_stack-impedance*b_stack;
                variance_data=sum(abs(data_residual).^2,2)/...
                (datapoint_stack-2);
               % calculate the variance of the impedance data
                variance_temp1=variance_data(1,1)*inv(b_stack*b_reference_stack' ...
                    *inv( b_reference_stack'* b_reference_stack)* ...
                     b_reference_stack*b_stack');
                variance_temp2=variance_data(2,1)*inv(b_stack*b_reference_stack' ...
                    *inv( b_reference_stack'* b_reference_stack)* ...
                     b_reference_stack*b_stack');
                variance_impedance=[variance_temp1(1,1),variance_temp1(2,2); ...
                   variance_temp2(1,1),variance_temp2(2,2) ];
            end
            impedance_error= variance_impedance;
        
        end
        
        
        
        function [impedance,impedance_error]=robust_method(obj,e_stack,b_stack, ...
                datapoint_stack,theta,a)
          % this function is used to calculate the impedance using robust
          % method iteratively, the objective function can be general, e.g.
          % robust method, huber method and also least square method in
          % different way flag_method is used to indicate which method are
          % used such as 1 for  robust method and 2 for huber method.
          % flag_sheme is used to identify the iteration scheme used to
          % calculate the impedance
          flag_scheme=1;
          flag_method=2;
          % solve the least square function as the initial solution of the
          % iteration 
          impedance=-a*tan(theta)/2*e_stack*b_stack'*inv(b_stack*b_stack');
% % %           bb_inv=inv(b_stack*b_stack');
% % %            impedance_temp1=bb_inv*...
% % %                     b_stack*e_stack(1,:)';
% % %            impedance_temp2=bb_inv* ...
% % %                     b_stack*e_stack(2,:)';
% % %            impedance=[impedance_temp1';impedance_temp2']; 
           %initialize det_residual
           det_residual=1;
           % calculate the deviation and residual for least sqaure
           % impedance as intial eastimate
           % % %            data_residual=e_stack-impedance*b_stack;
           data_residual=-a*tan(theta)/2*e_stack-impedance*b_stack;
           % data standard deviation calculation
           deviation_data=sqrt(sum(abs(data_residual).^2,2)/...
                  (datapoint_stack-2));
           % data variance calculation 
           variance_data=deviation_data.^2;
           new_residual=sum(sum(abs(data_residual).^2));
           % used to keep old electrical filed for cleanning
           e_old1=e_stack';
% % %            e_old2=e_stack(2,:)';
           % set the convergence criterial
           convergence_set=0.000001;
           while (det_residual>convergence_set)
           
              % calculate the weights, flag_method=1 for robust method
              % flag_method=2 for huber method
% % %               [e1_weights,e2_weights] =obj.weight_calculation( ...
% % %                 data_residual,deviation_data,flag_method);
         [e1_weights] =obj.weight_calculation( ...
                data_residual,deviation_data,flag_method);
              disp('weight calculation done...');
              % select the interation scheme for calcualtion
              if flag_scheme==1 
                % compute the impedance using using scheme 1 
                %weighted least square method
                % mutiplication of weight with magentic components 
% % %                 w_b1=[b_stack(1,:).*e1_weights;b_stack(2,:).*e1_weights];
                w_b1=b_stack.*e1_weights;
% % %           w_b2=[b_stack(1,:).*e2_weights;b_stack(2,:).*e2_weights];
                % impedance estimation using weighted least square method
%                 impedance_temp1=inv(w_b1 ...
%                   *b_stack')*w_b1*e_stack(1,:)';
                impedance=-a*tan(theta)/2*e_stack*w_b1'*inv(b_stack*w_b1');
% % %                 impedance_temp2=inv(w_b2 ...
% % %                    *b_stack')*w_b2*e_stack(2,:)';
                % organize the impedance components into impedance matrix
% % %                 impedance=[impedance_temp1';impedance_temp2'];
              elseif flag_scheme==2
              
                % mutiplication of weight with magentic components 
                w_b1=[b_stack(1,:).*e1_weights;b_stack(2,:).*e1_weights];
                w_b2=[b_stack(1,:).*e2_weights;b_stack(2,:).*e2_weights];
                % mutiplication of weight with electric components 
                w_e1=e1_weights'.*e_old1;
                w_e2=e2_weights'.*e_old2;
                % mutiplication of 1-weight with magnetic components 
                w1_b1=[(1-e1_weights)'.*b_stack(1,:)', ...
                    (1-e1_weights)'.*b_stack(2,:)'];
                w1_b2=[(1-e2_weights)'.*b_stack(1,:)',...
                    (1-e2_weights)'.*b_stack(2,:)'];
                % clean the data  to exclude poor data points for 
                % electric components
                e_new1=w_e1+w1_b1* ...
                    impedance_temp1;
                e_new2=w_e2+w1_b2* ...
                    impedance_temp2;
                % calculate the impedance using least square method
                % using cleaned data
                impedance_temp1=inv(b_stack*b_stack')*b_stack*e_new1;
                impedance_temp2=inv(b_stack*b_stack')*b_stack*e_new2;
                % organize impedance components into 2*2 matrix
                impedance=[impedance_temp1';impedance_temp2'];
              end
        
             % keep old residual for comparison
             old_residual=new_residual;
             % update the new data residual
             data_residual=-a*tan(theta)/2*e_stack-impedance*b_stack  % e_stack-impedance*b_stack;
             % calculate the standard deviation of the data
             deviation_data=sqrt(sum(abs(data_residual).^2,2)/...
                  (datapoint_stack-2));
             % data variance estimation
             variance_data=deviation_data.^2;
             % calculate the new data residual
             new_residual=sum(sum(abs(data_residual).^2));
             disp(['old_residuall:',num2str(old_residual)]);
              % update the det_residual
             det_residual=(new_residual-old_residual)/ ...
                  ((new_residual+old_residual)/2);
              disp(['det_residual :',num2str(det_residual)]);
           end
           
           [e1_weights,e1_weights_derivative]...
                =obj.weight_calculation(data_residual,deviation_data,flag_method);
            
           data_residual=(-a*tan(theta)/2*e_stack-impedance*b_stack)/(a*tan(theta)/2);
           
           if (flag_method==1)
           % error estimate based on matrix operation
            disp('linearized impedance error estimation for ruboust method...');
            
            % multiply weights with magnetic components
             w_b1=[b_stack(1,:).*e1_weights;b_stack(2,:).*e1_weights];
% % %              w_b2=[b_stack(1,:).*e2_weights;b_stack(2,:).*e2_weights];
             a_matrix1=inv(w_b1 ...
                   *b_stack')*w_b1;
            % variance for [zxx;zxy]
            variance_temp1=variance_data(1,1)*a_matrix1*a_matrix1';
            a_matrix2=inv(w_b2 ...
                   *b_stack')*w_b2;
            % variance for [zyx;zyy]
            variance_temp2=variance_data(1,1)*a_matrix2*a_matrix2';    
            impedance_error=[variance_temp1(1,1),variance_temp1(2,2); ...
                   variance_temp2(1,1),variance_temp2(2,2) ];
            % calculate variance for impedance based on M-estimator error 
            % function  from M-estimator error calculation
           elseif (flag_method==2)
                disp('using huber error equation... ');
                %compute covariance for zxi components
                weighted_residual1=datapoint_stack^2*sum(e1_weights'.* ...
                abs(data_residual(1,:)').^2);
                derivative_sum1=(datapoint_stack-2)*sum(e1_weights_derivative)^2;
                covariance_zx= (-a*tan(theta)/2)^2*weighted_residual1/(derivative_sum1)* ...
                    inv(b_stack*b_stack');
                % compute covariance for zyi components
% % %                 weighted_residual2=datapoint_stack^2*sum(e2_weights'.* ...
% % %                 abs(data_residual(2,:)').^2);
% % %                 derivative_sum2=(datapoint_stack-2)*sum(e2_weights_derivative)^2;
% % %                 covariance_zy= weighted_residual2/(derivative_sum2)* ...
% % %                     inv(b_stack*b_stack');
                  impedance_error=covariance_zx;
% % %
% impedance_error=[covariance_zx(1,1),covariance_zx(2,2); ...
% % %                    covariance_zy(1,1),covariance_zy(2,2) ];
            end
        end 
        
        function [e1_weights,e1_weights_derivative] =weight_calculation(obj,...
                data_residual,deviation_data,flag_method)
%        function [e1_weights,e2_weights,e1_weights_derivative, ...
%                 e2_weights_derivative] =weight_calculation(obj,...
%                 data_residual,deviation_data,flag_method)
            % this function is used to calculate weights for weighted least
            % square estimation deviation_data is normalization term,standard
            % deviation from residual of last iteration
            % data_residual is data residual 
            % flag_method is indicator for different =1 for which object function 
            % is used 1 for L1 function; 2 for huber function
            % for L1 function 
            if flag_method==1 
                % normalized weight calcualtion
                e1_weights=1./abs(data_residual(1,:)/deviation_data(1,1));
% % %           e2_weights=1./abs(data_residual(2,:)/deviation_data(2,1));
                % initialized the weight derivative
                e1_weights_derivative=e1_weights;
% % %           e2_weights_derivative=e2_weights;
                % caculate the weight derivative;
                e1_weights_derivative=0;
% % %           e2_weights_derivative=0;
            elseif flag_method==2
                % calculate the weights, flag_method=1 for robust method
                % flag_method=2 for huber method
                e1_residual=abs(data_residual(1,:)/deviation_data(1,1));
% % %                 e2_residual=abs(data_residual(2,:)/deviation_data(2,1)); 
                zero_residual=1.5;
                % initialized the weight derivative
                e1_weights_derivative=e1_residual;
% % %                 e2_weights_derivative=e2_residual;
                % apply the huber function
                % weights for e1_residual less r0 then
                e1_residual(e1_residual<zero_residual)=1;
% % %                 e2_residual(e2_residual<zero_residual)=1;
                % weights for e1_residual greater and equal r0 then
                e1_residual(e1_residual>zero_residual)=zero_residual./ ...
                   e1_residual(e1_residual>=zero_residual)  ;
% % %                 e2_residual(e2_residual>zero_residual)=zero_residual./ ...
% % %                    e2_residual(e2_residual>=zero_residual)  ;
               % calculate the derivative for residual less than
               % zero_residual
                e1_weights_derivative(e1_weights_derivative<zero_residual)=1;
% % %                 e2_weights_derivative(e2_weights_derivative<zero_residual)=1;
                % calculate the derivative for residual more and equal than
                % zero_residual
                e1_weights_derivative(e1_weights_derivative>=zero_residual)=0;
% % %                 e2_weights_derivative(e2_weights_derivative>=zero_residual)=0;
                % pass the value to e1_weights and e2_weights
                e1_weights=e1_residual;
% % %                 e2_weights=e2_residual; 
                % calculate the weight derivative;
                
            end 
        end
        
        function [apparent_resistivity,phase,resistivity_error,phase_error] ...
                 =resistivity_phase(obj,impedance,impedance_error,frequency)
                 % this function is used to calculate apparent resistivity
                 % and phase (corresponding errors) from impedance and its
                 % errors;
                 % impedance is N*2 matrix(ie Zxy,Zyx),impedance_error also N*2 matrix
                 % for impedance error, N is bandnunmber,frequency is
                 % N*1 vector
                 %¡¡apparent_resistivity and phase are N*2 matrix.
                 %¡¡resistivity_error adn phase_error are N*2 error matrix.
                 % estimate the apparent resistivity and phase standard
                 % deviation based on linearized calculation
                 
                 % calcualte the absolute impedance and impedance error.
                 absolute_impedance=abs(impedance);
                 absolute_error=sqrt(abs(impedance_error));
                 % calculate the apparent resisitivity without multiply
                 % period.
                 rho_temp=1/5*absolute_impedance.^2;
                 % calculate the period from frequency
                 T=1./frequency;
                 % calculate the apparent resisitivity by multiply period
                 apparent_resistivity=[rho_temp(:,1).*T, rho_temp(:,2).*T];
                 % calculate the phase for Zxy,Zyx in radians
                 phase_temp1=atan(imag(impedance(:,1)./real(impedance(:,1))));
                 phase_temp2=atan(imag(impedance(:,2)./real(impedance(:,2))));
                 % transform to degree from radians and combine them into
                 % phase matrix
                 phase_temp1=phase_temp1*(180/pi);
                 phase_temp2=phase_temp2*(180/pi);
                 phase=[phase_temp1,phase_temp2];
                 % calculate the linzed errors for phase and apparent resistivity
                 rho_error=[0.4*absolute_impedance(:,1).*absolute_error(:,1),...
                     0.4*absolute_impedance(:,2).*absolute_error(:,2)];
                 rho_error=[rho_error(:,1).*T,rho_error(:,2).*T];
                 resistivity_error=rho_error;
                 phi_error=[1./absolute_impedance(:,1).*absolute_error(:,1)/sqrt(2), ...
                        1./absolute_impedance(:,2).*absolute_error(:,2)/sqrt(2)];
                 % calculate the phase errror using alan chave's equation
                 phi_error1=[asin(absolute_error(:,1)./absolute_impedance(:,1)),...
                 asin(absolute_error(:,2)./absolute_impedance(:,2))];
                 % the error of phase in radius are changed to degree
                 phi_error=phi_error*180/pi;  
                 phi_error1=phi_error1*180/pi;
                 phase_error=phi_error;
                  
        end     
        
        
        
       function [impedance,impedance_error]=robust_reference(obj,e_stack,b_stack, ...
                b_reference_stack,datapoint_stack)
          % this function is used to calculate the impedance using robust
          % method iteratively, the objective function can be general, e.g.
          % robust method, huber method and also least square method in
          % different way flag_method is used to indicate which method are
          % used such as 1 for  robust method and 2 for huber method.
          % flag_sheme is used to identify the iteration scheme used to
          % calculate the impedance
          flag_method=2;
          % solve the least square function as the initial solution of the
          % iteration 
          impedance=e_stack*b_reference_stack'* ...
                   inv(b_stack*b_reference_stack');
           %initialize det_residual
           det_residual=1;
           % calculate the deviation and residual for least sqaure
           % impedance as intial eastimate
           data_residual=e_stack-impedance*b_stack;
           % data standard deviation calculation
           deviation_data=sqrt(sum(abs(data_residual).^2,2)/...
                  (datapoint_stack-2));
           % data variance calculation 
           variance_data=deviation_data.^2;
           new_residual=sum(sum(abs(data_residual).^2));
           % used to keep old electrical filed for cleanning
           e_old1=e_stack(1,:)';
           e_old2=e_stack(2,:)';
           % set the convergence criterial
           convergence_set=0.000001;
           while (det_residual>convergence_set)
           
              % calculate the weights, flag_method=1 for robust method
              % flag_method=2 for huber method
              [e1_weights,e2_weights] =obj.weight_calculation( ...
                data_residual,deviation_data,flag_method);
              disp('weight calculation done...');
              % select the interation scheme for calcualtion
              
              % compute the impedance using using scheme 1 
              %weighted least square method
              % mutiplication of weight with magentic components 
              w_b1=[b_reference_stack(1,:).*e1_weights;b_reference_stack(2,:).*e1_weights];
              w_b2=[b_reference_stack(1,:).*e2_weights;b_reference_stack(2,:).*e2_weights];
              % impedance estimation using weighted least square method
              impedance_temp1=inv(w_b1 ...
                  *b_stack')*w_b1*e_stack(1,:)';
              impedance_temp2=inv(w_b2 ...
                  *b_stack')*w_b2*e_stack(2,:)';
              % organize the impedance components into impedance matrix
              impedance=[impedance_temp1';impedance_temp2'];
        
             % keep old residual for comparison
             old_residual=new_residual;
             % update the new data residual
             data_residual=e_stack-impedance*b_stack;
             % calculate the standard deviation of the data
             deviation_data=sqrt(sum(abs(data_residual).^2,2)/...
                  (datapoint_stack-2));
             % data variance estimation
             variance_data=deviation_data.^2;
             % calculate the new data residual
             new_residual=sum(sum(abs(data_residual).^2));
             disp(['old_residuall:',num2str(old_residual)]);
              % update the det_residual
             det_residual=abs(new_residual-old_residual)/ ...
                  ((new_residual+old_residual)/2);
              disp(['det_residual :',num2str(det_residual)]);
           end
           % error estimate based on matrix operation
            disp('linearized impedance error estimation for ruboust method...');
            [e1_weights,e2_weights,e1_weights_derivative,e2_weights_derivative]...
                =obj.weight_calculation(data_residual,deviation_data,flag_method);
            % multiply weights with magnetic components
             w_b1=[b_reference_stack(1,:).*e1_weights;b_reference_stack(2,:).*e1_weights];
             w_b2=[b_reference_stack(1,:).*e2_weights;b_reference_stack(2,:).*e2_weights];
             a_matrix1=inv(w_b1 ...
                   *b_stack')*w_b1;
            % variance for [zxx;zxy]
            variance_temp1=variance_data(1,1)*a_matrix1*a_matrix1';
            a_matrix2=inv(w_b2 ...
                   *b_stack')*w_b2;
            % variance for [zyx;zyy]
            variance_temp2=variance_data(1,1)*a_matrix2*a_matrix2';    
            impedance_error=[variance_temp1(1,1),variance_temp1(2,2); ...
                   variance_temp2(1,1),variance_temp2(2,2) ];
      end 
                
        

    function  plot_rho_phi(obj,apparent_resistivity,phase,frequency, ...
        resistivity_error,phase_error)
       % this function is used to plot the apparent resistivity and phase for 
       % different frequencies with error bar,apparent_resistivity,phase, 
       % resistivity_error and phase_error are
       % n*2 matrix (ie rhoxy,rhoyx, phixy,phiyx)
       n_frequency=size(frequency,1);
       period=1./frequency;
       % error bar calcualtion for apparent resistivity
       Da1=[apparent_resistivity(:,1)-1*resistivity_error(:,1),...
           apparent_resistivity(:,2)-1*resistivity_error(:,2)];
       % aviod Da1 equal negative value and zero value
       Da1(Da1<0)=0.001;
       Da2=[apparent_resistivity(:,1)+1*resistivity_error(:,1),...
           apparent_resistivity(:,2)+1*resistivity_error(:,2)];
       % error bar calcualtion for phase.
       Dp1=[phase(:,1)-1*phase_error(:,1),...
           phase(:,2)-1*phase_error(:,2)];
       % aviod Dp1 equal negative value 
       Dp1(Dp1<0)=0;
        Dp2=[phase(:,1)+1*phase_error(:,1),...
           phase(:,2)+1*phase_error(:,2)];
       % start to plot
       % plot the apparent resistivity
       subplot(2,1,1);
       loglog(period,apparent_resistivity(:,1),'or','markersize',5,'linewidth',1);
       hold on;
       loglog(period,apparent_resistivity(:,2),'ob','markersize',5,'linewidth',1);
       legend('\rho_x_y','\rho_y_x');
       %plot the error bar for zxy
       for i=1:n_frequency
         DT1(i)=10^((log10(period(i))-0.04));
         DT2(i)=10^((log10(period(i))+0.04));
         loglog([DT1(i),DT2(i)],[Da1(i,1),Da1(i,1)],'r','linewidth',1);
         loglog([DT1(i),DT2(i)],[Da2(i,1),Da2(i,1)],'r','linewidth',1);
         loglog([period(i),period(i)],[Da1(i,1),Da2(i,1)],'r','linewidth',1)
       end 
       
       % plot the error bar zyx
       for i=1:n_frequency
         DT1(i)=10^((log10(period(i))-0.04));
         DT2(i)=10^((log10(period(i))+0.04));
         loglog([DT1(i),DT2(i)],[Da1(i,2),Da1(i,2)],'b','linewidth',1);
         loglog([DT1(i),DT2(i)],[Da2(i,2),Da2(i,2)],'b','linewidth',1);
         loglog([period(i),period(i)],[Da1(i,2),Da2(i,2)],'k','linewidth',1);
       end 
       ylabel('Apparent resistivity (\Omega\cdotm)');
       set(gca,'ticklength',[0.015 0.03],'xticklabel',[]);
    
       % find min and max of apparent resistivity for setting axis length 
       ymin=min(min(apparent_resistivity))-min(min(apparent_resistivity))/10;
       ymin=max(ymin,0.0001);
       ymax=max(max(apparent_resistivity))+...
            max(max(apparent_resistivity))/10;
       set(gca,'ylim',[ymin,100]);
       hold off;
       % plot the phase
       subplot(2,1,2);
       semilogx(period,phase(:,1),'or','markersize',5,'linewidth',1);
       hold on;
       semilogx(period,phase(:,2),'ob','markersize',5,'linewidth',1);
       legend('\Phi_x_y','\Phi_y_x');
       % plot the error bar for phase error bar
       for i=1:n_frequency
         DT1(i)=10^((log10(period(i))-0.04));
         DT2(i)=10^((log10(period(i))+0.04));
         semilogx([DT1(i),DT2(i)],[Dp1(i,1),Dp1(i,1)],'r','linewidth',1);
         semilogx([DT1(i),DT2(i)],[Dp2(i,1),Dp2(i,1)],'r','linewidth',1);
         semilogx([period(i),period(i)],[Dp1(i,1),Dp2(i,1)],'r','linewidth',1);
       end 
      
      % plot the error bar for phase error bar
       for i=1:n_frequency
         DT1(i)=10^((log10(period(i))-0.04));
         DT2(i)=10^((log10(period(i))+0.04));
         semilogx([DT1(i),DT2(i)],[Dp1(i,2),Dp1(i,2)],'b','linewidth',1);
         semilogx([DT1(i),DT2(i)],[Dp2(i,2),Dp2(i,2)],'b','linewidth',1);
         semilogx([period(i),period(i)],[Dp1(i,2),Dp2(i,2)],'k','linewidth',1);
       end 
        % find min and max of apparent resistivity for setting axis length 
        ymin=min(min(phase))-min(min(phase))/10;
        ymax=max(max(phase))+...
            max(max(phase))/10;
        ymin=max(ymin,0.0001);
       set(gca,'ylim',[0,80]);
       ylabel('Phase (\circ)');
       xlabel('Period (s)');
    
       hold off;


                 
    end
  end
end