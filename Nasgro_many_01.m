function num_err=Nasgro_many_01(stress)
%Initialize NASGRO variables and material properties
R=0.1;
C1=1.61e-10;delK_eff=1.8;p=0.35;alpha=2.73;m=2.77;Kc=50;% NASGRO values for AF LPBF Ti-6Al-4V
D=7.1176/1000; %Sample diameter
flow_stress=1150;
del_sigma_0=450;
del_Kth_lc=2.7;Yw=0.65;
sqarea_0=(1./pi).*((del_Kth_lc./(Yw.*del_sigma_0)).^2); %El-haddad parameter expressed with murakami's area
sqarea_0eff=(1./pi).*((delK_eff./(Yw.*del_sigma_0)).^2);
li=sqarea_0-sqarea_0eff; % crack extension for closure to build-up 
sig_0=1100;% 

%Newmans Closure formulation
A0=(0.8613-0.3387.*alpha+0.0465.*(alpha.^2)).*((cos((pi/2).*0.2)).^(1/alpha));
A1=(1.047-0.233.*alpha).*0.2;
A3=2*A0+A1-1;
A2=1-A0-A1-A3;
f=max(R,A0+A1*R+A2*(R^2)+A3.*(R.^3));
U=(1-f)./(1-R);


% Populate defect locations in polar coordinates and sizes 
theta=[2.544400909	
];
% Radial locations
radial=[3.394247128	
];
sizes=[0.055997925	
]; %defect sizes 

% Creat defect combined matrix for the algoriwth, note entries in row vectors

a_ext=zeros(1,size(sizes,2));% set all initial crack extensions to 0 
Y_place=zeros(1,size(sizes,2)); %placeholder for magnification factors

ai_mat=[sizes./1000;theta;radial./1000;Y_place;a_ext]; %Build matrix out of vectors

%Initialise loop variables
N_sum=0; %cycle counter
i=1;
N=100; %May need to reduce this

    

 while i<100000000000
    
% multicoalation3 merges what needs to be merged and sends back full according to surface code
%ai_mat(1,:)= sqrt(((ai_mat(1,:)).^2)./pi); %Convert sqrt(area) -> radii for coealescence (circle)
[merged_population]=multi_coalition_version3(ai_mat);
ai_mat=merged_population; %set ai to crack size from merged_population
%ai_mat(1,:)= sqrt(pi.*(ai_mat(1,:).^2)); %Convert radii -> sqrt(area) for NASGRO
       
del_K=stress.*ai_mat(4,:).*sqrt(pi.*ai_mat(1,:));
[b,ind]=max(ai_mat(1,:));%tracking the maximum size defect for checks
        
        
% Keep count of cycles
if i~=1
    N_sum=N_sum+N;
else
    N_sum=N;
end
        
        
if ai_mat(1,ind)> D %Size constraint on largest crack
    disp('Crack propogates through');
    break;
end
         
%Plasticity Correction term
rp=(1/6).*pi.*(((del_K./(1-R))./sig_0).^2);
       
%Independant growth important in the near-threshold growth
% Calculates a growing stress value to be used in the formulation
if ai_mat(1,ind) <(D/2) && max(del_K)<20 % net stress formulation limited 
    growing_stress=net_stress(ai_mat(1,ind),stress,D);
end
  
% Maierhof closure formula
delK_R_th=(delK_eff+(del_Kth_lc-delK_eff).*(1-exp(-ai_mat(5,:)/li)));
F1=1-(1-(U.^m)).*(1-exp(-ai_mat(5,:)/li));


% da/dN for NASGRO equation Maierhof formulation
% vector of crack growth rate for all defects
da_dN=(C1.*F1.*((growing_stress.*ai_mat(4,:).*sqrt(pi.*(ai_mat(1,:)+rp))).^(m-p)).*((growing_stress.*ai_mat(4,:).*sqrt(pi.*(ai_mat(1,:)+rp)))-delK_R_th).^(p));


%On first iteration check if cracks grow at all has to be above
%del_K effth

conv_factor=1;
        % Note: initially set to 1, for convergence possible may need to relax this condition
        % when delk_eff -> delk_lc as becomes more difficult to find stress
        % within tolerances for multiple cracks


if i==1
    if  max(del_K)<= conv_factor*delK_eff  % Check if crack initially grows
        disp('Nothing grows');
        N_sum=0;
        break;
    end
end
        
        
        % Crack growth rate slowing down due to build up of crack closure
        
       for ii=1:size(ai_mat(1,:),2)
            if del_K(ii) >= conv_factor*delK_eff %Something happens otherwise remains the same
                if real(da_dN(ii))>0 && imag(da_dN(ii))==0 % Grows on its own
                    ai_mat(1,ii)=ai_mat(1,ii)+da_dN(ii)*N; % Crack size vector
                    ai_mat(5,ii)=ai_mat(5,ii)+da_dN(ii)*N; % Crack extension vector
                elseif real(da_dN(ii))<0 || imag(da_dN(ii))~=0 %Above threshold but still weird
                    da_dN(ii)=1e-12; %Chosen suficiently small
                    ai_mat(1,ii)=ai_mat(1,ii)+da_dN(ii)*N; % Crack size vector
                    ai_mat(5,ii)=ai_mat(5,ii)+da_dN(ii)*N; % Crack extension vector
                end
            end
            
        end
        

        
 
% Failure condition 1:
%Monitor net_section stress caused by largest defect does not exceed flow stress
if ai_mat(1,ind) <(D/2) 
    net_section=net_stress(ai_mat(1,ind),stress,D);
    if net_section>flow_stress
        disp('Net section exceeds flow- stress too high');
        break;
    end
end
   
%Failure condition 2:
      
      % Calculate Kmax
      del_K_max=stress.*ai_mat(4,ind).*sqrt(pi.*ai_mat(1,ind));

        %Checks for stress too high, stress intensity and net_section
        if (del_K_max/(1-R))> Kc  %close to Kc and breaks
            disp('Stress exceeds plain strain fracture- Stress too high');
            disp(del_K_max);
            disp(N_sum);
            break;
        end
     
%Failure condition 3:
%Stress too low, needs to oscillate around convergence (1000 cycles used as a tolerance here)

if  N_sum > (10^7)+1000
    N_sum=0;
    disp('Too many cycles- Stress too low');
    break;
end

i=i+1; %increment counter


 end

num_err=abs(10^7-N_sum);

end
