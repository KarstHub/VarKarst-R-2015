%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VarKarst_R routine for karst recharge
% according to Hartmann et al. (2016, GMD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [sim]=VarKarst_R(input,parameter_x0)

%[Rech_avg]=VarKarst_Med_Model(Pre,Epot,a,Epi_MRT_mean,Soi_max_mean,Epi_max_mean)

    n=15; %number of model compartments
    
% Decompose input

    Pre=input{1};
    Epot=input{2};

%Variable Initiation for high speed calculations
    Soi=zeros(length(Pre),n); %soil
    Epi=zeros(length(Pre),n); %epikarst
    Inf=zeros(length(Pre),n); %infiltration
    Eact=zeros(length(Pre),n); %actual evaporation
    Exc_Soi=zeros(length(Pre),n); %excess soil water
    Exc_Epi=zeros(length(Pre),n); %saturation excess
    Rech=zeros(length(Pre),n); %recharge
    InSurf=zeros(length(Pre),n); %sufrace runoff to next unsaturated compartment
    RunOff=zeros(length(Pre),n);%sufrace runoff
    
%Parameter preparation    
    a=parameter_x0(1);
    Epi_MRT_mean=parameter_x0(2);
    Soi_max_mean=parameter_x0(3);
    Epi_max_mean=parameter_x0(4);

    MRT_dummy=1/Epi_MRT_mean*(a+1);
    Epi_MRT=1./min(1,MRT_dummy*((1:n)/n).^a);
    
    Soi_dummy=Soi_max_mean*2^(a/(a+1)); 
    Soi_max=Soi_dummy*((1:n)/n).^a;
    Epi_dummy=Epi_max_mean*2^(a/(a+1));
    Epi_max=Epi_dummy*((1:n)/n).^a;   
    
        
%Time loop
    for t=1:size(Pre,1)
        if t==1
            Soi_0=zeros(1,n);%Soi_Ini;
            Epi_0=zeros(1,n);%Epi_Ini;
        else
            Soi_0=Soi(t-1,:);
            Epi_0=Epi(t-1,:);
        end

        Inf(t,:)=Pre(t,1)+InSurf(t,:);
        Eact_dummy=Epot(t,1).*min(max((Soi_0)./(Soi_max),0),1);

        Soi_dummy=max(Soi_0+Inf(t,:)-Eact_dummy,0);
        Eact(t,:)=Eact_dummy+min(Soi_0+Inf(t,:)-Eact_dummy,0);
        Soi(t,:)=min(Soi_max,Soi_dummy);
        Exc_Soi(t,:)=max(0,Soi_dummy-Soi_max);           
        Rech_dummy=Epi_0./Epi_MRT;        
        Epi_dummy=max(Epi_0+Exc_Soi(t,:)-Rech_dummy,0);
        Rech(t,:)=min(Epi_0+Exc_Soi(t,:),Rech_dummy);
        Epi(t,:)=min(Epi_max,Epi_dummy);
        Exc_Epi(t,:)=max(0,Epi_dummy-Epi_max);
                
        if t<size(Pre,1) & find(Exc_Epi(t,:)>0, 1, 'last' )<n
            InSurf(t+1,find(Exc_Epi(t,:)>0, 1, 'last' )+1)=sum(Exc_Epi(t,1:n-1),2);
        elseif t<size(Pre,1) & find(Exc_Epi(t,:)>0, 1, 'last' )==n
            RunOff(t+1,n)=sum(Exc_Epi(t,1:n),2);
        end
        if Rech(t,:)<0
            error('Rech<0!');
        end
    end
    
   sim=[mean(Rech,2) mean(Eact,2)];