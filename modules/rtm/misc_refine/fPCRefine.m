function [PC_best] = fPCRefine(PC_start,PatternIn,rotmat,PatternInfo,Refine,SettingsXCF,screen_int,RTM_info,RTM_setup)
%% fPCRefine - McAuliffe PC Refinement using gradient descent
%
% Inputs - 
% PC_start = [PCx, PCy, PCz]
%
% Refine = structure
% Refine.ss=0.08; %initial probe volume
% Refine.p=2; %order of polynomial to fit to tetrahedron
% Refine.n_its=10; %number of interations
%
% EBSD_geom = geometry of the detector

%rotation matricies

phbest=0;

%create the descent kernel
xs=[1,3,2,2,2,2,2];
ys=[2,2,1,3,2,2,2];
zs=[2,2,2,2,1,3,2];
PC_best=PC_start;

%initalise a counter
skipgradcalc=0;

%count for a number of iterations
for n=1:Refine.n_its
    
    phs=zeros(3,3,3);

    for i=1:7
        x=xs(i);
        y=ys(i);
        z=zs(i);
        
        x_itval=(x-2).*Refine.ss;
        y_itval=(y-2).*Refine.ss;
        z_itval=(z-2).*Refine.ss;
        
        PC_refined=[PC_best(1)+x_itval,PC_best(2)+y_itval,PC_best(3)+z_itval];
        [ EBSD_geom ] = EBSP_Gnom( PatternInfo,PC_refined );
        
        %Optimise the fit, allowing for small orientation updates
        [~,regout] = refine5(PatternIn,EBSD_geom,PC_start,rotmat,SettingsXCF,screen_int,RTM_info.isHex,RTM_setup);
        
        %extract the PH from the refinement
        phs(x,y,z)=regout(4);
    end
    
    %fit to iterpolants and find gradients    
    xcol=reshape(phs(:,2,2),1,3);
    ycol=reshape(phs(2,:,2),1,3);
    zcol=reshape(phs(2,2,:),1,3);
    sv=[-Refine.ss,0,Refine.ss];
    
    %sum peakheights across the axes and fit a polynomial
    [fitx]=polyfit(sv,xcol,Refine.p); %xaxis
    [fity]=polyfit(sv,ycol,Refine.p); %yaxis
    [fitz]=polyfit(sv,zcol,Refine.p); %zaxis
    
    %evaluate the gradient at x=0
    grad_x=fitx(Refine.p);
    grad_y=fity(Refine.p);
    grad_z=fitz(Refine.p);
    
    %descend
    grad=0.1*[grad_x.*Refine.ss,grad_y.*Refine.ss,grad_z.*Refine.ss];
    
    Refine.PC_it(n,:)=PC_start;
    Refine.PH_it(n)=phs(2,2,2);
    
    phtrial=phs(2,2,2);
    
    % decide whether to continue
    if n<Refine.n_its && phtrial>phbest
        
        %update the current best PC and peak height
        %move the trial PC (for the next iteration)
        PC_best=PC_start;
        phbest=phtrial;
        
        if skipgradcalc==0
            %set the next PC to follow the gradient
            PC_start=[PC_best(1)+grad(1),PC_best(2)+grad(2),PC_best(3)+grad(3)];
        else
        end
        
        skipgradcalc=0;
        counter=0;
        
    elseif n<Refine.n_its %else reduce the step size
        Refine.ss=0.9*Refine.ss;
        %and move slightly randomly
        reduction=(5-mod(counter,5))/5;
        
        grad=reduction.*0.01.*[(rand-0.5),(rand-0.5),(rand-0.5)].*PC_start; % -0.25 to 0.25 percent of PC_start
        PC_best=[PC_best(1)+grad(1),PC_best(2)+grad(2),PC_best(3)+grad(3)];
        skipgradcalc=1;
        
        counter=counter+1;
        
    end
    
end

end


