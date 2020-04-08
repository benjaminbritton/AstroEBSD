function ipf_cols = Map_IPF(ipf_dir,MapQ)

ipf_cols=zeros([MapQ.ypts, MapQ.xpts,3]);

%set up some prefactors - useful for colouring the ipfs
rfac=1E6;
fudge=0.0001;
pow=0.8;

for y_co=1:MapQ.ypts
    for x_co=1:MapQ.xpts
        Q_s=MapQ.Q_map_sym{y_co,x_co};
        
        %rotate the IPF direction on the PF
        ipf_v=Q_VecRot(Q_s,ipf_dir(:)); %find where the 'ipf_dir' vector is pointing for orientation Q_s, so rotate ipf_dir by Q_s
        ipf_v(3,ipf_v(3,:)==0)=eps*1000; % deal with a numerical error
        
        switch MapQ.Crystal(y_co,x_co)
            case 1 %cubic
                %defines the R G B settings for the corners of the IPF
                col_grid=[0 0 1;[1 0 1]/sqrt(2);[1 1 1]/sqrt(3)];
                
                %find the point which fits in this corner
                n2=find(round(ipf_v(2,:)*rfac)>=-eps*1000 & round(ipf_v(1,:)*rfac)<=round(abs(ipf_v(3,:))*rfac) & round(ipf_v(1,:)*rfac)>=round(ipf_v(2,:)*rfac));
%          the vector in the fundamental zone has positive y, and has x>y>0 (in in the lower 45deg positive pie
%          slice of the IPF). Not sure what happens to x relating to z?

                 
                %make it positive (improper inversion)
                ipf_v(3,n2(1))=abs(ipf_v(3,n2(1)));
                
                %create the colors
                p_rgb_set=ipf_v(:,n2(1)).'/col_grid;
                
                %fudge the colors very slightly
                p_rgb_set=(p_rgb_set+fudge).^pow;
                
                %bring the max value to white
                cols_temp=p_rgb_set./max(p_rgb_set);
            case 2
                %calculate the pole figure position
                ipf_p=ipf_v(1:2,:).*repmat(sign(ipf_v(3,:))./(ipf_v(3,:)+sign(ipf_v(3,:))),2,1);
                %find the value that is in the correct triangle (HCP)
                n2=find(p(1,:)>=0 & p(2,:)>=0 & (p(1,:)*tand(30)-p(2,:))>=-eps*1000); %this is for hkl, what about bruker?
                p1=ipf_v(:,n2(1));
                
                theta=acos(p1(3));
                if theta>pi/2
                    theta=pi-theta;
                end
                rho=atan(p1(2)/p1(1));
                
                %TSL
                %cols_temp=[1-2*theta/pi,(2*theta/pi)*(1-(6*rho/pi)),(2*theta/pi)*(6*rho/pi)];
                
                %Bruker
                cols_temp=[1-2*theta/pi,(2*theta/pi)*(6*rho/pi),(2*theta/pi)*(1-(6*rho/pi))];
                cols_temp=(cols_temp+fudge).^pow;
        end
        
        ipf_cols(y_co,x_co,:)=cols_temp;
    end
end
