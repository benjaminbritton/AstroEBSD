function [Settings_Cor,Settings_Rad,Settings_PCin]=Astro_EBSPset(EBSP_one,Settings_Cor,Settings_Rad,Settings_PCin,InputUser)
% Assess how to index the pattern

phases_list=phasefolder_read(InputUser);

%use one background correct to populate settings_cor
if exist('Settings_Cor') ~= 1
    Settings_Cor=struct;
end
[~,Settings_Cor ] = EBSP_BGCor( EBSP_one.PatternIn,Settings_Cor );


%read and set up the PC start
if isfield(Settings_PCin,'start')
    PC_start=Settings_PCin.start;
else
    PC_start=[0.5 0.5 0.5];
end

%read and set up the PC start
if isfield(Settings_PCin,'range')
    PC_range=Settings_PCin.range;
else
    PC_range=[0.1 0.1 0.1];
end

%  Create and then hide the GUI as it is being constructed.
f = figure('Visible','off','Position',[300,300,800,600]);


%figure box positioning
xsep=80; xstart=50; xwid=70; ystart=100; ysep=-18; yhig=15;

h_close=uicontrol('style','pushbutton','string','Close',...
    'position',[xstart+xsep*8 ystart+5*ysep xwid yhig]);
h_close.Units='normalized';
h_close.Callback=@closefun;


% Create the data to plot.
h_raw = axes('Units','Pixels','Position',[50,400,200,185],'Parent',f);

h_indexed = axes('Units','Pixels','Position',[550 400 200 185],'Parent',f,'Visible','off');

h_cor = axes('Units','Pixels','Position',[300,400,200,185],'Parent',f);
h_cor.Units = 'normalized';
axis(h_cor,'equal');

h_rad = axes('Units','Pixels','Position',[50,150,200,185],'Parent',f,'title','Radon Transform');
h_rad.Units = 'normalized';

h_norms = axes('Units','Pixels','Position',[300,150,200,185],'Parent',f);
h_norms.Units = 'normalized';

EBSP_single.PatternIn=EBSP_one.PatternIn;
[EBSP_single.ebsp_cor,EBSP_single.PatternInfo ] = EBSP_BGCor( EBSP_single.PatternIn,Settings_Cor );
[ EBSP_single.Peak_Centre,EBSP_single.Single.Peak_Set_All,EBSP_single.Peak_Set_All,...
    EBSP_single.R_EBSP,EBSP_single.R_Edge,EBSP_single.R_rho,EBSP_single.R_theta ] ...
    = EBSP_RadHunt( EBSP_single.ebsp_cor,Settings_Rad);
%generate the geometry
% PC_start=get(h_pc_xe,'UserData');
[ EBSD_geometry ] = EBSP_Gnom( EBSP_single.PatternInfo,PC_start );
[ nhat_gnom] = EBSP_NormConv( EBSP_single.Peak_Centre,EBSP_single.PatternInfo.size,PC_start);

i_raw=imagesc(EBSP_single.PatternIn,'Parent',h_raw);
i_cor=imagesc(EBSP_single.ebsp_cor,'Parent',h_cor);
i_rad=imagesc(EBSP_single.R_theta,EBSP_single.R_rho,EBSP_single.R_EBSP,'Parent',h_rad);
% i_test=imagesc(EBSP_One.ebsp_cor,'Parent',h_indexed);  i_test.Visible='on';

hold(h_rad,'on'); s_rad=scatter(EBSP_single.Peak_Centre(:,1),EBSP_single.Peak_Centre(:,2),'r','Parent',h_rad);hold(h_rad,'off');

h_raw.Units='pix'; axis(h_raw,'equal'); axis(h_raw,'tight','xy'); h_raw.Units = 'normalized';
h_cor.Units='pix'; axis(h_cor,'equal'); axis(h_cor,'tight','xy'); h_cor.Units = 'normalized';
h_rad.Units='pix'; axis(h_rad,'equal'); axis(h_rad,'tight','xy'); h_rad.Units = 'normalized';
axis(h_rad,'square') %make the radon square looking

h_indexed.Units='pix'; axis(h_indexed,'equal'); axis(h_indexed,'tight','xy'); h_indexed.Units = 'normalized';

set(i_raw,'UserData',EBSP_single);
[i_norms,p_norms]=DrawNorms(EBSP_single.ebsp_cor,nhat_gnom,EBSD_geometry,h_norms);
colormap('gray');

%% BG Controls
%hot pixel
h_check_hot=uicontrol('style','checkbox','string', ...
    'Hot Pixel', 'value', Settings_Cor.hotpixel,...
    'position',[xstart ystart xwid yhig]);
h_check_hot.Units = 'normalized';

h_check_hotthres=uicontrol('style','edit','string', ...
    num2str(Settings_Cor.hot_thresh),...
    'position',[xstart ystart+ysep xwid yhig]);
h_check_hotthres.Units = 'normalized';

%gaussian
h_check_gauss=uicontrol('style','checkbox','string', ...
    'Gauss Filter', 'value', Settings_Cor.gfilt,...
    'position',[xstart+xsep ystart xwid yhig]);
h_check_gauss.Units = 'normalized';

h_check_gfilt_s=uicontrol('style','edit','string', ...
    num2str(Settings_Cor.gfilt_s),...
    'position',[xstart+xsep ystart+ysep xwid yhig]);
h_check_gfilt_s.Units = 'normalized';

%resize
h_check_resize=uicontrol('style','checkbox','string', ...
    'Resize', 'value', Settings_Cor.resize,...
    'position',[xstart+2*xsep ystart xwid yhig]);
h_check_resize.Units = 'normalized';


h_check_resize_t=uicontrol('style','edit','string', ...
    num2str(Settings_Cor.size),...
    'position',[xstart+2*xsep ystart+ysep xwid yhig]);
h_check_resize_t.Units = 'normalized';

%real bg
h_check_realbg=uicontrol('style','checkbox','string', ...
    'Real Background', 'value', Settings_Cor.RealBG,...
    'position',[xstart+4*xsep ystart xwid yhig]);
h_check_realbg.Units = 'normalized';

%radius
h_check_radius=uicontrol('style','checkbox','string', ...
    'Radius','value', Settings_Cor.radius,...
    'position',[xstart+3*xsep ystart xwid yhig]);
h_check_radius.Units = 'normalized';

h_check_radius_t=uicontrol('style','edit','string', ...
    num2str(round(Settings_Cor.radius_frac*100,1)),...
    'position',[xstart+3*xsep ystart+ysep xwid yhig]);
h_check_radius_t.Units = 'normalized';

h_bg_update=uicontrol('style','pushbutton','string','Update BG',...
    'position',[xstart ystart+5*ysep xwid yhig]);
h_bg_update.Units='normalized';
h_bg_update.Callback=@Update_BG;
h_bg_update.UserData=Settings_Cor;

h_rad_update=uicontrol('style','pushbutton','string','Update Radon',...
    'position',[xstart+xsep ystart+5*ysep xwid yhig]);
h_rad_update.Units='normalized';
h_rad_update.Callback=@Update_Radon;
h_rad_update.UserData=Settings_Rad;

h_phases=uicontrol('style','popupmenu','String',phases_list,...
     'position',[xstart+2*xsep ystart+5*ysep xwid*2 yhig]);
h_phases.Units='normalized';

try
    phase_val=find(logical(1-cellfun('isempty',strfind(phases_list,InputUser.Phase_Input{1}))) == true);
    if isempty(phase_val)
        phase_val=1;
    end
catch
    phase_val=1;
end

h_phases.Value=phase_val(1);

% Initialize the GUI.
% Change units to normalized so components resize
% automatically.
f.Units = 'normalized';

%Create a plot in the axes.

% Assign the GUI a name to appear in the window title.
f.Name = 'AstroEBSD - Pattern Analysis Settings';
% Move the GUI to the center of the screen.


        
%% Radon Controls
h_theta_range=uicontrol('style','text','string', ...
    'Theta Step',...
    'position',[xstart ystart+2*ysep xwid yhig]);
h_theta_range.Units = 'normalized';

% h_theta_range1=uicontrol('style','edit','string', ...
%     num2str(Settings_Rad.theta_range(1)),...
%     'position',[xstart ystart+3*ysep xwid/3 yhig]);
% h_theta_range1.Units = 'normalized';
%
% h_theta_range2=uicontrol('style','edit','string', ...
%     num2str(Settings_Rad.theta_range(2)),...
%     'position',[xstart+xwid/3 ystart+3*ysep xwid/3 yhig]);
% h_theta_range2.Units = 'normalized';

h_theta_range3=uicontrol('style','edit','string', ...
    num2str(Settings_Rad.theta_range(3)),...
    'position',[xstart ystart+3*ysep xwid yhig]);
h_theta_range3.Units = 'normalized';


h_peaks=uicontrol('style','text','string', ...
    'Pks [M R ID]',...
    'position',[xstart+xsep ystart+2*ysep xwid yhig]);
h_peaks.Units = 'normalized';

h_peaks_max=uicontrol('style','edit','string', ...
    num2str(Settings_Rad.max_peaks),...
    'position',[xstart+xsep ystart+3*ysep xwid/3 yhig]);
h_peaks_max.Units = 'normalized';

h_peaks_report=uicontrol('style','edit','string', ...
    num2str(Settings_Rad.num_peak),...
    'position',[xstart+xsep+xwid/3 ystart+3*ysep xwid/3 yhig]);
h_peaks_report.Units = 'normalized';

h_peaks_found=uicontrol('style','text','string', ...
    num2str(-1),... %NEEDS FIXING
    'position',[xstart+xsep+2*xwid/3 ystart+3*ysep xwid/3 yhig]);
h_peaks_found.Units = 'normalized';

h_rho_per=uicontrol('style','text','string', ...
    'Rho Sep %',...
    'position',[xstart+xsep*2 ystart+2*ysep xwid yhig]);
h_rho_per.Units = 'normalized';

h_rho_pere=uicontrol('style','edit','string', ...
    num2str(round(Settings_Rad.rho_search_per*100,2)),...
    'position',[xstart+2*xsep ystart+3*ysep xwid yhig]);
h_rho_pere.Units = 'normalized';

h_rho_min=uicontrol('style','text','string', ...
    'min Rho %',...
    'position',[xstart+xsep*3 ystart+2*ysep xwid yhig]);
h_rho_min.Units = 'normalized';

h_rho_mine=uicontrol('style','edit','string', ...
    num2str(round(Settings_Rad.min_peak_width*100,3)),...
    'position',[xstart+3*xsep ystart+3*ysep xwid yhig]);
h_rho_mine.Units = 'normalized';

h_theta_sep=uicontrol('style','text','string', ...
    'Theta sep',...
    'position',[xstart+xsep*4 ystart+2*ysep xwid yhig]);
h_theta_sep.Units = 'normalized';

h_theta_sepe=uicontrol('style','edit','string', ...
    num2str(Settings_Rad.theta_search_pix),...
    'position',[xstart+4*xsep ystart+3*ysep xwid yhig]);
h_theta_sepe.Units = 'normalized';

set(h_peaks_found,'string',num2str(size(p_norms,2)));

%% Pattern centre boxed

h_pc_x=uicontrol('style','text','string', ...
    'PCx',...
    'position',[xstart+xsep*6 ystart xwid yhig]);
h_pc_x.Units = 'normalized';

h_pc_y=uicontrol('style','text','string', ...
    'PCy',...
    'position',[xstart+xsep*7 ystart xwid yhig]);
h_pc_y.Units = 'normalized';

h_pc_z=uicontrol('style','text','string', ...
    'PCz',...
    'position',[xstart+xsep*8 ystart xwid yhig]);
h_pc_z.Units = 'normalized';

h_pc_xe=uicontrol('style','edit','string', ...
    num2str(round(PC_start(1),3)),...
    'position',[xstart+xsep*6 ystart+ysep xwid yhig],'Callback',@Update_PC);
h_pc_xe.Units = 'normalized';
h_pc_xe.UserData=PC_start;

h_pc_ye=uicontrol('style','edit','string', ...
    num2str(round(PC_start(2),3)),...
    'position',[xstart+xsep*7 ystart+ysep xwid yhig],'Callback',@Update_PC);
h_pc_ye.Units = 'normalized';

h_pc_ze=uicontrol('style','edit','string', ...
    num2str(round(PC_start(3),3)),...
    'position',[xstart+xsep*8 ystart+ysep xwid yhig],'Callback',@Update_PC);
h_pc_ze.Units = 'normalized';

h_dpc_x=uicontrol('style','text','string', ...
    'dPCx',...
    'position',[xstart+xsep*6 ystart+2*ysep xwid yhig],'Callback',@Update_PC);
h_dpc_x.Units = 'normalized';

h_dpc_y=uicontrol('style','text','string', ...
    'dPCy',...
    'position',[xstart+xsep*7 ystart+2*ysep xwid yhig],'Callback',@Update_PC);
h_dpc_y.Units = 'normalized';

h_dpc_z=uicontrol('style','text','string', ...
    'PCz',...
    'position',[xstart+xsep*8 ystart+2*ysep xwid yhig],'Callback',@Update_PC);
h_dpc_z.Units = 'normalized';

h_dpc_xe=uicontrol('style','edit','string', ...
    num2str(round(PC_range(1),3)),...
    'position',[xstart+xsep*6 ystart+3*ysep xwid yhig],'Callback',@Update_PC);
h_dpc_xe.Units = 'normalized';

h_dpc_ye=uicontrol('style','edit','string', ...
    num2str(round(PC_range(2),3)),...
    'position',[xstart+xsep*7 ystart+3*ysep xwid yhig],'Callback',@Update_PC);
h_dpc_ye.Units = 'normalized';

h_dpc_ze=uicontrol('style','edit','string', ...
    num2str(round(PC_range(3),3)),...
    'position',[xstart+xsep*8 ystart+3*ysep xwid yhig],'Callback',@Update_PC);
h_dpc_ze.Units = 'normalized';


h_pc_search=uicontrol('style','pushbutton','string','Search PC',...
    'position',[xstart+xsep*6 ystart+5*ysep xwid yhig]);
h_pc_search.Units='normalized';
h_pc_search.Callback=@Find_PC;
h_pc_search.UserData=PC_range;

h_pc_index=uicontrol('style','pushbutton','string','Index',...
    'position',[xstart+xsep*7 ystart+5*ysep xwid yhig]);
h_pc_index.Units='normalized';
h_pc_index.Callback=@Index_Pat;


%%
%PC
xsep=80; xstart=550; xwid=150; ystart=300; ysep=-18; yhig=15;
ebox_pc=uicontrol('style','edit','string',' ','position',[xstart ystart xwid yhig],'Visible','off');
ebox_mae=uicontrol('style','edit','string',' ','position',[xstart ystart+ysep xwid yhig],'Visible','off');
ebox_indband=uicontrol('style','edit','string',' ','position',[xstart ystart+2*ysep xwid yhig],'Visible','off');
ebox_repband=uicontrol('style','edit','string',' ','position',[xstart ystart+3*ysep xwid yhig],'Visible','off');
ebox_maxband=uicontrol('style','edit','string',' ','position',[xstart ystart+4*ysep xwid yhig],'Visible','off');

ebox_pc.Units = 'normalized';
ebox_mae.Units = 'normalized';
ebox_indband.Units = 'normalized';
ebox_repband.Units = 'normalized';
ebox_maxband.Units = 'normalized';
        
Update_PC
Update_BG
Update_Radon
        
%% Make the GUI visible.
movegui(f,'center')
drawnow;
f.Visible = 'on';

disp('AstroEBSD GUI is active')
disp('Close the window or use CRTL + C to return to the command line');
uiwait(f,300);
%% sub functions

    function Update_BG(~,eventdata)
        Settings_Cor.hotpixel=get(h_check_hot,'value');
        Settings_Cor.hot_thresh=str2double(get(h_check_hotthres,'string'));
        
        Settings_Cor.gfilt=get(h_check_gauss,'value');
        Settings_Cor.gfilt_s=str2double(get(h_check_gfilt_s,'string'));
        
        %extract the size dta
        Settings_Cor.resize=get(h_check_resize,'value');
        size_str=get(h_check_resize_t,'string');
        size_str_sp=strfind(size_str,' ');
        Settings_Cor.size=[str2double(size_str(1:size_str_sp(1)-1)),str2double(size_str(size_str_sp(2)+1:end))];
        
        Settings_Cor.RealBG=get(h_check_realbg,'value');
        
        Settings_Cor.radius=get(h_check_radius,'value');
        Settings_Cor.radius_frac=round(str2double(get(h_check_radius_t,'string')))/100;
        
        set(h_bg_update,'Userdata',Settings_Cor);
        
        EBSP_single=get(i_raw,'UserData');
        
        [EBSP_single.ebsp_cor,EBSP_single.PatternInfo ] = EBSP_BGCor( EBSP_single.PatternIn,Settings_Cor );
        
        [ EBSP_single.Peak_Centre,EBSP_single.Single.Peak_Set_All,EBSP_single.Peak_Set_All,...
            EBSP_single.R_EBSP,EBSP_single.R_Edge,EBSP_single.R_rho,EBSP_single.R_theta ] ...
            = EBSP_RadHunt( EBSP_single.ebsp_cor,Settings_Rad);
        %generate the geometry
        PC_start=get(h_pc_xe,'UserData');
        [ EBSD_geometry ] = EBSP_Gnom( EBSP_single.PatternInfo,PC_start );
        [ nhat_gnom] = EBSP_NormConv( EBSP_single.Peak_Centre,EBSP_single.PatternInfo.size,PC_start);
        
        
        i_raw=imagesc(EBSP_single.PatternIn,'Parent',h_raw);
        i_cor=imagesc(EBSP_single.ebsp_cor,'Parent',h_cor);
        i_rad=imagesc(EBSP_single.R_theta,EBSP_single.R_rho,EBSP_single.R_EBSP,'Parent',h_rad);
        
        h_raw.Units='pix'; axis(h_raw,'equal'); axis(h_raw,'tight','xy'); h_raw.Units = 'normalized';
        h_cor.Units='pix'; axis(h_cor,'equal'); axis(h_cor,'tight','xy'); h_cor.Units = 'normalized';
        h_rad.Units='pix'; axis(h_rad,'equal'); axis(h_rad,'tight','xy'); h_rad.Units = 'normalized';
        caxis(h_rad,[nanmin(EBSP_single.R_EBSP(:)) nanmax(EBSP_single.R_EBSP(:))])

        hold(h_rad,'on');
        
        delete(s_rad);
        s_rad=scatter(EBSP_single.Peak_Centre(:,1),EBSP_single.Peak_Centre(:,2),'r','Parent',h_rad);
        
        set(i_raw,'UserData',EBSP_single);
        
        for n=1:size(p_norms,2)
            if ishandle(p_norms(n)) && p_norms(n)~= 0
                delete(p_norms(n));
            end
        end
        
        [i_norms,p_norms]=DrawNorms(EBSP_single.ebsp_cor,nhat_gnom,EBSD_geometry,h_norms);
        set(h_peaks_found,'string',num2str(size(p_norms,2)));
                
        h_raw.Units='pix'; axis(h_raw,'equal'); axis(h_raw,'tight'); h_raw.Units = 'normalized';
        h_cor.Units='pix'; axis(h_cor,'equal'); axis(h_cor,'tight'); h_cor.Units = 'normalized';
        h_rad.Units='pix'; axis(h_rad,'equal'); axis(h_rad,'tight'); h_rad.Units = 'normalized';
        ylim(h_rad,[EBSP_single.R_rho(1) EBSP_single.R_rho(end)]);
        axis(h_rad,'square') %make the radon square looking

        guidata(f);
    end

    function Update_PC(source,eventdata)
        PC_start(1)=round(str2double(get(h_pc_xe,'string')),3);
        PC_start(2)=round(str2double(get(h_pc_ye,'string')),3);
        PC_start(3)=round(str2double(get(h_pc_ze,'string')),3);
        h_pc_xe.UserData=PC_start;
    end

    function closefun(source,eventdata)
        %run the update functions
        Update_PC
        Update_BG
        Update_Radon
        
        Settings_PCin.start=get(h_pc_xe,'UserData');
        Settings_PCin.range=get(h_pc_search,'UserData');
        
        assignin('base', 'Settings_Cor', Settings_Cor);
        assignin('base', 'EBSP_single', EBSP_single);
        assignin('base', 'Settings_Rad', Settings_Rad);
        assignin('base', 'Settings_PCin', Settings_PCin);
        assignin('base', 'InputUser', InputUser);
        close(gcf)
    end
function Index_Pat(source,eventdata)
        
        Update_PC
        Update_BG
        Update_Radon
        
        set(h_pc_index,'Enable','off');
        drawnow();
        
        EBSP_single=get(i_raw,'UserData');
        PC_start=get(h_pc_xe,'UserData');
        [ EBSP_single.nhat_gnom] = EBSP_NormConv( EBSP_single.Peak_Centre,EBSP_single.PatternInfo.size,PC_start);
        
        %get the current phase
        phase_list=get(h_phases,'string');
        phase_num=get(h_phases,'value');
        %build the phase
        [ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num ] = Phase_Builder( phase_list(phase_num),InputUser.Phase_Folder );
        %verver the pattern
        [EBSP_single.rotdata,banddata_one]=EBSP_Index(nhat_gnom,Crystal_LUT{1},Settings_LUT{1}.thresh_trig,Crystal_UCell{1},eye(3)); %#ok<PFBNS>
        EBSP_single.PatternGeometry= EBSP_Gnom( EBSP_single.PatternInfo,PC_start );

        EBSP_single.PC=PC_start;
        EBSP_single.Phase= phase_list(phase_num);
        set(i_raw,'UserData',EBSP_single);
        drawnow();
        indexed_info;

        Plot_EBSPAnnotated( EBSP_single.ebsp_cor,EBSP_single.PatternGeometry,EBSP_single.nhat_gnom,EBSP_single.rotdata.detector,Crystal_UCell{1},Crystal_LUT{1}.family_norm_list,h_indexed);
        
        drawnow();
        
        set(h_indexed,'Units','pix'); axis(h_indexed,'equal'); axis(h_indexed,'tight'); set(h_indexed,'Units','normalized');
        set(h_indexed,'XLim',[EBSP_single.PatternGeometry.x_gn_min EBSP_single.PatternGeometry.x_gn_max]);
        set(h_indexed,'YLim',[EBSP_single.PatternGeometry.y_gn_min EBSP_single.PatternGeometry.y_gn_max]);
        
        set(h_indexed,'Visible','on');
        
        drawnow();
        
        set(h_pc_index,'Enable','on');

       
end

    function indexed_info
        EBSP_single=get(i_raw,'UserData');
        
        xsep=80; xstart=550; xwid=150; ystart=300; ysep=-18; yhig=15;
        
        %PC
        PC_str=['PC = [' num2str(round(EBSP_single.PC(1),3)) ',' num2str(round(EBSP_single.PC(2),3)) ',' num2str(round(EBSP_single.PC(3),3)) ']' ];
        error_str=['MAE = ' num2str(round(EBSP_single.rotdata.error*180/pi,4)) ' deg'];
        indband_str=['Indexed Bands = ' num2str(EBSP_single.rotdata.maxok)];
        repband_str=['Bands Found = ' get(h_peaks_max,'string')];
        maxband_str=['Max Bands = ' get(h_peaks_found,'string')];
        
        set(ebox_pc,'String',PC_str);
        set(ebox_mae,'String',error_str);
        set(ebox_indband,'String',indband_str);
        set(ebox_repband,'String',repband_str);
        set(ebox_maxband,'String',maxband_str);
        
        set(ebox_pc,'Visible','on');
        set(ebox_mae,'Visible','on');
        set(ebox_indband,'Visible','on');
        set(ebox_repband,'Visible','on');
        set(ebox_maxband,'Visible','on');
        
        drawnow()
        
    end

    function Find_PC(source,eventdata)
        
        set(h_pc_search,'Enable','off');
        drawnow();
        PC_range(1)=round(str2double(get(h_dpc_xe,'string')),3);
        PC_range(2)=round(str2double(get(h_dpc_ye,'string')),3);
        PC_range(3)=round(str2double(get(h_dpc_ze,'string')),3);
        h_pc_search.UserData=PC_range;
        PC_start=get(h_pc_xe,'UserData');
  
        %set up the search
        PC_GA_options = optimoptions('ga');
        PC_GA_options.FunctionTolerance=1E-3;
        PC_GA_options.UseParallel=false;
        PC_GA_options.MaxGenerations=15;
        PC_GA_options.PopulationSize=30;
        PC_GA_options.MaxStallGenerations=20;
        PC_GA_options.Display='iter';
        
        PC_GA_ub=PC_start+PC_range;
        PC_GA_lb=PC_start-PC_range;
        
        
        EBSP_single=get(i_raw,'UserData');

        %get the current phase
        phase_list=get(h_phases,'string');
        phase_num=get(h_phases,'value');
        %build the phase
        [ Crystal_UCell,Crystal_Family,Crystal_LUT,Settings_LUT,Phase_Num ] = Phase_Builder( phase_list(phase_num),InputUser.Phase_Folder );
   
        
        FitFunc = @(PC_test) PC_GAOpt( PC_test,EBSP_single.Peak_Centre,EBSP_single.PatternInfo.size,Crystal_LUT,Crystal_UCell,1);
        [PC_out, EBSP_single.PC_err] = ga(FitFunc, 3, [], [], [], [],PC_GA_lb, PC_GA_ub,[],PC_GA_options);
        
        set(h_pc_xe,'UserData',PC_out);        
        set(h_pc_xe,'string',num2str(round(PC_out(1),3)));
        set(h_pc_ye,'string',num2str(round(PC_out(2),3)));
        set(h_pc_ze,'string',num2str(round(PC_out(3),3)));                
        
        %need to redraw the gnom
        [ EBSD_geometry ] = EBSP_Gnom( EBSP_single.PatternInfo,PC_out );
        [ nhat_gnom] = EBSP_NormConv( EBSP_single.Peak_Centre,EBSP_single.PatternInfo.size,PC_out);
        
        set(i_raw,'UserData',EBSP_single);
        [i_norms,p_norms]=DrawNorms(EBSP_single.ebsp_cor,nhat_gnom,EBSD_geometry,h_norms);
        
        set(h_pc_search,'Enable','on');
        drawnow();
        
        Index_Pat;
        drawnow();
        
    end

    function Update_Radon(source,eventdata)
        
        Settings_Cor=get(h_bg_update,'Userdata');
        
        %get variables
        Settings_Rad.max_peaks=str2double(get(h_peaks_max,'string'));
        Settings_Rad.num_peak=str2double(get(h_peaks_report,'string'));
        Settings_Rad.rho_search_per=str2double(get(h_rho_pere,'string'))/100;
        Settings_Rad.min_peak_width=str2double(get(h_rho_mine,'string'))/100;
        Settings_Rad.theta_search_pix=str2double(get(h_theta_sepe,'string'));
        %         Settings_Rad.theta_range(1)=str2double(get(h_theta_range1,'string'));
        %         Settings_Rad.theta_range(2)=str2double(get(h_theta_range2,'string'));
        Settings_Rad.theta_range(3)=str2double(get(h_theta_range3,'string'));
               
        %update the analysis
        EBSP_single=get(i_raw,'UserData');
        
        %update the radon settings
        [EBSP_single.ebsp_cor,EBSP_single.PatternInfo ] = EBSP_BGCor( EBSP_single.PatternIn,Settings_Cor );
        [ EBSP_single.Peak_Centre,EBSP_single.Single.Peak_Set_All,EBSP_single.Peak_Set_All,...
            EBSP_single.R_EBSP,EBSP_single.R_Edge,EBSP_single.R_rho,EBSP_single.R_theta ] ...
            = EBSP_RadHunt( EBSP_single.ebsp_cor,Settings_Rad);
        
        %generate the geometry
        [ EBSD_geometry ] = EBSP_Gnom( EBSP_single.PatternInfo,PC_start );
        [ nhat_gnom] = EBSP_NormConv( EBSP_single.Peak_Centre,EBSP_single.PatternInfo.size,PC_start);
        
        i_rad=imagesc(EBSP_single.R_theta,EBSP_single.R_rho,EBSP_single.R_EBSP,'Parent',h_rad);
        ylim(h_rad,[EBSP_single.R_rho(1) EBSP_single.R_rho(end)]);
        caxis(h_rad,[nanmin(EBSP_single.R_EBSP(:)) nanmax(EBSP_single.R_EBSP(:))])
        hold(h_rad,'on');
        
        delete(s_rad);
        s_rad=scatter(EBSP_single.Peak_Centre(:,1),EBSP_single.Peak_Centre(:,2),'r','Parent',h_rad);
        
        set(i_raw,'UserData',EBSP_single);
        
        for n=1:size(p_norms,2)
            if ishandle(p_norms(n)) && p_norms(n)~= 0
                delete(p_norms(n));
            end
        end
        [i_norms,p_norms]=DrawNorms(EBSP_single.ebsp_cor,nhat_gnom,EBSD_geometry,h_norms);
        set(h_peaks_found,'string',num2str(size(p_norms,2)));
                
        h_raw.Units='pix'; axis(h_raw,'equal'); axis(h_raw,'tight'); h_raw.Units = 'normalized';
        h_cor.Units='pix'; axis(h_cor,'equal'); axis(h_cor,'tight'); h_cor.Units = 'normalized';
        h_rad.Units='pix'; axis(h_rad,'equal'); axis(h_rad,'tight'); h_rad.Units = 'normalized';
        axis(h_rad,'square') %make the radon square looking

        guidata(f);
    end

    function [i1,p1]=DrawNorms(EBSD_Pattern_norm,nhat_gnom,EBSD_Geometry_in,s1)
        x2_gnom=zeros(size(nhat_gnom,1),2);
        y2_gnom=zeros(size(nhat_gnom,1),2);
        
        for n=1:size(nhat_gnom,1)
            x2_gnom(n,:)=[EBSD_Geometry_in.x_screen(1) EBSD_Geometry_in.x_screen(end)];
            y2_gnom(n,:)=(-x2_gnom(n,:).*nhat_gnom(n,1)-nhat_gnom(n,3))./nhat_gnom(n,2);
        end
        
        i1=imagesc(EBSD_Geometry_in.x_screen,EBSD_Geometry_in.y_screen,EBSD_Pattern_norm,'Parent',s1);
        p1=[];
        for n=1:size(nhat_gnom)
            hold on;
            %     plot(x1_gnom(n,:),y1_gnom(n,:),'r');
            p1(n)=plot(x2_gnom(n,:),y2_gnom(n,:),'-','color','r','LineWidth',2,'Parent',s1);
            % plot(x2_gnom(n,:),y2_gnom(n,:),'.','color','w','LineWidth',2,'Parent',s1);
        end
        
        %plot the PC
        p1(end+1)=scatter(0,0,15,'w','filled','Parent',s1);
        p1(end+2)=scatter(0,0,20,'m','Parent',s1);
        
        caxis(s1,[nanmin(EBSD_Pattern_norm(:)) nanmax(EBSD_Pattern_norm(:))]);
        s1.Units='pix'; axis(s1,'equal','xy'); s1.Units = 'normalized';
        s1.XLim=[EBSD_Geometry_in.x_gn_min EBSD_Geometry_in.x_gn_max];
        s1.YLim=[EBSD_Geometry_in.y_gn_min EBSD_Geometry_in.y_gn_max];

    end

    function pha_names=phasefolder_read(InputUser)
        if exist(InputUser.Phase_Folder,'dir') ~=7
            error(['The folder ' InputUser.Phase_Folder ' does not exist']);
        end
        
        folder_dir = dir(InputUser.Phase_Folder);
        [~,idx] = sort([folder_dir.datenum]);
        folder_dir=folder_dir(idx);
        folder_files={folder_dir.name};
       
        %find the pha names
        pha_names=folder_files(logical(1-cellfun('isempty',strfind(folder_files,'pha'))));
        for n=1:size(pha_names,2)
            pha_names{n}=pha_names{n}(1:end-4);
        end
 
    end

end
