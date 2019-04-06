%% Simulate 3D SNN reservoir in Liquid State Machine (LSM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulates reservoir in a LSM and is useful in :                         
%   * tuning the synaptic weights to get various dynamics(oscillatory,
%   decaying, saturation)
%   * observe effect of changing individual network parameter to find
%   sensitivity of parameters
%
% Author - Ajinkya Gorad 
% LSMReservoir@2019 (version-0.8)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function LSMReservoirApp
%% Reset random number generator
rng(0);
%% Default configuration settings (could be modified in GUI)
STDP = 0;       % Spike Timing Dependent Plasticity
SAVE_VIDEO=0;   % Save the simulation frame by frame for each param setting
NOISE = 0;      % noise in membrane potential of neurons
TWO_SPIKES = 0; % Only give two input spikes
POISSON =0;     % Poisson spike inputs
%% Configuration parameters
eta = 0;        % noise level (NOISE only)
lambdaIn = 10;  % Poisson input spike rate (POISSON only)
DeltaT = 0;     % Separation between two spike input (TWO_SPIKES only)
%% Reservoir parameters
W0 = [3 6;-2 -2];   % Initial weights in the reservoir (E->E, E->I; I->E, I->I) Excitatory(E), Inhibitory(I)
K0 = [0.45 0.3;0.6 0.15];   % Probability coefficients for reservoir connections
Wee  = W0(1,1); Wei= W0(1,2); Wie= W0(2,1); Wii = W0(2,2); % assign individual weights
Wm =50;              % Maximum limit for reservoir weight synaptic scaling (SLIDER)
alpha_G = 1;            % default synaptic scaling parameter
resSize= [5 5 5];    % 3D reservoir size
r0 = 2;              % effective connectivity distance parameter for gaussian profile
f_inhibit=0.2;       % fraction of inhibitory neurons in reservoir
d_syn = 1E-3;        % Synaptic Delay in reservoir (FIXED)
%%
if(SAVE_VIDEO)
    writerObj = VideoWriter('media\STOCH_RES2.mp4');
    open(writerObj);
end


%% Create 3D reservoir network
% createNetworkM creates an multidimensional graph on N-D grid with
% dimensions(M) = lenght(resSize) and neuron at each ND-point (n1,...,nM)
% X : source node of network graph
% Xn : destination node of network graph
% T : synaptic delay for each neuron (FIXED)
% W : weight of each synapse
% R : coordinates of each neuron
% E : value for choosing excitatory(1)/inhibitory(-1) neuron

[X,Xn,T,W,R,E] = createNetworkM(resSize,W0,r0,K0,f_inhibit,d_syn); % gaussian probablistically connected network

%% Get secondary parameters
Nres = length(E);   % number of reservoir neurons
V = zeros(Nres,1);  % Membrane potential of reservoir Neurons
v = zeros(Nres,3);  % Synaptic traces
I = zeros(Nres,3);  % Neuron input current
spiked = logical(zeros(Nres,1)); % to store spiked neurons at an instant
spikedLog =logical(sparse([])); % to store spiked neurons in all simulated time <K
G = sparse(X,Xn,W); % Sparse matrix is better to store loosely connected networks
G0=G;    % save default configuration for network weights and connectivity

% STDP definitions
dGp = -1; dGn = 1; % weight/conductance (G) update rate
A = sparse(zeros(size(G)));
A(E>0,E>0) = 1;A(E>0,E<0) = 1; % Define coefficients for STDP
A(E<0,E>0) = 1;A(E<0,E<0) = 1; %  -;;-
A = A.*logical(G);  % get the graph connectivity matrix A

%% Input mapping of spikes
Gin=ones(Nres,1)*8;     % Defines input weights for each neuron
Iapp = zeros(Nres,1);   % Externally applied current for each neuron
inFrac=0.9;             % Fraction of reservoir  neurons which take input
%% Values from the Literature : Zhang et al. 2015
dt = 1E-3; % Define time step of simulation
K = 2000;  % Define time of simulation
alphaV = (1-dt/20E-3); % membrane decay constant (discrete)
alpha1E = (1-dt/40E-3);alpha2E =(1-dt/6E-3); % synaptic trace decay constants
alpha3stdp = (1-dt/10E-3); % used in STDP
alpha1I = (1-dt/30E-3); alpha2I = (1-dt/4E-3);
% create and assign matrices
alpha1 = zeros(Nres,1); alpha2 = alpha1; alpha3 = alpha1;
Kt0e=0;Kt0i=0;Kt0=zeros(Nres,1);
alpha1(E>0) = alpha1E; alpha2(E>0) = alpha2E; alpha3(:) = alpha3stdp;
alpha1(E<0) = alpha1I; alpha2(E<0) = alpha2I;
Vth = 20;       % Threshold voltage (mV)

%% Update synaptic scaling
%updateKt0();

%% Create figure handles and GUI
fig_handle = figure('name','simulation'); % Main figure window
set(fig_handle, 'Position', [0, 0, 1280, 600]); % Defines the size of main window
text_subplot = axes('Position',[0 0 0.2 1],'Visible','off');infoText=text(0.05,0.65,' ','FontSize',10); % Information region
raster_subplot = axes('Position',[.2 0.1 .55 .35]);     % reservoir spikes
input_subplot = axes('Position',[.2 0.6 .55 .35]);      % input Spikes (to each neuron)
psp_subplot = axes('Position',[.8 0.6 .15 .35]);        % Synaptic model/Reservoir distribution
activity_subplot = axes('Position',[.8 0.1 .15 0.35]);  % Spike activity with time
glog_subplot = axes('Position',[0.025 0.8 0.15 0.2],'Visible','off'); % Weight update (STDP only)
set(gcf,'Color','w');
%% Create sliders
axes(text_subplot);
slidebarW=200; % width in pixel for slider
slhan=uicontrol('style','slider','String','alpha_G','TooltipString','alpha_G','position',[0 0 slidebarW 10],'min',0,'max',1);
slhan2=uicontrol('style','slider','String','inFrac','TooltipString','inFrac','position',[0 10 slidebarW 10],'min',0,'max',1);
slhan3=uicontrol('style','slider','String','alpha1E','TooltipString','alpha1E','position',[0 20 slidebarW 10],'min',0,'max',1);
slhan4=uicontrol('style','slider','String','alpha2E','TooltipString','alpha2E','position',[0 30 slidebarW 10],'min',0,'max',1);
slhan5=uicontrol('style','slider','String','alphaV','TooltipString','alphaV','position',[0 40 slidebarW 10],'min',0,'max',1);
slhan6=uicontrol('style','slider','String','alpha1I','TooltipString','alpha1I','position',[0 50 slidebarW 10],'min',0,'max',1);
slhan7=uicontrol('style','slider','String','alpha2I','TooltipString','alpha2I','position',[0 60 slidebarW 10],'min',0,'max',1);
slhan8=uicontrol('style','slider','String','Vth','TooltipString','Vth','position',[0 70 slidebarW 10],'min',0,'max',1);
slhan9=uicontrol('style','slider','String','Wee','TooltipString','Wee','position',[0 80 slidebarW 10],'min',0,'max',1,'SliderStep',[1/10/50 2/50],'BackgroundColor',[0.0 0.5 0.2],'Value',Wee/Wm);
slhan10=uicontrol('style','slider','String','Wei','TooltipString','Wei','position',[0 90 slidebarW 10],'min',0,'max',1,'SliderStep',[1/10/50 2/50],'BackgroundColor',[0.2 0.5 0.2],'Value',Wei/Wm);
slhan11=uicontrol('style','slider','String','Wie','TooltipString','Wie','position',[0 100 slidebarW 10],'min',0,'max',1,'SliderStep',[1/10/50 2/50],'BackgroundColor',[0.5 0.2 0.2],'Value',-Wie/Wm);
slhan12=uicontrol('style','slider','String','Wii','TooltipString','Wii','position',[0 110 slidebarW 10],'min',0,'max',1,'SliderStep',[1/10/50 2/50],'BackgroundColor',[0.5 0.0 0.2],'Value',-Wii/Wm);
slhan13=uicontrol('style','slider','String','Noise','TooltipString','Noise','position',[0 120 slidebarW 10],'min',0,'max',1);
slhan14=uicontrol('style','slider','String','lambdaIn','TooltipString','lambdaIn','position',[0 130 slidebarW 10],'min',0,'max',1);
slhan15=uicontrol('style','slider','String','DeltaT','TooltipString','DeltaT','position',[0 140 slidebarW 10],'min',0,'max',1, 'SliderStep',[1/800 50/800]);

%% Create checkboxes
endsliderY = 140; % width in pixel for checkbox
chkbx1 = uicontrol('style','checkbox','String','STDP','position',[0 endsliderY+10 80 10],'Value',STDP,'callback',@ifSTDP);
chkbx2 = uicontrol('style','checkbox','String','PoissonTrain','position',[0 endsliderY+20 80 10],'Value',POISSON,'callback',@ifPoisson);
chkbx2 = uicontrol('style','checkbox','String','Noise','position',[0 endsliderY+30 80 10],'Value',NOISE,'callback',@ifNoise);
chkbx2 = uicontrol('style','checkbox','String','2Spike','position',[0 endsliderY+40 80 10],'Value',TWO_SPIKES,'callback',@if2Spikes);
pb = uicontrol('Style','pushbutton','Callback',@pb_callback,'position',[0 endsliderY+50 80 10]);
%% Text information string
infoStringTmpl ={'\alpha_G = <alpha_G>'
    'frac_{in} = <inFrac>'
    '\tau_{1e} = <alpha1E>'
    '\tau_{2e} = <alpha2E>'
    '\tau_{V} = <alphaV>'
    '\tau_{1i} = <alpha1I>'
    '\tau_{2i} = <alpha2I>'
    'V_{th} = <Vth>'
    '{\eta} = <eta>'
    '{\lambda_{in}} = <lambdaIn>'
    '{\DeltaT} = <DeltaT>'
    '-------------------'
    'W_{EE} = <Wee>'
    'W_{EI} = <Wei>'
    'W_{IE} = <Wie>'
    'W_{II} = <Wii>'};
% slhan6=uicontrol('style','slider','String','wE->E','position',[0 0 100 10],'min',0,'max',1);
% slhan7=uicontrol('style','slider','String','inFrac','position',[100 0 100 10],'min',0,'max',1);
% slhan8=uicontrol('style','slider','String','alpha1','position',[0 10 100 10],'min',0,'max',1);
% slhan9=uicontrol('style','slider','String','alpha2','position',[100 10 100 10],'min',0,'max',1);
%% Create listener to push buttons
hSliderListener = addlistener(slhan, 'Value','PostSet',@updatealpha1);
hSliderListener2 = addlistener(slhan2, 'Value','PostSet',@updatealpha2);
hSliderListener3 = addlistener(slhan3, 'Value','PostSet',@updatealpha3);
hSliderListener4 = addlistener(slhan4, 'Value','PostSet',@updatet4);
hSliderListener5 = addlistener(slhan5, 'Value','PostSet',@updatet5);
hSliderListener6 = addlistener(slhan6, 'Value','PostSet',@updatet6);
hSliderListener7 = addlistener(slhan7, 'Value','PostSet',@updatet7);
hSliderListener8 = addlistener(slhan8, 'Value','PostSet',@updatet8);
hSliderListener9 = addlistener(slhan9, 'Value','PostSet',@updatet9);
hSliderListener10 = addlistener(slhan10, 'Value','PostSet',@update10);
hSliderListener11 = addlistener(slhan11, 'Value','PostSet',@update11);
hSliderListener12 = addlistener(slhan12, 'Value','PostSet',@update12);
hSliderListener13 = addlistener(slhan13, 'Value','PostSet',@update13);
hSliderListener14 = addlistener(slhan14, 'Value','PostSet',@update14);
hSliderListener15 = addlistener(slhan15, 'Value','PostSet',@update15);
%hSliderListener3 = addlistener(slhan3, 'Value','PostSet',@updatealpha3);

%% Run main code
updateKt0();   % updatepspplot();
updateplot();

%% Helper Functions
% Below contains the callback and simulation functions 
% They are evaluated based on the GUI response
%% Generate New Network
    function pb_callback(source,eventdata)
        % Generate new network with different random seed
        rng('shuffle');
        [X,Xn,T,W,R,E] = createNetworkM([5,5,5],W0,2,K0,0.2,1E-3); % spatially connected network
        %[X,Xn,Tau,W,R,E,N] = p2N5etwork(200,0.2,Wee,Wie);
        Nres = length(E);
        V = zeros(Nres,1);
        v = zeros(Nres,3);
        I = zeros(Nres,3);
        spiked = logical(zeros(Nres,1));
        G = sparse(X,Xn,W);G0=G;
        
        % STDP definitions
        dGp = -1; dGn = 1;
        A = sparse(zeros(size(G)));
        A(E>0,E>0) = 1;A(E>0,E<0) = 1;
        A(E<0,E>0) = 1;A(E<0,E<0) = 1;
        A = A.*logical(G);
        updateplot();
    end
%% Update config values
    function ifSTDP(source,eventdata)
        STDP = source.Value;
        updateplot();
    end
    function ifPoisson(source,eventdata)
        POISSON = source.Value;
        updateplot();
    end
    function ifNoise(source,eventdata)
        NOISE = source.Value;
        updateplot();
    end
    function if2Spikes(source,eventdata)
        TWO_SPIKES = source.Value;
        updateplot();
    end
    function updatealpha1(hObject,eventdata)
        alpha_G = get(eventdata.AffectedObject, 'Value');
        updateplot();
    end
    function updatealpha2(source,eventdata)
        inFrac = get(eventdata.AffectedObject, 'Value');
        updateplot();
    end
    function updatealpha3(source,eventdata)
        alpha1E = get(eventdata.AffectedObject, 'Value');
        updateKt0();    updatepspplot();    updateplot();
    end
    function updateKt0()
        Kt0e = (1-alpha1E)*(1-alpha2E)/(alpha1E-alpha2E)/dt; % normalize by area
        Kt0i = (1-alpha1I)*(1-alpha2I)/(alpha1I-alpha2I)/dt;
        Kt0(E>0) = Kt0e;Kt0(E<0) = Kt0i;
    end
    function updatet4(source,eventdata)
        alpha2E = get(eventdata.AffectedObject, 'Value');
        updateKt0();    updatepspplot();    updateplot();
    end
    function updatet5(source,eventdata)
        alphaV = get(eventdata.AffectedObject, 'Value');
        updatepspplot();
        updateplot();
    end

    function updatet6(source,eventdata)
        alpha1I = get(eventdata.AffectedObject, 'Value');
        updateKt0();    updatepspplot();    updateplot();
    end
    function updatet7(source,eventdata)
        alpha2I = get(eventdata.AffectedObject, 'Value');
        updateKt0();    updatepspplot();    updateplot();
    end
    function updatet8(source,eventdata)
        Vth = 20*get(eventdata.AffectedObject, 'Value');
        updateKt0();    updateplot();
    end
    function updatet9(source,eventdata)
        Wee = Wm*get(eventdata.AffectedObject, 'Value');
        G0(E>0,E>0) = Wee*(A(E>0,E>0)~=0);  W0(1,1) = Wee;
        updateplot();
    end
    function update10(source,eventdata)
        Wei =  Wm*get(eventdata.AffectedObject, 'Value');
        G0(E>0,E<0) = Wei*(A(E>0,E<0)~=0);  W0(1,2) = Wei;
        updateplot();
    end
    function update11(source,eventdata)
        Wie = -Wm*get(eventdata.AffectedObject, 'Value');
        G0(E<0,E>0) = Wie*(A(E<0,E>0)~=0);  W0(2,1) = Wie;
        updateplot();
    end
    function update12(source,eventdata)
        Wii = -Wm*get(eventdata.AffectedObject, 'Value');
        G0(E<0,E<0) = Wii*(A(E<0,E<0)~=0);  W0(2,2) = Wii;
        updateplot();
    end
    function update13(source,eventdata)
        eta = Vth*get(eventdata.AffectedObject, 'Value');
        updateplot();
    end
    function update14(source,eventdata)
        lambdaIn = 1/5/dt*get(eventdata.AffectedObject, 'Value');
        updateplot();
    end
    function update15(source,eventdata)
        DeltaT = 800*get(eventdata.AffectedObject, 'Value'); % in ms
        DeltaT = ceil(DeltaT-400); % make it integer
        updateplot();
    end

    function updateInfoString()
        infoString = infoStringTmpl;
        s = sprintf('%.4f' ,alpha_G); infoString = strrep(infoString,'<alpha_G>',s);
        s = sprintf('%.4f' ,inFrac); infoString = strrep(infoString,'<inFrac>',s);
        s = sprintf('%.4f' ,alpha1E); infoString = strrep(infoString,'<alpha1E>',s);
        s = sprintf('%.4f' ,alpha2E); infoString = strrep(infoString,'<alpha2E>',s);
        s = sprintf('%.4f' ,alphaV); infoString = strrep(infoString,'<alphaV>',s);
        s = sprintf('%.4f' ,alpha1I); infoString = strrep(infoString,'<alpha1I>',s);
        s = sprintf('%.4f' ,alpha2I); infoString = strrep(infoString,'<alpha2I>',s);
        s = sprintf('%.4f' ,Vth); infoString = strrep(infoString,'<Vth>',s);
        s = sprintf('%.4f' ,Wee); infoString = strrep(infoString,'<Wee>',s);
        s = sprintf('%.4f' ,Wei); infoString = strrep(infoString,'<Wei>',s);
        s = sprintf('%.4f' ,Wie); infoString = strrep(infoString,'<Wie>',s);
        s = sprintf('%.4f' ,Wii); infoString = strrep(infoString,'<Wii>',s);
        s = sprintf('%.4f' ,eta); infoString = strrep(infoString,'<eta>',s);
        s = sprintf('%.4f' ,lambdaIn); infoString = strrep(infoString,'<lambdaIn>',s);
        s = sprintf('%.4f' ,DeltaT); infoString = strrep(infoString,'<DeltaT>',s);
        infoText.String = infoString;
        
    end

%% Show synaptic model
    function updatepspplot()
        ve = zeros(1,2);
        vi = zeros(1,2);
        v = zeros(1);
        for k = 1:100
            ve = ve.*[alpha1E alpha2E]+double(sum(k==[1]));
            vi = vi.*[alpha1I alpha2I]+double(sum(k==[1]));
            v = v.*alphaV+double(sum(k==[1]));
            vloge(k) =  Kt0e*ve*[1 -1]';
            vlogi(k) =  Kt0i*vi*[1 -1]';
            vlog(k) = (1-alphaV)/alphaV/dt* v;
        end
        plot(vloge,'b'); hold on;
        plot(vlogi,'r');
        plot(vlog,'g');hold off;
        drawnow;title('PSP & LIF Model(Normalized)');
        xlabel('time(ms'); ylabel('V');
        legend('Excitatory PSP','Inhibitory PSP','LIF Model');
    end
%% Show reservoir spike distribution by neuron
    function showDist()
        % show spike distribution in reservoir
        axes(psp_subplot);
        resDist=sort(sum(1.0*spikedLog,2))/(K*dt);
        h0=fill([1 1:length(resDist) length(resDist)],[0; resDist; 0],'g','FaceAlpha',0.4,'EdgeAlpha',0.1); drawnow;
        xlabel('Sorted Neuron #'); ylabel('Spike Rate(s^{-1}/neuron)');
        title('Reservoir Spike Distribution');
    end
%% Update main plot
% simulates the dynamics and plots the output
    function updateplot()
        figure(fig_handle);
        s = sprintf("Kx(alpha1E,alpha2E):%.3f(%1.3f,%1.3f) Kx(alpha1I,alpha2I):%.3f(%1.3f,%1.3f) alpha_G:%.5f alphaV:%.3f Vth:%.3f",Kt0e,alpha1E,alpha2E,Kt0i,alpha1I,alpha2I,alpha_G,alphaV,Vth);
        axes(input_subplot);hold off;
        s2 = sprintf("Wee : %.2f Wei : %.2f Wie : %.2f Wii : %.2f\r\n",Wee,Wei,Wie,Wii);
        %text(0,Nres,char(s),'Color','r','FontSize',10);
        fprintf("Plotupdated %s\r\n%s\r\n",s,s2);
        updateInfoString();
        VLog =[];
        GsumLog = zeros(0,2,2);
        
        IappLog=sparse([]);
        V = zeros(Nres,1); % reset valuess
        v = zeros(Nres,3);
        I = zeros(Nres,3);
        G = G0;
        spikedLog(:)=0;
        k_spiked = -Inf*ones(Nres,1);
        inhN = find(E<0);
        rng(0);
        for k =1:K
            Iapp(:)=0;
            %if(k<500)Iapp(1:floor(Nres*inFrac)) = rand(floor(Nres*inFrac),1)<0.05; end;
            %if(k==10 ||(k==250 || k==260)|| (k==500 || k==510))Iapp(1:floor(Nres*inFrac)) =100000; end;%impulse current (all excitatory spike)
            %if((k>=10 && k<200) || (k>=400 && k<600))
            if(TWO_SPIKES)
                if(k==K/2 || k==K/2+DeltaT)
                    Iapp(1:floor(Nres*inFrac)) = 10000*(ones(1,floor(Nres*inFrac)));
                end
            end
            if(POISSON)
                if(k>=10 && k<K/4)
                    Iapp(1:floor(Nres*inFrac)) = 10000*(rand(1,floor(Nres*inFrac))<lambdaIn*dt);
                end
            else if(~TWO_SPIKES)
                    if(k==10 || k==110||k==210||k==310)
                        
                        Iapp(1:floor(Nres*inFrac)) = 10000*(rand(1,floor(Nres*inFrac))>0);
                        if(k==400) rng(0);end
                    end
                end
            end
            Iapp(inhN) = 0; IappLog(:,end+1) = Iapp; % no current to inhibitory neurons, log applied current
            V = alphaV*V+alpha_G*G'*I(:,mod(k-1,2)+1)*dt+Gin.*Iapp*dt+NOISE*eta*randn(Nres,1)*sqrt(dt);
            spikedOld = spiked; spiked = V>Vth;             V(spiked) = 0;V(V<0) = 0;
            spiked((k-k_spiked)<=2)=0;
            k_spiked(spiked)=k; % refractory period
            
            v(:,1:2) = [alpha1 alpha2].*v(:,1:2)+spiked;
            I(:,mod(k,3)+1) = Kt0.*v(:,1:2)*[1 -1]'; % second order synaptic model
            % Update reservoir weights
            if(STDP)
                v(:,3) = alpha3.*v(:,3)+spikedOld;
                Ag = spiked.*v(:,3)';
                G = G+(dGp*Ag+dGn*Ag').*A;
                GsumLog(k,1,1) =sum(sum(G(E>0,E>0)));GsumLog(k,1,2) =sum(sum(G(E>0,E<0)));
                GsumLog(k,2,1) =sum(sum(G(E<0,E>0)));GsumLog(k,2,2) =sum(sum(G(E<0,E<0)));
            end
            VLog(:,end+1)=V;
            spikedLog(:,k)=spiked;
        end
        if(STDP)
            axes(glog_subplot);
            plot(GsumLog(:,1,1),'b'); hold on;plot(GsumLog(:,1,2),'c');
            plot(GsumLog(:,2,1),'m');plot(GsumLog(:,2,2),'r');hold off;
        else
            plot(0,0);
        end
        axes(raster_subplot);
        imagesc(VLog); colorbar;hold on;
        [s,ts]=find(spikedLog.*(E>0));   plot(ts,s,'.','MarkerSize',0.01,'Color','g');
        [s,ts]=find(spikedLog.*(E<0));    plot(ts,s,'.','MarkerSize',0.01,'Color','r'); hold off;
        spikesK = full(sum(spikedLog,1));
        spikesKE= full(sum(spikedLog(E>0,:),1));
        spikesKI =full(sum(spikedLog(E<0,:),1));
        spikesK1=spikesK;
        spikesK = conv(spikesK,ones(1,20)/20); % plot spikes over 2ms window (=refractory period)
        spikesKE = conv(spikesKE,ones(1,20)/20); % plot spikes over 2ms window (=refractory period)
        spikesKI = conv(spikesKI,ones(1,20)/20); % plot spikes over 2ms window (=refractory period)
        xlim([0 k]); ylim([0 Nres]); title('Output Raster and Neuron Potential');
        xlabel('time(ms)'); ylabel('Neuron ID');
        axes(input_subplot);
        spikedLog=logical(spikedLog);
        [s,ts]=find(IappLog);   plot(ts,s,'.','MarkerSize',0.1,'Color','k'); hold on;
        [s,ts]=find(spikedLog.*(E>0));   plot(ts,s,'.','MarkerSize',0.01,'Color','b');
        [s,ts]=find(spikedLog.*(E<0));    plot(ts,s,'.','MarkerSize',0.01,'Color','r'); hold off;
        hold off; xlabel('time(ms)'); ylabel('Neuron ID');
        xlim([0 k]); ylim([0 Nres]);title('Input Spikes & Output Raster');
        axes(activity_subplot);
        plot(spikesKE(1:end),'-o','Color','b','MarkerSize',2); hold on;
        plot(spikesKI(1:end),'-o','Color','r','MarkerSize',2);
        plot(spikesK(1:end),'-o','Color','g','MarkerSize',2); hold off; xlim([1 K]);
        title('Network Activity'); legend('Excitatory','Inhibitory','Total');
        xlabel('time(ms)'); ylabel('# of firing neurons');
        showDist();
        set(gcf,'Color','w');set(findobj(gcf,'type','axes'),'FontName','Consolas','FontSize',10,'FontWeight','Bold', 'LineWidth', 1);
        if(SAVE_VIDEO)
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
        end
    end
end
