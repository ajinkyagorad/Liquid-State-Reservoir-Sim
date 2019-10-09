%% Create reservoir network for Liquid State Machine (LSM)
% according to Zhang's paper for LSM
% pure random/spatial random/ pure spatial
% resSize (reservoir size) : 1x3 dimensional array with entries as k l m
% create a multidimensional network
% w (weight) : is 2x2 matrix with fixed weight for entries {E,I}->{E,I} [E->E E->I;I->E I->I];
function [X,Xn,T,W,R,E]= createNetworkM(resSize,w,r0,k0,f_inhibit,tau,show,UC)
if(nargin<1) resSize = [3 3 15]; end
if(nargin<2) w = [3 6 ;-2 -2]; end
if(nargin<3) r0 = 2; end
if(nargin<4) k0 = [0.45 0.3;0.6 0.15]; end
if(nargin<5) f_inhibit = 0.2; end
if(nargin<6) tau = 1E-3; end
if(nargin<7) show = 0; end
if(nargin<8) UC = 0; end
REMOVE_REFLECTED_LOOPS = 0;
%% generate Nd coordinates of neurons
Ndim = length(resSize);
N = prod(resSize); % number of neurons
R = [1:resSize(1)]';
for i_dim = 2:Ndim
    dim = resSize(i_dim);
    R = repmat(R,dim,1);
    r =[1:dim]'*ones(1,prod(resSize(1:i_dim-1)));
    r = reshape(r',length(r(:)),1);
    R(:,i_dim) = r;
end
%% Assign excitatory/inhibitory behaviour
if(UC(1)==0) % randomly defined E/I neurons
    E = ones(N,1); % E : excitatory (+1)
E(rand(N,1)<f_inhibit) = -1; %  inhibitory(-1)

else % Unit cell based regular E/I positions
    sizeUC=size(UC); if(numel(sizeUC)<3); sizeUC(3)=1;end
    rpt=resSize./sizeUC;rpt=ceil(rpt);
    E=repmat(UC,rpt(1),rpt(2),rpt(3));E=E(1:resSize(1),1:resSize(2),1:resSize(3));
    E=E(:);
end
%% Assign synapses to network
% Pure random/ spatial random/ pure spatial
%XXX : doesn't work in old matlab versions
[~,ver] = version;
if(posixtime(datetime(ver))<posixtime(datetime('March 27, 2017')))
    D = (repmat(R(:,1),1,N)-repmat(R(:,1),1,N)').^2;
    for i = 2:Ndim
        D = D + (repmat(R(:,i),1,N)-repmat(R(:,i),1,N)').^2;
    end
else
    
    D = (R(:,1)-R(:,1)').^2;
    for i = 2:Ndim
        D = D + (R(:,i)-R(:,i)').^2;
    end
end

D = sqrt(D);

%%
nE = sum(E>0);
nI = sum(E<0);
% assign on basis of E/I information
% sort distances based on E
[~,sort_E] =sort(E);
[~,sort_back] = sort(sort_E);
D = D(sort_E,sort_E);
% XXX : (optimise following for performance & memory)
ConnProb = [k0(2,2)*exp(-D(E<0,E<0).^2/r0^2) k0(2,1)*exp(-D(E<0,E>0).^2/r0^2);...
    k0(1,2)*exp(-D(E>0,E<0).^2/r0^2) k0(1,1)*exp(-D(E>0,E>0).^2/r0^2)];
D(ConnProb<rand(N,N)) =0;
clear ConnProb; % XXX : delete all NxN matrices if not needed
D = D(sort_back,sort_back);
%% remove reflected coefficients by choosing one random connection
% DR = D;%for distance  histogram
if(REMOVE_REFLECTED_LOOPS)
DL = logical(D);
[r,c] = find(triu(DL.*DL'));
reflected_loops = length(r); % reflected loops in the network
if(reflected_loops>0)
    loop_id = randperm(reflected_loops);
    partition_id = ceil(reflected_loops/2); % create halfway partition for lower and upper triangular matrix
    upper_loops = loop_id(1:partition_id);
    lower_loops = loop_id(partition_id+1:end);
    
    for i=upper_loops
        D(r(i),c(i)) = 0;
        DR(r(i),c(i))=0;
    end
    for i = lower_loops
        D(c(i),r(i)) = 0;
        D(c(i),r(i)) = 0;
    end
end
end
if(find(D.*D'))
    s = sprintf('Network has %i reflected loops',length(find(D.*D'~=0)));
    warning(s);
end
%D(logical(tril(D.*D'))) = 0; % remove all remaining reflected coefficients
%% Assign weights
%DR(D==0) = []; for distance histogram
X = logical(D);

%% Assign synaptic time delay
if(tau~=0)
    T = tau*logical(X);
else
    T = D*1E-3;
end
T = reshape(T,[1 numel(T)])';
T = T(T~=0);
clear D;
fan_in = sum(X,1); % fan in of neuron
fan_out = sum(X,2); % fan out of neuron
W = sparse(size(X));
W(E>0,E>0) = w(1,1); % changed
W(E>0,E<0) = w(1,2);
W(E<0,E>0) = w(2,1);
W(E<0,E<0) = w(2,2);
W(X==0) = 0;
w_in = sum(W,1); % net input weight to neuron
w_out = sum(W,2); % net output weight of neuron
%% Create arrays out of 2D matrices and reduce them
W = reshape(W,[1 numel(W)])';
W = W(W~=0);
n = 1:N; % id of neuron with  order (((1...K)...L)...M)
% Connections formed are of the form X->Xn
Xn = X*diag(n);
Xn = reshape(Xn,[1 numel(X)])';
Xn = Xn(Xn~=0); % destination neuron

X = diag(n)*X;
X = reshape(X,[1 numel(X)])';
X = X(X~=0);  % input neurons for destination neuron

%% Display out information
fprintf('Reservoir with size %s created %i/%i \r\n',mat2str(resSize),length(X),length(E));

if(nargout<1 || show >=1)
    % Color coding
    Color = zeros(N,3); % color coding of neurons
    Color(find(E>0),:) = ones(length(find(E>0)),1).*[0 0 1]; % excitatory neuronsosn coded as blue
    Color(find(E<0),:) = ones(length(find(E<0)),1).*[1 0 0]; % inhibitory neurons coded red
    % Display Network
    figure('name','Neuron Positions');
    subplot(3,3,[1 2 4 5]);
    scatter3(R(:,1),R(:,2),R(:,3),10,Color); hold on;% R is for color
    if(length(X)<10000)
        h=plot3([R(Xn,1)'; R(X,1)'],[R(Xn,2)' ;R(X,2)'],[R(Xn,3)'; R(X,3)'],'-');
    else
        fprintf('Too many synapses to display\r\n');
    end
    subplot(3,3,3);
    hist(T,1000);
    title('Delays'); ylabel('#Synapses'); xlabel('{\tau_{delay}}(s)');
    subplot(3,3,7);     hist(fan_in,1000);
    ylabel('#Neurons'); xlabel('Fan-in');
    subplot(3,3,8);    hist(fan_out,1000);
    ylabel('#Neurons'); xlabel('Fan-out');
    subplot(3,3,9);
    [N_W,W_edge] = histcounts(w_in,500);
    bar(W_edge(2:end),log10(N_W));
    title('Mean W_{in}'); ylabel('#Neurons'); xlabel('Weight');
    subplot(3,3,6);
    [N_W,W_edge] = histcounts(w_out,500);
    bar(W_edge(2:end),log10(N_W));
    title('Mean W_{out}'); ylabel('#Neurons'); xlabel('Weight');
    drawnow;
    hold off;
end

end