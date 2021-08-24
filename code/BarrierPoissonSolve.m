
show_figs = cpf.figs;

try
    load(poisson_mat_out);
catch ME
    error('Parameterization failed, poisson.mat does not exist');
end

centers = (V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))/3;
nf=length(F);

normals = cross(V(F(:,2),:)-V(F(:,1),:), V(F(:,3),:)-V(F(:,1),:),2);
normals=normals./MESH.normv(normals);
if show_figs
    figure; title('gradients of parameterization functions');
    MESH_VIS.mesh(cpf.M); hold on
    quiver3(centers(:,1), centers(:,2), centers(:,3), ...
        normals(:, 1),normals(:, 2),normals(:, 3));
    hold on
end

% % Solve parameterization in Matlab

B1 = V(F(:,2),:)-V(F(:,1),:);
B1=B1./MESH.normv(B1);
B2 = cross(normals, B1);
B2=B2./MESH.normv(B2);

%creating 3D->2D gradient
IReduc = repmat((1:2*N*nf)',1,3);
JReduc = repmat(reshape((1:3*N*nf)',3,N*nf)', 1,2);
JReduc = reshape(JReduc', 3, 2*N*nf)';
B1N = reshape(repmat(B1', N, 1), 3, N*nf);
B2N = reshape(repmat(B2', N, 1), 3, N*nf);
SReduc = reshape([B1N;B2N], 3, 2*N*nf)';

reducMat = fast_sparse(IReduc, JReduc, SReduc, 2*N*nf, 3*N*nf);
G2 = reducMat*G;


rawFieldGx0 = G*x0;
rawFieldGx0=reshape(rawFieldGx0, 3*N, length(F))';
funcValues = x2CornerMat*x0;
funcValues = reshape(funcValues, N, length(V))';

if show_figs
    for i=0:N-1
        axis equal
        cameratoolbar
        quiver3(centers(:,1), centers(:,2), centers(:,3), ...
            rawFieldGx0(:, i*3+1),rawFieldGx0(:, i*3+2),...
            rawFieldGx0(:, i*3+3));
    end
end

%testing G2
rawFieldGx02 = G2*x0;
rawFieldGx02=reshape(rawFieldGx02, 2*N, length(F))';
rawFieldGx0From2 = zeros(length(F), 3 * N);
for i=0:N-1
    rawFieldGx0From2(:,i*3+1:i*3+3) = rawFieldGx02(:,2*i+1).*B1 +rawFieldGx02(:,2*i+2).*B2;
end
G2Error = max(max(abs(rawFieldGx0From2-rawFieldGx0)));

funcValues = x2CornerMat*x0;
funcValues = reshape(funcValues, N, length(V))';

if show_figs
    for i=0:N-1
        axis equal
        cameratoolbar
        quiver3(centers(:,1), centers(:,2), centers(:,3), rawFieldGx0(:, i*3+1),rawFieldGx0(:, i*3+2),rawFieldGx0(:, i*3+3));
    end
end

% if show_figs
% %figure
% %hold on
% %patch('vertices', funcValues(:,1:2), 'faces', F, 'faceColor', 'w', 'edgeColor','b');
% end

fixedIndices=[];
fixedValues=[];

leftIndices=integerIndices;
wobj=1;
wconst=10e5;
wbarrier=0.0001;
wclose = 0.01;
s=0.1;
xcurr=x0;
fraction=1;

IImagField=repmat((1:N*nf)',1,4);
JImagField=IImagField;
for i=0:nf-1
    varOffset=i*2*N+1;
    barOffset=i*N+1;
    JImagField(barOffset:barOffset+N-1,1:2)=reshape(varOffset:varOffset+2*N-1, 2, N)';
    JImagField(barOffset:barOffset+N-1,3:4)=JImagField([barOffset+1:barOffset+N-1,barOffset],1:2);
end

origFieldVolumes = zeros(nf,N);
field = G2*x0;
for i=0:nf-1
    varOffset=i*2*N+1;
    barOffset=i*N+1;
    faceField = reshape(field(varOffset:varOffset+2*N-1),2,N)';
    faceFieldNext = faceField([2:N,1],:);
    
    origFieldVolumes(i+1,:) = (faceField(:,1).*faceFieldNext(:, 2) - faceField(:,2).*faceFieldNext(:, 1));
end
i=1;
packetSize=1;
roundCoeff=1;
success=1;
while (i<=length(integerIndices))
    % i
    xprev=xcurr;
    roundDiffs= abs(fraction*xcurr(leftIndices)-roundCoeff*round(fraction*xcurr(leftIndices)/roundCoeff));
    [~,minRoundIndices]=sort(roundDiffs);
    topIndex=min([packetSize;length(minRoundIndices)]);
    minRoundIndices=minRoundIndices(1:topIndex);
    origValues = xcurr(leftIndices(minRoundIndices));
    roundValues = roundCoeff*round(fraction*xcurr(leftIndices(minRoundIndices))/roundCoeff)/fraction;
    prevFixedIndices = fixedIndices;
    prevFixedValues=fixedValues;
    prevLeftIndices=leftIndices;
    fixedIndices=[fixedIndices;leftIndices(minRoundIndices)];
    fixedValues=[fixedValues;roundValues];
    leftIndices(minRoundIndices)=[];
    options=optimoptions('lsqnonlin','Algorithm','levenberg-marquardt',...
        'display', 'off','MaxIterations', 250, ...
        'SpecifyObjectiveGradient',true,...
        'FiniteDifferenceType', 'central');%,'CheckGradients', true); 
    objfun = @(xcurr)objective(xcurr,xprev,A,b,C, G2, N, FN, fixedIndices, fixedValues, s, wobj,wclose, wconst, wbarrier, IImagField, JImagField,origFieldVolumes);
    
    %for checking gradients only!
    %x0=rand(size(x0));
    [xcurr ,~,~,exitflag,~]= lsqnonlin(objfun,xcurr,[],[],options);
    roundDiffs= abs(fraction*xcurr(fixedIndices)-roundCoeff*round(fraction*xcurr(fixedIndices)/roundCoeff));
    if (max(roundDiffs)>10e-7) || (size(C,1)~=0 && max(abs(C*xcurr))>10e-7)   %did not converge, try again with different integers
        if (packetSize>1)
            packetSize=floor(packetSize/2)
        else
            fraction=fraction*2
            warning('BarrierPoissonSolve: raising fraction');
            if (fraction==2)
                success=-1;
                break;
            end
            packetSize=1
        end
        fixedIndices=prevFixedIndices;
        fixedValues=prevFixedValues;
        leftIndices=prevLeftIndices;
        xcurr=xprev;
    else
        i=i+packetSize;
    end
end

xcurr=xcurr*fraction; %so that we do get full integers at the end
xcurr(integerIndices)=round(xcurr(integerIndices));

rawFieldGx0Curr = G*xcurr;
funcValuesCurr = x2CornerMat*xcurr;
funcValuesCurr = reshape(funcValuesCurr, N, length(V))';

%checking bijectivity
minBijectivity=zeros(length(FN),1);
for i=0:length(FN)-1
    varOffset=i*3*N+1;
    barOffset=i*N+1;
    faceField = reshape(rawFieldGx0Curr(varOffset:varOffset+3*N-1),3,N)';
    faceField=faceField./repmat(MESH.normv(faceField), 1,3);
    tripleProducts = dot(repmat(FN(i+1,:),N,1), cross(faceField, faceField([2:end,1],:)),2);
    minBijectivity(i+1)=min(tripleProducts);
end

if show_figs
%     figure; title('param (using first 2 coordinates), blue = input, red = output');
%     dockfig;
%     hold on
%     patch('vertices', funcValues(:,1:2), 'faces', F, 'faceColor', 'none', 'edgeColor','b');
%     patch('vertices', funcValuesCurr(:,1:2), 'faces', F, 'faceColor', 'none', 'edgeColor','r');

    figure; title('min bijectivity on output');
    set(gcf, 'WindowStyle', 'docked'); colorbar;
    hold on
    patch('vertices', funcValuesCurr(:,1:2), 'faces', F, 'faceColor', 'flat', 'CData', minBijectivity);

%     for i=0:N-1
%         figure; 
%         title(sprintf('input param function %d and gradient',i+1));
%         dockfig;
%         hold on
%         patch('vertices', V, 'faces', F, 'faceColor', 'interp', 'CData', funcValues(:,i+1));
%         axis equal
%         cameratoolbar
%         quiver3(centers(:,1), centers(:,2), centers(:,3), rawFieldGx0(:, i*3+1),rawFieldGx0(:, i*3+2),rawFieldGx0(:, i*3+3));
%     end
end


function [f,g]=objective(xcurr,xprev,A,b,C,G2,  N, FN, fixedIndices, fixedValues, s, wobj,wclose, wconst, wbarrier, IImagField, JImagField, origFieldVolumes)

fAb = A*xcurr-b;
gAb = A;
fClose = xcurr-xprev;
gClose=speye(length(xcurr));
fLinConst=C*xcurr;
gLinConst=C;
fConst = xcurr(fixedIndices)-fixedValues;
gConst=fast_sparse((1:length(fixedIndices))', fixedIndices, ones(length(fixedIndices),1), length(fixedIndices),length(xcurr));
nf = length(FN);
field = G2*xcurr;
fBarrier = zeros(N*nf,1);
fBarrierDerivative= zeros(N*nf,1);

SImagField=IImagField;

for i=0:nf-1
    varOffset=i*2*N+1;
    barOffset=i*N+1;
    faceField = reshape(field(varOffset:varOffset+2*N-1),2,N)';
    faceFieldNext = faceField([2:N,1],:);
    
    imagProduct = (faceField(:,1).*faceFieldNext(:, 2) - faceField(:,2).*faceFieldNext(:, 1))./origFieldVolumes(i+1,:)';
    barResult = (imagProduct/s).^3 - 3*(imagProduct/s).^2 + 3*(imagProduct/s);
    barResult2 = 1./barResult -1;
    barResult2(imagProduct<=0)=Inf;
    barResult2(imagProduct>=s)=0;
    fBarrier(barOffset:barOffset+N-1)=barResult2;
    
    barDerivative=(3*(imagProduct.^2/s^3) -6*(imagProduct/s^2) + 3/s)./origFieldVolumes(i+1,:)';
    barDerivative(imagProduct<=0)=Inf;
    barDerivative(imagProduct>=s)=0;
    fBarrierDerivative(barOffset:barOffset+N-1) =  barDerivative;
    
    SImagField(barOffset:barOffset+N-1,:)=[faceFieldNext(:, 2), -faceFieldNext(:, 1), -faceField(:,2),faceField(:,1)];
end

gImagField=fast_sparse(IImagField, JImagField, SImagField, N*nf, 2*N*nf);
barDerVec=-fBarrierDerivative./(fBarrier.^2);
barDerVec(isinf(barDerVec))=0;
barDerVec(isnan(barDerVec))=0;
gBarrierFunc = spdiags(barDerVec, 0, length(fBarrier), length(fBarrier));
gBarrier=gBarrierFunc*gImagField*G2;
f=[wobj*fAb;wclose*fClose;wconst*fConst;wconst*fLinConst;wbarrier*fBarrier];
g=[wobj*gAb;wclose*gClose;wconst*gConst;wconst*gLinConst;wbarrier*gBarrier];
end





