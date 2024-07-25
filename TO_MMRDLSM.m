%**************************************************************************
% MeshGeneration
[ndof,noe,NodeCoord,EDofMat,ElNodes,NodeRep,Phi,IC,nelx,nely] = ...
    MeshGeneration(GP.nor,GP.lx,GP.ly,GP.els,GP.X0,GP.Y0,GP.noph);
%**************************************************************************
% Boundary conditions implementation
[F,Dis,FixedNodes,ukdis,kdis] = BoundaryConditionsImplementation(ndof,...
    GP.EBC,GP.NBC,NodeCoord,GP.els);
%**************************************************************************
% Indexes which is used in stiffness assembly
iK = reshape(kron(EDofMat,ones(8,1))',8*8*noe,1);
jK = reshape(kron(EDofMat,ones(1,8))',8*8*noe,1);
%**************************************************************************
% Element stiffness calculation
[Ke,Vole,Bmat,Emat] = PlaneElementStiffnessCalculation(GP.ngp,GP.els,MP.v);
%**************************************************************************
% Calculation of T1 and T2 for each element
[T1e,T2e] = RD_T1T2Generator(GP.ngp,GP.els);
% T1 and T2 assembly
iT = reshape(kron(ElNodes,ones(4,1))',4*4*noe,1);
jT = reshape(kron(ElNodes,ones(1,4))',4*4*noe,1);
st1 = reshape(T1e(:)*ones(1,noe),4*4*noe,1);
st2 = reshape(T2e(:)*ones(1,noe),4*4*noe,1);
T1 = sparse(iT,jT,st1);
T2 = sparse(iT,jT,st2);
%% Optimization Loop
% Parameters
% opc = Optimization counter
opc = 1;
% [r,s] = Grid on each element for Phi interpolation
[r,s] = meshgrid(linspace(-1,1,20));
% ns is the number of points in each element for Phi interpolation
ns = numel(s);
% tmpPhi = Interpolated Phi values on element for volume fraction
% calculation
tmpPhi = zeros(numel(s),noe,GP.noph);
% vlfe = The percentage of materials from each phase in the element
vlfe = zeros(GP.noph,noe);
% VC is volume constraint
VC = zeros(GP.noph,OP.NOI);
% Compliance is objective function
Compliance = zeros(OP.NOI,1);
% derivative of volume fraction wrt phi
vlfe_phi = zeros(GP.noph,noe,GP.noph);
% Lambda in augmented Lagrngian
Lambda = zeros(GP.noph,1);
%%
while opc<=OP.NOI
    % Volume fraction calculation for each element
    for i=1:GP.noph
        tmpPhi(:,:,i) = 0.25*((1-r(:)).*(1-s(:))*Phi(ElNodes(:,1),:,i)'+...
                              (1+r(:)).*(1-s(:))*Phi(ElNodes(:,2),:,i)'+...
                              (1-r(:)).*(1+s(:))*Phi(ElNodes(:,3),:,i)'+...
                              (1+r(:)).*(1+s(:))*Phi(ElNodes(:,4),:,i)');
    end
    % hphie = heavisided phi in elements
    hphie = tmpPhi>0;
%**************************************************************************
    for i=1:GP.noph
%**************************************************************************
% Calculation of vlfe
        if i<GP.noph
            cond = sum(hphie(:,:,1:i),3)+(~hphie(:,:,i+1));
            vlfe(i,:) = sum(cond==(i+1))/ns;
            for j=1:i
                if (i==1 && j==1)
                    vlfe_phi(1,:,1) = sum(1-hphie(:,:,2))/ns;
                else
                    ind = [setdiff(1:i,j)];
                    cond = sum(hphie(:,:,ind),3)+(~hphie(:,:,i+1));
                    vlfe_phi(i,:,j) = sum(cond==(i-1))/ns;
                end
            end
            cond = sum(hphie(:,:,1:i),3);
            vlfe_phi(i,:,i+1) =-sum(cond==i)/ns;
        else
            cond = sum(hphie(:,:,1:i),3);
            vlfe(i,:) = sum(cond==i)/ns;
            for j=1:i
                ind = [setdiff(1:i,j)];
                cond = sum(hphie(:,:,ind),3);
                vlfe_phi(i,:,j) = sum(cond==(i-1))/ns;
            end
        end
        VC(i,opc) = sum(vlfe(i,:))/noe;
%**************************************************************************
    % Lambda calcultion in augmented Lagrngian
        if opc<OP.nRelax(i)
            Lambda(i)=OP.Mu(i)*...
                (VC(i,opc)-VC(i,1)+(VC(i,1)-OP.VC0(i))*opc/OP.nRelax(i));
        else
            Lambda(i) =Lambda(i)+OP.Gamma(i)*(VC(i,opc)-OP.VC0(i));
            OP.Gamma(i) = min(OP.Gamma(i)+OP.dGamma(i),OP.maxGamma(i));
        end
    end
%**************************************************************************
    % ErModel=Ersatz material model and stiffness calculation
    ErModel = sum((vlfe.^MP.P).*MP.Emax')+(1-sum(vlfe))*MP.Emin;
%**************************************************************************
    % Assembling of stiffness matrix
    sK = reshape(Ke(:)*(ErModel.*ones(1,noe)),8*8*noe,1);
    K = sparse(iK,jK,sK); 
    K = (K+K')/2;
%**************************************************************************
% Solving system of equlibrium equations
    Dis(ukdis)=K(ukdis,ukdis)\(F(ukdis)-K(ukdis,kdis)*Dis(kdis));
%**************************************************************************
    % ElemComp = Elemet compliance
    ElemComp=sum(0.5*(Ke*Dis(EDofMat)').*(Dis(EDofMat)'.*ErModel));
    Compliance(opc) = sum(ElemComp);
    %% Sensitivity and velocity calculation
    ErModel_phi = MP.P*sum(vlfe_phi.*(vlfe.^(MP.P-1)).*MP.Emax');
    VC_phi = sum(Lambda.*vlfe_phi);
%**************************************************************************
    % Loop over materials for updating design domain   
    for i=1:GP.noph
%**************************************************************************
    % Comp_phi = objective function sensitivity wrt Phi
        ElemComp_phi=-sum(0.5*(Ke*Dis(EDofMat)').*(Dis(EDofMat)'.*ErModel_phi(:,:,i)));
        Comp_phi = sparse(ElNodes,ones(noe,4),0.25*ElemComp_phi'.*ones(noe,4));
%**************************************************************************
    % A_phi = area function sensitivity wrt Phi
        A_phi = sparse(ElNodes,ones(noe,4),0.25*VC_phi(:,:,i)'.*ones(noe,4));
%**************************************************************************
    % V = Boundary velocity in LS method
        V = Comp_phi/mean(abs(Comp_phi))+A_phi;
%**************************************************************************
 % Update scheme
        T = (T1/OP.delta_T+OP.Tho*T2);
        Yy = (T1*(Phi(:,:,i)/OP.delta_T-V));
        Phi(:,:,i)=T\Yy;
        Phi(:,:,i) = min(max(Phi(:,:,i),-1),1);
    end
    if OP.PlotResult==1
        hPhi = Phi>0;
        Xdef = NodeCoord(:,1);
        Ydef = NodeCoord(:,2);
        hold on
        for i=1:GP.noph
            for j=1:GP.nor
                x = reshape(Xdef(IC{j}),(nelx(j)+1),(nely(j)+1))';
                y = reshape(Ydef(IC{j}),(nelx(j)+1),(nely(j)+1))';
                z = reshape(Phi(IC{j},:,i),(nelx(j)+1),(nely(j)+1))';
                contourf(x,y,z,[0,0],"FaceColor",MP.Colors(i,:));
            end
        end
        for j=1:GP.nor
            x = reshape(Xdef(IC{j}),(nelx(j)+1),(nely(j)+1))';
            y = reshape(Ydef(IC{j}),(nelx(j)+1),(nely(j)+1))';
            z = reshape(-Phi(IC{j},:,1),(nelx(j)+1),(nely(j)+1))';
            contourf(x,y,z,[0,0],"FaceColor",[1,1,1]);
        end
        axis('equal')
        title('Optimum Shape');
        hold off
        drawnow
        writeVideo(v,getframe);
    end
    [opc,VC(:,opc)']
    opc=opc+1;
end
%% Mesh generation function
function [ndof,noe,NodeCoord,EDofMat,ElNodes,NodeRep,Phi,IC,nelx,...
    nely] = MeshGeneration(nor,lx,ly,els,X0,Y0,noph)
%==========================================================================
    % Inputs
    % nor = Number of rectangles
    % lx = array which has x_length of rectangles
    % ly = array which has y_length of rectangles
    % els = Element size at each direction
    % noph = Number of phases
%==========================================================================
    % Preallocation arrays for the 
    % nummber of elements for each rectangle at each direction
    nelx = zeros(nor,1); % nelx(i) = int64(lx(i)/els); 
    nely = zeros(nor,1); % nely(i) = int64(ly(i)/els); 
%==========================================================================
    % coordinates{i} = [x,y,z] where x,y and z are node coordinates for
    % i'th rectangle
    coordinates = cell(nor,1);
%==========================================================================
    % NodeCoord = [coordinates{1};...,coordinates{nor}] by removing
    % duplicate nodes
    NodeCoord = [];
%==========================================================================
    % Loop over rectangles
    for i=1:nor
        nelx(i) = int64(lx(i)/els); 
        nely(i) = int64(ly(i)/els); 
        nonx = nelx(i)+1; 
        nony = nely(i)+1; 
%==========================================================================
        % X0 = X-coordinate of left-bottom-back for each rectangle
        % Y0 = Y-coordinate of left-bottom-back for each rectangle
        % nonx and nony are number of nodes for i'th rectangle for
        % x-y direction, respectively which can be calculated by
        x_node = linspace(X0(i),X0(i)+lx(i),nonx);
        y_node = linspace(Y0(i),Y0(i)+ly(i),nony);
%==========================================================================
        x = repmat(x_node',nony,1);
        y = repmat(kron(y_node',ones(nonx,1)),1);
        coordinates{i} = [x,y];
        NodeCoord = [NodeCoord;coordinates{i}];
    end
%==========================================================================
    % removing duplicate nodes
    NodeCoord = unique(NodeCoord,'stable','rows');
%==========================================================================
    % Phi specifies if material exist at nodes. size(Phi)=number of rows in
    % NodeCoord
    % Phi>=0 => material exist
    % Phi = 0=>boundary representation
    % Phi<0 => no material
    Phi = 0.1*ones(size(NodeCoord,1),1,noph);
%==========================================================================
    % ElNodes = global indexes for nodes of each element 
    ElNodes = [];
%==========================================================================
    % IC variable for prohibitting duplicate indexing is used, because 
    % nodes with the same coordinates in different rectangles must have 
    % unique global indexes     
    IC = cell(nor,1);
%==========================================================================
    % Loop over rectangles
    for i=1:nor
        nodes = (1:size(coordinates{i},1))';
        [~,ic,iC] = intersect(coordinates{i},NodeCoord,'stable','rows');
        nodes(ic) = iC;IC{i}=iC;
        a = repmat([0,1,nelx(i)+1,nelx(i)+2],nelx(i),1);
        b = repmat(a,nely(i),1)+kron((0:nely(i)-1)',ones(nelx(i),1));
        eleNode = b+(1:nelx(i)*nely(i))';
        if i>1
            eleNode = nodes(eleNode);
        end
        ElNodes = [ElNodes;eleNode]; %ok
    end
%==========================================================================
    % EDofMat = Degrees of freedom for each element
    EDofMat = kron(ElNodes,[2,2])+repmat([-1,0],1,4);
%==========================================================================
    % noe = Number of elements
    % ndof = Number of degrees of freedom
    % nonodes = Number of nodes
    noe = sum(nelx.*nely);
    nonodes = max(ElNodes,[],'all');
    ndof = 2*nonodes;
%==========================================================================
    % NodeRep = reppetition of each node in all element;
    NodeRep=groupcounts(ElNodes(:));
end
%% Boundary conditions implementation function
function [F,Dis,FixedNodes,ukdis,kdis] = ...
    BoundaryConditionsImplementation(ndof,EBC,NBC,NodeCoord,els)
%==========================================================================
    % Dis = Displacement vector
    % F = Force vector
    Dis = nan(ndof,1);
    F = zeros(ndof,1);
%==========================================================================
    for s=1:2
        if s==1
            BC = EBC;
            ND = Dis;
        else
            BC = NBC;
            ND = F;
        end
        for i=1:size(BC,1)
            Cond = NodeCoord(:,1)>=BC(i,1) & ...
                   NodeCoord(:,1)<=BC(i,2) & ...
                   NodeCoord(:,2)>=BC(i,3) & ...
                   NodeCoord(:,2)<=BC(i,4);
            Cond = find(Cond);
            if ~isempty(Cond)
                ND(Cond*2-1) = BC(i,5);
                ND(Cond*2) = BC(i,6);
            else
                distance = sum((NodeCoord-[BC(i,1),BC(i,3)]).^2,2);
                Cond = find(distance<=els);
                ND(Cond*2-1) = BC(i,5)/numel(Cond);
                ND(Cond*2) = BC(i,6)/numel(Cond);
            end
        end
        if s==1
            Dis = ND;
        else
            F = ND;
%==========================================================================
            % FixedNodes = Nodes should remain in design 
            % domain in optimization loop, for example NBC nodes
            FixedNodes = Cond;
        end
    end
%==========================================================================
    % ukdis = Unknown displacement index
    % kdis = Known displacement index
    ukdis = isnan(Dis);
    kdis = ~isnan(Dis);
end
%% Plane element stiffness calculation function
function [Ke,Vole,Bmat,Emat] = PlaneElementStiffnessCalculation(ngp,els,v)

% Notes!!!
% This function is usable for structres with elements with equal size at 
% each direction
%==========================================================================
% ngp = Number of gauss points for each direction
% els = Element size at each direction
%==========================================================================
    % Global coordinates for each element
    % because the element size dose not change in the domain of structre,
    % this global coordinate is used for all element.
    eXcor = [0,els,0,els]';
    eYcor = [0,0,els,els]';
%     eXcor = [0,els,0,els]';
%     eYcor = [0,0,-els,-els]';
%==========================================================================
 % Elasticity matrix by cosidering 1 as modulus of elasticity
Emat = 1/(1-v^2)*[1,v,0;...
    v,1,0;...
    0,0,(1-v)/2];
%==========================================================================
    % Strain-displacement matrix
    Bmat = zeros(3,8);
%==========================================================================
% Volume at gussian points
    Volg = zeros(ngp);
%==========================================================================
    % Vole = Volume for each element
    Vole = 0;
%==========================================================================
    % Ke = stiffness matrix for element
    Ke = zeros(8);
%==========================================================================
% gauss points and their wheights
    [gp,wgp]=makegaussianpoint(ngp);
%==========================================================================
   % Loop over gauss points for Stiffness calcualtion
    for j=1:ngp
        s = gp(j);
        for k=1:ngp
            r = gp(k);
%==========================================================================
% 4-node shape functions for plane element
%                 N=0.25*[(1-r)*(1-s),(1+r)*(1-s),...
%                          (1-r)*(1+s),(1+r)*(1+s)];
%                 N=0.25*[(1-r)*(1+s),(1+r)*(1+s),...
%                          (1-r)*(1-s),(1+r)*(1-s)];
%==========================================================================
            % Calculation of shape function derivative wrt local
            % coordiantes
            N_r = 0.25*[-(1-s),(1-s),-(1+s),(1+s)];

            N_s = 0.25*[-(1-r),-(1+r),(1-r),(1+r)];
%             N_r = 0.25*[-(1+s),(1+s),-(1-s),(1-s)];
% 
%             N_s = 0.25*[(1-r),(1+r),-(1-r),-(1+r)];
%==========================================================================
            % Calcualtion derivitave of global coordiantes 
            % wrt local coordiantes
            X_r = N_r*eXcor; Y_r = N_r*eYcor;
            X_s = N_s*eXcor; Y_s = N_s*eYcor;
%==========================================================================
            % Jacobian matrix calculation
            J = [X_r,Y_r;X_s,Y_s];
            N_X_Y = J\[N_r;N_s];
%==========================================================================
            % For Normal strain at x direction
                Bmat(1,1:2:end)=N_X_Y(1,:);

            % For Normal strain at y direction
                Bmat(2,2:2:end)=N_X_Y(2,:);

            % For Shear strain at x and y direction
                Bmat(3,1:2:end)=N_X_Y(2,:);
                Bmat(3,2:2:end)=N_X_Y(1,:);
%==========================================================================
            % Volume at gussian points
                Volg(k,j) = det(J)*wgp(k)*wgp(j);

            % Volume for each element = sum(Arg)
                Vole = Vole+Volg(k,j);
%==========================================================================
            % Stiffness for each element = sum(Stiffness at each
            % gussian point)
                Ke = Ke+Bmat'*Emat*Bmat*Volg(k,j);
%==========================================================================
        end
    end
end
%% T1 and T2 Generator function
function [T1e,T2e] = RD_T1T2Generator(ngp,els)

% Notes!!!
% This function is usable for structres with elements with equal size at 
% each direction
%==========================================================================
% ngp = Number of gauss points for each direction
% els = Element size at each direction
%==========================================================================
    % Global coordinates for each element
    % because the element size dose not change in the domain of structre,
    % this global coordinate is used for all element.
    eXcor = [0,els,0,els]';
    eYcor = [0,0,els,els]';
%==========================================================================
% Volume at gussian points
    Volg = zeros(ngp,ngp);
%==========================================================================
   % Loop over gauss points for Stiffness calcualtion
    % In RDLS 
    % T1 = integral(N'*N)
    % T2 = integral(grad(N')*grad(N))
    % where N is shape function vector with the size of 1X8 and grad(N) and
    % grad(N) is gradient of shape function vector wrt x,y and z with the
    % size of 3X8
    T1e = zeros(4);
    T2e = zeros(4);
%==========================================================================
% gauss points and their wheights
    [gp,wgp]=makegaussianpoint(ngp);
%==========================================================================
% Loop over gauss points for Stiffness calcualtion
    for j=1:ngp
        s = gp(j);
        for k=1:ngp
            r = gp(k);
%==========================================================================
% 4-node shape functions for plane element
                N=0.25*[(1-r)*(1-s),(1+r)*(1-s),...
                         (1-r)*(1+s),(1+r)*(1+s)];
%==========================================================================
            % Calculation of shape function derivative wrt local
            % coordiantes
            N_r = 0.25*[-(1-s),(1-s),-(1+s),(1+s)];

            N_s = 0.25*[-(1-r),-(1+r),(1-r),(1+r)];
%==========================================================================
            % Calcualtion derivitave of global coordiantes 
            % wrt local coordiantes
            X_r = N_r*eXcor; Y_r = N_r*eYcor;
            X_s = N_s*eXcor; Y_s = N_s*eYcor;
%==========================================================================
            % Jacobian matrix calculation
            J = [X_r,Y_r;X_s,Y_s];
            N_X_Y = J\[N_r;N_s];
%==========================================================================
            % Volume at gussian points
                Volg(k,j) = det(J)*wgp(k)*wgp(j);
%==========================================================================
            % T1 and T2 for each element = sum(Stiffness at each
            % gussian point)
                T1e = T1e+N'*N*Volg(k,j);
                T2e = T2e+N_X_Y'*N_X_Y*Volg(k,j);
%==========================================================================
        end
    end
end