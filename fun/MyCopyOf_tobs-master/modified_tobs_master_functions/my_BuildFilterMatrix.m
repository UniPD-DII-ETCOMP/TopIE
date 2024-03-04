%% Build filter matrix
function [Hfilt] = my_BuildFilterMatrix(Ind,mesh,radius)

    disp(['...Building filter matrix.'])
    
    R=[];
    ifilt=[];
    jfilt=[];
    
    % Compute baricenters
    BT=(mesh.P(mesh.T(:,1),:)+mesh.P(mesh.T(:,2),:)+mesh.P(mesh.T(:,3),:))./3; %barycenters
    % Retrieve design domain
    ind_design=[];
    Nmat=size(Ind,2);
    for ii=1:Nmat
        if strcmp(Ind(ii).tag,'design_domain')
            ind_design=[ind_design;Ind(ii).ind];
        end
    end

    % Selecting coordinates of design nodes
    X_den = BT(ind_design,1);
    Y_den = BT(ind_design,2);

    % Auxiliary counter
    k = 0;

    % Looping at the elements inside the design domain
    for i = 1:length(ind_design)

        % Element
        el = ind_design(i);

        % Coordinates of the centroid of the element
        x = BT(el,1);
        y = BT(el,2);

        % Weight function (linear with radial distance)
        r = max(0,1-(((X_den-x).^2+(Y_den-y).^2).^(1/2)/radius));

        % Auxiliary counter update
        k = k+nnz(r);

        % Index vectors building
        ifilt(k-nnz(r)+1:k) = i*ones(nnz(r),1);
        jfilt(k-nnz(r)+1:k) = find(r);

        % Weights update
        R(k-nnz(r)+1:k) = r(r>0)/sum(r);

    end

    % H filter matrix
    Hfilt = sparse(ifilt(1:k),jfilt(1:k),R(1:k),length(ind_design),length(ind_design));

    disp(['...Done.'])
end % end BuildFilterMatrix
        