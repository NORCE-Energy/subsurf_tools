function ekofisk = EkofiskGeostat_perm3D



% This function is used to specify all geostatistical parameters used to

% generate an initial ensemble. Mean values for permeability, porosity,

% net-to-gross, region multipliers, oil-water contacts, and fault

% multipliers, are loaded or copied from the Statoil benchmark case. 

%

% For more information we refer to the paper: 

%

% Lorentzen, R., Luo, X., Bhakta, T., Valestrand, R.: "History matching

% ekofisk reservoir and petroelastic models using seismic impedance with

% correlated noise" Submitted to SPE Journal.

% 

% We also ask that the above paper is cited in publications aided by this

% code.

%

% This program is free software: you can redistribute it and/or modify

% it under the terms of the GNU General Public License as published by

% the Free Software Foundation, either version 3 of the License, or

% (at your option) any later version.

%

% This program is distributed in the hope that it will be useful,

% but WITHOUT ANY WARRANTY; without even the implied warranty of

% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the

% GNU General Public License for more details.

%

% You should have received a copy of the GNU General Public License

% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%

% Copyright (C): IRIS (International Research Institute of Stavanger), 2017. 

% Contact: Rolf.Lorentzen@iris.no

% actnum
load('inputData.mat', 'options');
actnum = options.actnum;
ekofisk.actnum=options.actnum;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ekofisk.dim=options.dim;
dim = options.dim;
ldim=dim(1)*dim(2); % layer dimension

% %%%%% Get porosity from another project
% load('/media/tubh/IntData/Bhakta/NIOR Project/Work in progress/MatlabCodes/Ens_Rerevoir_Para_est/Ekofisk_Field/data_backup/Ensembles_06.10.2017/Ekofisk_Data_LOFS2.mat', 'PORE')
% 
%  Eko_poro_std= std(PORE,[],2); %%% I try to 
%  Eko_poro_std= Eko_poro_std';
%  Eko_poro_std = reshape(Eko_poro_std,options.dim);
%  Eko_poro_std = fliplr(Eko_poro_std);
%  Eko_poro_std(actnum==0) = [];

%Eko_poro_std = 0.01*ones(sum(actnum),1);

%%%% Also create permeabililty field as a vector same length of porosity
%%%% std
Eko_permx_std = 0.8.*ones(sum(actnum),1);
Eko_permy_std = 0.8.*ones(sum(actnum),1);
Eko_permz_std = 0.5.*ones(sum(actnum),1);
%Eko_permx_std = Eko_poro_std / mean(Eko_poro_std);
%Eko_permx_std(Eko_permx_std>1)=1;

% porosity
% p = read_petrel_datafile('PETW_DEC14_RESV_V4_PROP_PORO.GRDECL','Oct14_PORO');
% p = p.data;
% p(actnum==0) = [];

% for nr=1:dim(3)
%     index=ldim*(nr-1)+1:1:ldim*nr;
%     values=p(sum(actnum(1:index(1)-1))+1:sum(actnum(1:index(end))));   
%     meanv(nr)=mean(values); %%#ok<*AGROW>
%     stdv(nr)=std(values);
% end
% 
% ekofisk.poroMean = p;
% ekofisk.poroLayerMean=meanv';
% ekofisk.poroLayerStd=stdv';
% ekofisk.poroStd= Eko_poro_std; %% 0.01;
% ekofisk.poroLB=0.001;
% ekofisk.poroUB=0.5;
% ekofisk.poroRange=26;

% permeability
k = read_petrel_datafile('PETW_DEC14_RESV_V4_PROP_PERMX.GRDECL','Oct14_PERMX');
k = k.data;
k(actnum==0) = [];
k = log(k); % Convert permeability field to log-perm.

for nr=1:dim(3)
    index=ldim*(nr-1)+1:1:ldim*nr;
    values=k(sum(actnum(1:index(1)-1))+1:sum(actnum(1:index(end))));    
    meanv(nr)=mean(values);
    stdv(nr)=std(values);
end

ekofisk.permxLogMean = k;
ekofisk.permxLayerLnMean=meanv';
ekofisk.permxLayerStd=stdv';
ekofisk.permxStd= Eko_permx_std ; 
ekofisk.permyStd= Eko_permy_std ; 
ekofisk.permzStd= Eko_permz_std ; 
ekofisk.permxLB=-5;
ekofisk.permxUB=10;
ekofisk.permxRange=26;

% ekofisk.permzRange=4;

% % correlation between layers
% 
% for nr=1:dim(3)-1
%     index=ldim*(nr-1)+1:1:ldim*nr;
%     index2=ldim*(nr)+1:1:ldim*(nr+1);
%     actlayer1=actnum(index);    
%     actlayer2=actnum(index2);
%     active=actlayer1.*actlayer2;
%     values1=[k(sum(actnum(1:index(1)-1))+1:sum(actnum(1:index(end)))) ;...
%         p(sum(actnum(1:index(1)-1))+1:sum(actnum(1:index(end))))];
%     values2=[k(sum(actnum(1:index2(1)-1))+1:sum(actnum(1:index2(end)))) ;...
%         p(sum(actnum(1:index2(1)-1))+1:sum(actnum(1:index2(end))))];
%     v1=[actlayer1;actlayer1];
%     v1(v1==1) = values1;
%     v2=[actlayer2;actlayer2];
%     v2(v2==1) = values2;
%     co=corrcoef(v1(active==1), v2(active==1));
%     corrWithNextLayer(nr)=co(1,2);
% 
% end

% ekofisk.corrWithNextLayer=corrWithNextLayer';

% correlation between porosity and permeability
ekofisk.poroPermxCorr=0.0;

% correlation between permx and permy 
ekofisk.PermxPermyCorr=0.7;

% correlation between permx and permz 
ekofisk.PermxPermzCorr=0.3;













