/* ------------------------ BC interface ---------------------------*/
vector sf_origin;
vector sf_zDir;
const word sf_dictName;
const word vc_dictName;

/* ------------------------ Inflow generation ---------------------------*/
scalar sf_cutOff;
scalarField sf_seeds;
vectorField sf_tProfile;
vectorField sf_uProfile;
vectorField sf_sProfile;
vectorField sf_lProfile;
vectorField sf_p;
vectorField sf_q;
vectorField sf_k;

/* ------------------------ Inflow correction ---------------------------*/
scalar vc_corr;

vector vc_t1;
vector vc_t2;

scalarField vc_LdOp;
vectorField vc_LoOp;
scalarField vc_UdOp;
vectorField vc_UoOp;

scalarField vc_LUpc;
scalarField vc_LUpr;

vectorField vc_nodes;
scalarField vc_faces;
vectorField vc_faceCentres;
vectorField vc_grad1Op;
vectorField vc_grad2Op;
vectorField vc_bcEdges;
vectorField vc_bcNorm;