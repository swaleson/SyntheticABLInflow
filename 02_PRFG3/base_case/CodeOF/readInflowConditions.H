// ------------------------ Inflow genereation ---------------------------//
sf_cutOff = sf_dict.get<scalar>("cutOff");
sf_seeds = sf_dict.get<scalarField>("seeds");
sf_tProfile = sf_dict.get<vectorField>("tProfile");
sf_uProfile = sf_dict.get<vectorField>("uProfile");
sf_sProfile = sf_dict.get<vectorField>("sProfile");
sf_lProfile = sf_dict.get<vectorField>("lProfile");
sf_p = sf_dict.get<vectorField>("pMainSeed");
sf_q = sf_dict.get<vectorField>("qMainSeed");
sf_k = sf_dict.get<vectorField>("kMainSeed");

// ------------------------ Inflow correction ---------------------------//
if (vc_dictName!="_notPresent"){
	vc_corr = vc_dict.get<scalar>("corr");
	vc_t1= vc_dict.get<vector>("t1");
	vc_t2 = vc_dict.get<vector>("t2");
	vc_LdOp = vc_dict.get<scalarField>("LdOp");
	vc_LoOp = vc_dict.get<vectorField>("LoOp");
	vc_UdOp = vc_dict.get<scalarField>("UdOp");
	vc_UoOp = vc_dict.get<vectorField>("UoOp");
	vc_LUpc = vc_dict.get<scalarField>("LUpc");
	vc_LUpr = vc_dict.get<scalarField>("LUpr");
	vc_nodes = vc_dict.get<vectorField>("nodes");
	vc_faces = vc_dict.get<scalarField>("faces");
	vc_faceCentres = vc_dict.get<vectorField>("faceCentres");
	vc_grad1Op = vc_dict.get<vectorField>("grad1Op");
	vc_grad2Op = vc_dict.get<vectorField>("grad2Op");
	vc_bcEdges = vc_dict.get<vectorField>("bcEdges");
	vc_bcNorm = vc_dict.get<vectorField>("bcNorm");
}
