/*--------------------- Calculating fluctuations ----------------------------------------------- */
scalar pi = 3.14159265359;
scalar k_dot_x = 0;
scalar coskx = 0.0;
scalar sinkx = 0.0;
scalar phase = 0;
vector uvw_local = vector::zero;
int iSeed = 0;
scalar lScale = 0;

// Assembling seed1
iSeed = int(sf_seedsWeights[iN][0]+intOff);
lScale = max(sf_seedsProps[iSeed][2],sf_lSmall);

forAll(sf_k,iK)
{
	// If fluctuations cannot be supported by the mesh and fall below the cutOff they are not assembled
	if ((mag(sf_k[iK])*sf_dxScale/lScale)>sf_cutOff){
         continue;
		}
	phase = 1.0E2*iSeed+1.0E1*iK;
    k_dot_x = ((sf_k[iK] & sf_nodesTrans[iN]) - sf_k[iK][0]*sf_seedsProps[iSeed][0]*time + phase)/lScale;
    coskx = cos(2.0*pi*k_dot_x);
    sinkx = sin(2.0*pi*k_dot_x);

    uvw_local += sf_p[iK]*coskx * sf_seedsProps[iSeed][1]*sqrt(sf_seedsWeights[iN][2]);
    uvw_local += sf_q[iK]*sinkx * sf_seedsProps[iSeed][1]*sqrt(sf_seedsWeights[iN][2]);
}


// Assembling seed2 if required
if((sf_seedsWeights[iN][1]+intOff)>=0){
	iSeed = int(sf_seedsWeights[iN][1]+intOff);
	lScale = max(sf_seedsProps[iSeed][2],sf_lSmall);
	forAll(sf_k,iK)
	{
	    // If fluctuations cannot be supported by the mesh and fall below the cutOff they are not assembled
		if ((mag(sf_k[iK])*sf_dxScale/lScale)>sf_cutOff){
			 continue;
			}
		phase = 1.0E2*iSeed+1.0E1*iK;
		k_dot_x = ((sf_k[iK] & sf_nodesTrans[iN]) - sf_k[iK][0]*sf_seedsProps[iSeed][0]*time + phase)/lScale;
		coskx = cos(2.0*pi*k_dot_x);
		sinkx = sin(2.0*pi*k_dot_x);

		uvw_local += sf_p[iK]*coskx * sf_seedsProps[iSeed][1]*sqrt(1.0-sf_seedsWeights[iN][2]);
		uvw_local += sf_q[iK]*sinkx * sf_seedsProps[iSeed][1]*sqrt(1.0-sf_seedsWeights[iN][2]);
	}
}

// Limiting value to sf_uSmall
uvw_local[0] = max(uvw_local[0],sf_uSmall-sf_nodesProps[iN][0]);

// Transforming to global ref
uvw_flu = uvw_local[0]*sf_xDir + uvw_local[1]*sf_yDir + uvw_local[2]*sf_zDir;

