/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "synInflowFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolateXY.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::synInflowFvPatchVectorField::
synInflowFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    sf_origin(vector::zero),
    sf_zDir(vector(0.0,0.0,1.0)),
    sf_dictName("inflowGenerationDict"),
    vc_dictName("inflowCorrectionDict")
{}


Foam::synInflowFvPatchVectorField::
synInflowFvPatchVectorField
(
    const synInflowFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    sf_origin(ptf.sf_origin),
    sf_zDir(ptf.sf_zDir),
    sf_dictName(ptf.sf_dictName),
    vc_dictName(ptf.vc_dictName)
{}


Foam::synInflowFvPatchVectorField::
synInflowFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    sf_origin(dict.getOrDefault<vector>("origin", vector::zero)),
    sf_zDir(dict.getOrDefault<vector>("zDir", vector(0.0,0.0,1.0))),
    sf_dictName(dict.getOrDefault<word>("sf_dictName", "inflowGenerationDict")),
    vc_dictName(dict.getOrDefault<word>("vc_dictName", "_notPresent"))
{

const fvMesh& mesh=p.boundaryMesh().mesh();

   IOdictionary sf_dict
    (
    		IOobject
			(
			sf_dictName,
			mesh.time().system(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::NO_WRITE
			)
    );

   IOdictionary vc_dict
    (
    		IOobject
			(
			vc_dictName,
			mesh.time().system(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::NO_WRITE
			)
    );

    #include "readInflowConditions.H"

    // Preparing to use solver
    ls_Ld = vc_LdOp;
    ls_Lo = vc_LoOp;
    ls_Ud = vc_UdOp;
    ls_Uo = vc_UoOp;
    ls_LUpc = vc_LUpc;
    ls_LUpr = vc_LUpr;

    // Getting face data
    sf_faceC = p.Cf();
    sf_faceAreas = p.magSf();
    // Setting local ref system (xDir taken as average)
    sf_zDir = sf_zDir/mag(sf_zDir);
	sf_xDir = gSum(-p.Sf()/p.magSf())/gSum(p.magSf()*0.0+1);
    sf_yDir = - sf_xDir ^ sf_zDir;
    Info << "SynFlow: xDir="  << sf_xDir << "    yDir:=" << sf_yDir << "    zDir:=" << sf_zDir << "." << endl;

    // Initializing local variables
    sf_facesIds = scalarField(sf_faceAreas.size(),0.0);
    sf_seedsProps = vectorField(sf_seeds.size(),vector::zero);
    sf_nodesTrans = vectorField(vc_nodes.size(),vector::zero);
    sf_uTot_c = vectorField(sf_faceC.size(),vector::zero);
    sf_nodesProps = vectorField(vc_nodes.size(),vector::zero);
    sf_nodalArea = scalarField(vc_nodes.size(),0.0);
    vectorField sf_uAvg_n = vectorField(vc_nodes.size(),vector::zero);
    vectorField sf_uFlu_n = vectorField(vc_nodes.size(),vector::zero);

	vector uvw_avg = vector::zero;
	vector uvw_flu = vector::zero;
    vector xyz = vector::zero;
    scalar time = this->db().time().value();
    scalar deltaT = this->db().time().deltaTValue();
    time = max(time,deltaT);
    sf_curTimeIndex = -1;
    sf_dxScale = pow(gSum(p.magSf())/gSum(p.magSf()*0.0+1),0.5);
    sf_lSmall = 10*sf_dxScale;
    scalar tMult = 0.0;
    scalar intOff = 0.01;
    sf_calculateCorrections = true;

    if (vc_dictName=="_notPresent"){
		sf_calculateCorrections = false;
	}

    // Check on mesh
    if (sf_calculateCorrections && (abs(gSum(p.magSf()*0.0+1)-vc_faceCentres.size())>0)){
		Info << "SynFlow: ATTENTION InflowCorrectionDict appears to be not corresponding to the present mesh." << endl;
	}

    // Taking care in case vcDict is not present, calculating velocities directly at cell centres
    if (!sf_calculateCorrections){
        vc_nodes = sf_faceC;                                                              // ATTENTION! NODES ARE SUBSTITUTED WITH FACE CENTRES
        sf_nodesTrans = vectorField(vc_nodes.size(),vector::zero);                        // ATTENTION! NODES ARE SUBSTITUTED WITH FACE CENTRES
        sf_nodesProps = vectorField(vc_nodes.size(),vector::zero);                        // ATTENTION! NODES ARE SUBSTITUTED WITH FACE CENTRES
        scalarField uu = scalarField(sf_uProfile.component(1)*sf_uProfile.component(2));
        sf_uSmall = 0.05*sum(uu)/sum(uu*0.0+1);
	}
    // Interpolating u, s, l at seeds
    forAll(sf_seeds, iS)
    {
        sf_seedsProps[iS][0] = interpolateXY(sf_seeds[iS],scalarField(sf_uProfile.component(0)),scalarField(sf_uProfile.component(1)*sf_uProfile.component(2)));
        sf_seedsProps[iS][1] = interpolateXY(sf_seeds[iS],scalarField(sf_sProfile.component(0)),scalarField(sf_sProfile.component(1)*sf_sProfile.component(2)));
        sf_seedsProps[iS][2] = interpolateXY(sf_seeds[iS],scalarField(sf_lProfile.component(0)),scalarField(sf_lProfile.component(1)*sf_lProfile.component(2)));
    }
    // Transfroming nodes coordinates and interpolating u, s and l at nodes
    forAll(vc_nodes, iN)
    {
		xyz = vc_nodes[iN];

        sf_nodesTrans[iN][0] = (xyz - sf_origin) & sf_xDir;
        sf_nodesTrans[iN][1] = (xyz - sf_origin) & sf_yDir;
        sf_nodesTrans[iN][2] = (xyz - sf_origin) & sf_zDir;

        sf_nodesProps[iN][0] = interpolateXY(sf_nodesTrans[iN][2],scalarField(sf_uProfile.component(0)),scalarField(sf_uProfile.component(1)*sf_uProfile.component(2)));
        sf_nodesProps[iN][1] = interpolateXY(sf_nodesTrans[iN][2],scalarField(sf_sProfile.component(0)),scalarField(sf_sProfile.component(1)*sf_sProfile.component(2)));
        sf_nodesProps[iN][2] = interpolateXY(sf_nodesTrans[iN][2],scalarField(sf_lProfile.component(0)),scalarField(sf_lProfile.component(1)*sf_lProfile.component(2)));
    }
    // Attributing nodes to seeds and calculating weights [seed1, seed2, w1] with w2 = 1-w1 (-1 no data)
    sf_seedsWeights = vectorField(vc_nodes.size(),vector::zero);
    if (sf_seeds.size()<=1){
		sf_seedsWeights = sf_seedsWeights*0.0 + vector(0.0,0.0,1.0);
	}else{
		scalar zz = 0.0;
		scalar ww;
		scalar nSig = 3.0;
		scalar cVal = 1/(erf(nSig/sqrt(2.0)));

		forAll(vc_nodes, iN)
		{
			xyz = sf_nodesTrans[iN];
			zz = xyz[2];

			if (zz<=sf_seeds[0]){
				 sf_seedsWeights[iN] = vector(0.0,-1.0,1.0);
			}else if (zz>=sf_seeds[sf_seeds.size()-1]) {
				 sf_seedsWeights[iN] = vector(sf_seeds.size()-1,-1.0,1.0);
			} else {
				forAll(sf_seeds,iS){
					if (sf_seeds[iS]>zz)
					{
					ww = 0.5 - cVal/2*erf((zz-0.5*(sf_seeds[iS-1]+sf_seeds[iS]))/((sf_seeds[iS]-sf_seeds[iS-1])*sqrt(2.0)/(2.0*nSig)));
					sf_seedsWeights[iN] = vector(iS-1,iS,ww);
					break;
					}
				}
			}
		}
	}
    // If corrections are not available no need to go further
    if (!sf_calculateCorrections){
        return;
	}

    // Calculating velocities at nodes
    forAll(vc_nodes, iN)
    {
		uvw_avg = uvw_avg*0.0;
        uvw_flu = uvw_flu*0.0;
        #include "averageField.H"
        #include "fluctuatingField.H"
        sf_uAvg_n[iN] = uvw_avg;
        sf_uFlu_n[iN] = uvw_flu;
    }

    // Getting current processor faces ids based on face centre position
    int thisId = 0;
    scalar thisDist = 0.0;
    scalar smallDist = 0.0;
    scalar maxDist = 0.0;
    forAll(sf_faceC, iC)
    {
        xyz = sf_faceC[iC];

        thisId = -1;
        smallDist = 1.0E12;
        forAll(vc_faceCentres, iFc)
        {
			thisDist = mag(xyz-vc_faceCentres[iFc]);
			if (thisDist<smallDist)
			{
                thisId = iFc;
                smallDist = thisDist;
            }
        }
        sf_facesIds[iC] = thisId;
        maxDist = max(maxDist,smallDist);
    }
    reduce(maxDist,maxOp<scalar>());

    // Calculating nodal area
	int id1 = 0;
	int id2 = 0;
	int id3 = 0;
	int id4 = 0;
	int iCf = 0;
	forAll(sf_faceAreas, iC)
	{
	    iCf = int(sf_facesIds[iC]);
	    // Getting indices
	    id1 = int(vc_faces[iCf*4]+intOff);
	    id2 = int(vc_faces[iCf*4+1]+intOff);
	    id3 = int(vc_faces[iCf*4+2]+intOff);
	    id4 = int(vc_faces[iCf*4+3]+intOff);
        // Adding area
	    sf_nodalArea[id1] += sf_faceAreas[iC]/4.0;
	    sf_nodalArea[id2] += sf_faceAreas[iC]/4.0;
	    sf_nodalArea[id3] += sf_faceAreas[iC]/4.0;
	    sf_nodalArea[id4] += sf_faceAreas[iC]/4.0;
    }
    reduce(sf_nodalArea,sumOp<scalarField>());

    // Little checks on mesh
    if (sf_calculateCorrections && (abs(gSum(p.magSf()*0.0+1)-vc_faceCentres.size())>0) | (maxDist>1.0E-2*sf_dxScale)){
		Info << "SynFlow: ATTENTION InflowCorrectionDict appears to be not corresponding to the present mesh." << endl;
	}
	else{
	    Info << "SynFlow: Checks on InflowCorrectionDict and mesh correspondence passed." << endl;
	}

    scalar avgFlux = sum((sf_xDir & sf_uAvg_n)*sf_nodalArea);
    scalar totFlux = sum((sf_xDir & (sf_uAvg_n + sf_uFlu_n))*sf_nodalArea);
    // Scaling flux
    sf_uAvg_n = sf_uAvg_n * avgFlux/totFlux;
    sf_uFlu_n = sf_uFlu_n * avgFlux/totFlux;

    // Setting small values
    sf_uSmall = 0.05*avgFlux/sum(sf_nodalArea);

    // Getting tProfile value
    tMult = max(0.01,interpolateXY(time,scalarField(sf_tProfile.component(0)),scalarField(sf_tProfile.component(1)*sf_tProfile.component(2))));
    // Calculating sf_uTot_c from nodal values (to be used as uStar at cell centres)
    forAll(sf_uTot_c, iC)
    {
	   iCf = int(sf_facesIds[iC]);
	   // Getting indices
	   id1 = int(vc_faces[iCf*4]+intOff);
	   id2 = int(vc_faces[iCf*4+1]+intOff);
	   id3 = int(vc_faces[iCf*4+2]+intOff);
	   id4 = int(vc_faces[iCf*4+3]+intOff);

	   sf_uTot_c[iC] = tMult*0.25*(sf_uAvg_n[id1] + sf_uFlu_n[id1] + sf_uAvg_n[id2] + sf_uFlu_n[id2] + sf_uAvg_n[id3] + sf_uFlu_n[id3] + sf_uAvg_n[id4] + sf_uFlu_n[id4]);
    }
    // Actually calculating field
    updateCoeffs();

}


Foam::synInflowFvPatchVectorField::
synInflowFvPatchVectorField
(
    const synInflowFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    sf_origin(ptf.sf_origin),
    sf_zDir(ptf.sf_zDir),
    sf_dictName(ptf.sf_dictName),
    vc_dictName(ptf.vc_dictName)
{}


Foam::synInflowFvPatchVectorField::
synInflowFvPatchVectorField
(
    const synInflowFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    sf_origin(ptf.sf_origin),
    sf_zDir(ptf.sf_zDir),
    sf_dictName(ptf.sf_dictName),
    vc_dictName(ptf.vc_dictName)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::synInflowFvPatchVectorField::updateCoeffs()
{
    // Already updated
	if ((sf_curTimeIndex==this->db().time().timeIndex()) | (this->db().time().timeIndex()==-1))
	{
	    return;
    }
    vectorField sf_uAvg_n = vectorField(vc_nodes.size(),vector::zero);
    vectorField sf_uFlu_n = vectorField(vc_nodes.size(),vector::zero);
    vector uvw_avg(vector::zero);
    vector uvw_flu(vector::zero);
	scalar intOff = 0.01;
	scalar tMult = 0.0;
	scalar time = this->db().time().value();
	scalar avgFlux = 0.0;
	scalar totFlux = 0.0;

    // Getting tProfile value
    tMult = max(0.01,interpolateXY(time,scalarField(sf_tProfile.component(0)),scalarField(sf_tProfile.component(1)*sf_tProfile.component(2))));
    // Part called only when corrections are not calculated
    if (!sf_calculateCorrections){
		forAll(vc_nodes, iN)
		   {
			uvw_avg = uvw_avg*0.0;
			uvw_flu = uvw_flu*0.0;
			#include "averageField.H"
			#include "fluctuatingField.H"
			uvw_flu = uvw_flu;
			uvw_avg = uvw_avg;
			sf_uAvg_n[iN] = uvw_avg;
			sf_uFlu_n[iN] = uvw_flu;
		   }

        // Caclulating overall fluxes (when corrections are not available values are face ones, not nodal)
        avgFlux = gSum((sf_xDir & sf_uAvg_n)*sf_faceAreas);
        totFlux = gSum((sf_xDir & (sf_uAvg_n + sf_uFlu_n))*sf_faceAreas);
        // Assigning BC
        fvPatchField<vector>::operator==(tMult*(sf_uAvg_n+sf_uFlu_n)*avgFlux/totFlux);
 		fixedValueFvPatchField<vector>::updateCoeffs();
 		Info << "SynFlow: Applied velocity BC without corrections with tProfile value " << tMult << endl;
 		sf_curTimeIndex = this->db().time().timeIndex();
 		return;
    }
    // Starting part called when corrections are present
    // Initializations
	vector valVec(vector::zero);
	scalar val = 0.0;
	int id1 = 0;
	int id2 = 0;
	int id3 = 0;
	int id4 = 0;
	int iCf = 0;
	scalar deltaT = this->db().time().deltaTValue();
	scalarField vc_bfluxes(vc_nodes.size(),0.0);
	scalarField vc_du1sdx1(vc_faceCentres.size(),0.0);
	scalarField vc_du2sdx2(vc_faceCentres.size(),0.0);
	scalarField vc_R(vc_nodes.size(),0.0);
	scalarField vc_lambda(vc_nodes.size(),0.0);
	vectorField vc_uvwCorr(vc_faceCentres.size(),vector::zero);
    vectorField uTot_c = vectorField(sf_faceC.size(),vector::zero);
	scalar uStar = 0.0;

	// For solver
	scalarField ls_b = vc_LdOp*0.0;

    // Preparing for parallelization
    int thisProcN = Pstream::myProcNo();
    int nNodes = vc_nodes.size();
    int nNodesPerProc = int(nNodes/Pstream::nProcs());
    int thisProcnNodes = nNodesPerProc;
    int iN = 0;

    if (thisProcN==(Pstream::nProcs()-1))
    {
         thisProcnNodes = nNodes - (nNodesPerProc*(Pstream::nProcs()-1));
	}
    // Calculating new synthetic nodal velocities
    for(int iNd=0;iNd<thisProcnNodes;iNd++)
       {
        iN = nNodesPerProc*thisProcN + iNd;
        uvw_avg = uvw_avg*0.0;
        uvw_flu = uvw_flu*0.0;
        #include "averageField.H"
        #include "fluctuatingField.H"
        uvw_flu = uvw_flu;
        uvw_avg = uvw_avg;
        sf_uAvg_n[iN] = uvw_avg;
        sf_uFlu_n[iN] = uvw_flu;
       }

	reduce(sf_uFlu_n,sumOp<vectorField>());
	reduce(sf_uAvg_n,sumOp<vectorField>());

	avgFlux = sum((sf_xDir & sf_uAvg_n)*sf_nodalArea);
	totFlux = sum((sf_xDir & (sf_uAvg_n + sf_uFlu_n))*sf_nodalArea);

	// Scaling flux and applying time multiplier
	sf_uAvg_n = sf_uAvg_n*tMult*avgFlux/totFlux;
	sf_uFlu_n = sf_uFlu_n*tMult*avgFlux/totFlux;

	// Computing fluxes from BCs
	forAll(vc_bcEdges, iC)
	{
		id1 = int(vc_bcEdges[iC][0]+intOff);
		id2 = int(vc_bcEdges[iC][1]+intOff);

		uvw_flu = 0.5*(sf_uFlu_n[id1] + sf_uFlu_n[id2]);

		val = vc_bcEdges[iC][2]*(uvw_flu & (vc_bcNorm[iC][0]*vc_t1 + vc_bcNorm[iC][1]*vc_t2));

		vc_bfluxes[id1] += 0.5*val;
		vc_bfluxes[id2] += 0.5*val;
	}
	// Calculating R
	// Computing gradients
	forAll(vc_grad1Op, iC)
	{
		id1 = int(vc_grad1Op[iC][0]+intOff);
		id2 = int(vc_grad1Op[iC][1]+intOff);
		vc_du1sdx1[id1] += vc_grad1Op[iC][2]*(sf_uFlu_n[id2] & vc_t1);
	}

	forAll(vc_grad2Op, iC)
	{
		id1 = int(vc_grad2Op[iC][0]+intOff);
		id2 = int(vc_grad2Op[iC][1]+intOff);
		vc_du2sdx2[id1] += vc_grad2Op[iC][2]*(sf_uFlu_n[id2] & vc_t2);
	}

	// Assembling R
	forAll(sf_faceC, iC)
	{
	   iCf = int(sf_facesIds[iC]);
	   // Getting indices
	   id1 = int(vc_faces[iCf*4]+intOff);
	   id2 = int(vc_faces[iCf*4+1]+intOff);
	   id3 = int(vc_faces[iCf*4+2]+intOff);
	   id4 = int(vc_faces[iCf*4+3]+intOff);
	   // Adding u related term
	   val = (0.25*(sf_uAvg_n[id1] + sf_uAvg_n[id2] + sf_uAvg_n[id3] + sf_uAvg_n[id4]) & sf_xDir);
	   val = val + (0.25*(sf_uFlu_n[id1] + sf_uFlu_n[id2] + sf_uFlu_n[id3] + sf_uFlu_n[id4]) & sf_xDir);
	   // Calculating convective velocity as average between old and new
	   uStar = max(0.5*(val + (sf_uTot_c[iC] & sf_xDir)),sf_uSmall);
	   val =  -((val - (sf_uTot_c[iC] & sf_xDir))/deltaT)/uStar + vc_du1sdx1[iCf] + vc_du2sdx2[iCf];

	   // Integrating and giving 1/4 to each node
	   vc_R[id1] += 0.25*val*sf_faceAreas[iC];
	   vc_R[id2] += 0.25*val*sf_faceAreas[iC];
	   vc_R[id3] += 0.25*val*sf_faceAreas[iC];
	   vc_R[id4] += 0.25*val*sf_faceAreas[iC];

	}

    // Reducing R
    reduce(vc_R,sumOp<scalarField>());

    // Solving for correction potential
    ls_b = vc_R - vc_bfluxes;
    ls_b = ls_b - sum(ls_b)/ls_b.size();
    // Fixing level at node 0 (depends on how the matrix is here provided)
    ls_b[0] = 0.0;
    // Solving
    #include "LUSolver.H"
    vc_lambda = ls_xxn;

    // Calculating corrections for v and w (u1 and u2 in local system)
    vc_uvwCorr = vc_uvwCorr*0.0;
    forAll(vc_grad1Op, iC)
    {
		id1 = int(vc_grad1Op[iC][0]+intOff);
		id2 = int(vc_grad1Op[iC][1]+intOff);
		vc_uvwCorr[id1] += vc_grad1Op[iC][2]*vc_lambda[id2]*vc_t1;
    }

    forAll(vc_grad2Op, iC)
    {
		id1 = int(vc_grad2Op[iC][0]+intOff);
		id2 = int(vc_grad2Op[iC][1]+intOff);
		vc_uvwCorr[id1] += vc_grad2Op[iC][2]*vc_lambda[id2]*vc_t2;
    }

    // Calculating corrected velocities
    sf_uTot_c = sf_uTot_c*0.0;
    forAll(sf_uTot_c, iC)
    {
	   iCf = int(sf_facesIds[iC]);
	   // Getting indices
	   id1 = int(vc_faces[iCf*4]+intOff);
	   id2 = int(vc_faces[iCf*4+1]+intOff);
	   id3 = int(vc_faces[iCf*4+2]+intOff);
	   id4 = int(vc_faces[iCf*4+3]+intOff);

	   // Adding average
	   valVec = 0.25*(sf_uAvg_n[id1] + sf_uAvg_n[id2] + sf_uAvg_n[id3] + sf_uAvg_n[id4]);
	   // Interpolation of fluctuating synthetic
	   valVec += 0.25*(sf_uFlu_n[id1] + sf_uFlu_n[id2] + sf_uFlu_n[id3] + sf_uFlu_n[id4]);
	   // Applying corrections
	   valVec += vc_corr*vc_uvwCorr[iCf];
	   // Setting minimal convective velocity, just to be sure
	   valVec = valVec - (valVec & sf_xDir)*sf_xDir + max(valVec & sf_xDir,sf_uSmall)*sf_xDir;
	   // Storing new BC
	   sf_uTot_c[iC] = valVec;
	   uTot_c[iC] = valVec;
    }
    // Applying field
    fvPatchField<vector>::operator==(uTot_c);
    fixedValueFvPatchField<vector>::updateCoeffs();
    Info << "SynFlow: Applied velocity BC with tProfile value " << tMult << endl;
    sf_curTimeIndex = this->db().time().timeIndex();

}


void Foam::synInflowFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntry("origin", sf_origin);
    os.writeEntry("zDir", sf_zDir);
    os.writeEntry("sf_dictName", sf_dictName);
    os.writeEntry("vc_dictName", vc_dictName);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       synInflowFvPatchVectorField
   );
}


// ************************************************************************* //
