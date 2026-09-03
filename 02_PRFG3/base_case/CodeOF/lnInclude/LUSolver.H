/* Implements the solution by using LU decomposition

Assumes the following inputs
ls_xxn = scalarField for solution
ls_Ld = scalarField for L diagonal
ls_Lo = vectorField containing L matrix in sparse format ordered for forward substitution
ls_Ud = scalarField for U diagonal
ls_Uo = vectorField containing U matrix in sparse format ordered for backward substitution
ls_LUpc = scalarField permutation for columns
ls_LUpr = scalarField permutation for raws
ls_b = scalarField for known term
*/

scalar ls_intOff = 0.01;
scalarField ls_xxn = ls_b*0.0;
scalarField ls_xx = ls_b*0.0;
scalarField ls_y = ls_b*0.0;

// Calculating intermediate vector y
forAll(ls_b,iN)
{
	id1 = int(ls_LUpc[iN]+ls_intOff);
    ls_y[id1] = ls_b[iN]/ls_Ld[iN];
}

forAll(ls_Lo,iN)
{
    id1 = int(ls_Lo[iN][0]+ls_intOff);
    id2 = int(ls_Lo[iN][1]+ls_intOff);
    ls_y[id1] -= ls_Lo[iN][2]*ls_y[id2]/ls_Ld[id1];
}

// Calculating solution
forAll(ls_y,iN)
{
    ls_xxn[iN] = ls_y[iN]/ls_Ud[iN];
}

forAll(ls_Uo,iN)
{
    id1 = int(ls_Uo[iN][0]+ls_intOff);
    id2 = int(ls_Uo[iN][1]+ls_intOff);
    ls_xxn[id1] -= ls_Uo[iN][2]*ls_xxn[id2]/ls_Ud[id1];
}

// Reordering
ls_xx = ls_xxn;
forAll(ls_xx,iN)
{
	id1 = int(ls_LUpr[iN]+ls_intOff);
    ls_xxn[iN] = ls_xx[id1];
}
