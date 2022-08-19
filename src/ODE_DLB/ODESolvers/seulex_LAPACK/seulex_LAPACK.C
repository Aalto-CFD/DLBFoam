/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | DLBFoam: Dynamic Load Balancing 
   \\    /   O peration     | for fast reactive simulations
    \\  /    A nd           | 
     \\/     M anipulation  | 2020, Aalto University, Finland
-------------------------------------------------------------------------------
License
    This file is part of DLBFoam library, derived from OpenFOAM.

    https://github.com/Aalto-CFD/DLBFoam

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

#include "seulex_LAPACK.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(seulex_LAPACK, 0);

    addToRunTimeSelectionTable(ODESolver, seulex_LAPACK, dictionary);

    const scalar
        seulex_LAPACK::stepFactor1_ = 0.8, // Safety factor for step control, Hairer default 0.8, OF default 0.6 
        seulex_LAPACK::stepFactor2_ = 0.93, //Safety factor for step control HNEW=H*sf2_*(sf1_*TOL/ERR)^(1/(J-1)) 
        seulex_LAPACK::stepFactor3_ = 0.1, // Step size selection: (sf3_*(1/(J-1)))/sF4_ <= hnew(J)/hold <= 1/(sf3_*(1/(J-1)))
        seulex_LAPACK::stepFactor4_ = 4.0, // Step size selection, see above
        seulex_LAPACK::stepFactor5_ = 0.5, // Step size hnew = sf5_*hold, whenever err_j >= err_j-1 and j>3
        seulex_LAPACK::kFactor1_ = 0.7, // Order is decreased if W(K-1) <= W(K)*kFactor1_
        seulex_LAPACK::kFactor2_ = 0.9; // Order is increased if W(K) <= W(K-1)*kFactor2_
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::seulex_LAPACK::seulex_LAPACK(const ODESystem& ode, const dictionary& dict)
:
    ODESolver(ode, dict),
    jacRedo_(min(1e-4, min(relTol_))),
    nSeq_(iMaxx_),
    cpu_(iMaxx_),
    coeff_(iMaxx_, iMaxx_),
    theta_(2*jacRedo_),
    table_(kMaxx_, n_),
    dfdx_(n_),
    dfdy_(n_),
    a_(n_),
    pivotIndices_(n_),
    dxOpt_(iMaxx_),
    temp_(iMaxx_),
    y0_(n_),
    ySequence_(n_),
    scale_(n_),
    dy_(n_),
    yTemp_(n_),
    dydx_(n_),
  	NASAP_mintemp(200),
	NASAP_maxtemp(5000)

{
    // The CPU time factors for the major parts of the algorithm
    const scalar cpuFunc = 1, cpuJac = 5, cpuLU = 1, cpuSolve = 1;

    nSeq_[0] = 2;
    nSeq_[1] = 3;

    for (int i=2; i<iMaxx_; i++)
    {
        nSeq_[i] = 2*nSeq_[i-2];
    }
    cpu_[0] = cpuJac + cpuLU + nSeq_[0]*(cpuFunc + cpuSolve);

    for (int k=0; k<kMaxx_; k++)
    {
        cpu_[k+1] = cpu_[k] + (nSeq_[k+1]-1)*(cpuFunc + cpuSolve) + cpuLU;
    }

    // Set the extrapolation coefficients array
    for (int k=0; k<iMaxx_; k++)
    {
        for (int l=0; l<k; l++)
        {
            scalar ratio = scalar(nSeq_[k])/nSeq_[l];
            coeff_(k, l) = 1/(ratio - 1);
        }
    }    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::seulex_LAPACK::seul
(
    const scalar x0,
    const scalarField& y0,
    const label li,
    const scalar dxTot,
    const label k,
    scalarField& y,
    const scalarField& scale
) const
{

    label nSteps = nSeq_[k];
    scalar dx = dxTot/nSteps;
    // LAPACK variables:
    lapack_int     npow2 = n_*n_;
    char    TRANS = 'N';
    lapack_int INFO = n_;
    lapack_int LDA  = n_;
    lapack_int LDB  = n_;
    lapack_int N = n_;
    lapack_int NRHS = 1;
    lapack_int IPIV[n_];
    // System matrix for LAPACK
    double A[npow2];
    memset( A, 0, (npow2)*sizeof(double) );
    // RHS vector for LAPACK
	double b[n_];
    memset( b, 0, (n_)*sizeof(double) );

    for (label j=0; j<n_; j++)
    {
        for (label i=0; i<n_; i++)
        {
            A[j*n_ + i] = -dfdy_[i][j];
        }
        A[j*n_ + j] += 1/dx;
    }

    dgetrf_(&N,&N,A,&LDA,IPIV,&INFO);
    if(INFO){std::cout << "an error occured in seul LAPACK LU phase: "<< INFO << endl << endl;}

    scalar xnew = x0 + dx;
    odes_.derivatives(xnew, y0,li,dy_);

    /* solve the linear system */
    for (label i=0; i<n_; i++)
    {
        b[i] = dy_[i];
    }
    dgetrs_(&TRANS,&N,&NRHS,A,&LDA,IPIV,b,&LDB,&INFO);
	if(INFO){std::cout << "an error occured in seul LAPACK BS 1: "<< INFO << endl;}
	for (label i=0; i<n_; i++)
	{
		dy_[i] = b[i];
	}
    /* solving linear system ends */

    yTemp_ = y0;



    for (label nn=1; nn<nSteps; nn++)
    {
        //Temp_ += dy_;
        for (label i=0; i<n_; i++)
	    {
	        yTemp_[i] = yTemp_[i] + dy_[i];
	    }
        
        xnew += dx;

        /*****************************************************/
        /* Stability check begins, original Hairer's         */ 
        /* implementation on fortran  lines 1067-1126 and    */ 
        /* idea  from book's p. 140 eq. 9.30                 */
        /*****************************************************/
        // check that temperature remains in range           */
   	    if (yTemp_[0] <= NASAP_mintemp || yTemp_[0] > NASAP_maxtemp)
        {
		    return false;
        }

        if (nn == 1 && k<=1)
        {
            scalar dy1 = 0.0;
            for (label i=0; i<n_; i++)
            {
                dy1 += sqr(dy_[i]/scale[i]);
            }
            dy1 = sqrt(dy1);

            odes_.derivatives(x0 + dx, yTemp_,li, dydx_);
            for (label i=0; i<n_; i++)
            {
                dy_[i] = dydx_[i] - dy_[i]/dx;
            }

            /* solve the linear system */
            for (label i=0; i<n_; i++)
            {
                b[i] = dy_[i];
            }
            dgetrs_(&TRANS,&N,&NRHS,A,&LDA,IPIV,b,&LDB,&INFO);
            if(INFO){std::cout << "an error occured in seul LAPACK BS 2: "<< INFO << endl;}
            for (label i=0; i<n_; i++)
            {
                dy_[i] = b[i];
            }
            /* solving linear system ends */

            scalar dy2 = 0.0;

	        for (label i=0; i<n_; i++)
	        {
               // check magnitude of change in solution vector, compared to VGREAT
	           if (mag(dy_[i]) > sqrt(vGreat) ){ dy2 = vGreat; break;}
		       dy2 += sqr(dy_[i]/scale[i]);
	        }

	    	dy2 = sqrt(dy2);
            // check convergence of Newton method: apply monotonicity check
		    theta_ = dy2/max(1.0, dy1 + small);
		    if (theta_ > 1.0 ) //dy2 > dy1 and convergence is not going to happen, return false
	    	{
		        return false;
	    	}

            /* The following monotonicity check differs from OF-2.4.x but follows Hairer's original */
            /*
            const scalar denom = max(1, dy1);
            scalar dy2 = 0;
            for (label i=0; i<n_; i++)
            {
                // Test of dy_[i] to avoid overflow
                if (mag(dy_[i]) > scale[i]*denom)
                {
                    theta_ = 1;
                    return false;
                }

                dy2 += sqr(dy_[i]/scale[i]);
            }
            dy2 = sqrt(dy2);
            theta_ = dy2/denom;

		    if (theta_ > 1.0 ) //dy2 > dy1 and convergence is not going to happen, return false
	    	{
		        return false;
	    	}
            */
        }

        odes_.derivatives(xnew, yTemp_,li, dy_);

        /* solve the linear system */
        for (label i=0; i<n_; i++)
        {
            b[i] = dy_[i];
        }
        dgetrs_(&TRANS,&N,&NRHS,A,&LDA,IPIV,b,&LDB,&INFO);
        if(INFO){std::cout << "an error occured in seul LAPACK BS 3: "<< INFO << endl;}
        for (label i=0; i<n_; i++)
        {
            dy_[i] = b[i];
        }
        /* solving linear system ends */
    }

    for (label i=0; i<n_; i++)
    {
        y[i] = yTemp_[i] + dy_[i];
    }

    return true;
}

/* Polyomial extrapolation, originally on lines 1235-> in fortran code */
/* and eq. 9.26 in the book                                            */
void Foam::seulex_LAPACK::extrapolate
(
    const label k,
    scalarRectangularMatrix& table,
    scalarField& y
) const
{
	for (int j=k-1; j>0; j--)
	{
		for (label i=0; i<n_; i++)
		{
			table[j-1][i] = table[j][i] + coeff_[k][j]*(table[j][i] - table[j-1][i]);
		}
	}

	for (int i=0; i<n_; i++)
	{
		y[i] = table[0][i] + coeff_[k][0]*(table[0][i] - y[i]);
	}
}


bool Foam::seulex_LAPACK::resize()
{
    if (ODESolver::resize())
    {
        table_.shallowResize(kMaxx_, n_);
        resizeField(dfdx_);
        resizeMatrix(dfdy_);
        resizeMatrix(a_);
        resizeField(pivotIndices_);
        resizeField(y0_);
        resizeField(ySequence_);
        resizeField(scale_);
        resizeField(dy_);
        resizeField(yTemp_);
        resizeField(dydx_);

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::seulex_LAPACK::solve
(
    scalar& x,
    scalarField& y,
    const label li,
    stepState& step
) const
{
    temp_[0] = GREAT;
    scalar dx = step.dxTry;
    y0_ = y;
    dxOpt_[0] = mag(0.1*dx);

    if (step.first || step.prevReject)
    {
        theta_ = 2*jacRedo_;
    }

    if (step.first)
    {
        // NOTE: the first element of relTol_ and absTol_ are used here.
        scalar logTol = -log10(relTol_[0] + absTol_[0])*0.6 + 0.5;
        kTarg_ = max(1, min(kMaxx_ - 1, int(logTol)));
    }

    forAll(scale_, i) // defined similarly in fortran on line 874 by Hairer
    {
        scale_[i] = absTol_[i] + relTol_[i]*mag(y[i]);
    }

    bool jacUpdated = false;

    if (theta_ > jacRedo_)
    {
        odes_.jacobian(x, y, li, dfdx_, dfdy_);
        jacUpdated = true;
    }

    int k;
    scalar dxNew = mag(dx);
    bool firstk = true;

    while (firstk || step.reject)
    {
        dx = step.forward ? dxNew : -dxNew;
        firstk = false;
        step.reject = false;

        if (mag(dx) <= mag(x)*sqr(SMALL))
        {
             WarningInFunction
                    << "step size underflow :"  << dx << endl;
        }

        scalar errOld = 0;


		for (k=0; k<=kTarg_+1; k++)
		{


            bool success = seul(x, y0_, li, dx, k, ySequence_, scale_);

            if (!success)
			{
                if(debug) Info << "seul function returned false --> dx is halved" << endl;
				step.reject = true;
				dxNew = mag(dx)*stepFactor5_;
				break;
			}
			if (k == 0)
			{
				forAll(y,i)
				{
					y[i] = ySequence_[i];
				}
			}
			else
			{
				forAll(ySequence_, i)
				{
					table_[k-1][i] = ySequence_[i];
				}
			}
			if (k != 0)
			{
				extrapolate(k, table_, y);
                
                /* Compute the optimal step sizes   */
                /* Corresponds fortran lines 1144-> */
				scalar err = 0.0;
				forAll(y, i)
				{
					scale_[i] = absTol_[i] + relTol_[i]*mag(y0_[i]);
					err += sqr((y[i] - table_[0][i])/scale_[i]);
				}
				err = sqrt(err/n_);
				if (err > 1.0/SMALL || (k > 1 && err >= errOld))
				{
					step.reject = true;
					dxNew = mag(dx)*stepFactor5_;
					break;
				}
				// compute optimal timestep estimation
				errOld = max(4.0*err, 1);
				scalar expo = 1.0/(k + 1);
				scalar facmin = pow(stepFactor3_, expo);
				scalar fac;
				//Despite the peculiar use of if/else, this is equivalent
                //to the original fortran implementation on line 1155
				if (err == 0.0)
				{
					fac = 1.0/facmin;
				}
				else
				{
					fac = stepFactor2_/pow(err/stepFactor1_, expo);
					fac = max(facmin/stepFactor4_, min(1.0/facmin, fac));
				}
				dxOpt_[k] = mag(dx*fac);
				temp_[k] = cpu_[k]/dxOpt_[k];

                /* Computing the optimal order... */
                /* not yet sure the fortran lines */
				if ((step.first || step.last) && err <= 1.0)
				{
					break;
				}
				if(
						k == kTarg_ - 1
						&& !step.prevReject
						&& !step.first && !step.last
				  )
				{
					if (err <= 1.0)
					{
						break;
					}
					else if (err > nSeq_[kTarg_]*nSeq_[kTarg_ + 1]*4.0)
					{
						step.reject = true;
						kTarg_ = k;
						if (kTarg_>1 && temp_[k-1] < kFactor1_*temp_[k])
						{
							kTarg_--;
						}
						dxNew = dxOpt_[kTarg_];
						break;
					}
				}
				if (k == kTarg_)
				{
					if (err <= 1.0)
					{
						break;
					}
					else if (err > nSeq_[k + 1]*2.0)
					{
						step.reject = true;
						if (kTarg_>1 && temp_[k-1] < kFactor1_*temp_[k])
						{
							kTarg_--;
						}
						dxNew = dxOpt_[kTarg_];
						break;
					}
				}
				if (k == kTarg_+1)
				{
					if (err > 1.0)
					{
						step.reject = true;
						if
							(
								kTarg_ > 1
								&& temp_[kTarg_-1] < kFactor1_*temp_[kTarg_]
								)
						{
							kTarg_--;
						}
						dxNew = dxOpt_[kTarg_];
					}
					break;
				}
			}
		}
		if (step.reject)
		{
			step.prevReject = true;
			if (!jacUpdated)
			{
				theta_ = 2.0*jacRedo_;
				if (theta_ > jacRedo_ && !jacUpdated)
				{

   					odes_.jacobian(x, y, li, dfdx_, dfdy_);
					jacUpdated = true;
	                

		    }
			}
		}
	}
	jacUpdated = false;
	step.dxDid = dx;
	x += dx;
	label kopt;
    // Compute the optimal order kopt
	if (k == 1)
	{
		kopt = 2;
	}
	else if (k <= kTarg_)
	{
		kopt=k;
		if (temp_[k-1] < kFactor1_*temp_[k])
		{
			kopt = k - 1;
		}
		else if (temp_[k] < kFactor2_*temp_[k - 1]) //line 951
		{
			kopt = min(k + 1, kMaxx_ - 1); //see line 952
		}
	}
	else
	{
		kopt = k - 1;
		if (k > 2 && temp_[k-2] < kFactor1_*temp_[k - 1])
		{
			kopt = k - 2;
		}
		if (temp_[k] < kFactor2_*temp_[kopt])
		{
			kopt = min(k, kMaxx_ - 1);
		}
	}
	if (step.prevReject)
	{
		kTarg_ = min(kopt, k);
		dxNew = min(mag(dx), dxOpt_[kTarg_]);
		step.prevReject = false;
	}
	else
	{
		if (kopt <= k)
		{
			dxNew = dxOpt_[kopt];
		}
		else
		{
			if (k < kTarg_ && temp_[k] < kFactor2_*temp_[k - 1])
			{
				dxNew = dxOpt_[k]*cpu_[kopt + 1]/cpu_[k];
			}
			else
			{
				dxNew = dxOpt_[k]*cpu_[kopt]/cpu_[k];
			}
		}
		kTarg_ = kopt;
	}

	step.dxTry = step.forward ? dxNew : -dxNew; //if you would force here a too large dxNew, convergence deteriorates
 
}



// ************************************************************************* //
