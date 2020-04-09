// When using the reference cell mapping approach and load balancing, it is possible to end up an situation
// that your number of chemically active cells changes suddenly due to something, e.g. restart or anomaly 
// in CFD solution. Then it is possible that load balancing gets stuck in gatherAndSendGuestData.H at
// "if(cells2SendCounter==nGuestCells)" because nGuestCells is much larger than any value of cells2SendCounter.
// We have to check before hand that nActiveCells on the current step is possibly over nGuestCells.

label nActiveCellsTot = 0;

if(this->useReferenceCell_)
{
    forAll(rho, celli)
    {
        for (label i=0; i<this->nSpecie_; i++)
        {
            this->c_[i]  = this->Y_[i][celli]; 
            c0[i] = this->c_[i];
        }
        
        //////////////////////////////////
        // Define error w.r.t. ref cell //
        //////////////////////////////////
        scalar c2cref_L2 = 0.0;
        scalar cref_L2   = 0.0;
        for (label i=0; i<nRefSp; i++)
        {
             c2cref_L2 += sqr(this->c_[refSpInds[i]] - refSpYs[i]);
             cref_L2   += sqr( refSpYs[i] );
        }
        c2cref_L2 = sqrt(c2cref_L2/cref_L2);
        ////////////////////////////////////
     
        //Reference cell check
        if(c2cref_L2 > 1e-4) // 1e-4 limit chosen by testing. e.g. 1e-5 is exceeded due processor boundary interpolations,
        {  
            nActiveCellsTot += 1;
        }
    }
}
else
{
    nActiveCellsTot += 1;
}

Info<<"nActiveCellsTot is : "<<nActiveCellsTot<<endl;
