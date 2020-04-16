// When using the reference cell mapping approach and load balancing, it is possible to end up an situation
// that your number of chemically active cells changes suddenly due to something, e.g. restart or anomaly
// in CFD solution. Then it is possible that load balancing gets stuck in gatherAndSendGuestData.H at
// "if(cells2SendCounter==nGuestCells)" because nGuestCells is much larger than any value of cells2SendCounter.
// We have to check before hand that nActiveCells on the current step is possibly over nGuestCells.

label nActiveCellsTot = 0;

if(refcell_mapper_->active())
{
    forAll(rho, celli)
    {
        if(!refcell_mapper_-> applyMapping(Y_,celli))
        {
            nActiveCellsTot +=1;
        }
    }
}
else
{
    nActiveCellsTot = this->mesh().cells().size();
}

reduce(nActiveCellsTot, sumOp<label>());
Info<<"nActiveCellsTot is : "<<nActiveCellsTot<<endl;
