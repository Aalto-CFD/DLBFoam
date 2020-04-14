    //- Initialize the time progress    
    scalar timeLeft = deltaT[celli];

    // Calculate the chemical source terms
    while (timeLeft > small)
    {
        scalar dt = timeLeft;
       
        this->solve(this->c_, Ti, pi, dt, this->deltaTChem_[celli]);
        timeLeft -= dt;
    }

    deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

    //RR_ should be weighted by density, due to NS formulation. rr = Y0-Y1 / dt, but RR_ = rho(Y0-Y1)/dt , the effect at Sh() function too!
    for (label i=0; i<this->nSpecie_; i++)
    {
       
        this->RR_[i][celli] = rhoi*(this->c_[i] - c0[i])/deltaT[celli];
        
    }

