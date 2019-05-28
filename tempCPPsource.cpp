//----------------------------------------------
dim1 nodesOrder;
for(int i = 0; i < N; i++)
    nodesOrder.push_back(i);

double lastPushMeanYPrime, OldOmega, NewOmega;

for(int step = 0; step < num_steps; ++step){
    std::random_shuffle(nodesOrder.begin(), nodesOrder.end());
    sumAcceptanceRewirig = 0;
    lastPushMeanYPrime = Mean(dydt(y, Cij));
    MeanYPrime.push_back(lastPushMeanYPrime);
    /// local strategy
    for(int j=0; j<N; j++){
        OldOmega = Mean(MeanYPrime, true);
        NewCij = rewiring(nodesOrder[j], nodesOrder);
        NewY = runge_kutta4_integrator(y, NewCij);
        lastPushMeanYPrime = MeanYPrime.back();
        MeanYPrime.pop_back();
        MeanYPrime.push_back(Mean(dydt(NewY, NewCij)));
        NewOmega = Mean(MeanYPrime, true);
        if(abs(Omega[nodesOrder[j]]-OldOmega) > abs(Omega[nodesOrder[j]]-NewOmega) ){
            Cij = NewCij;
            y = NewY;
            sumAcceptanceRewirig++;
        }
        else{
            MeanYPrime.pop_back();
            MeanYPrime.push_back(lastPushMeanYPrime);
        }
    }

    //global strategy

    for(int j=0; j<N; j++){
  RGlobalBeforeRewiring = order_parameter(y);
  NewCij = rewiring(nodesOrder[j], nodesOrder);
  NewY = runge_kutta4_integrator(y, NewCij);
  RGlobalAfterRewiring = order_parameter(NewY);
  if(RGlobalAfterRewiring[0] > RGlobalBeforeRewiring[0]){
    Cij = NewCij;
    sumAcceptanceRewirig++;
  }
}

    AcceptanceRateRewiring.push_back(sumAcceptanceRewirig);
    //std::cout<<"step="<<step<<"\n";
    r1 = order_parameter(y);
    //std::cout<<"\n"<<"r1="<<r1[0];
    r2 = order_parameter_k(y);
    Order1.push_back(r1[0]);
    Psi1.push_back(r1[1]);
    Order2.push_back(r2[0]);
    Psi2.push_back(r2[1]);

    y = runge_kutta4_integrator(y, Cij);
}
//---------------------------------------------------




/*------------------------------------------------------------*/
dim2 ODE::rewiring(int indexFocusNode, dim1 nodesOrder){
    std::random_shuffle(nodesOrder.begin(), nodesOrder.end());
    bool zeroToOne, OneToZero;
    zeroToOne = 0;
    OneToZero = 0;
    for(int i=0; i<N && (!zeroToOne) && (!OneToZero); i++){
        if(Cij[indexFocusNode][nodesOrder[i]]==1 && !OneToZero){
            Cij[indexFocusNode][nodesOrder[i]] = 0;
            OneToZero = 1;
        }else if(Cij[indexFocusNode][nodesOrder[i]]==0 && !zeroToOne){
                 Cij[indexFocusNode][nodesOrder[i]] = 1;
                 zeroToOne = 1;
        }
    }
    return Cij;
}
/*------------------------------------------------------------*/
