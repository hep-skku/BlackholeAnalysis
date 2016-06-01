void dfactor()
{
  const double nu[9] = {0.579,0.517,0.473,0.742,0.558,-0.339,0.711,0.666,0.533}; 
  const double pi = TMath::Pi();
  const double kL = 11.3*pi;
  double phase[9];
  for(int i=0; i < 9; ++i) 
  {
    phase[i] = 1-2*nu[i];
    double dFactor = exp(-2*kL)*phase[i]*4.*((exp((1.+phase[i]/2.)*kL)-1)**2)/pi/(exp(kL*phase[i])-1)/((2+phase[i])**2);
    cout << pi*dFactor << endl;
  }
}
