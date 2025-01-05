public void GenerateVectorEphi()
{
    E_phi = new GlobalVector[A_phi.Length];   
    for (int i = 0; i < E_phi.Length; i++)
    {
        if (i == 0 || i == 1)
            E_phi[i] = new GlobalVector(A_phi[i].Size);
        else
        {
            double ti = Time[i];
            double ti_1 = Time[i - 1];
            double ti_2 = Time[i - 2];
            E_phi[i] = -1.0D / (ti_1 - ti_2) * A_phi[i - 2] + (ti - ti_2) / ((ti_1 - ti_2) * (ti - ti_1)) * A_phi[i - 1] - 
            (2 * ti - ti_1 - ti_2) / ((ti - ti_2) * (ti - ti_1)) * A_phi[i];
        }
    }
}