private static void Add(LocalMatrix lm, ref GlobalMatrix gm, List<int> elem)
{
    int ii = 0;
    foreach (var i in elem)
    {
        int jj = 0;
        foreach (var j in elem)
        {
            int ind = 0;
            double val = 0.0D;
            switch(i - j)
            {
                case 0:
                    val = lm[ii, jj];
                    gm._diag[i] += lm[ii, jj];
                    break;
                case < 0:
                    ind = gm._ig[j];
                    for (; ind <= gm._ig[j + 1] - 1; ind++)
                        if (gm._jg[ind] == i) break;
                    val = lm[ii, jj];
                    gm._au[ind] += lm[ii, jj];
                    break;
                case > 0:
                    ind = gm._ig[i];
                    for (; ind <= gm._ig[i + 1] - 1; ind++)
                        if (gm._jg[ind] == j) break;
                    val = lm[ii, jj];
                    gm._al[ind] += lm[ii, jj];
                    break;
            }
            jj++;
        }
        ii++;
    }
}

private static void Add(LocalVector lv, ref GlobalVector gv, List<int> elems)
{
    for (int i = 0; i < elems.Count; i++)
        gv[elems[i]] += lv[i];
}