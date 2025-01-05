public static void ConsiderBoundaryConditions(ref GlobalMatrix gm, ref GlobalVector gv, ArrayOfPoints2D arrp, ArrayOfBorders arrBd, double t)
{
    foreach (var border in arrBd)
    {
        switch (border[0])
        {
            // КУ - I-го рода
            case 1:
            for (int i = 2; i < 4; i++)
            {
                for (int j = gm._ig[border[i]]; j < gm._ig[border[i] + 1]; j++)
                    gm._al[j] = 0.0D;
                gm._diag[border[i]] = 1.0D;
                for (int j = 0; j < gm._jg.Count; j++)
                    if (gm._jg[j] == border[i])
                        gm._au[j] = 0.0D;
            }
        
            // Гаусс для симметрии:
            if (border[0] == 1)
            {
                for (int i = 2; i < 4; i ++)
                {
                    for (int j = gm._ig[border[i]]; j < gm._ig[border[i] + 1]; j++)
                    {
                        var k = gm._au[j];
                        var f = -1.0D * k * gv[border[i]];
                        gv[gm._jg[j]] += f;
                        gm._au[j] = 0.0D;
                    }
                    for (int j = 0; j < gm._jg.Count; j++)
                        if (gm._jg[j] == border[i])
                        {
                            var k = gm._al[j];
                            var f = -1.0D * k * gv[border[i]];
                            gv[gm._jg[j]] += f;
                            gm._al[j] = 0.0D;
                        }
                }
            }
        
            switch (border[1])
            {
                case 1:
                for (int i = 2; i < 4; i++)
                    gv[border[i]] = U1_1(arrp[border[i]], t);
                break;

                case 2:
                for (int i = 2; i < 4; i++)
                    gv[border[i]] = U1_2(arrp[border[i]], t);
                break;
                case 3:
                for (int i = 2; i < 4; i++)
                    gv[border[i]] = U1_3(arrp[border[i]], t);
                break;

                case 4:
                for (int i = 2; i < 4; i++)
                    gv[border[i]] = U1_4(arrp[border[i]], t);
                break;

                default:
                throw new Exception("No such bord");
            }
            break;
           // КУ - II-го рода
           case 2:
               for (int i = 2; i < 4; i++)
                   gv[border[i]] += 0.0D;
               break;
           // КУ - III-го рода
           case 3: throw new ArgumentException("Пока нет возможности учитывать КУ III-го рода");
        }
    }
}