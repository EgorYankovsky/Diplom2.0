namespace Functions;

public static class BasisFunctions2D
{
    /*
      q3               q4
        *-------------*
        |             |
        |             |
        |             |
        |             |
        *-------------*
      q1               q2
        ^z
        |
        |
        ----->r

    */
    public static double GetValue(double q1, double q2, double q3, double q4,
                                  double r0, double r1, double z0, double z1,
                                  double ri, double zi) 
        => q1 * R1(r1, r0, ri) * Z1(z1, z0, zi) + q2 * R2(r1, r0, ri) * Z1(z1, z0, zi) +
           q3 * R1(r1, r0, ri) * Z2(z1, z0, zi) + q4 * R2(r1, r0, ri) * Z2(z1, z0, zi);

    public static double R1(double r1, double r0, double r) => (r1 - r) / (r1 - r0);

    public static double R2(double r1, double r0, double r) => (r - r0) / (r1 - r0);

    public static double Z1(double z1, double z0, double z) => (z1 - z) / (z1 - z0);

    public static double Z2(double z1, double z0, double z) => (z - z0) / (z1 - z0);

    public static double dR1(double r1, double r0, double r) => -1.0D / (r1 - r0);

    public static double dR2(double r1, double r0, double r) => 1.0D / (r1 - r0);

    public static double dZ1(double z1, double z0, double z) => -1.0D / (z1 - z0);

    public static double dZ2(double z1, double z0, double z) => 1.0D / (z1 - z0);
}


public static class BasisFunctions3D
{
    /*
          q7                q8
            *-------------*
           /|            /|
          / |           / |
         /  |          /  |
        *-------------* q6|
      q5|   *-------- | --*
        |  / q3       |  / q4
        | /           | /
        |/            |/
        *-------------*
      q1               q2
    
    z
    ^  ^y
    | / 
    |/
    |---->x
    
    */
    public static double GetValue(double q1, double q2, double q3, double q4,
                                  double q5, double q6, double q7, double q8,
                                  double x0, double x1, double y0, double y1,
                                  double z0, double z1, double xi, double yi,
                                  double zi)
        => q1 * X1(x1, x0, xi) * Y1(y1, y0, yi) * Z1(z1, z0, zi) + q2 * X2(x1, x0, xi) * Y1(y1, y0, yi) * Z1(z1, z0, zi) +
               q3 * X1(x1, x0, xi) * Y2(y1, y0, yi) * Z1(z1, z0, zi) + q4 * X2(x1, x0, xi) * Y2(y1, y0, yi) * Z1(z1, z0, zi) +
               q5 * X1(x1, x0, xi) * Y1(y1, y0, yi) * Z2(z1, z0, zi) + q6 * X2(x1, x0, xi) * Y1(y1, y0, yi) * Z2(z1, z0, zi) +
               q7 * X1(x1, x0, xi) * Y2(y1, y0, yi) * Z2(z1, z0, zi) + q8 * X2(x1, x0, xi) * Y2(y1, y0, yi) * Z2(z1, z0, zi);

    public static double X1(double x1, double x0, double x) => (x1 - x) / (x1 - x0);
    public static double X2(double x1, double x0, double x) => (x - x0) / (x1 - x0);
    public static double Y1(double y1, double y0, double y) => (y1 - y) / (y1 - y0);
    public static double Y2(double y1, double y0, double y) => (y - y0) / (y1 - y0);
    public static double Z1(double z1, double z0, double z) => (z1 - z) / (z1 - z0);
    public static double Z2(double z1, double z0, double z) => (z - z0) / (z1 - z0);
}

public static class BasisFunctions3DVec
{

    public static (double, double, double) GetValue(double eps, double nu, double khi, double[] q)
    {
        if (q.Length != 12) throw new ArgumentException("q doesn't contain 12 values");
        return (BF1(q[0], eps, nu, khi) + BF2(q[1], eps, nu, khi) + BF3(q[2], eps, nu, khi) + BF4(q[3], eps, nu, khi),
                BF5(q[4], eps, nu, khi) + BF6(q[5], eps, nu, khi) + BF7(q[6], eps, nu, khi) + BF8(q[7], eps, nu, khi),
                BF9(q[8], eps, nu, khi) + BF10(q[9], eps, nu, khi) + BF11(q[10], eps, nu, khi) + BF12(q[11], eps, nu, khi));
    }


    //                                          x           y          z
    private static double BF1(double q, double eps, double nu, double khi) => q * LF1(nu) * LF1(khi);
    private static double BF2(double q, double eps, double nu, double khi) => q * LF2(nu) * LF1(khi);
    private static double BF3(double q, double eps, double nu, double khi) => q * LF1(nu) * LF2(khi);
    private static double BF4(double q, double eps, double nu, double khi) => q * LF2(nu) * LF2(khi);


    private static double BF5(double q, double eps, double nu, double khi) => q * LF1(eps) * LF1(khi);
    private static double BF6(double q, double eps, double nu, double khi) => q * LF2(eps) * LF1(khi);
    private static double BF7(double q, double eps, double nu, double khi) => q * LF1(eps) * LF2(khi);
    private static double BF8(double q, double eps, double nu, double khi) => q * LF2(eps) * LF2(khi);


    private static double BF9(double q, double eps, double nu, double khi) => q * LF1(eps) * LF1(nu);
    private static double BF10(double q, double eps, double nu, double khi) => q * LF2(eps) * LF1(nu);
    private static double BF11(double q, double eps, double nu, double khi) => q * LF1(eps) * LF2(nu);
    private static double BF12(double q, double eps, double nu, double khi) => q * LF2(eps) * LF2(nu);

    private static double LF1(double v) => 1 - v;
    private static double LF2(double v) => v;
}