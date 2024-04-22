using DataStructs;
using System.Numerics;

namespace Functions;

public static class Function
{
    // Функция правой части.
    public static double F(double r, double z, double t) => -8.0D * r - 6.0D * z + z * z * z / (r * r);

    public static (double, double, double) F(double x, double y, double z, double t) => (2.0D,
                                                                                         2.0D,
                                                                                         2.0D);


    // Тестируемая функция.
    public static double U(double r, double z, double t) => r * r * r + z * z * z;

    public static (double, double, double) A(double x, double y, double z, double t) => (2.0D,
                                                                                         2.0D,
                                                                                         2.0D);

/*
                    3
        ________________________
        |                      |
        |                      |
      4 |                      | 2
        |                      |
        |                      |
        |______________________|
                    1
*/



    #region Первые краевые условия.

    public static double U1_1(Point p, double t) => p.R * p.R * p.R + 1.0D;

    public static double U1_2(Point p, double t) => 8.0D + p.Z * p.Z * p.Z;

    public static double U1_3(Point p, double t) => p.R * p.R * p.R + 8.0D;

    public static double U1_4(Point p, double t) => 1.0D + p.Z * p.Z * p.Z;

    #endregion


    #region Вторые краевые условия.
    
    public static double dUdn1_1(Point p) => 0.0D;

    public static double dUdn1_2(Point p) => 0.0D;

    public static double dUdn1_3(Point p) => 0.0D;

    public static double dUdn1_4(Point p) => 0.0D;

    #endregion


    #region Третьи краевые условия.

    public static double U3dUdn3_1(Point p) => throw new Exception("Can't commit III bc.");

    public static double U3dUdn3_2(Point p) => throw new Exception("Can't commit III bc.");

    public static double U3dUdn3_3(Point p) => throw new Exception("Can't commit III bc.");

    public static double U3dUdn3_4(Point p) => throw new Exception("Can't commit III bc.");

    #endregion
}