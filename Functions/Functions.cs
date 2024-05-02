using DataStructs;
using System.Numerics;

namespace Functions;

public static class Function
{
    // Функция правой части.
    public static double F(double r, double z, double t) => Math.Pow(10.0D, 2);

    public static (double, double, double) F(double x, double y, double z, double t) => (y * y - 2.0D,
                                                                                         z * z - 2.0D,
                                                                                         x * x - 2.0D);


    // Тестируемая функция.
    public static double U(double r, double z, double t) => r * r * r + z * z * z;

    public static (double, double, double) A(double x, double y, double z, double t) => (y * y,
                                                                                         z * z,
                                                                                         x * x);

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

    public static double U1_1(Point2D p, double t) => 0.0D;

    public static double U1_2(Point2D p, double t) => 0.0D;

    public static double U1_3(Point2D p, double t) => 0.0D;

    public static double U1_4(Point2D p, double t) => 0.0D;

    #endregion


    #region Вторые краевые условия.
    
    public static double dUdn1_1(Point2D p) => 0.0D;

    public static double dUdn1_2(Point2D p) => 0.0D;

    public static double dUdn1_3(Point2D p) => 0.0D;

    public static double dUdn1_4(Point2D p) => 0.0D;

    #endregion


    #region Третьи краевые условия.

    public static double U3dUdn3_1(Point2D p) => throw new Exception("Can't commit III bc.");

    public static double U3dUdn3_2(Point2D p) => throw new Exception("Can't commit III bc.");

    public static double U3dUdn3_3(Point2D p) => throw new Exception("Can't commit III bc.");

    public static double U3dUdn3_4(Point2D p) => throw new Exception("Can't commit III bc.");

    #endregion
}