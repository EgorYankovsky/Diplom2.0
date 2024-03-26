using MathObjects;
/*
  | --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- |    | --- |      | ---- |
  | 1.0	0.0			0.0	0.0										                          |    | 1.0 |      |  1.0 |
  | 0.0	1.0	0.0		0.0	0.0	0.0									                      |    | 2.0 |      |  2.0 |
  |     0.0	1.0	0.0		0.0	0.0	0.0								                    |    | 3.0 |      |  3.0 |
  |         0.0	1.0			0.0	0.0								                      |    | 4.0 |      |  4.0 |
  | 0.4	0.3			3.0	0.2			0.3	0.4						                      |    | 5.0 |      | 22.4 |
  | 0.2	0.1	0.1		0.2	2.0	0.2		0.3	0.2	0.1					                |    | 6.0 |      | 19.5 |
  |     0.3	0.4	0.3		0.2	1.0	0.2		0.2	0.3	0.4				              |    | 7.0 |      | 18.0 |
  |         0.1	0.1			0.2	2.0			0.2	0.1				                  |    | 8.0 | ---- | 19.8 |
  |                 0.3	0.4			3.0	0.1			0.2	0.3		              | *  | 8.0 | ---- | 30.3 |
  |                 0.3	0.2	0.1		0.1	4.0	0.3		0.4	0.3	0.2	        |    | 7.0 |      | 36.9 |
  |                     0.2	0.3	0.4		0.3	3.0	0.1		0.1	0.2	0.3     |    | 6.0 |      | 28.1 |
  |                         0.2	0.1			0.1	2.0			0.3	0.2         |    | 5.0 |      | 13.6 |
  |                                 0.2	0.3			1.0	0.1		          |    | 4.0 |      |  8.0 |
  |                                 0.4	0.3	0.2		0.1	2.0	0.4	      |    | 3.0 |      | 13.7 |
  |                                     0.1	0.2	0.3		0.4	3.0	0.1   |    | 2.0 |      | 10.7 |
  |                                         0.3	0.2			0.1	4.0     |    | 1.0 |      |  7.0 |
  | --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- |    | --- |      | ---- |
                                A                                          x            b

*/

public static class TestClass
{
    public static int[] ig = {0, 0, 1, 2, 3, 5, 9, 13, 16, 18, 22, 26, 29, 31, 35, 39, 42};

    public static List<int> jg = new()
                    {  0,  1,  2,  0,  1,  0,  1,
                       2,  4,  1,  2,  3,  5,  2,
                       3,  6,  4,  5,  4,  5,  6,
                       8,  5,  6,  7,  9,  6,  7,
                      10,  8,  9,  8,  9, 10, 12,
                       9, 10, 11, 13, 10, 11, 14};

    public static double[] diag = {1.0, 1.0, 1.0, 1.0,
                                   3.0, 2.0, 1.0, 2.0,
                                   3.0, 4.0, 3.0, 2.0,
                                   1.0, 2.0, 3.0, 4.0};

    public static double[] al = {0.0, 0.0, 0.0, 0.4, 0.3, 0.2, 0.1,
                                 0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1,
                                 0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1,
                                 0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1,
                                 0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1,
                                 0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1};

    public static double[] au = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                 0.0, 0.2, 0.0, 0.0, 0.0, 0.2, 0.0,
                                 0.0, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1,
                                 0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1,
                                 0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1,
                                 0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1};

    public static GlobalMatrix A = new (ig, jg, diag, al, au);

    public static GlobalVector b = new (new double[] { 1.0,  2.0,  3.0,  4.0,
                                                      22.4, 19.5, 18.0, 19.8,
                                                      30.3, 36.9, 28.1, 13.6,
                                                       8.0, 13.7, 10.7,  7.0});
}