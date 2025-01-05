using Project;
using System.Globalization;
using MathObjects;
using DataStructs;
using Solver;
CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;

const string CalculationArea = @"D:\CodeRepos\CS\Diplom\Data\Input\Info.dat";
const string BordersInfo = @"D:\CodeRepos\CS\Diplom\Data\Input\Borders.dat";
const string AnswerPath = @"D:\CodeRepos\CS\Diplom\Data\Output\Answer.dat";


FEM myFEM = new();
myFEM.ReadData(CalculationArea, BordersInfo);
myFEM.ConstructMesh();
myFEM.BuildMatrixAndVector();
myFEM.SetSolver(new MCG());
myFEM.Solve();
myFEM.WriteData(AnswerPath);