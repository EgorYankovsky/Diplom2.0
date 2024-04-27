using Project;
using System.Globalization;
using Dumpify;
using Solver;
using Manager;
using System.Numerics;
using Processor;
using Solution;
using MathObjects;
using Functions;
using Grid;
using DataStructs;
using static Grid.MeshReader;
using static Grid.MeshGenerator;
using static Manager.FolderManager;

CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;
DumpConfig.Default.ColorConfig.ColumnNameColor = DumpColor.FromHexString("#FF0000");
DumpConfig.Default.ColorConfig.PropertyValueColor = DumpColor.FromHexString("#00FF00");
DumpConfig.Default.ColorConfig.NullValueColor = DumpColor.FromHexString("#0000FF");

string CalculationArea = Path.GetFullPath("../../../../Data/Input/WholeMesh.txt");
string LayersArea = Path.GetFullPath("../../../../Data/Input/Fields.txt");
string BordersInfo = Path.GetFullPath("../../../../Data/Input/Borders.txt");
string TimePath = Path.GetFullPath("../../../../Data/Input/Time.txt");
string SubtotalsPath = Path.GetFullPath("../../../../Data/Subtotals/");
string AnswerPath = Path.GetFullPath("../../../../Data/Output/");
string PicturesPath = Path.GetFullPath("../../../../Drawer/Pictures/");

//var A = TestClass.A;
//var v = TestClass.b;
//var Solver = new LU_LOS();
//var x = Solver.Solve(A, v);
//
//for (int i = 0; i < x.Item1.Size; i++)
//    Console.WriteLine(x.Item1[i]);
//
//return 0;
bool isSolving2DimTask = Checker(AnswerPath);

// Pre-processor. Clearing output folders.

if (isSolving2DimTask)
    ClearFolders(new List<string> {SubtotalsPath + "/2_dim/",
                                   PicturesPath + "/E_phi/",
                                   PicturesPath + "/A_phi/",
                                   AnswerPath + "/E_phi/Discrepancy/",
                                   AnswerPath + "/E_phi/Answer/",
                                   AnswerPath + "/A_phi/Discrepancy/",
                                   AnswerPath + "/A_phi/Answer/"});


/*
int nx = 4;
int ny = 4;
int nz = 3;

List<double> nodesX = [0.0, 1.0, 2.0, 3.0];
List<double> nodesY = [0.0, 1.0, 2.0, 3.0];
List<double> nodesZ = [0.0, 1.0, 2.0];

double xMin = nodesX[0];
double xMax = nodesX[^1];
double yMin = nodesY[0];
double yMax = nodesY[^1];
double zMin = nodesZ[0];
double zMax = nodesZ[^1];

ArrayOfPoints arrpnt = new(nx * ny * nz);

foreach (var Z in nodesZ)
    foreach (var Y in nodesY)
        foreach (var X in nodesX)
            arrpnt.Append(new Point(X, Y, Z));

foreach (var pnt in arrpnt)
{
    if (pnt.Z == zMax)
        pnt.Type = Location.BoundaryII;
    else if (pnt.Z == zMin || pnt.X == xMin ||pnt.X == xMax || pnt.Y == yMin || pnt.Y == yMax)
        pnt.Type = Location.BoundaryI;
    else
        pnt.Type = Location.Inside;
}

ArrayOfRibs arr = new(3 * nx * ny * nz - nx * ny - nx * nz - ny * nz);

int nxny = nx * ny;

for (int k = 0; k < nz; k++)
{
    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx - 1; i++) 
            arr.Add(new Rib(arrpnt[k * nxny + nx * j + i], arrpnt[k * nxny + nx * j + i + 1]));
        if (j != ny - 1)
            for (int i = 0; i < nx; i++)
                arr.Add(new Rib(arrpnt[k * nxny + nx * j + i], arrpnt[k * nxny + nx * (j + 1) + i]));
    }
    if (k != nz - 1)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++)
                arr.Add(new Rib(arrpnt[k * nxny + nx * j + i], arrpnt[(k + 1) * nxny + nx * j + i]));
}

//return 0;
*/

// Genereting mesh.
ReadMesh(CalculationArea, BordersInfo);
ReadTimeMesh(TimePath);

Mesh3Dim mesh3D = new(NodesX, InfoAboutX, NodesY, InfoAboutY,
                      NodesZ, InfoAboutZ, Elems, Borders);

Mesh2Dim mesh2D = new(NodesR, InfoAboutR, NodesZ, InfoAboutZ,
                      Elems, Math.Sqrt(Math.Pow(mesh3D.nodesX[^1], 2) + Math.Pow(mesh3D.nodesY[^1], 2)));

var timeMesh = GenerateTimeMesh(Time.Item1, Time.Item2, tn, tk);

// Main process of 2-dim task.
mesh2D.SetBorders(mesh3D.borders);
ConstructMesh(ref mesh2D);
FEM2D myFEM2D = new(mesh2D, timeMesh);
if (isSolving2DimTask)
{
    myFEM2D.SetSolver(new LU_LOS());
    myFEM2D.Solve();
    myFEM2D.GenerateVectorEphi();
    myFEM2D.WriteData(AnswerPath);
    myFEM2D.WriteDiscrepancy(AnswerPath);

    // Post-processor of zero-layer. Drawing A_phi and E_phi.
    int a = Postprocessor.DrawA_phi();
    int b = Postprocessor.DrawE_phi();
    Console.WriteLine($"Drawing A_phi finished with code: {a}\n" +
                      $"Drawing E_phi finished with code: {b}");
}
else
    myFEM2D.ReadAnswer(AnswerPath);

ConstructMesh(ref mesh3D);
FEM3D myFEM3D = new(mesh3D,  timeMesh);
myFEM3D.ConvertResultTo3Dim(myFEM2D);
myFEM3D.Layers = ReadFields(LayersArea);
//myFEM3D.CommitFields();
myFEM3D.ConstructMatrixAndVector();
myFEM3D.SetSolver(new LOS());


static bool Checker(string Answer)
{
    if (CountFilesAmount(Answer + "A_phi/Answer/") == CountFilesAmount(Answer + "E_phi/Answer/"))
        while (true)
        {
            string? ans;
            Console.WriteLine("Detected answer for 2-dim task. Would you like to solve 2-dim task again? [Y/n]");
            ans = Console.ReadLine()?.ToLower();
            if (ans == "y" || ans == "yes" || ans == "да" || ans == "д")
                return true;
            else if (ans == "n" || ans == "no" || ans == "нет" || ans == "н")
                return false;
            Console.WriteLine("Unexpected answer!");
        }
    return true;
}