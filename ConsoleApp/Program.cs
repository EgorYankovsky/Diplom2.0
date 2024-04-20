using Project;
using System.Globalization;
using Dumpify;
using Solver;
using Manager;
using Processor;
using Solution;
using MathObjects;
using Grid;
using DataStructs;
using static Grid.MeshReader;
using static Grid.MeshGenerator;

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

bool isUsing2dimAnswer = Checker();

// Pre-processor. Clearing output folders.
if (!isUsing2dimAnswer)
    FolderManager.ClearFolders(new List<string> {SubtotalsPath,
                                                 PicturesPath + "/E_phi/",
                                                 PicturesPath + "/A_phi/",
                                                 AnswerPath + "/E_phi/Discrepancy/",
                                                 AnswerPath + "/E_phi/Answer/",
                                                 AnswerPath + "/A_phi/Discrepancy/",
                                                 AnswerPath + "/A_phi/Answer/"});


// Genereting mesh.
ReadMesh(CalculationArea, BordersInfo);
ReadTimeMesh(TimePath);

Mesh3Dim mesh3D = new(NodesX, InfoAboutX, NodesY, InfoAboutY,
                      NodesZ, InfoAboutZ, Elems, Borders);

Mesh2Dim mesh2D = new(NodesR, InfoAboutR, NodesZ, InfoAboutZ,
                      Elems, Math.Sqrt(Math.Pow(mesh3D.nodesX[^1], 2) + Math.Pow(mesh3D.nodesY[^1], 2)));
mesh2D.SetBorders(mesh3D.borders);

ConstructMesh(ref mesh2D);
var timeMesh = GenerateTimeMesh(Time.Item1, Time.Item2, tn, tk);


// Main process of 2-dim task.
FEM2D myFEM2D = new(mesh2D, timeMesh);
myFEM2D.SetSolver(new LOS());
myFEM2D.Solve();
myFEM2D.GenerateVectorEphi();
myFEM2D.WriteData(AnswerPath);

// Post-processor of zero-layer. Drawing A_phi and E_phi.
int a = Postprocessor.DrawA_phi();
int b = Postprocessor.DrawE_phi();
Console.WriteLine($"Drawing A_phi finished with code: {a}\n" +
                  $"Drawing E_phi finished with code: {b}");


// Converting 2-dim results into 3-dim form
return 0;
MeshGenerator.ConstructMesh(ref mesh3D);

FEM3D myFEM3D = new(myFEM2D);
myFEM3D.ConstructMesh(myFEM2D);
myFEM3D.ConvertResultTo3Dim(myFEM2D);
myFEM3D.GenerateArrays();
myFEM3D.AddField(MeshReader.ReadField(LayersArea));
myFEM3D.CommitFields();
myFEM3D.ConstructMatrixAndVector();

myFEM3D.SetSolver(new LOS());


static bool Checker()
{
    while (true)
    {
        var ans = string.Empty;
        Console.WriteLine("Use result of 2-dim task? [Y/n]");
        ans = Console.ReadLine()?.ToLower();
        if (ans == "y" || ans == "yes" || ans == "да" || ans == "д")
            return true;
        else if (ans == "n" || ans == "no" || ans == "нет" || ans == "н")
            return false;
        Console.WriteLine("Unexpected answer!");
    }
}