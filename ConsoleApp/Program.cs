using Project;
using System.Globalization;
using Solver;
using Processor;
using Grid;
using static Grid.MeshReader;
using static Grid.MeshGenerator;
using static Manager.FolderManager;
using DataStructs;

CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;

string InputDirectory = Path.GetFullPath("../../../../Data/Input/");
string SubtotalsDirectory = Path.GetFullPath("../../../../Data/Subtotals/");
string OutputDirectory = Path.GetFullPath("../../../../Data/Output/");
string PicturesDirectory = Path.GetFullPath("../../../../Drawer/Pictures/");

bool isSolving2DimTask = Checker(OutputDirectory);

// Pre-processor. Clearing output folders.
if (isSolving2DimTask)
    ClearFolders(new List<string> {SubtotalsDirectory + "/2_dim/",
                                   PicturesDirectory + "/E_phi/",
                                   PicturesDirectory + "/A_phi/",
                                   OutputDirectory});

ClearFolders(new List<string> {SubtotalsDirectory + "3_dim\\",
                               OutputDirectory + "A_phi\\Answer3D\\"});

// Reading mesh.
ReadMesh(InputDirectory + "WholeMesh.txt");
ReadTimeMesh(InputDirectory + "Time.txt");

// Set recivers
List<Point3D> recivers = [new(2500.0, 0.0, -100.0), new(2500.0, 0.0, -200.0),
                          new(10.0, 0.0, -700.0), new(1000.0, 0.0, -1250.0)];

Mesh3Dim mesh3D = new(NodesX, InfoAboutX, NodesY, InfoAboutY,
                      NodesZ, InfoAboutZ, Elems, Borders);

Mesh2Dim mesh2D = new(NodesR, InfoAboutR, NodesZ, InfoAboutZ,
                      Elems, Math.Sqrt(Math.Pow(mesh3D.nodesX[^1], 2) + Math.Pow(mesh3D.nodesY[^1], 2)));
mesh2D.SetBorders(mesh3D.borders);

var timeMesh = GenerateTimeMesh(Time.Item1, Time.Item2, tn, tk);

// Main process of 2-dim task.
ConstructMesh(ref mesh2D);
FEM2D myFEM2D = new(mesh2D, timeMesh);
if (isSolving2DimTask)
{
    myFEM2D.SetSolver(new LOS());
    myFEM2D.Solve();
    myFEM2D.GenerateVectorEphi();
    myFEM2D.WriteData(OutputDirectory);
    myFEM2D.WriteDiscrepancy(OutputDirectory);
    myFEM2D.WritePointsToDraw(OutputDirectory + "ToDraw\\2_dim\\Aphi\\",
                              OutputDirectory + "ToDraw\\2_dim\\Ephi\\");
    myFEM2D.MeasureValuesOnReceivers(recivers, OutputDirectory + "ToDraw\\2_dim\\Receivers\\");

    // Post-processor of normal layer. Drawing A_phi, E_phi and graphics.
    int a = Postprocessor.DrawA_phi();
    int b = Postprocessor.DrawE_phi();
    int c = Postprocessor.DrawGraphics2D();
    Console.WriteLine($"Drawing A_phi finished with code: {a}\n" +
                      $"Drawing E_phi finished with code: {b}\n" + 
                      $"Drawing graphics finished with code: {c}\n");
}
else
    myFEM2D.ReadAnswer(OutputDirectory);

ConstructMesh(ref mesh3D);
FEM3D myFEM3D = new(mesh3D,  timeMesh);
myFEM3D.ConvertResultTo3Dim(myFEM2D);
//myFEM3D.CheckSolution(recivers);
myFEM3D.WriteData(OutputDirectory + "A_phi\\Answer3D\\ConvertedTo3D\\");

// Solving first layer: groundwater.
ReadAnomaly(InputDirectory + "Anomalies\\Anomaly1.txt");

Mesh3Dim mesh3D_a1 = new(NodesX, InfoAboutX, NodesY, InfoAboutY,
                         NodesZ, InfoAboutZ, Elems, Borders);
mesh3D_a1.CommitAnomalyBorders(FieldBorders);
ConstructMeshAnomaly(ref mesh3D_a1, SubtotalsDirectory + "3_dim\\Anomaly0\\");
FEM3D fem3D_a1 = new(mesh3D_a1, timeMesh, myFEM3D, 0);
fem3D_a1.SetSolver(new LU_LOS());
fem3D_a1.Solve();

return 0;
myFEM3D.AddSolution(fem3D_a1);
myFEM3D.GenerateVectorE();
myFEM3D.WriteData(OutputDirectory + $"A_phi\\Answer3D\\AfterField1\\");



return 0;

// Main process of 3D task.
List<Layer> Layers = [new Layer(0.0, 0.0, 0.0, 0.0)];

ConstructMeshAnomaly(ref mesh3D_a1, SubtotalsDirectory + "3_dim\\Field0\\");
fem3D_a1.SetSolver(new LU_LOS());
fem3D_a1.Solve();

myFEM3D.AddSolution(fem3D_a1);
myFEM3D.GenerateVectorE();
myFEM3D.WriteData(OutputDirectory + $"A_phi\\Answer3D\\AfterField1\\");
return 0;

// Solving first layer from -2000 < z < -1500
//ReadMesh(CalculationArea, BordersInfo);
Mesh3Dim mesh3D_l1 = new(NodesX, InfoAboutX, NodesY, InfoAboutY,
                         NodesZ, InfoAboutZ, Elems, Borders);

ConstructMesh(ref mesh3D_l1, Layers[0], 0);
FEM3D fem3D_l1 = new(mesh3D_l1, timeMesh, myFEM3D, 0);
fem3D_l1.ConstructMatrixes();
fem3D_l1.SetSolver(new LU_LOS());
fem3D_l1.Solve();
myFEM3D.AddSolution(fem3D_l1);
myFEM3D.GenerateVectorE();
myFEM3D.WriteData(OutputDirectory + $"A_phi\\Answer3D\\AfterField1\\");



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