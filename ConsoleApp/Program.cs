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
string OutputDirectory = Path.GetFullPath("..\\..\\..\\..\\Data\\Output\\");
string PicturesDirectory = Path.GetFullPath("../../../../Drawer/Pictures/");

//var am = TestClass.A;
//var bm = TestClass.b;
//var solver = new LU_LOS();
//var solver1 = new LOS();
//var solver2 = new BCG();
//var q = solver.Solve(am, bm);
//Console.WriteLine();
//var q1 = solver1.Solve(am, bm);
//Console.WriteLine();
//var q2 = solver2.Solve(am, bm);
//for (int i = 0; i < q.Item1.Size; i++)
//   Console.WriteLine(q.Item1[i]);
//Console.WriteLine();
//for (int i = 0; i < q1.Item1.Size; i++)
//   Console.WriteLine(q1.Item1[i]);
//Console.WriteLine();
//for (int i = 0; i < q2.Item1.Size; i++)
//   Console.WriteLine(q2.Item1[i]);
//return 0;

bool isSolving2DimTask = Checker(OutputDirectory);

// Pre-processor. Clearing output folders.
if (isSolving2DimTask)
    ClearFolders(new List<string> {SubtotalsDirectory + "/2_dim/",
                                   PicturesDirectory + "/E_phi/",
                                   PicturesDirectory + "/A_phi/",
                                   OutputDirectory});

//ClearFolders(new List<string> {SubtotalsDirectory + "3_dim\\",
//                               OutputDirectory + "A_phi\\Answer3D\\"});

// Reading mesh.
ReadMesh(InputDirectory + "WholeMesh.txt");
ReadTimeMesh(InputDirectory + "Time.txt");

// Set recivers
List<Point3D> recivers = [new(1050.0, 0.0, 0.0), new(950.0, 0.0, 0.0),
                          new(707.10, 707.10, 0.0), new(800.0, 0.0, -8000.0)];

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
    myFEM2D.SetSolver(new LU_LOS());
    myFEM2D.Solve();
    myFEM2D.GenerateVectorEphi();
    myFEM2D.WriteData(OutputDirectory);
    myFEM2D.WriteDiscrepancy(OutputDirectory);
    myFEM2D.WritePointsToDraw(OutputDirectory + "ToDraw\\2_dim\\Aphi\\",
                              OutputDirectory + "ToDraw\\2_dim\\Ephi\\");
    myFEM2D.MeasureValuesOnReceivers(recivers, OutputDirectory + "ToDraw\\2_dim\\Receivers\\");

    // Post-processor of normal layer. Drawing A_phi, E_phi and graphics.
    //int a = Postprocessor.DrawA_phi();
    //int b = Postprocessor.DrawE_phi();
    //int c = Postprocessor.DrawGraphics2D();
    //Console.WriteLine($"Drawing A_phi finished with code: {a}\n" +
    //                  $"Drawing E_phi finished with code: {b}\n" + 
    //                  $"Drawing graphics finished with code: {c}\n");
}
else
{
    myFEM2D.ReadAnswer(OutputDirectory);
    Console.WriteLine("2D answer read");
}

myFEM2D.GenerateVectorEphi();
myFEM2D.MeasureValuesOnReceivers(recivers, OutputDirectory + "ToDraw\\2_dim\\Receivers\\");

return 0;

ConstructMesh(ref mesh3D);
Console.WriteLine("3D mesh constructed");
FEM3D myFEM3D = new(mesh3D,  timeMesh);
myFEM3D.ConvertResultTo3Dim(myFEM2D);
Console.WriteLine("2D answer converted to 3D");
//myFEM3D.CheckSolution(recivers);
//myFEM3D.WriteData(OutputDirectory + "A_phi\\Answer3D\\ConvertedTo3D\\");

// Solving first layer: groundwater.
ReadAnomaly(InputDirectory + "Anomalies\\Anomaly1.txt");
Console.WriteLine("Anomaly read");

Mesh3Dim mesh3D_a1 = new(NodesX, InfoAboutX, NodesY, InfoAboutY,
                         NodesZ, InfoAboutZ, Elems, Borders);
mesh3D_a1.CommitAnomalyBorders(FieldBorders);
ConstructMeshAnomaly(ref mesh3D_a1, SubtotalsDirectory + "3_dim\\Anomaly0\\");
Console.WriteLine("Anomaly mesh built");
FEM3D fem3D_a1 = new(mesh3D_a1, timeMesh, myFEM3D, 0);
//fem3D_a1.ReadData(OutputDirectory + "A_phi/Answer3D/AfterField1/");
fem3D_a1.SetSolver(new LOS());
Console.WriteLine("Solving begun");
fem3D_a1.Solve();
Console.WriteLine("Solved");
fem3D_a1.GenerateVectorE();
Console.WriteLine("E generated");
myFEM3D.AddSolution(fem3D_a1);
myFEM3D.WriteDataToDraw2DimSolution(OutputDirectory + "ToDraw\\3_dim\\");
myFEM3D.MeasureValuesOnReceivers(recivers, OutputDirectory + "ToDraw\\3_dim\\Receivers\\");
//myFEM3D.MeasureValuesOnReceivers(recivers, OutputDirectory + "ToDraw\\3_dim\\Receivers\\");
//fem3D_a1.WriteDataAtLine(OutputDirectory + "ToDraw\\3_dim\\Graphics\\");

//fem3D_a1.WriteDataToDraw2DimSolution(OutputDirectory + "ToDraw\\3_dim\\");
//fem3D_a1.WriteData(OutputDirectory + $"A_phi\\Answer3D\\AfterField1\\");
//int a = Postprocessor.DrawA3D(OutputDirectory + "ToDraw/3_dim/A/", PicturesDirectory + "A3d\\");
int e = Postprocessor.DrawE3D(OutputDirectory + "ToDraw/3_dim/E/", PicturesDirectory + "E3d\\");
Console.WriteLine($"{e}");
//myFEM3D.AddSolution(fem3D_a1);
//myFEM3D.WriteDataToDraw2DimSolution(OutputDirectory + "ToDraw\\3_dim\\");
//Postprocessor.DrawA3D(OutputDirectory + "ToDraw/3_dim/A/", PicturesDirectory + "A3d\\");
//Postprocessor.DrawE3D(OutputDirectory + "ToDraw/3_dim/E/", PicturesDirectory + "E3d\\");

return 0;
// Solving second layer: oil.
ReadAnomaly(InputDirectory + "Anomalies\\Anomaly2.txt");
Mesh3Dim mesh3D_a2 = new(NodesX, InfoAboutX, NodesY, InfoAboutY,
                         NodesZ, InfoAboutZ, Elems, Borders);
mesh3D_a2.CommitAnomalyBorders(FieldBorders);
ConstructMeshAnomaly(ref mesh3D_a2, SubtotalsDirectory + "3_dim\\Anomaly1\\");
FEM3D fem3D_a2 = new(mesh3D_a2, timeMesh, myFEM3D, 1);
fem3D_a2.ReadData(OutputDirectory + "A_phi/Answer3D/AfterField2/");
//fem3D_a2.SetSolver(new LOS());
//fem3D_a2.Solve();
fem3D_a2.GenerateVectorE();
myFEM3D.AddSolution(fem3D_a2);
//myFEM3D.WriteDataToDraw2DimSolution(OutputDirectory + "ToDraw\\3_dim\\");
myFEM3D.MeasureValuesOnReceivers(recivers, OutputDirectory + "ToDraw\\3_dim\\Receivers\\");
//myFEM3D.WriteDataToDraw2DimSolution(OutputDirectory + "ToDraw\\3_dim\\");
//fem3D_a2.WriteDataAtLine(OutputDirectory + "ToDraw\\3_dim\\Graphics\\");

//fem3D_a2.WriteDataToDraw2DimSolution(OutputDirectory + "ToDraw\\3_dim\\");
//fem3D_a2.WriteData(OutputDirectory + $"A_phi\\Answer3D\\AfterField2\\");


//myFEM3D.AddSolution(fem3D_a1);

return 0;
//myFEM3D.AddSolution(fem3D_a1);
myFEM3D.GenerateVectorE();
myFEM3D.WriteData(OutputDirectory + $"A_phi\\Answer3D\\AfterField1\\");



return 0;

// Main process of 3D task.
List<Layer> Layers = [new Layer(0.0, 0.0, 0.0, 0.0)];

//ConstructMeshAnomaly(ref mesh3D_a1, SubtotalsDirectory + "3_dim\\Field0\\");
//fem3D_a1.SetSolver(new LU_LOS());
//fem3D_a1.Solve();
//
//myFEM3D.AddSolution(fem3D_a1);
//myFEM3D.GenerateVectorE();
//myFEM3D.WriteData(OutputDirectory + $"A_phi\\Answer3D\\AfterField1\\");
//return 0;

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